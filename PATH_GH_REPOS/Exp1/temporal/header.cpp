#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "header.h"

int DefaultParser(OptionsParser &args,const char *mesh_file,int order,
                   bool static_cond,bool pa, const char *device_config,
                   bool visualization,bool algebraic_ceed){
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                   " isoparametric space.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
#ifdef MFEM_USE_CEED
   args.AddOption(&algebraic_ceed, "-a", "--algebraic", "-no-a", "--no-algebraic",
                  "Use algebraic Ceed solver");
#endif
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
      cout << "args are not good" << '\n';
   }
   args.PrintOptions(cout);
return 0;
}


void Refine_mesh(Mesh &mesh,int refine_param){
    int dim = mesh.Dimension();
    int ref_levels =
        (int)floor(log(refine_param/mesh.GetNE())/log(2.)/dim);
    for (int l = 0; l < ref_levels; l++)
    {
         mesh.UniformRefinement();
    }
}




void get_laplacian(const GridFunction &X,GridFunction &LX,
                   FiniteElementSpace &VectFes);

void PerformExp(int order,int refinment,
                int num_iterations,double tol){


//Set parameters of solution
    flowParameters FP;
    FP.refine_param = refinment;
    FP.order = order;
    double vel = FP.speed;
    double radius = FP.radius;
    double side = FP.side;

// 1. Parse command-line options.
   const char *mesh_file = "mesh.msh";
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;
   bool algebraic_ceed = false;

   Device device(device_config);
   device.Print();




   // 3. Read the mesh from the given mesh file
   Mesh mesh(mesh_file);
   int dim = mesh.Dimension();

   // 4. Refine the mesh

   Refine_mesh(mesh,FP.refine_param);
   cout << "Order:\t" << order << "\n";
   cout << "Refine Parameter:\t" << FP.refine_param << "\n";

   double hmin,hmax,kappamin,kappamax;
   mesh.GetCharacteristics(hmin,hmax,kappamin,kappamax);
   cout << "Mesh parameters" <<"\n";
   cout << "h_min\t" << hmin << '\n';
   cout << "h_max\t" << hmax << '\n';
   cout << "kappa_min\t" << kappamin << '\n';
   cout << "kappa_max\t" << kappamax << '\n';

// 5. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   bool delete_fec;
   fec = new H1_FECollection(order, dim);
   delete_fec = true;

   FiniteElementSpace fespace(&mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace.GetTrueVSize() << endl;


   //True degrees of freedom
   Array<int> ess_tdof_list;

   cout << "Number of attributes\t"<< mesh.bdr_attributes.Max() <<"\n";

   Array<int> nbc_bdr2(mesh.bdr_attributes.Max());
   Array<int> dbc_bdr3(mesh.bdr_attributes.Max());
   Array<int> nbc_bdr_others(mesh.bdr_attributes.Max());

   nbc_bdr2 = 0; nbc_bdr2[2] = 1;
   dbc_bdr3 = 0; dbc_bdr3[3] = 1; //dbc_bdr3[2] = 1;
   nbc_bdr_others = 1; nbc_bdr_others[2] = 0; nbc_bdr_others[3] = 0;
   fespace.GetEssentialTrueDofs(dbc_bdr3, ess_tdof_list);
   // 7. Set up the linear form b(.)
   LinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(0.0);
   ConstantCoefficient nbc_coeff2(-vel);
   //b.AddDomainIntegrator(new DomainLFIntegrator(zero));
   b.AddBoundaryIntegrator(new BoundaryLFIntegrator(nbc_coeff2), nbc_bdr2);
   //b.AddBoundaryIntegrator(new BoundaryLFIntegrator(zero), nbc_bdr_others);
   b.Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace);
   x = 0.0;
   x.ProjectBdrCoefficient(zero,dbc_bdr3);
   //dbc_bdr3 = 0; dbc_bdr3[2] = 1;
   //x.ProjectCoefficient(one,dbc_bdr3);

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator(one));


   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   a.Assemble();

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   cout << "Size of linear system: " << A->Height() << endl;

   // 11. Solve the linear system A X = B.
#ifndef MFEM_USE_SUITESPARSE
      // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
      GSSmoother M((SparseMatrix&)(*A));
      PCG(*A, M, B, X, 1,num_iterations,tol, 0.0);
#else
      // If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
      UMFPackSolver umf_solver;
      umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      umf_solver.SetOperator(*A);
      umf_solver.Mult(B, X);
#endif


   // 12. Recover the solution as a finite element grid function.
   a.RecoverFEMSolution(X, b, x);

   string nameFile = to_string(FP.refine_param) +"-order-"+ to_string(order);
   string solFile = nameFile + ".gf";
   string meshFile = nameFile + "refined.mesh";
   // 13. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs(meshFile);
   mesh_ofs.precision(8);
   mesh.Print(mesh_ofs);
   ofstream sol_ofs(solFile);
   sol_ofs.precision(8);
   x.Save(sol_ofs);

   // 14. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << x << flush;
   }


//Create Vectorial Space and get the Laplacian of x
   H1_FECollection vectorfec = H1_FECollection(order, mesh.Dimension());
   FiniteElementSpace vectorfes = FiniteElementSpace(&mesh,&vectorfec,mesh.Dimension());

   GridFunction Lx(&fespace),Ly(&fespace);
   Lx = 0.0;
   get_laplacian(x,Lx,vectorfes);

   //Compute error in laplacian
   double area = side*side - M_PI*radius*radius;
   cout << "|Laplace_x |_{L^2}/Area = " << Lx.ComputeL2Error(zero)/area << '\n';



   // 15. Save data in the ParaView format
   ParaViewDataCollection paraview_dc(nameFile, &mesh);
   paraview_dc.SetPrefixPath("/media/sf_Paraview_samples/23-05-2022-Exps/Exp1");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("x",&x);
   paraview_dc.RegisterField("Laplace_x",&Lx);
   paraview_dc.Save();

   // 15. Free the used memory.
   delete fec;

}


void get_laplacian(const GridFunction &X,
                   GridFunction &LX,
                   FiniteElementSpace &VectFes){
   GridFunction V(&VectFes);
   V = 0.0;
   GradientGridFunctionCoefficient grad(&X);
   V.ProjectCoefficient(grad);
   DivergenceGridFunctionCoefficient Laplace(&V);
   LX.ProjectCoefficient(Laplace);
}
