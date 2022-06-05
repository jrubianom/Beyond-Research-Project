#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "header.h"


double dcoef = 0.0;

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

//Refine the mesh a number of iterations equal to ref_levels
int Refine_mesh(Mesh &mesh,int ref_levels){

    for (int l = 0; l < ref_levels; l++)
    {
         mesh.UniformRefinement();
    }
    return mesh.GetNE();
}



//Get the gradient given X
GridFunction get_gradient(const GridFunction &X,
                   FiniteElementSpace &VectFes);

double f_exact(const Vector &x);

void grad_exact(const Vector &X, Vector &Grad);

flowParameters PP;
double vel = PP.speed;
double radius = PP.radius;
double side = PP.side;

void PerformExp(flowParameters &FP,int order,
                int num_iterations,double tol){

   FP.order = order;

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

   int NE = Refine_mesh(mesh,FP.refine_param); //Number of elements


   FiniteElementCollection *fec;

   fec = new H1_FECollection(order, dim);

   FiniteElementSpace fespace(&mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace.GetTrueVSize() << endl;

   //Set the ess_bdr attributes
   Array<int> ess_tdof_list;
   Array<int> ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 1;
   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   cout << "Number of attributes:\t" << mesh.bdr_attributes.Max()<<'\n';

   double hmin,hmax,kappamin,kappamax;
   mesh.GetCharacteristics(hmin,hmax,kappamin,kappamax);
   cout << "Mesh parameters" <<"\n";
   cout << "h_min\t" << "h_max\t" << "kappa_min\t" <<  "kappa_max" << '\n';
   cout << hmin << '\t' << hmax << '\t' << kappamin << '\t' << kappamax << '\n';
   cout << "Number of elements\t" << NE << '\n';
   FP.hmax = hmax;

   //Define and set the linear form
   LinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   ConstantCoefficient zero(0.0);
   FunctionCoefficient f_sol(f_exact);
   //b.AddDomainIntegrator(new DomainLFIntegrator(coef));
   b.Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace),exact_sol(&fespace);
   x = 0.0;exact_sol = 0.0;
   x.ProjectBdrCoefficient(f_sol,ess_bdr);
   exact_sol.ProjectCoefficient(f_sol);

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new DiffusionIntegrator(one));
   a.Assemble();

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   cout << "Size of linear system: " << A->Height() << endl;

   // 11. Solve the linear system A X = B.

#ifndef MFEM_USE_SUITESPARSE
      // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
      GSSmoother M((SparseMatrix&)(*A));
      PCG(*A, M, B, X, 1, num_iterations, tol, 0.0);
#endif

   // 12. Recover the solution as a finite element grid function.
   a.RecoverFEMSolution(X, b, x);

   // 13. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
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

   GridFunction Gradx = get_gradient(x,vectorfes);
   GridFunction Grad_exact_sol = get_gradient(exact_sol,vectorfes);

   VectorFunctionCoefficient vec_grad_exact(mesh.Dimension(),grad_exact);


//14.5 Compute L2 error
   const IntegrationRule *irs[Geometry::NumGeom];
   int order_quad = max(2, 2*order+1);
    for (int ii=0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));
    cout << "|Grad F_h -F Grad |_{L^2} = " << Gradx.ComputeL2Error(vec_grad_exact,irs) << '\n';

   // 15. Save data in the ParaView format
   ParaViewDataCollection paraview_dc(nameFile, &mesh);
   paraview_dc.SetPrefixPath("/media/sf_Paraview_samples/23-05-2022-Exps/Exp2");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("x",&x);
   paraview_dc.RegisterField("exact_sol",&exact_sol);
   paraview_dc.RegisterField("Grad_x",&Gradx);
   paraview_dc.RegisterField("Grad_exact_sol",&Grad_exact_sol);
   paraview_dc.Save();


   // 15. Free the used memory.

   delete fec;

}


GridFunction get_gradient(const GridFunction &X,
                   FiniteElementSpace &VectFes){
   GridFunction Gradx(&VectFes);
   Gradx = 0.0;
   GradientGridFunctionCoefficient grad(&X);
   Gradx.ProjectCoefficient(grad);
   return Gradx;
}


double f_exact(const Vector &x){
   double r2 = x(0)*x(0) + x(1)*x(1);
   return vel*(1+radius*radius/r2)*x(0);
}

void grad_exact(const Vector &X, Vector &Grad){
   double r2 = X(0)*X(0) + X(1)*X(1);
   Grad(0) = vel*(1 + pow(radius/r2,2)*(X(1)*X(1)-X(0)*X(0)));
   Grad(1) = -2*vel*pow(radius/r2,2)*X(0)*X(1);
}
