#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "header.h"

using namespace std;
using namespace mfem;


double vel = 1.0;
double radius = 2.0;
double dcoef = 0.0;
double f_exact(const Vector &x);
void get_laplacian(const GridFunction &X,GridFunction &LX,
                   FiniteElementSpace &VectFes);

int main(int argc, char *argv[])
{
   double paras,h_min = radius/10;
// 1. Parse command-line options.
   const char *mesh_file = "../mesh.msh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;
   bool algebraic_ceed = false;

   OptionsParser args(argc, argv);
   int dp = DefaultParser(args,mesh_file,order,static_cond,pa,device_config,
                          visualization,algebraic_ceed);

   Device device(device_config);
   device.Print();


   //Set parameters of solution
   flowParameters FP;
   FP.refine_param = 100001;

   // 3. Read the mesh from the given mesh file
   Mesh mesh(mesh_file);
   int dim = mesh.Dimension();

   // 4. Refine the mesh

   Refine_mesh(mesh,FP.refine_param);
   cout << "Refine param\t"<< FP.refine_param << '\n';

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


   //Define and set the linear form
   LinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   ConstantCoefficient coef(dcoef);
   FunctionCoefficient f_sol(f_exact);
   b.AddDomainIntegrator(new DomainLFIntegrator(coef));
   b.Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace),y(&fespace);
   x = 0.0;y = 0.0;
   x.ProjectBdrCoefficient(f_sol,ess_bdr);
   y.ProjectCoefficient(f_sol);

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
      PCG(*A, M, B, X, 1, 1000, 1e-9, 0.0);
#endif

   // 12. Recover the solution as a finite element grid function.
   a.RecoverFEMSolution(X, b, x);

//12.5 Compute L2 error
   const IntegrationRule *irs[Geometry::NumGeom];
   int order_quad = max(2, 2*order+1);
    for (int ii=0; ii < Geometry::NumGeom; ++ii)
        irs[ii] = &(IntRules.Get(ii, order_quad));
   cout << "|F_h -F|_{L^2}/|F| = " << x.ComputeL2Error(f_sol)/ComputeLpNorm(2,f_sol,mesh,irs) << '\n';

   // 13. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh.Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
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
   Lx = 0.0; Ly = 0.0;
   get_laplacian(x,Lx,vectorfes);
   get_laplacian(y,Ly,vectorfes);


   // 15. Save data in the ParaView format
   string paraViewFile = to_string(FP.refine_param);
   ParaViewDataCollection paraview_dc(paraViewFile, &mesh);
   paraview_dc.SetPrefixPath("/media/sf_Paraview_samples/Exp2");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("x",&x);
   paraview_dc.RegisterField("exact_sol",&y);
   paraview_dc.RegisterField("Laplace_x",&Lx);
   paraview_dc.RegisterField("Laplace_y",&Ly);
   paraview_dc.Save();

   // 15. Free the used memory.

   delete fec;


   return 0;
}

double f_exact(const Vector &x){
   double r2 = x(0)*x(0) + x(1)*x(1);
   return vel*(1+radius*radius/r2)*x(0);//-dcoef/2.0*x(0)*x(0)+x(1)*x(1);//vel*(1+a*a/r2)*x(0);
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
