#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "header.h"

using namespace std;
using namespace mfem;


double vel = 1;
double radius = 2.0,side = 30.0;
double f_exact(const Vector &x);
double neum_exact(const Vector &x);

int main(int argc, char *argv[])
{
   double paras,h_min = radius/10;
// 1. Parse command-line options.
   const char *mesh_file = "mesh.msh";
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
   FP.refine_param = 950;

   // 3. Read the mesh from the given mesh file
   Mesh mesh(mesh_file);
   int dim = mesh.Dimension();

   // 4. Refine the mesh

   Refine_mesh(mesh,FP.refine_param);
   cout << "Refine Parameter:\t" << FP.refine_param << "\n";


   // 5. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec;
   bool delete_fec;
   if (order > 0)
   {
      fec = new H1_FECollection(order, dim);
      delete_fec = true;
   }
   else if (mesh.GetNodes())
   {
      fec = mesh.GetNodes()->OwnFEC();
      delete_fec = false;
      cout << "Using isoparametric FEs: " << fec->Name() << endl;
   }
   else
   {
      fec = new H1_FECollection(order = 1, dim);
      delete_fec = true;
   }
   FiniteElementSpace fespace(&mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace.GetTrueVSize() << endl;
   //mesh.GetCharacteristics(h_min,paras,paras,paras);
   //cout << "Mesh parameter" << h_min << '\n';
   // 6. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   Array<int> ess_tdof_list;

   cout << "Number of attributes\t"<< mesh.bdr_attributes.Max() <<"\n";

   Array<int> dbc_bdr(mesh.bdr_attributes.Max());
   Array<int> nbc_bdr(mesh.bdr_attributes.Max());

   dbc_bdr = 1; dbc_bdr[4] = 0;
   nbc_bdr = 0; nbc_bdr[4] = 1;
   fespace.GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);

   // 7. Set up the linear form b(.)
   LinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   FunctionCoefficient nbd_coeff(neum_exact);
   FunctionCoefficient dbc_coeff(f_exact);
   FunctionCoefficient f_sol(f_exact);
   //b.AddDomainIntegrator(new DomainLFIntegrator(zero));
   b.AddBoundaryIntegrator(new BoundaryLFIntegrator(nbd_coeff), nbc_bdr);
   b.Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace),y(&fespace);
   x = 0.0; y = 0.0;
   x.ProjectBdrCoefficient(dbc_coeff,dbc_bdr);
   //x.ProjectCoefficient(zero,dbc_bdr3);
   //dbc_bdr3 = 0; dbc_bdr3[2] = 1;
   //x.ProjectCoefficient(one,dbc_bdr3);
   y.ProjectCoefficient(f_sol);

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm a(&fespace);
   if (pa) { a.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   a.AddDomainIntegrator(new DiffusionIntegrator(one));


   // 10. Assemble the bilinear form and the corresponding linear system,
   //     applying any necessary transformations such as: eliminating boundary
   //     conditions, applying conforming constraints for non-conforming AMR,
   //     static condensation, etc.
   if (static_cond) { a.EnableStaticCondensation(); }
   a.Assemble();

   OperatorPtr A;
   Vector B, X;
   a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

   cout << "Size of linear system: " << A->Height() << endl;

   // 11. Solve the linear system A X = B.
   if (!pa)
   {
#ifndef MFEM_USE_SUITESPARSE
      // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
      GSSmoother M((SparseMatrix&)(*A));
      PCG(*A, M, B, X, 1, 1000, 1e-9, 0.0);
#else
      // If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
      UMFPackSolver umf_solver;
      umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
      umf_solver.SetOperator(*A);
      umf_solver.Mult(B, X);
#endif
   }
   else
   {
      if (UsesTensorBasis(fespace))
      {
         if (algebraic_ceed)
         {
            ceed::AlgebraicSolver M(a, ess_tdof_list);
            PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
         }
         else
         {
            OperatorJacobiSmoother M(a, ess_tdof_list);
            PCG(*A, M, B, X, 1, 400, 1e-12, 0.0);
         }
      }
      else
      {
         CG(*A, B, X, 1, 400, 1e-12, 0.0);
      }
   }


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


   // 15. Save data in the ParaView format
   string paraViewFile = to_string(FP.refine_param);
   ParaViewDataCollection paraview_dc(paraViewFile, &mesh);
   paraview_dc.SetPrefixPath("/media/sf_Paraview_samples/Exp1");
   paraview_dc.SetLevelsOfDetail(order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("x",&x);
   paraview_dc.RegisterField("exact_sol",&y);
   paraview_dc.Save();

   // 15. Free the used memory.
   if (delete_fec)
   {
      delete fec;
   }
return 0;
}

double f_exact(const Vector &x){
   double r2 = x(0)*x(0) + x(1)*x(1);
   return vel*(1+radius*radius/r2)*x(0);
}

double neum_exact(const Vector &x){
   double tol = 1e-12;
   return 0;
}
