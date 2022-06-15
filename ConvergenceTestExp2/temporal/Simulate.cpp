#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "Simulate.h"

using namespace std;
using namespace mfem;


void PerformExp(flowParameters &FP,int order,
                int num_iterations,double tol){
   FP.order = order;
   FP.num_iterations = num_iterations;
   FP.tol = tol;

// Parse command-line options.
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;
   bool algebraic_ceed = false;
   Device device(device_config);
   device.Print();


   //Read the mesh from the given mesh file
   int dim = FP.mesh->Dimension();

   //Refine the mesh

   int NE = Refine_mesh(FP.mesh,FP.refine_param); //Number of elements
   FP.NE = NE;
   FiniteElementCollection *fec;

   fec = new H1_FECollection(order, dim);

   FiniteElementSpace fespace(FP.mesh, fec);
   cout << "Number of finite element unknowns: "
        << fespace.GetTrueVSize() << endl;

   //Set the ess_bdr attributes
   Array<int> ess_tdof_list;
   Array<int> ess_bdr(FP.mesh->bdr_attributes.Max());
   ess_bdr = 1;
   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   cout << "Number of attributes:\t" << FP.mesh->bdr_attributes.Max()<<'\n';

   double hmin,hmax,kappamin,kappamax;
   FP.mesh->GetCharacteristics(hmin,hmax,kappamin,kappamax);
   FP.hmin = hmin; FP.hmax = hmax;
   FP.kappamin = kappamin; FP.kappamax = kappamax;

   //Define and set the linear form
   LinearForm b(&fespace);
   AssembleLinearForm(fespace,b);

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace),exact_sol(&fespace);
   x = 0.0;exact_sol = 0.0;
   FunctionCoefficient f_sol(f_exact);
   x.ProjectBdrCoefficient(f_sol,ess_bdr);
   exact_sol.ProjectCoefficient(f_sol);

   // 9. Set up the bilinear form a(.,.)
   BilinearForm a(&fespace);
   AssembleBilinearForm(fespace,a);

   // Solve the linear system and recover the solution
   FP.sizeSys = SolveSystem(a,b,ess_tdof_list,x,FP);



//Create Vectorial Space and get the Laplacian of x
   H1_FECollection vectorfec = H1_FECollection(order, FP.mesh->Dimension());
   FiniteElementSpace vectorfes = FiniteElementSpace(FP.mesh,&vectorfec,FP.mesh->Dimension());
   GridFunction Gradx = get_gradient(x,vectorfes);
   GridFunction Grad_exact_sol = get_gradient(exact_sol,vectorfes);
   VectorFunctionCoefficient vec_grad_exact(FP.mesh->Dimension(),grad_exact);


//14.5 Compute L2 grad error
   FP.grad_error = get_error(Gradx,vec_grad_exact,FP.order);
   FP.SaveConvergenceInfo();

   // Save the refined mesh and the solution.
    SaveParams(FP,x,FP.mesh);
    //Send the solution by socket to a GLVis server.
    Visualize(visualization,FP.mesh,x);
    //  Save data in the ParaView format
    SaveInParaView(FP,FP.mesh,x,exact_sol,Gradx,Grad_exact_sol);


   // 15. Free the used memory.

   delete fec;

}
