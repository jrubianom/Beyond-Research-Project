#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "Assemble.h"

using namespace std;
using namespace mfem;

void AssembleLinearForm(FiniteElementSpace &fespace,LinearForm &b){
    b.Assemble();
}

void AssembleBilinearForm(FiniteElementSpace &fespace,BilinearForm &a){
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    a.Assemble();
}

int SolveSystem(BilinearForm &a,LinearForm &b,Array<int> ess_tdof_list,
                 GridFunction &x,flowParameters &FP){
    OperatorPtr A;
    Vector B,X;
    a.FormLinearSystem(ess_tdof_list,x,b,A,X,B);
#ifndef MFEM_USE_SUITESPARSE
      // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
      GSSmoother M((SparseMatrix&)(*A));
      PCG(*A, M, B, X, 1, FP.num_iterations, FP.tol, 0.0);
#endif

   //  Recover the solution as a finite element grid function.
      a.RecoverFEMSolution(X, b, x);

      return A->Height();
}

//Refine the mesh a number of iterations equal to ref_levels
int Refine_mesh(Mesh &mesh,int ref_levels){

    for (int l = 0; l < ref_levels; l++)
    {
         mesh.UniformRefinement();
    }
    return mesh.GetNE();
}
