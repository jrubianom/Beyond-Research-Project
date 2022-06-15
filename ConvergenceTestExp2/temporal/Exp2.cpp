#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "Simulate.h"

using namespace std;
using namespace mfem;


int main(int argc, char *argv[])
{
   const char *mesh_file = "mesh.msh";
   Mesh mesh(mesh_file);
   int ref_levels = 5;
   int max_order = 5;
   flowParameters FP;
   FP.refine_param = 1;
   FP.num_iterations = 7000;
   FP.tol = 1e-12;
   FP.ptrmesh(&mesh);
   for(int order = 1; order < max_order; order++){
      FP.order = order;
      for(int l = 0; l < ref_levels; l++){
         PerformExp(FP,order,7000,1e-12);
         FP.start = false;
      }
   }
   delete FP.mesh;
return 0;
}
