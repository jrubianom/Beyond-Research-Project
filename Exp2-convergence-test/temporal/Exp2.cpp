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
   int ref_levels = 20;
   flowParameters FP;
   FP.refine_param = 1;
   FP.num_iterations = 7000;
   FP.tol = 1e-12;
   for(int order = 1; order < 5; order++){
      FP.order = order;
      for(int l = 0; l < ref_levels; l++){
         FP.refine_param = 1;
         PerformExp(FP,order,7000,1e-12);
         FP.start = false;
      }
   }
return 0;
}
