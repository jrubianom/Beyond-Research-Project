#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "header.h"

using namespace std;
using namespace mfem;



int main(int argc, char *argv[])
{
   const char *mesh_file = "mesh.msh";
   //Mesh mesh(mesh_file);
   int ref_levels = 1;
   flowParameters FP;
   FP.refine_param = 1;
   for(int order = 1; order < 2; order++){
      for(int l = 0; l < ref_levels; l++){
         FP.refine_param = l;
         PerformExp(FP,order,7000,1e-12);
      }
   }
return 0;
}
