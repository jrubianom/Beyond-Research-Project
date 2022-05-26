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
   int k = 1000;
   int refinments[] = {k,10*k,100*k,150*k,500*k};
   int len_Ref = sizeof(refinments)/sizeof(refinments[0]);
   for(int order = 1; order < 3; order++){
      for(int l = 0; l < len_Ref; l++){
         PerformExp(order,refinments[l],7000,1e-12);
      }
   }
return 0;
}
