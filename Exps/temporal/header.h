#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct flowParameters{
    double speed = 100;
    double radius = 2.0;
    int refine_param = 70000;

};

int DefaultParser(OptionsParser &args,const char *mesh_file,int order,
                   bool static_cond,bool pa, const char *device_config,
                  bool visualization,bool algebraic_ceed);

void Refine_mesh(Mesh &mesh,int refine_param);
