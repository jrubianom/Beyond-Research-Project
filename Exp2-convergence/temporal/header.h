#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"

using namespace std;
using namespace mfem;

struct flowParameters{
    double speed = 1.0;
    double radius = 2.0;
    int refine_param = 10;
    int order = 1;
    double side = 60;
    double hmax;
};

int DefaultParser(OptionsParser &args,const char *mesh_file,int order,
                   bool static_cond,bool pa, const char *device_config,
                  bool visualization,bool algebraic_ceed);

int Refine_mesh(Mesh &mesh, int ref_levels);

void PerformExp(flowParameters &FP,int order = 1,int num_iterations = 1000,
                double tol = 1e-9);
