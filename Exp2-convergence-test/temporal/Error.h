#ifndef _ERROR_H
#define _ERROR_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"
#include "header.h"

using namespace std;
using namespace mfem;

//Get the gradient given X as a Gridfunction
GridFunction get_gradient(const GridFunction &X,
                          FiniteElementSpace &VectFes);

//Function that gives the grad of the exact solution
void grad_exact(const Vector &X, Vector &Grad);

//exact solution
double f_exact(const Vector &x);

double get_error(GridFunction &GF,VectorFunctionCoefficient &VF,int order);

#endif
