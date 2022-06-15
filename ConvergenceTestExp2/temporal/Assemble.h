#ifndef _ASSEMBLE_H
#define _ASSEMBLE_H

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "header.h"


void AssembleLinearForm(FiniteElementSpace &fespace,LinearForm &b);

void AssembleBilinearForm(FiniteElementSpace &fespace,BilinearForm &a);

int SolveSystem(BilinearForm &a,LinearForm &b,Array<int> ess_tdof_list,
                 GridFunction &x,flowParameters &FP);

int Refine_mesh(Mesh *mesh,int ref_levels);

#endif
