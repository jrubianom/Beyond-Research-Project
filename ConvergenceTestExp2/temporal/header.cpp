#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include "header.h"


double dcoef = 0.0;

int DefaultParser(OptionsParser &args,const char *mesh_file,int order,
                   bool static_cond,bool pa, const char *device_config,
                   bool visualization,bool algebraic_ceed){
    args.AddOption(&mesh_file, "-m", "--mesh",
                   "Mesh file to use.");
    args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                   " isoparametric space.");
    args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
    args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
    args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
#ifdef MFEM_USE_CEED
   args.AddOption(&algebraic_ceed, "-a", "--algebraic", "-no-a", "--no-algebraic",
                  "Use algebraic Ceed solver");
#endif
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
      cout << "args are not good" << '\n';
   }
   args.PrintOptions(cout);
return 0;
}

//Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g file.gf".
void SaveParams(flowParameters &FP, GridFunction &x,Mesh *mesh){
   ofstream mesh_ofs(FP.meshFile());
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   ofstream sol_ofs(FP.SolFile());
   sol_ofs.precision(8);
   x.Save(sol_ofs);
}

void Visualize(bool visualization,Mesh *mesh, GridFunction &x){
   // 14. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << mesh << x << flush;
   }
}

//Save data in the ParaView format
void SaveInParaView(flowParameters &FP,Mesh *mesh,
                    GridFunction &x,GridFunction &exact_sol,
                    GridFunction &Gradx,GridFunction &Grad_exact_sol){

   ParaViewDataCollection paraview_dc(FP.nameFile(), mesh);
   paraview_dc.SetPrefixPath("DiferentOrdersPV");
   paraview_dc.SetLevelsOfDetail(FP.order);
   paraview_dc.SetCycle(0);
   paraview_dc.SetDataFormat(VTKFormat::BINARY);
   paraview_dc.SetHighOrderOutput(true);
   paraview_dc.SetTime(0.0); // set the time
   paraview_dc.RegisterField("x",&x);
   paraview_dc.RegisterField("exact_sol",&exact_sol);
   paraview_dc.RegisterField("Grad_x",&Gradx);
   paraview_dc.RegisterField("Grad_exact_sol",&Grad_exact_sol);
   paraview_dc.Save();

}



void flowParameters::SaveConvergenceInfo(){
   int sw = 16,precision = 10;
   ofstream results;
   results.precision(precision);
   results.open("ResultsConvergence.txt", ios::app);
   if(start){
      results << left << setw(sw) << "order" << setw(sw) << "Size"<< setw(sw) <<
      "NE" << setw(sw) << "GradError" << setw(sw) <<
      "hmin" << setw(sw) <<  " hmax " << setw(sw) <<
      " kmin" << setw(sw) << " kmax \n";
   }
   results << left <<  setw(sw) << order << setw(sw) << sizeSys << setw(sw) <<
      NE << setw(sw) << grad_error << setw(sw) << hmin << setw(sw) <<
      hmax << setw(sw) << kappamin << setw(sw) << kappamax << setw(sw) << '\n';
   results.close();
}
