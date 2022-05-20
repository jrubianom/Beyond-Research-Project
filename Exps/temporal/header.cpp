#include "header.h"

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


void Refine_mesh(Mesh &mesh,int refine_param){
    int dim = mesh.Dimension();
    int ref_levels =
        (int)floor(log(refine_param/mesh.GetNE())/log(2.)/dim);
    for (int l = 0; l < ref_levels; l++)
    {
         mesh.UniformRefinement();
    }
}
