#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

double func(const Vector &X);
double Divfunc(const Vector &X);
double Chi(const Vector &X){
   double r2 = X(0)*X(0) + X(1)*X(1) + X(2)*X(2);
   if(abs(r2) < 4){
      return 0;
   }
   else{
      return 1;
   }
}


class FF: public Coefficient
{
   private:
      double a;
      double (*Function)(const Vector &);
   public:

      FF(double (*F) (const Vector &))
         :Function(F){};
      virtual ~FF(){};
      virtual double Eval( ElementTransformation &T,
                        const IntegrationPoint &ip)
      {
         double transip[3];
         Vector x(transip, 3);
         T.Transform(ip, x);
         return Function(x);
      }
};

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   const char *mesh_file = "mesh.msh";
   int order = 1;
   bool static_cond = false;
   bool pa = false;
   const char *device_config = "cpu";
   bool visualization = true;
   bool algebraic_ceed = false;
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();
   {
      int ref_levels = 0;
      cout << "ref_levels\t" << ref_levels << endl;
      for (int l = 0; l < ref_levels; l++)
      {
         mesh.UniformRefinement();
      }
   }
   cout << "Finalize ref\n";
   cout << "Number of attr\t" << mesh.bdr_attributes.Max() <<endl;
    cout << "Number of elements\t" << mesh.GetNE() <<endl;
   cout << "Dim\t" << dim << endl;

   FiniteElementCollection *fec;
   bool delete_fec;
   fec = new H1_FECollection(order, dim);
   FiniteElementSpace fespace(&mesh, fec);

   //Characteristic function for sphere bdr
   Vector Xi(mesh.bdr_attributes.Max());
   Xi = 0.0;
   Xi[4] = 1.0;
   //Characteristic function for boundary 4
   PWConstCoefficient XiCoeff = PWConstCoefficient(Xi);

   Array<int> ess_tdof_list, ess_bdr(mesh.bdr_attributes.Max());
   ess_bdr = 1;

   //fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);


   GridFunction x(&fespace);
   x = 0.0;
   ConstantCoefficient one(1.0);
   FunctionCoefficient funcCoeff(func);
   x.ProjectBdrCoefficient(XiCoeff,ess_bdr);

   FF normalF(func);
   LinearForm lf_bdr(&fespace);
   lf_bdr.AddBoundaryIntegrator(new BoundaryLFIntegrator(funcCoeff));
   lf_bdr.Assemble();
   double bdr_int = lf_bdr(x);

   double radius = 2;
   //cout << "Error\t" << abs(bdr_int- 2*M_PI*radius) << endl;
   cout << "Error\t" << abs(bdr_int- 5*M_PI*pow(radius,3)) << endl;

   delete fec;
   //delete XiCoeff;

   return 0;
}


double func(const Vector &X){
   return 2*X(0)*X(0) + 3*X(1)*X(1);
}

double Divfunc(const Vector &X){
   double r2 = X(0)*X(0) + X(1)*X(1) + X(2)*X(2);
   return 3*X(0)/sqrt(r2);
}
