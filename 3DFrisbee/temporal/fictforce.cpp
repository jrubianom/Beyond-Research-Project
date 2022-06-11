#include "fictforce.h"




//Initialice the frisbee
void FictitiousForce::InitFrisbee(Config &params_init){
    Frisbee.Init(params_init.mass0,params_init.I0(0),
                 params_init.I0(1),params_init.I0(2),
                 params_init.angles0,params_init.w0,params_init.Accel0,
                 params_init.Velocity0);
}


void FictitiousForce::UpdateFrisbee(double dt, Vector3d Torques,
                                    Vector3d Forces){
    Frisbee.UpdateAll(dt,Torques,Forces);
}

//Auxiliar Function that given the velocity of thhe fluid U(x), get the
//Ficticious Force and store it in the vector FF
void GetFictForce(Vector X,Vector U,Vector &FF){
    Vector3d Force_aux;
    Vector3d U_aux = {U(0),U(1),U(2)};
    Vector3d X_aux = {X(0),X(1),X(2)};
    Force_aux = 2*Frisbee.Omega.cross(U_aux) +
        Frisbee.Omega.cross(Frisbee.Omega.cross(X_aux)) +
        Frisbee.Accel + Frisbee.OmegaDot.cross(X_aux);
    FF(0) = Force_aux(0); FF(1) = Force_aux(1); FF(2) = Force_aux(2);
}

void FictitiousForce::Eval(Vector &FictForce ,ElementTransformation &T,
                           const IntegrationPoint &ip){
    Vector X;
    Vector U;

    U.SetSize(GetVDim()); U(0) = 0; U(1) = 0; U(2) = 0;
    FictForce(GetVDim());
    FictForce(0) = 0; FictForce(1) = 0; FictForce(2) = 0;

    T.Transform(ip,X);
    u->GetVectorValue(T,ip,U);
    GetFictForce(X,U, FictForce);
}

//Init Characteristic Coefficient Xi and linearForm
void FictitiousForce::Init_Integrals(Config &params_init,ParMesh *pmesh,
                                     ParFiniteElementSpace *fespace){
    //Characteristic function for bdr
   Vector Xi(pmesh->bdr_attributes.Max());
   Xi = 0.0;
   Xi[4] = 1.0;
   XiCoeff = PWConstCoefficient(Xi);

   Array<int> ess_tdof_list, ess_bdr(pmesh->bdr_attributes.Max());
   ess_bdr = 1;
   lf_bdr = new ParLinearForm(fespace);
   lf_bdr->AddBoundaryIntegrator(new BoundaryLFIntegrator(XiCoeff));
   lf_bdr->Assemble();
   //I nedd to think about get the vector Fictitiuos force
}


/* void FictitiousForce::GetForce_and_Torque(...){
        //code using the  lf_bdr ...
        Frisbee.GetTorque(...);
        Frisbee.GetForce(...);
        }
*/
