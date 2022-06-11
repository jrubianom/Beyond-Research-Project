#ifndef _FICTFORCE_H
#define _FICTFORCE_H


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "mfem.hpp"
#include "EngineRigidBody.h"

using namespace mfem;
using namespace std;


struct Config{
    double mass0;
    // attributes related with rotational movement
    Vector3d I0 = {1.5, 1.5,0.5}; //Inertia vector
    Vector3d w0 = {1,0,9}; //omega in Non inertial frame
    Vector3d angles0 = {M_PI/6, -M_PI/2, 0}; // Euler angles
    Vector3d Accel0 = {0,0,0}; //Initial linear acceleration
    Vector3d Velocity0 = {0,0,0}; //Initial linear Velocity

    //
    int order = 1;
};


class FictitiousForce:public VectorCoefficient{
    protected:
        RigigidBody Frisbee;
        GridFunction *u = nullptr; // Velocity in Non Inertial Frame
        ParLinearForm *lf_bdr = nullptr;
        PWConstCoefficient XiCoeff;
        double t;

    public:
        FictitiousForce(int dim):VectorCoefficient(dim){};
        virtual ~FictitiousForce(){};
        void InitFrisbee(Config &params_init);
        void Init_Integrals(Config &params_init,ParMesh *pmesh,
                            ParFiniteElementSpace *fespace);

//Una vez sabiendo coom sacar integrales
        void SetVelocity(GridFunction* GfVel){u = GFVel};
        void SetIntegrals(LinearForm* LF){lf_bdr = LF};
        void SetTime(double time){t = time};
        /*void GetForce_and_Torque(){
        //code ...
        //Frisbee.GetTorque(...);
        //Frisbee.GetForce(...);
        }*/
        //void GetForce_and_Torque();
//Update all
        void UpdateFrisbee(double dt, Vector3d Torques, Vector3d Forces);

//Auxiliar function that given the velocity of thhe fluid U(x), get the
//Ficticious Force and store it in the vector FF
        void GetFictForce(Vector X,Vector U,Vector &FF);
        virtual void Eval(Vector &FictForce ,ElementTransformation &T,
                          const IntegrationPoint &ip);
};

#endif
