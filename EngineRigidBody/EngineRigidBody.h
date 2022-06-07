//Engine of Rigid Body Motion.
//Given a Torque get the angular velocity in each iteration
#ifndef _ENGINERIGIDBODY_H
#define _ENGINERIGIDBODY_H

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;
class RigidBody{
        private:

                double mass;
                // attributes related with rotational movement
                VectorXd I = VectorXd(3); //Inertia vector
                VectorXd q = VectorXd(4); //quaternions
                VectorXd w = VectorXd(3); //omega in Non inertial frame
                VectorXd Omega = VectorXd(3); // omega in innertial frame
                VectorXd OmegaDot = VectorXd(3); // D Omega / Dt in inertial frame
                VectorXd T = VectorXd(3); //Torque
                MatrixXd A = MatrixXd(3,3); //Change of basis from IF to NIF
                MatrixXd Ainv  = MatrixXd(3,3); //Change of basis from NIF to I

                // attributes related with translational movement
                VectorXd Force = VectorXd(3); //Force mmesured in Inertial Frame
                VectorXd Accel = VectorXd(3); // Aceleration mmeaasured in IF
                VectorXd Velocity = VectorXd(3); // Velocity measured in IF



        public:
                //angles are euler angles
                void Init(double mass0, double Ix,double Iy, double Iz,VectorXd angles0,
                          VectorXd w0);
                void UpdateA();
                void UpdateAinv();
                void Rotate(double dt);
                VectorXd GetVectorNIF(VectorXd &X); // get vector in the non inertial frame
                VectorXd GetVectorIF(VectorXd &X);//get vector in inertial frame
                void GetTorqueNIF(VectorXd &Tif); ///get the torque in non-inertial-frame
                void GetForceIF(VectorXd Fif); //pass the net force in IF
                void UpdateRotational(double dt);
                void UpdateTranslational(double dt);
                void UpdateAll(double dt,VectorXd &Tif,VectorXd Fif);
                VectorXd getW();
                VectorXd getOmega();
                VectorXd getOmegaDot();
                VectorXd getAccel();
                VectorXd getVelocity();
                double getTheta();
};

#endif
