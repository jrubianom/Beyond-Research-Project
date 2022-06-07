//Engine of Rigid Body Motion.
//Given a Torque in the non inertial frame
///get the angular velocity in each iteratio in the inertial frame

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "EngineRigidBody.h"

using namespace std;
using namespace Eigen;

void RigidBody::Init(double mass0, double Ix,double Iy, double Iz,VectorXd angles0,
                     VectorXd w0){
  mass = mass0;
  I(0) = Ix; I(1) = Iy; I(2) = Iz;
  w(0) = w0(0); w(1) = w0(1); w(2) = w0(2);
  double theta0 = angles0(0),phi0 = angles0(1), psi0 = angles0(2);
  q(0) = cos(0.5*theta0)*cos(0.5*(phi0+psi0));
  q(1) = sin(0.5*theta0)*cos(0.5*(phi0-psi0));
  q(2) = sin(0.5*theta0)*sin(0.5*(phi0-psi0));
  q(3) = cos(0.5*theta0)*sin(0.5*(phi0+psi0));
}

void RigidBody::UpdateA(){
  MatrixXd B = MatrixXd(3,3);
  B << q(0)*q(0) + q(1)*q(1)-q(2)*q(2)-q(3)*q(3) , 2*(q(1)*q(2) + q(0)*q(3)), 2*(q(1)*q(3) - q(0)*q(2)),
    2*(q(1)*q(2) - q(0)*q(3)), q(0)*q(0) - q(1)*q(1) + q(2)*q(2) - q(3)*q(3), 2*(q(2)*q(3) + q(0)*q(1)),
    2*(q(1)*q(3) + q(0)*q(2)), 2*(q(2)*q(3) - q(0)*q(1)), q(0)*q(0) - q(1)*q(1) - q(2)*q(2) + q(3)*q(3);
    A = B;
}


void RigidBody::UpdateAinv(){
    Ainv = A.inverse();
}

VectorXd RigidBody::GetVectorNIF(VectorXd &X){
    return A * X;
}

VectorXd RigidBody::GetVectorIF(VectorXd &X){
    return Ainv * X;
}

void RigidBody::GetTorqueNIF(VectorXd &Tif){
    T = GetVectorNIF(Tif);
}

void RigidBody::UpdateRotational(double dt){
  VectorXd q_before = VectorXd(4);
  q_before << q(0),q(1),q(2),q(3);
  VectorXd w_before = VectorXd(3);
  w_before << w(0),w(1),w(2);
  VectorXd wDot = VectorXd(3);
  //Matrix for updating quaternions
  MatrixXd Q_w = MatrixXd(4,4);
  Q_w << 0,-w_before(0),-w_before(1),-w_before(2),
    w_before(0),0,w_before(2),-w_before(1),
    w_before(1),-w_before(2),0,w_before(0),
    w_before(2),w_before(1),-w_before(0),0;

  //update quaternions
  q += dt*0.5* Q_w * q_before;

  //Update omega  dot with euler eqs
  wDot(0) = (T(0)/I(0)) + w_before(1)*w_before(2)*(I(1)-I(2))/I(0);
  wDot(1) = (T(1)/I(1)) + w_before(0)*w_before(2)*(I(2)-I(0))/I(1);
  wDot(2) = (T(2)/I(2)) + w_before(0)*w_before(1)*(I(0)-I(1))/I(2);

  //Update omega with Euler Eqs
  w += dt*wDot;

  //update vector Omega and OmegaDot = DDomega/Dt
  //mesured by the Inertial Frame
  UpdateA();
  UpdateAinv();
  Omega = Ainv * w;
  OmegaDot = Ainv * wDot;
}

void RigidBody::GetForceIF(VectorXd Fif){
  Force =  Fif;
}

void RigidBody::UpdateTranslational(double dt){
  Accel = Force/mass;
  Velocity += dt*Accel;
}


void RigidBody::UpdateAll(double dt,VectorXd &Tif,VectorXd Fif){
  GetForceIF(Fif);
  UpdateTranslational(dt);
  GetTorqueNIF(Tif);
  UpdateRotational(dt);
}

//Get functions

VectorXd RigidBody::getW(){
    return w;
}

VectorXd RigidBody::getOmega(){
  return Omega;
}

VectorXd RigidBody::getOmegaDot(){
  return OmegaDot;
}

VectorXd RigidBody::getAccel(){
  return Accel;
}

VectorXd RigidBody::getVelocity(){
  return Velocity;
}

double RigidBody::getTheta(){
    return acos((q(0)*q(0)-q(1)*q(1)-q(2)*q(2)+q(3)*q(3))/(q(0)*q(0)+q(1)*q(1)+q(2)*q(2)+q(3)*q(3)));
}
