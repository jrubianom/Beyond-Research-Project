#include "navier_solver.hpp"
#include <fstream>
#include <string>
#include "fictforce.h"

using namespace mfem;
using namespace navier;

//Configuration of initial and other params

// Initial and boundary functions
void Initial_Velocity(const Vector &x, double t, Vector &u);
void Vel_Boundary_Condition(const Vector &x, double t, Vector &u);
double Press_Boundary_Condition(const Vector &x, double t);
void Gravity(const Vector &x, double t, Vector &a);

//Main function
int main(int argc, char *argv[])
{   
    //Init MPI
    MPI_Session mpi(argc, argv);

    //Parameters
    double t = 0.0;
    int vis_print = 0;
    bool last_step = false;

    //Load Mesh (Pointer to Delete it After Parallel Mesh is Created)
    Mesh *mesh = new Mesh("mesh.msh");
    mesh->EnsureNodes();
    int dim = mesh->Dimension();

    //Refine Serial Mesh
    for (int i = 0; i < Parameters.serial_refinements; ++i)
        mesh->UniformRefinement();

    //Make Parallel Mesh
    ParMesh pmesh = ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    
    //Refine Parallel Mesh
    for (int ii = 0; ii < Parameters.parallel_refinements; ii++)
        pmesh.UniformRefinement();

    //Create H1 Element Collections
    H1_FECollection vfec = H1_FECollection(Parameters.order,
                                           pmesh.Dimension());
    H1_FECollection pfec = H1_FECollection(Parameters.order);

    //Create Finite Element Spaces
    ParFiniteElementSpace vfes = ParFiniteElementSpace(&pmesh, &vfec,
                                                       pmesh.Dimension()); //Vector Space
    ParFiniteElementSpace pfes = ParFiniteElementSpace(&pmesh, &pfec);                   //Scalar Space

    //Define Navier Solver
    NavierSolver flowsolver(&pmesh, Parameters.order, Parameters.kinvis);
    flowsolver.EnablePA(true);
    flowsolver.EnableNI(true);
    flowsolver.EnableVerbose(true);

    //Set Velocity Initial Conditions
    ParGridFunction *u_0 = flowsolver.GetCurrentVelocity();
    VectorFunctionCoefficient u_bdr(pmesh.Dimension(), Initial_Velocity);
    u_0->ProjectCoefficient(u_bdr);

    //Set Pressure Initial Conditions
    ParGridFunction *p_0 = flowsolver.GetCurrentPressure();
    ConstantCoefficient  p_bdr(Parameters.atm_pressure);
    p_0->ProjectCoefficient(p_bdr);

    //Set Velocity Boundary Conditions
    Array<int> bdr_attributes_vel(pmesh.bdr_attributes.Max());
    //code ...
    flowsolver.AddVelDirichletBC(Vel_Boundary_Condition, bdr_attributes_vel);

    //Set  Pressure Boundary Conditions
    Array<int> bdr_attributes_p(pmesh.bdr_attributes.Max());
    //code
    flowsolver.AddPresDirichletBC(Press_Boundary_Condition, bdr_attributes_p);

    //Define Solution Pointers 
    ParGridFunction *u = flowsolver.GetCurrentVelocity(); //Velocity Solution
    ParGridFunction *p = flowsolver.GetCurrentPressure(); //Pressure Solution

    //get Vorticity
    ParGridFunction w = ParGridFunction(&vfes);
    CurlGridFunctionCoefficient curl_u(u);
    w.ProjectCoefficient(curl_u);

    //Create Domain Attr Array
    Array<int> domain_attr(pmesh.bdr_attributes.Max()); 
    domain_attr=1;

    Add Gravity Acceleration Term
    flowsolver.AddAccelTerm(Gravity,domain_attr);

    //Set the RHS due to fictitious forces and set the frisbee conditions
    FictitiousForce ForceRHS(pmesh.Dimension());
    ForceRHS.InitFrisbee(Parameters);
    ForceRHS.Init_Integrals(Parameters,&pmesh,&pfec);
    ForceRHS.SetVelocity(flowsolver.GetCurrentVelocity());
    ForceRHS.SetTime(t);
    flowsolver.AddAccelTerm(&ForceRHS,domain_attr);

    //Set Up Solver
    flowsolver.Setup(Parameters.dt);

    //Paraview Visualization
    ParaViewDataCollection paraview_out = ParaViewDataCollection("results/graph", &pmesh);
    paraview_out.SetLevelsOfDetail(Parameters.order);
    paraview_out.SetDataFormat(VTKFormat::BINARY);
    paraview_out.RegisterField("Pressure", p);
    paraview_out.RegisterField("Velocity", u);
    paraview_out.RegisterField("Vorticity", &w);
    paraview_out.SetCycle(vis_print);
    paraview_out.SetTime(t);
    paraview_out.Save();

    //Perform Time Integration
    for (int step = 0; !last_step; ++step)
    {
        //Check if Last Step
        if (t + Parameters.dt >= Parameters.t_final - Parameters.dt / 2)
            last_step = true;

        //Time Step
        flowsolver.Step(t, Parameters.dt, step);

        //Update Time of The Velocity Boundary Condition
        u_bdr.SetTime(t);

        //Update RHS, i.e, update the current state of frisbee
        //and therfore the Fictituos force
        ForceRHS.SetVelolicty(flowsolver.GetCurrentVelocity());
        ForceRHS.SetTime(t);
        ForceRHS.GetForce_and_Torque(...);
        ForceRHS.UpdateFrisbee(...);

        Compute CFL Condition, Must be <= 1
        double cfl = flowsolver.ComputeCFL(*u, Parameters.dt);

        if(mpi.Root()){
            std::cout << "step" << "\t" << "t" << "\t" << "dt" << "\t" << "cfl" << "\n";
            std::cout << step << "\t" << t << "\t" << Parameters.dt << "\t" << cfl << "\n";
        }

        //Print Data for Visualization
        if (step%Parameters.vis_freq==0)
        {   
            vis_print++;
            CurlGridFunctionCoefficient u_curl(u);
            w.ProjectCoefficient(u_curl);
            paraview_out.SetCycle(vis_print);
            paraview_out.SetTime(t);
            paraview_out.Save();
        }
    }

    flowsolver.PrintTimingData();

    return 0;
}
