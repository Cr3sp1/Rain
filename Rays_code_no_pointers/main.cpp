#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "body.h"
#include "ray.h"
# include "VectorStat.h"


using namespace std;


int main (int argc, char *argv[]){

    // Declare stuff
    vector<double> rain_vel(3);
    vector<double> rel_vel(3);
    vector<double> box(3); 
    double body_vel, dx;
    int nstep;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> box[0] >> box[1] >> box[2];
    ReadInput >> dx;
    ReadInput >> body_vel;
    ReadInput >> nstep;
    ReadInput >> rain_vel[0] >> rain_vel[1] >> rain_vel[2];
    rel_vel = rain_vel;
    rel_vel[0] -= body_vel;
    ReadInput.close();


    // Outputs settings and stuff
    cout << "relative velocity = [ " << rel_vel[0] << ", " << rel_vel[1] << ", " << rel_vel[2] << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", nstep = " << nstep << ", dx = " << dx << endl;
    // cout << "Time = " << tot_time << " s" << endl;

    // Check many-body and motion
    vector<vector<double>> Rot({{1,1,1}, {1,1,2}});
    vector<vector<double>> Tran{{0,0,0.5}};
    double W = 0.5*M_PI;

    ManyBody Mb;
    Mb.AddBody( Sphere((double)0.2*box, (double)0.2*box[0], "S1", {}, 0, {} ) );
    Mb.AddBody(Pippo( box*(double)0.5, { {box[0]/4, 0, 0}, {0, box[1]/5, 0}, {0, 0, box[2]/6} }, "P1", {box*(double)0.5, box*(double)0.5+ vector<double>{0,1,0} }, W, {} ));
    Mb.AddBody( Capsule( (double)0.25*box, (double)0.75*box,  (double)0.1*box[0], "C1", {}, 0, {} ));
    Mb.PrintState("../data/Mb1.dat");
    Mb.Move(0.25);
    Mb.PrintState("../data/Mb2.dat");

    Mb.Move(0);
    Mb.Attach( "S1", "P1");
    Mb.Attach( "C1", "S1");
    Mb.Move(0.25);
    Mb.PrintState("../data/Mb3.dat");


     
    // // Builds objects
    // Sphere trialS( (double)0.7*box, (double)0.25*box[0] );
    // Pippo trialP(box*(double)0.5, { {box[0]/4, 0, 0}, {0, box[1]/4, 0}, {0, 0, box[2]/4} });
    // Capsule trialC( (double)0.25*box, (double)0.75*box,  (double)0.1*box[0] );
    // ManyBody trialM1( {trialS}, {trialP}, {trialC} );
    // trialM1.PrintState("../data/BodyState1.dat");
    // Sphere inner( (double)0.7*box, (double)0.15*box[0] );
    // ManyBody trialM2( {trialS, inner}, {trialP}, {trialC} );
    

    // // Draw shadow
    // ProjSurface Plz( box, rel_vel, dx );
    // Plz.PrintH("../data/H.dat");
    // Plz.PrintR("../data/RayOrigins.dat");
    // Plz.BodyProj(trialM1);
    // Plz.PrintR("../data/RayOriginsCap.dat");

    

    // // Simulate different body velocities
    // vector<vector<double>> resultsS = CompareAN( box, trialS, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsP = CompareAN( box, trialP, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsC = CompareAN( box, trialC, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsM = CompareBB( box, trialM1, trialM2, rain_vel, 2, 7, nstep, dx );


    // // output
    // Print("../data/CompareSphere.dat", resultsS);
    // Print("../data/ComparePippo.dat", resultsP);
    // Print("../data/CompareCapsule.dat", resultsC);
    // Print("../data/CompareManyBody.dat", resultsM);

    return 0;
}
