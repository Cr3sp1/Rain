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
    unsigned int nstep_v, nstep_t;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> box[0] >> box[1] >> box[2];
    ReadInput >> dx;
    ReadInput >> body_vel;
    ReadInput >> nstep_v;
    ReadInput >> nstep_t;
    ReadInput >> rain_vel[0] >> rain_vel[1] >> rain_vel[2];
    rel_vel = rain_vel;
    rel_vel[0] -= body_vel;
    ReadInput.close();


    // Outputs settings and stuff
    cout << "relative velocity = [ " << rel_vel[0] << ", " << rel_vel[1] << ", " << rel_vel[2] << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", nstep_v = " << nstep_v << ", dx = " << dx << endl;



    // // Walking man
    // ManyBody Walk("../Bodies/WalkinMan.in");
    // ProjSurface BodSurfW( box, rel_vel, dx );
    // BodSurfW.BodyProj( Walk );
    // BodSurfW.PrintRaysFlat("../data/Walk/WalkProjF.dat");
    // for( size_t i = 0; i < 11; i++ ){
    //     Walk.Move(0.1*i);
    //     Walk.PrintState( ("../data/Walk/Walk" + to_string(i) + ".dat"));
    // }

    // // Running man
    // ManyBody Run("../Bodies/RunningMan.in");
    // for( size_t i = 0; i < 11; i++ ){
    //     Run.Move(0.1*i);
    //     Run.PrintState( ("../data/Run/Run" + to_string(i) + ".dat"));
    // }

    // vector<vector<double>> resultsWalk = Simulate( box, Walk, rain_vel, 2, 7, nstep_v, dx, 0, 1, nstep_t );
    // Print("../data/Walk/WalkWet.dat", resultsWalk);


    // Temp
    ManyBody Temp("../Bodies/TempMan.in");
    Temp.Move(0.25);
    Temp.PrintState("../data/Temp/TempWalk.dat");

    // // Walk
    // ManyBody Walk("../Bodies/WalkingMan.in");
    // Walk.PrintState("../data/Temp/Walk.dat");
    
     
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
    // Plz.BodyProj(trialM1);
    // Plz.PrintR("../data/RayOriginsCap.dat");
    // Plz.PrintRaysFlat("../data/RayOriginsCapF.dat");

    

    // Simulate different body velocities
    // vector<vector<double>> resultsS = CompareAN( box, trialS, rain_vel, 2, 7, nstep_v, dx );
    // vector<vector<double>> resultsP = CompareAN( box, trialP, rain_vel, 2, 7, nstep_v, dx );
    // vector<vector<double>> resultsC = CompareAN( box, trialC, rain_vel, 2, 7, nstep_v, dx );
    // vector<vector<double>> resultsM = CompareBB( box, trialM1, trialM2, rain_vel, 2, 7, nstep_v, dx );


    // output
    // Print("../data/CompareSphere.dat", resultsS);
    // Print("../data/ComparePippo.dat", resultsP);
    // Print("../data/CompareCapsule.dat", resultsC);
    // Print("../data/CompareManyBody.dat", resultsM);

    return 0;
}

