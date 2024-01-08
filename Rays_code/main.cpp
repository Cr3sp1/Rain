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
    // ManyBody Walk("../Bodies/WalkingMan.in");
    // // ProjSurface BodSurfW( box, rel_vel, dx );
    // // BodSurfW.BodyProj( Walk );
    // //BodSurfW.PrintRaysFlat("../data/Walk/WalkProjF.dat");
    // for( size_t i = 0; i < 9; i++ ){
    //     double t = 0.125*i;
    //     Walk.Move(t);
    //     Walk.PrintState("../data/Walk/Walk"+to_string(t).substr(0,5)+".dat");
    // }

    // // Running man
    // ManyBody Run("../Bodies/RunningMan.in");
    // for( size_t i = 0; i < 9; i++ ){
    //     double t = 0.125*i;
    //     Run.Move(t);
    //     Run.PrintState("../data/Run/Run"+to_string(t).substr(0,5)+".dat");
    // }


    // ManyBody Temp("../Bodies/TempMan.in");
    // Temp.Move(0.25);
    // Temp.PrintState("../data/Temp/TempRun.dat");
    // vector<vector<double>> resultsWalk = Simulate( box, Walk, rain_vel, 2, 7, nstep_v, dx, 0, 1, nstep_t );
    // Print("../data/Walk/WalkWet.dat", resultsWalk);
    

    // Error analysis
    ManyBody TrialS("../Bodies/Sphere.in");
    vector<vector<double>> resultsS;
    // resultsS = SimErr( box, TrialS, rel_vel, body_vel, 200, 0.0001, 1 );
    // Print( "../data/Sphere/ErrorS.dat", resultsS, 15);
        
    // Error analysis
    ManyBody TrialP("../Bodies/Pippo.in");
    vector<vector<double>> resultsP;
    // resultsP = SimErr( box, TrialP, rel_vel, body_vel, 200, 0.0001, 1 );
    // Print( "../data/Pippo/ErrorP.dat", resultsP, 15);

    // Error analysis
    ManyBody TrialC("../Bodies/Capsule.in");
    vector<vector<double>> resultsC;
    // resultsC = SimErr( box, TrialC, rel_vel, body_vel, 200, 0.0001, 1 );
    // Print( "../data/Capsule/ErrorC.dat", resultsC, 15);
    
    // // Draw shadow
    // ProjSurface Plz( box, rel_vel, dx );
    // Plz.BodyProj(trialM1);
    // Plz.PrintR("../data/RayOriginsCap.dat");
    // Plz.PrintRaysFlat("../data/RayOriginsCapF.dat");


    // Simulate different body velocities
    // resultsS = CompareAN( box, *TrialS.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Sphere/CompareS.dat", resultsS, 15);
    // resultsP = CompareAN( box, *TrialP.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Pippo/CompareP.dat", resultsP, 15);
    // resultsC = CompareAN( box, *TrialC.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Capsule/CompareC.dat", resultsC, 15);


    // Simulation of two pippos compenetrating
    ManyBody Trial2P("../Bodies/DoublePippo.in");
    vector<double> dist;
    vector<double> wet2P;
    ProjSurface Surf2P( box, rel_vel, dx );
    for( size_t i = 0; i < nstep_t; i++ ) {
        Trial2P.Move( asin( (double)i/(nstep_t-1))/(2*M_PI));   // it just works ;)
        vector<double> cent1 = dynamic_cast<Pippo*>(Trial2P.Find("Still"))->GetCent();
        vector<double> cent2 = dynamic_cast<Pippo*>(Trial2P.Find("Moving"))->GetCent();
        dist.push_back( Norm( cent1 - cent2 ));
        Surf2P.reset();
        wet2P.push_back(Surf2P.BodyProj(Trial2P)*Norm(rel_vel)/body_vel);
        // Surf2P.PrintRaysFlat("../data/Pippo/Proj2P/dist" + to_string(dist[i]) + ".dat");
    }

    vector<vector<double>> results2P = { dist, wet2P};
    results2P = Transpose(results2P);
    Print("../data/Pippo/DoubleP.dat", results2P, 12 );
    

    return 0;
}

