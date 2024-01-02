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


    ManyBody Temp("../Bodies/TempMan.in");
    Temp.Move(0.25);
    Temp.PrintState("../data/Temp/TempRun.dat");
    // vector<vector<double>> resultsWalk = Simulate( box, Walk, rain_vel, 2, 7, nstep_v, dx, 0, 1, nstep_t );
    // Print("../data/Walk/WalkWet.dat", resultsWalk);
    

    // Error analysis
    ManyBody TrialS("../Bodies/Sphere.in");
    vector<vector<double>> resultsS = SimErr( box, TrialS, rel_vel, body_vel, 100, 0.01, 1 );
    Print( "../data/Sphere/ErrorS.dat", resultsS);
        
    // Error analysis
    ManyBody TrialP("../Bodies/Pippo.in");
    vector<vector<double>> resultsP = SimErr( box, TrialP, rel_vel, body_vel, 100, 0.01, 1 );
    Print( "../data/Pippo/ErrorP.dat", resultsP);

    // Error analysis
    ManyBody TrialC("../Bodies/Capsule.in");
    vector<vector<double>> resultsC = SimErr( box, TrialC, rel_vel, body_vel, 100, 0.01, 1 );
    Print( "../data/Capsule/ErrorC.dat", resultsC);
    
    // // Draw shadow
    // ProjSurface Plz( box, rel_vel, dx );
    // Plz.BodyProj(trialM1);
    // Plz.PrintR("../data/RayOriginsCap.dat");
    // Plz.PrintRaysFlat("../data/RayOriginsCapF.dat");

    

    // Simulate different body velocities
    resultsS = CompareAN( box, *TrialS.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    resultsP = CompareAN( box, *TrialP.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    resultsC = CompareAN( box, *TrialC.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    // vector<vector<double>> resultsM = CompareBB( box, trialM1, trialM2, rain_vel, 2, 7, nstep_v, dx );

    Print("../data/Sphere/CompareS.dat", resultsS);
    Print("../data/Pippo/CompareP.dat", resultsP);
    Print("../data/Capsule/CompareC.dat", resultsC);
    // Print("../data/CompareManyBody.dat", resultsM);

    return 0;
}

