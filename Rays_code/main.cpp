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



    // Walking man
    ManyBody Walk("../Bodies/WalkinMan.in");
    ProjSurface BodSurf( box, rel_vel, dx );
    BodSurf.BodyProj( Walk );
    BodSurf.PrintR("../data/WalkProj.dat");
    BodSurf.PrintRaysFlat("../data/WalkProjF.dat");
    for( size_t i = 0; i < 11; i++ ){
        Walk.Move(0.025*i);
        Walk.PrintState( ("../data/Walk" + to_string(i) + ".dat"));
    }
    
     
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
    // vector<vector<double>> resultsS = CompareAN( box, trialS, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsP = CompareAN( box, trialP, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsC = CompareAN( box, trialC, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsM = CompareBB( box, trialM1, trialM2, rain_vel, 2, 7, nstep, dx );


    // output
    // Print("../data/CompareSphere.dat", resultsS);
    // Print("../data/ComparePippo.dat", resultsP);
    // Print("../data/CompareCapsule.dat", resultsC);
    // Print("../data/CompareManyBody.dat", resultsM);

    return 0;
}

