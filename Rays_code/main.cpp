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
    vector<vector<double>> Rot({{0,0,0}, {0,0,1}});
    double W = 0.5*M_PI;

    Body body1( Rot, W, {} );
    Sphere S1( {1,0,0}, 1, Rot, W, {});
    
    // check pointsegdist
    // vector<double> point{2,1,0};
    // vector<double> l1{1,0,0};
    // vector<double> l2{7,0,0};
    // cout << "Distance = " << PointSegDist(point, l1, l2) << endl;
     
    // Builds objects
    // Sphere trialS( (double)0.5*box, (double)0.3*box[0]);
    // Pippo trialP(box*(double)0.5, { {box[0]/4, 0, 0}, {0, box[1]/4, 0}, {0, 0, box[2]/4} });
    Capsule trialC( (double)0.25*box, (double)0.75*box,  (double)0.2*box[0] );
    // ManyBody trialM1( {trialS}, {trialP}, {trialC} );
    // trialM1.PrintState("../data/BodyState1.dat");
    // Sphere inner((double)0.5*box, (double)0.2*box[0]);
    // ManyBody trialM2( {trialS, inner}, {trialP}, {trialC} );

    // Debugging capsule
    ProjSurface Plz( box, rel_vel, dx );
    Plz.PrintH("../data/H.dat");
    Plz.PrintR("../data/RayOrigins.dat");
    Plz.BodyProj(trialC);
    Plz.PrintR("../data/RayOriginsCap.dat");

    

    // Simulate different body velocities
    // vector<vector<double>> resultsS = CompareAN( box, trialS, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsP = CompareAN( box, trialP, rain_vel, 2, 7, nstep, dx );
    vector<vector<double>> resultsC = CompareAN( box, trialC, rain_vel, 2, 7, nstep, dx );
    // vector<vector<double>> resultsM = CompareBB( box, trialM1, trialM2, rain_vel, 2, 7, nstep, dx );


    // output
    // Print("../data/CompareSphere.dat", resultsS);
    // Print("../data/ComparePippo.dat", resultsP);
    Print("../data/CompareCapsule.dat", resultsC);
    // Print("../data/CompareManyBody.dat", resultsM);

    // ofstream outputFileS("../data/CompareSphere.dat");
    // for (size_t i = 0; i < resultsS[0].size(); ++i) {
    //     outputFileS << resultsS[0][i] << " " << resultsS[1][i] << " " << resultsS[2][i] << endl;
    // }
    // outputFileS.close();
    

    // ofstream outputFileP("../data/ComparePippo.dat");
    // for (size_t i = 0; i < resultsP[0].size(); ++i) {
    //     outputFileP << resultsP[0][i] << " " << resultsP[1][i] << " " << resultsP[2][i] << endl;
    // }
    // outputFileP.close();

    // ofstream outputFileC("../data/CompareCapsule.dat");
    // for (size_t i = 0; i < resultsS[0].size(); ++i) {
    //     outputFileC << resultsC[0][i] << " " << resultsS[1][i] << " " << resultsS[2][i] << endl;
    // }
    // outputFileC.close();

    // ofstream outputFileM("../data/CompareManyBody.dat");
    // for (size_t i = 0; i < resultsM[0].size(); ++i) {
    //     outputFileM << resultsM[0][i] << " " << resultsM[1][i] << " " << resultsM[2][i] << endl;
    // }
    // outputFileM.close();


    return 0;
}

