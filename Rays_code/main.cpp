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
    vector<long double> rain_vel(3);
    vector<long double> rel_vel(3);
    vector<long double> box(3); 
    long double body_vel, dx;
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

    // Builds Rays
    ProjSurface hex( box, rel_vel, dx );

    // Prints Rays
    hex.PrintR("../data/RayOrigins.dat");
    hex.PrintH("../data/H.dat");

     
    // Builds objects
    Sphere trialS( (long double)0.5*box, (long double)0.3*box[0]);
    Pippo trialP(box*(long double)0.5, { {box[0]/3, 0, 0}, {0, box[1]/3, 0}, {0, 0, box[2]/3} });
    Capsule trialC( (long double)0.3*box, (long double)0.7*box,  (long double)0.2*box[0] );
    ManyBody trialM1( {trialS}, {trialP}, {trialC} );
    Sphere inner((long double)0.5*box, (long double)0.2*box[0]);
    ManyBody trialM2( {trialS, inner}, {trialP}, {trialC} );

    // Simulate different body velocities
    vector<vector<long double>> resultsS = CompareAN( box, trialS, rain_vel, 2, 7, nstep, dx );
    vector<vector<long double>> resultsP = CompareAN( box, trialP, rain_vel, 2, 7, nstep, dx );
    vector<vector<long double>> resultsC = CompareAN( box, trialC, rain_vel, 2, 7, nstep, dx );
    vector<vector<long double>> resultsM = CompareBB( box, trialM1, trialM2, rain_vel, 2, 7, nstep, dx );



    // output
    ofstream outputFileS("../data/CompareSphere.dat");
    for (size_t i = 0; i < resultsS[0].size(); ++i) {
        outputFileS << resultsS[0][i] << " " << resultsS[1][i] << " " << resultsS[2][i] << endl;
    }
    outputFileS.close();

    ofstream outputFileP("../data/ComparePippo.dat");
    for (size_t i = 0; i < resultsP[0].size(); ++i) {
        outputFileP << resultsP[0][i] << " " << resultsP[1][i] << " " << resultsP[2][i] << endl;
    }
    outputFileP.close();

    ofstream outputFileC("../data/CompareCapsule.dat");
    for (size_t i = 0; i < resultsS[0].size(); ++i) {
        outputFileC << resultsS[0][i] << " " << resultsS[1][i] << " " << resultsS[2][i] << endl;
    }
    outputFileC.close();

    ofstream outputFileM("../data/CompareManyBody.dat");
    for (size_t i = 0; i < resultsM[0].size(); ++i) {
        outputFileM << resultsM[0][i] << " " << resultsM[1][i] << " " << resultsM[2][i] << endl;
    }
    outputFileM.close();

    return 0;
}

