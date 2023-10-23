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
    long double body_vel, tot_dist;
    int nray, nstep;
    vector<long double> vel_ratio;
    vector<long double> NumResults;
    vector<long double> AnalResults;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> box[0] >> box[1] >> box[2];
    ReadInput >> nray;
    ReadInput >> tot_dist;
    ReadInput >> body_vel;
    long double tot_time = (body_vel > 0) ? tot_dist/body_vel : 60 ;
    ReadInput >> nstep;
    long double dt = tot_time/nstep;
    ReadInput >> rain_vel[0] >> rain_vel[1] >> rain_vel[2];
    rel_vel = rain_vel;
    rel_vel[0] -= body_vel;
    ReadInput.close();


    // Outputs settings and stuff
    cout << "relative velocity = [ " << rel_vel[0] << ", " << rel_vel[1] << ", " << rel_vel[2] << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", dt = " << dt << ", nstep = " << nstep << ", nray = " << nray << endl;
    cout << "Time = " << tot_time << " s" << endl;

    // Builds Rays
    ProjSurface hex( box, rel_vel, nray );

    // Prints Rays
    hex.PrintR("../data/RayOrigins.dat");
    hex.PrintH("../data/H.dat");

    
    // Builds objects
    Sphere trial( box*(long double)0.5, box[0]/2 );


    // Simulating parallelepiped
    long double AreaS = hex.BodyProj(trial)/(M_PI*trial.GetRad()*trial.GetRad());

    // Output
    cout << "Area of the sphere is " << AreaS << "*PI*r^2" << endl;

    // Simulate different body velocities
    for( int i = 0; i < 20; i++ ){
        vel_ratio.push_back( (4*i+20)*0.03);
        body_vel= vel_ratio[i]*Norm(rain_vel);
        rel_vel = rain_vel;
        rel_vel[0] -= body_vel;

        AnalResults.push_back(trial.Anal(rel_vel, tot_dist, body_vel));
        ProjSurface temp( box, rel_vel, nray);
        NumResults.push_back(temp.BodyProj( trial )*Norm(rel_vel)*tot_dist/body_vel);
    }

    vector<vector<long double>> results = { vel_ratio, AnalResults, NumResults };
    ofstream outputFile("../data/CompareSphere.txt");

    for (size_t i = 0; i < results[0].size(); ++i) {
        outputFile << results[0][i] << " " << results[1][i] << " " << results[2][i] << endl;
    }

    outputFile.close();

    return 0;
}

