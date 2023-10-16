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
    hex.PrintR("RayOrigins.dat");
    hex.PrintH("H.dat");


    // Checking rays generation
    // ofstream Pout("RaysGen.dat");
    // for( int i = 0; i < 10000; i+=50 ){
    //     ProjSurface temp( box, rel_vel, i+1 );
    //     Pout << (double)temp.GetNRays()/(i+1) << endl;
    // }
    // Pout.close();



    // Builds objects


    // Simulating Sphere


    // Output


    return 0;
}

