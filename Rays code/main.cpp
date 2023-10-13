#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "body.h"
#include "ray.h"


using namespace std;


int main (int argc, char *argv[]){

    // Declare stuff
    vector<double> rain_vel(3);
    vector<double> rel_vel(3);
    vector<double> box(3); 
    double body_vel, tot_dist;
    int nray, nstep;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> box[0] >> box[1] >> box[2];
    ReadInput >> nray;
    ReadInput >> tot_dist;
    ReadInput >> body_vel;
    double tot_time = (body_vel > 0) ? tot_dist/body_vel : 60 ;
    ReadInput >> nstep;
    double dt = tot_time/nstep;
    ReadInput >> rain_vel[0] >> rain_vel[1] >> rain_vel[2];
    rel_vel = rain_vel;
    rel_vel[0] -= body_vel;
    ReadInput.close();


    // Outputs settings and stuff
    cout << "relative velocity = [ " << rel_vel[0] << ", " << rel_vel[1] << ", " << rel_vel[2] << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", dt = " << dt << ", nstep = " << nstep << ", nray = " << nray << endl;
    cout << "Time = " << tot_time << " s" << endl;

    // Builds Rays
    ProjSurface hex( box, rel_vel, 100 );


    // Builds objects


    // Simulating Sphere


    // Output


    return 0;
}

