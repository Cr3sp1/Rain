#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "body.h"


using namespace std;


int main (int argc, char *argv[]){

    // Declare stuff
    vector<double> wind(3);
    vector<double> box(3); 
    double body_vel, tot_dist;
    int nray, nstep;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> box[0] >> box[1] >> box[2];
    ReadInput >> tot_dist;
    ReadInput >> body_vel;
    double tot_time = (body_vel > 0) ? tot_dist/body_vel : 60 ;
    ReadInput >> nstep;
    double dt = tot_time/nstep;
    ReadInput >> wind[0] >> wind[1] >> wind[2];
    ReadInput.close();


    // Outputs settings and stuff
    cout << "Wind = [ " << wind[0] << ", " << wind[1] << ", " << wind[2] << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", dt = " << dt << ", nstep = " << nstep << ", nray = " << nray << endl;
    cout << "Time = " << tot_time << " s" << endl;

    // Builds Rain 



    // Builds objects
    Sphere sfera();
    Pippo cube();


    // Simulating Sphere


    // Simulating Cube


    // Output
    // cout << "Sphere wetness: " << sfera.Wetness()/1000000 << " l" << endl;
    // cout << "Cube wetness: " << cube.Wetness()/1000000 << " l" << endl;

    return 0;
}

