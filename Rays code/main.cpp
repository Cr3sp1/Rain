#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "body.h"


using namespace std;


int main (int argc, char *argv[]){

    Random rand("Primes", "seed.in");

    // Declare stuff
    vector<double> wind(3);
    double body_vel, tot_dist;
    int nrain, nstep;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> tot_dist;
    ReadInput >> body_vel;
    double tot_time = (body_vel > 0) ? tot_dist/body_vel : 600 ;
    ReadInput >> nstep;
    double dt = tot_time/nstep;
    ReadInput >> wind[0] >> wind[1] >> wind[2];

    ReadInput.close();


    // Outputs settings and stuff
    cout << "Wind = [ " << wind[0] << ", " << wind[1] << ", " << wind[2] << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", dt = " << dt << ", nstep = " << nstep << ", nrain = " << nrain << endl;
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

    rand.SaveSeed();

    return 0;
}

