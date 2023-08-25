#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <cmath>
#include "random.h"
#include "DataBlocking.h"
#include "rain.h"
#include "body.h"
#include "rainfun.h"



using namespace std;
using namespace arma;


int main (int argc, char *argv[]){

    Random rand("Primes", "seed.in");

    // Declare stuff
    vec3 wind, box;
    double body_vel, tot_dist;
    int nrain, nstep;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> box(0) >> box(1) >> box(2);
    ReadInput >> nrain;
    ReadInput >> tot_dist;
    ReadInput >> body_vel;
    double tot_time = (body_vel > 0) ? tot_dist/body_vel : 600 ;
    ReadInput >> nstep;
    double dt = tot_time/nstep;
    ReadInput >> wind(0) >> wind(1) >> wind(2);

    ReadInput.close();


    // Outputs settings and stuff
    cout << "Wind = [ " << wind(0) << ", " << wind(1) << ", " << wind(2) << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", dt = " << dt << ", nstep = " << nstep << ", nrain = " << nrain << endl;
    cout << "Time = " << tot_time << " s" << endl;

    // Builds Rain vector
    vector<Raindrop> rain = Rain( nrain, box, wind, body_vel, rand );


    // Builds objects
    Sphere sfera( box/2, box(0)/4 );
    Pippo cube( box/2, sqrt(M_PI)*box/4 );


    // Simulating Sphere
    for( int i = 0; i < nstep; i++){
        sfera.Move(dt);
        for( int j = 0; j < (int) rain.size(); j++ ) {
            rain[j].Move(dt, box);
            sfera.Check(rain[j]);
        }
    }

    // Simulating Cube
    for( int i = 0; i < nstep; i++){
        cube.Move(dt);
        for( int j = 0; j < (int) rain.size(); j++ ) {
            if(i == 0) rain[j].SetHit(true);
            rain[j].Move(dt, box);
            cube.Check(rain[j]);
        }
    }


    // Output
    cout << "Sphere wetness: " << sfera.Wetness()/1000000 << " l" << endl;
    cout << "Cube wetness: " << cube.Wetness()/1000000 << " l" << endl;

    rand.SaveSeed();

    return 0;
}

