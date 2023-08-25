#include <armadillo>
#include <fstream>
#include <cmath>
#include "random.h"
#include "rain.h"


using namespace std;
using namespace arma;


// Default constuctor
Raindrop::Raindrop() {}

// Complete constructor ( [wind_vel] = [body_vel] = [mm/s] )
Raindrop::Raindrop(double radius, vec3 position, vec3 wind_vel, double body_vel) {
    rad = radius;
    vol = M_PI*pow(radius, 3)*4./3.;
    hit = true;
    pos = position;
    vec3 drop_vel;
    //if( radius < )
    vel = wind_vel + drop_vel;
    vel(0) -= body_vel;
}

// Puts the drop in a state of having been hit
void Raindrop::SetHit(bool status) {
    hit = status; 
    // cout << "Hit status: " << hit << endl;
    }

// Returns the radius
double Raindrop::Rad() const { return rad; }

// Returns the volume
double Raindrop::Vol() const { return vol; }

// Returns the hit status
bool Raindrop::Hit() const { return hit; }

// Returns the position 
vec3 Raindrop::Pos() const { return pos; }

// Returns the velocity
vec3 Raindrop::Vel() const { return vel; }

// Time evolution without boundaries ([dt] = [s])
void Raindrop::Move( double dt ) {
    pos += vel*dt;
}

// Time evolution with periodic boundary conditions. pos(i) are bound to [0, max(i)] ([dt] = [s], [max] = [mm])
void Raindrop::Move( double dt, vec3 max ) {
    pos += vel*dt;

    // This enacts periodic boundary conditions
    for ( int i = 0; i < 3; i++ ) {
        while( pos(i) < 0 ){
            pos(i) += max(i);
            hit = false;
            // cout << "unhit" << endl;
        }
        while( pos(i) > max(i) ) {
            pos(i) -= max(i);
            hit = false;
            // cout << "unhit" << endl;
        }
    }
}
