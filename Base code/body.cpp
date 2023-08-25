#include <armadillo>
#include <fstream>
#include <cmath>
#include "body.h"
#include "random.h"
#include "rain.h"



using namespace std;
using namespace arma;


// Default constuctor
Body::Body() {wetness = 0;}

// Returns the wetness
double Body::Wetness() const { return wetness; }


// Complete Sphere constructor ( [center] = [mm], [radius] = [mm] )
Sphere::Sphere( vec3 center, double radius ): Body() {
    cent = center;
    rad = radius;
}

// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
void Sphere::Check( Raindrop& drop ) {
    if( drop.Hit() == false && arma::norm( drop.Pos() - cent ) <= rad + drop.Rad() ) {
        wetness += drop.Vol();
        drop.SetHit(true);
        // cout << "hit" << endl;
    }
}



// Complete Parallelepiped constructor ( [center] = [mm], [radius] = [mm] )
Pippo::Pippo( vec3 center, vec3 dimensions ): Body() {
    cent = center;
    dim = dimensions;
}

// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
void Pippo::Check( Raindrop& drop ) {
    double radius = drop.Rad();
    vec3 position = drop.Pos();
    if( drop.Hit() == false && all(cent - dim/2 - radius < position  &&  position < cent + dim/2 + radius ) ) {
        wetness += drop.Vol();
        drop.SetHit(true);
        // cout << "hit" << endl;
    }
}


