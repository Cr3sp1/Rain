#include "body.h"

using namespace std;


// Default constuctor
Body::Body() {wetness = 0;}

// Returns the wetness
double Body::Wetness() const { return wetness; }


// Complete Sphere constructor ( [center] = [mm], [radius] = [mm] )
Sphere::Sphere( vector<double> center, double radius ): Body() {
    cent = center;
    rad = radius;
}

// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
void Sphere::Check( ) {
    
}



// Complete Parallelepiped constructor ( [center] = [mm], [radius] = [mm] )
Pippo::Pippo( vector<double> center, vector<double> dimensions ): Body() {
    cent = center;
    dim = dimensions;
}

// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
void Pippo::Check( ) {
    
}


