#include "body.h"

using namespace std;


// Default constuctor
Body::Body() {}

// Complete Sphere constructor ( [center] = [mm], [radius] = [mm] )
Sphere::Sphere( vector<long double> center, long double radius ): Body() {
    cent = center;
    rad = radius;
}

// Checks if the body is making contact with a ray
// bool Sphere::Check( Ray ray ) {
    
// }



// Complete Parallelepiped constructor ( [center] = [mm], [radius] = [mm] )
Pippo::Pippo( vector<long double> center, vector<long double> dimensions ): Body() {
    cent = center;
    dim = dimensions;
}

// Checks if the body is making contact with a ray
// bool Pippo::Check( Ray ray ) {

// }


