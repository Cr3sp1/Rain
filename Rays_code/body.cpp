#include "body.h"
#include "ray.h"

using namespace std;


// Default constuctor
Body::Body() {}

// Virtual functions
bool Body::Check( Ray& ray ) {return 0;}
void Body::Move( long double dt ) {}


// Complete Sphere constructor ( [center] = [mm], [radius] = [mm] )
Sphere::Sphere( vector<long double> center, long double radius ): Body() {
    cent = center;
    rad = radius;
}

// Checks if the body is making contact with a ray
bool Sphere::Check( Ray& ray ) {
    if( ray.IsOn() != true ) return false;
    vector<long double> v = ray.GetV();
    vector<long double> R0 = ray.GetR0();
    long double a = v*v;
    long double b = 2*((R0-cent)*v);
    long double c = (v-R0)*(v-R0) - rad*rad;
    if( b*b >= 4*a*c ) return true;
    
    return false;
}




// Checks if the body is making contact with a ray
// bool Pippo::Check( Ray ray ) {

// }

