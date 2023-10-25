#include "body.h"
#include "ray.h"

using namespace std;


// Default constuctor
Body::Body() {}

// Virtual functions
bool Body::Check( Ray& ray ) {return 0;}
void Body::Move( long double dt ) {}
void Body::Prime( vector<long double> p, vector<long double> v  ) {}
long double Body::Anal( vector<long double> v, long double bodyvel  ) {return -1;}


// Complete Sphere constructor ( [center] = [mm], [radius] = [mm] )
Sphere::Sphere( vector<long double> center, long double radius ): Body() {
    cent = center;
    rad = radius;
}

// Primes the body to be checked (projects the center of the sphere onto the surface)
void Sphere::Prime( vector<long double> p, vector<long double> v  ) {
    Hcent = Project( cent, p, v);
}

// Checks if the Sphere is making contact with a ray
bool Sphere::Check( Ray& ray ) {
    if( ray.IsOn() != true ) return false;
    vector<long double> Xrel = ray.GetR0() - Hcent;
    if( Xrel*Xrel <= rad*rad ) {
        ray.Off();
        return true;
    }
    return false;
}

// Analytical solution of rain intercepted. v is relative velocity
long double Sphere::Anal( vector<long double> v, long double bodyvel ) {
    long double surface = M_PI*rad*rad;
    return Norm(v)*surface/bodyvel;
}
	



// Complete Parallelepiped constructor 
Pippo::Pippo( vector<long double> P, vector<vector<long double>> Side ): Body() {
    p = P;
    side = Side;
}

// Primes the body to be checked (finds hexagonal projection on the same plane as the origins of the rays)
void Pippo::Prime( vector<long double> P, vector<long double> V ) {
    H = FindHexProj( p, side, V, P );
}

// Checks if the body is making contact with a ray
bool Pippo::Check( Ray& ray ) {
    if( ray.IsOn() != true ) return false;
    if( PointIsInsideT( ray.GetR0(), H )) {
        ray.Off();
        return true;
    }
    return false;
}

// Analytical solution of rain intercepted
long double Pippo::Anal( vector<long double> v, long double bodyvel  ) {
    Prime({0,0,0}, v );
    long double surf = 0;
    for( int i = 1; i < 7; i+=2 ) {
        surf += Norm( CrossProduct( H[0]-H[i], H[i]-H[i+1] ) );
    }
    return Norm(v)*surf/bodyvel;
}




