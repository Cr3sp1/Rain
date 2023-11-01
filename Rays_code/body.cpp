#include "body.h"
#include "ray.h"

using namespace std;


// Default constuctor
Body::Body() { t = 0; }

// Virtual functions
bool Body::Check( Ray& ray ) {return 0;}
void Body::Move( long double T ) {}
void Body::Prime( vector<long double> p, vector<long double> v  ) {}
long double Body::Anal( vector<long double> v, long double bodyvel  ) {return -1;}




// Complete Sphere constructor 
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
    if( ray.IsOn() == false ) return false;
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
Pippo::Pippo( vector<long double> Cent, vector<vector<long double>> Side ): Body() {
    cent = Cent;
    side = Side;
}

// Primes the body to be checked (finds hexagonal projection on the same plane as the origins of the rays)
void Pippo::Prime( vector<long double> P, vector<long double> V ) {
    vector<long double> p = cent + (long double)0.5*(side[0] + side[1] + side[2]);
    H = FindHexProj( p, side, V, P );
}

// Checks if the body is making contact with a ray
bool Pippo::Check( Ray& ray ) {
    if( ray.IsOn() == false ) return false;
    if( PointIsInsideT( ray.GetR0(), H )) {
        ray.Off();
        return true;
    }
    return false;
}

// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
long double Pippo::Anal( vector<long double> v, long double bodyvel  ) {
    long double flux = 0;
    flux += abs(CrossProduct(side[0], side[1]) * v);
    flux += abs(CrossProduct(side[1], side[2]) * v);
    flux += abs(CrossProduct(side[2], side[0]) * v);
    return flux/bodyvel;
}




// Complete Capsule constructor 
Capsule::Capsule( vector<long double> L1, vector<long double> L2,  long double radius ): Body() {
    l1 = L1;
    l2 = L2;
    rad = radius;
}

// Primes the body to be checked (projects the center of the sphere onto the surface)
void Capsule::Prime( vector<long double> p, vector<long double> v  ) {
    H1 = Project( l1, p, v);
    H2 = Project( l2, p, v);
}

// Checks if the Capsule is making contact with a ray
bool Capsule::Check( Ray& ray ) {
    if( ray.IsOn() == false ) return false;

    if( PointSegDist( ray.GetR0(), H1, H2 ) <= rad ) {
        ray.Off();
        return true;
    }
    return false;
}

// Analytical solution of rain intercepted. v is relative velocity
long double Capsule::Anal( vector<long double> v, long double bodyvel ) {
    long double L = abs((l1 - l2)*v)/Norm(v);
    long double surface = M_PI*rad*rad + L*rad;
    return Norm(v)*surface/bodyvel;
}




// Complete ManyBody constructor 
ManyBody::ManyBody( vector<Sphere> Spheres, vector<Pippo> Pippos, vector<Capsule> Capsules ): Body() {
    spheres = Spheres;
    pippos = Pippos;
    capsules = Capsules;
}

// Primes the body to be checked. Primes each body
void ManyBody::Prime( vector<long double> p, vector<long double> v  ) {
    for( long unsigned int i = 0; i < spheres.size(); i++ ) spheres[i].Prime( p, v );
    for( long unsigned int i = 0; i < pippos.size(); i++ ) pippos[i].Prime( p, v );
    for( long unsigned int i = 0; i < capsules.size(); i++ ) capsules[i].Prime( p, v );
}

// Checks if the ManyBody is making contact with a ray
bool ManyBody::Check( Ray& ray ) {
    for( long unsigned int i = 0; i < spheres.size(); i++ ) if(spheres[i].Check( ray )) return true;
    for( long unsigned int i = 0; i < pippos.size(); i++ ) if(pippos[i].Check( ray )) return true;
    for( long unsigned int i = 0; i < capsules.size(); i++ ) if(capsules[i].Check( ray )) return true;
    return false;
}

// Gets stuff
vector<long double> ManyBody::GetSphCent( unsigned int index ){
    if( index >= spheres.size() ){
        cout << "No sphere with index " << index << endl;
        return {}; 
    }
    return spheres[index].GetCent();
}

long double ManyBody::GetSphRad( unsigned int index ){
    if( index >= spheres.size() ){
        cout << "No sphere with index " << index << endl;
        return 0; 
    }
    return spheres[index].GetRad();
}

vector<long double> ManyBody::GetPipCent( unsigned int index ){
    if( index >= pippos.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return pippos[index].GetCent();
}

vector<vector<long double>>  ManyBody::GetPipSide( unsigned int index ) {
    if( index >= pippos.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return pippos[index].GetSide();
}

vector<long double> ManyBody::GetCapL1( unsigned int index ){
    if( index >= capsules.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return capsules[index].GetL1();
}
	
vector<long double> ManyBody::GetCapL2( unsigned int index ){
    if( index >= capsules.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return capsules[index].GetL2();
}
	
long double ManyBody::GetCapRad( unsigned int index ){
    if( index >= capsules.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return 0; 
    }
    return capsules[index].GetRad();
}
