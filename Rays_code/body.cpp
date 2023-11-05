#include "body.h"
#include "ray.h"

using namespace std;



// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Body::Move( double T ) {
    if( T == t ) return;

    SuperBody->Move(T);

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( long unsigned int i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double> point : rot ){
        point += delta;
    }

    // Generate rotation matrix
    double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
    vector<vector<double>> rotmat = RotMat( rot[1]-rot[0], theta );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Body::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double> point : rot ) point += Delta;

    // Rotate
    for( vector<double> point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double> vec : trans ) Rotate( vec, Rot0, Rotmat );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
}

// Copy constructor
Body::Body(const Body& other) : t(other.t), rot(other.rot), w(other.w), trans(other.trans) {
    // Deep copy sub-bodies
    for (const Body* body : other.SubBodies) {
        SubBodies.push_back(new Body(*body));
    }

    // Deep copy super-body
    if (other.SuperBody) {
        SuperBody= new Body(*other.SuperBody);
    } else {
        SuperBody = nullptr;
    }
}


// Copy assignment operator
Body& Body::operator=(const Body& other) {
    if (this != &other) {
        // Clean up existing sub-bodies
        for (Body* body : SubBodies) {
            delete body;
        }
        SubBodies.clear();

        // Deep copy connected bodies
        for (const Body* body : other.SubBodies) {
            SubBodies.push_back(new Body(*body));
        }

        // Clean up existing linked body
        delete SuperBody;

        // Deep copy linked body
        if (other.SuperBody) {
            SuperBody = new Body(*other.SuperBody);
        } else {
            SuperBody = nullptr;
        }

        // Copy other stuff
        t = other.t;
        rot = other.rot;
        w = other.w;
        trans = other.trans;
    }
    return *this;
}


// Destructor
Body::~Body() {
    for (Body* body : SubBodies) {
        delete body;
    }
    delete SuperBody;
}
    



// Primes the body to be checked (projects the center of the sphere onto the surface)
void Sphere::Prime( vector<double> p, vector<double> v  ) {
    Hcent = Project( cent, p, v);
}

// Checks if the Sphere is making contact with a ray
bool Sphere::Check( Ray& ray ) {
    if( ray.IsOn() == false ) return false;
    vector<double> Xrel = ray.GetR0() - Hcent;
    if( Xrel*Xrel <= rad*rad ) {
        ray.Off();
        return true;
    }
    return false;
}

// Analytical solution of rain intercepted. v is relative velocity
double Sphere::Anal( vector<double> v, double bodyvel ) {
    double surface = M_PI*rad*rad;
    return Norm(v)*surface/bodyvel;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Sphere::Move( double T ) {
    if( T == t ) return;

    SuperBody->Move(T);

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( long unsigned int i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double> point : rot ) point += delta;
    cent += delta;

    // Generate rotation matrix
    double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
    vector<vector<double>> rotmat = RotMat( rot[1]-rot[0], theta );

    // Rotate
    Rotate( cent, rot[0], rotmat );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Sphere::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double> point : rot ) point += Delta;
    cent += Delta;

    // Rotate
    for( vector<double> point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double> vec : trans ) Rotate( vec, Rot0, Rotmat );
    Rotate( cent, Rot0, Rotmat );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
}
	



// Primes the body to be checked (finds hexagonal projection on the same plane as the origins of the rays)
void Pippo::Prime( vector<double> P, vector<double> V ) {
    vector<double> p = cent + (double)0.5*(side[0] + side[1] + side[2]);
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
double Pippo::Anal( vector<double> v, double bodyvel  ) {
    double flux = 0;
    flux += abs(CrossProduct(side[0], side[1]) * v);
    flux += abs(CrossProduct(side[1], side[2]) * v);
    flux += abs(CrossProduct(side[2], side[0]) * v);
    return flux/bodyvel;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Pippo::Move( double T ) {
    if( T == t ) return;

    SuperBody->Move(T);

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( long unsigned int i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double> point : rot ) point += delta;
    cent += delta;
    for( vector<double> point : side ) point += delta;

    // Generate rotation matrix
    double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
    vector<vector<double>> rotmat = RotMat( rot[1]-rot[0], theta );

    // Rotate
    Rotate( cent, rot[0], rotmat );
    for( vector<double> point : side ) Rotate( point, rot[0], rotmat );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Pippo::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double> point : rot ) point += Delta;
    cent += Delta;
    for( vector<double> point : side ) point += Delta;

    // Rotate
    for( vector<double> point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double> vec : trans ) Rotate( vec, Rot0, Rotmat );
    Rotate( cent, Rot0, Rotmat );
    for( vector<double> point : side ) Rotate( point, Rot0, Rotmat );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
}




// Primes the body to be checked (projects the center of the sphere onto the surface)
void Capsule::Prime( vector<double> p, vector<double> v  ) {
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
double Capsule::Anal( vector<double> v, double bodyvel ) {
    double L = abs((l1 - l2)*v)/Norm(v);
    double surface = M_PI*rad*rad + L*rad;
    return Norm(v)*surface/bodyvel;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Capsule::Move( double T ) {
    if( T == t ) return;

    SuperBody->Move(T);

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( long unsigned int i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double> point : rot )  point += delta;
    l1 += delta;
    l2 += delta;

    // Generate rotation matrix
    double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
    vector<vector<double>> rotmat = RotMat( rot[1]-rot[0], theta );

    // Rotate
    Rotate( l1, rot[0], rotmat);
    Rotate( l2, rot[0], rotmat);

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Capsule::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double> point : rot ) point += Delta;
    l1 += Delta;
    l2 += Delta;

    // Rotate
    for( vector<double> point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double> vec : trans ) Rotate( vec, Rot0, Rotmat );
    Rotate( l1, Rot0, Rotmat );
    Rotate( l2, Rot0, Rotmat );
    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
}

// Adds a sub-body
void Body::AddSubBody(Body* subbody) { 
    SubBodies.push_back(subbody); 
    subbody->SetSuperBody(this);
}




// Complete ManyBody constructor 
ManyBody::ManyBody( vector<Sphere> Spheres, vector<Pippo> Pippos, vector<Capsule> Capsules ): Body() {
    spheres = Spheres;
    pippos = Pippos;
    capsules = Capsules;
}

// Primes the body to be checked. Primes each body
void ManyBody::Prime( vector<double> p, vector<double> v  ) {
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

// Time evolution of the body
void ManyBody::Move( double T ) {
    if( T == t ) return;
    // Move parts
    for( Sphere sphere : spheres ) sphere.Move(T);
    for( Pippo pippo : pippos ) pippo.Move(T);
    for( Capsule capsule : capsules ) capsule.Move(T);
}

// Gets stuff
vector<double> ManyBody::GetSphCent( unsigned int index ){
    if( index >= spheres.size() ){
        cout << "No sphere with index " << index << endl;
        return {}; 
    }
    return spheres[index].GetCent();
}

double ManyBody::GetSphRad( unsigned int index ){
    if( index >= spheres.size() ){
        cout << "No sphere with index " << index << endl;
        return 0; 
    }
    return spheres[index].GetRad();
}

vector<double> ManyBody::GetPipCent( unsigned int index ){
    if( index >= pippos.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return pippos[index].GetCent();
}

vector<vector<double>>  ManyBody::GetPipSide( unsigned int index ) {
    if( index >= pippos.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return pippos[index].GetSide();
}

vector<double> ManyBody::GetCapL1( unsigned int index ){
    if( index >= capsules.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return capsules[index].GetL1();
}
	
vector<double> ManyBody::GetCapL2( unsigned int index ){
    if( index >= capsules.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return {}; 
    }
    return capsules[index].GetL2();
}
	
double ManyBody::GetCapRad( unsigned int index ){
    if( index >= capsules.size() ){
        cout << "No parallelepiped with index " << index << endl;
        return 0; 
    }
    return capsules[index].GetRad();
}
