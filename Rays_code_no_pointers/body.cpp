#include "body.h"
#include "ray.h"

using namespace std;


// Time evolution of the body in its own frame of reference
void Body::Move( double T ) {
    if( T == t ) return;

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double>& point : rot ){
        point += delta;
    }

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }

    t = T;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Body::Move( double T, ManyBody& many ) {
    if( T == t ) return;

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double>& point : rot ){
        point += delta;
    }

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( delta, rot[0], rotmat, many );
    }

    t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Body::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat, ManyBody& many ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) Rotate( vec, Rot0, Rotmat );

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( Delta, Rot0, Rotmat, many );
    }
}

// Attaches the Body to a SuperBody
void Body::AttachTo( Body& SupBody ){
    SuperBody = SupBody.GetName();
    SupBody.AddSubBody(*this);
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

// Time evolution of the body in its own frame of reference
void Sphere::Move( double T ) {
    if( T == t ) return;

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double>& point : rot ) point += delta;
    cent += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }

    // Rotate
    Rotate( cent, rot[0], rotmat );

    t = T;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Sphere::Move( double T, ManyBody& many ) {
    if( T == t ) return;

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double>& point : rot ) point += delta;
    cent += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }

    // Rotate
    Rotate( cent, rot[0], rotmat );

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( delta, rot[0], rotmat, many );
    }

    t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Sphere::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat, ManyBody& many ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;
    cent += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) Rotate( vec, Rot0, Rotmat );
    Rotate( cent, Rot0, Rotmat );

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( Delta, Rot0, Rotmat, many );
    }
}
	



// Primes the body to be checked (finds hexagonal projection on the same plane as the origins of the rays)
void Pippo::Prime( vector<double> P, vector<double> V ) {
    vector<double> p = cent - (double)0.5*(side[0] + side[1] + side[2]);
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

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double>& point : rot ) point += delta;
    cent += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }

    // Rotate
    for( vector<double>& point : side ) {
        point += cent;
        Rotate( point, rot[0], rotmat );
    }
    Rotate( cent, rot[0], rotmat );
    for( vector<double>& point : side ) point -= cent;

    t = T;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Pippo::Move( double T, ManyBody& many  ) {
    if( T == t ) return;

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double>& point : rot ) point += delta;
    cent += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }

    // Rotate
    for( vector<double>& point : side ) {
        point += cent;
        Rotate( point, rot[0], rotmat );
    }
    Rotate( cent, rot[0], rotmat );
    for( vector<double>& point : side ) point -= cent;

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( delta, rot[0], rotmat, many );
    }

    t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Pippo::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat, ManyBody& many ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;
    cent += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) Rotate( vec, Rot0, Rotmat );
    for( vector<double>& point : side ) {
        point += cent;
        Rotate( point, Rot0, Rotmat );
    }
    Rotate( cent, Rot0, Rotmat );
    for( vector<double>& point : side ) point -= cent;

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( Delta, Rot0, Rotmat, many );
    }
}

// Returns all 8 vertices of the parallelepiped
vector<vector<double>>  Pippo::GetVertices() {
    vector<vector<double>> vertices{};
    for( double i = -0.5; i <= 0.5; i ++ ) {
        for( double j = -0.5; j <= 0.5; j ++ ) {
            for( double k = -0.5; k <= 0.5; k ++ ){
                vertices.push_back( cent + i*side[0] + j*side[1] + k*side[2] );
            }
        }
    }
    return vertices;
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
    // double L = abs((l1 - l2)*v)/Norm(v);
    vector<double> axis = l1 - l2;
    axis -= v*axis*v/(v*v);
    double L = Norm(axis);
    double surface = M_PI*rad*rad + L*2*rad;
    return Norm(v)*surface/bodyvel;
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Capsule::Move( double T ) {
    if( T == t ) return;

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double>& point : rot )  point += delta;
    l1 += delta;
    l2 += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }
}

 // Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
void Capsule::Move( double T, ManyBody& many ) {
    if( T == t ) return;

    // Calculate the total translation for the step
    vector<double> delta({0,0,0});
    for( size_t i = 0; i < trans.size(); i++ ){
        delta += trans[i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double>& point : rot )  point += delta;
    l1 += delta;
    l2 += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w != 0 and rot.size() == 2 ) {
        double theta = w*( sin(T*2*M_PI) - sin(t*2*M_PI) );
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
        rot = {{0,0,0}};
    }

    // Rotate
    Rotate( l1, rot[0], rotmat);
    Rotate( l2, rot[0], rotmat);

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( delta, rot[0], rotmat, many );
    }
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Capsule::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat, ManyBody& many ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;
    l1 += Delta;
    l2 += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) Rotate( vec, Rot0, Rotmat );
    Rotate( l1, Rot0, Rotmat );
    Rotate( l2, Rot0, Rotmat );

    // Moves sub-bodies
    for( string& SubBody : SubBodies ){
        many.Find( SubBody)->BeMoved( Delta, Rot0, Rotmat, many );
    }
}




// Complete ManyBody constructor 
ManyBody::ManyBody( const vector<Sphere>& Spheres, const vector<Pippo>& Pippos, const vector<Capsule>& Capsules ): Body() {
    spheres = Spheres;
    pippos = Pippos;
    capsules = Capsules;
}

// Primes the body to be checked. Primes each body
void ManyBody::Prime( vector<double> p, vector<double> v  ) {
    for( size_t i = 0; i < spheres.size(); i++ ) spheres[i].Prime( p, v );
    for( size_t i = 0; i < pippos.size(); i++ ) pippos[i].Prime( p, v );
    for( size_t i = 0; i < capsules.size(); i++ ) capsules[i].Prime( p, v );
}

// Checks if the ManyBody is making contact with a ray
bool ManyBody::Check( Ray& ray ) {
    for( size_t i = 0; i < spheres.size(); i++ ) if(spheres[i].Check( ray )) return true;
    for( size_t i = 0; i < pippos.size(); i++ ) if(pippos[i].Check( ray )) return true;
    for( size_t i = 0; i < capsules.size(); i++ ) if(capsules[i].Check( ray )) return true;
    return false;
}

// Time evolution of the body
void ManyBody::Move( double T ) {
    if( T == t ) return;

    // Move parts
    for( Sphere& sphere : spheres ) sphere.Move(T, *this);
    for( Pippo& pippo : pippos ) pippo.Move(T, *this);
    for( Capsule& capsule : capsules ) capsule.Move(T, *this);
    t = T;
}

// Attaches the sub-body to the super-body
void ManyBody::Attach( Body SubBody, string SuperName ){
    for( Sphere& sphere : spheres ){
        if( sphere.GetName() == SuperName ) sphere.AddSubBody( SubBody );
    }
    for( Pippo& pippo : pippos ){
        if( pippo.GetName() == SuperName ) pippo.AddSubBody( SubBody );
    }
    for( Capsule& capsule : capsules ){
        if( capsule.GetName() == SuperName ) capsule.AddSubBody( SubBody );
    }
}

void ManyBody::Attach( string SubName, string SuperName ){
    for( Sphere& sphere : spheres ){
        if( sphere.GetName() == SubName ) Attach( sphere, SuperName );
    }
    for( Pippo& pippo : pippos ){
        if( pippo.GetName() == SubName ) Attach( pippo, SuperName );
    }
    for( Capsule& capsule : capsules ){
        if( capsule.GetName() == SubName ) Attach( capsule, SuperName );
    }
}

// Pointer to the body with that name
Body* ManyBody::Find( string name ) {
    for( size_t i = 0; i < spheres.size(); i++ ) {
        if( spheres[i].GetName() == name ) return &spheres[i];
    }
    for( size_t i = 0; i < pippos.size(); i++ ) {
        if( pippos[i].GetName() == name ) return &pippos[i];
    }
    for( size_t i = 0; i < capsules.size(); i++ ) {
        if( capsules[i].GetName() == name ) return &capsules[i];
    }
    return nullptr;
}

// Prints to file the state (all the bodies and their parameters)
void ManyBody::PrintState( ofstream &fout ) {
    for( Sphere& sphere : spheres ) {
        vector<double> cent = sphere.GetCent();
        double rad = sphere.GetRad();
        fout << "S," << cent[0] << "," << cent[1] << "," << cent[2] << "," << rad << endl;
    }

    for( Pippo& pippo : pippos ) {
        vector<vector<double>> Side = pippo.GetSide();
        vector<double> cent = pippo.GetCent();
        fout << "P," << cent[0] << "," << cent[1] << "," << cent[2] << ",";
        for( size_t i = 0; i < Side.size(); i++ ) {
            fout << Side[i][0] << "," << Side[i][1] << "," << Side[i][2];
            if( i+1 != Side.size() ) fout << ",";
        }
        fout << endl;
    }

    for( Capsule& capsule : capsules ) {
        vector<double> l1 = capsule.GetL1();
        vector<double> l2 = capsule.GetL2();
        double rad = capsule.GetRad();
        fout << "C," << l1[0] << "," << l1[1] << "," << l1[2] << "," << l2[0] << "," << l2[1] << "," << l2[2] << "," << rad << endl;
    }
    
}

void ManyBody::PrintState( string outfile ) {
    ofstream fout(outfile);
    PrintState( fout );
    fout.close();
} 