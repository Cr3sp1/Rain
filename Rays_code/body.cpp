#include "body.h"
#include "ray.h"

using namespace std;



// Virtual destructor
Body::~Body() { 
    // cout << "Destroying " << name << endl; 
}

// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
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

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );

    t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Body::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) vec = Rotmat*vec;

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
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

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );

    t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Sphere::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;
    cent += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) vec = Rotmat*vec;
    Rotate( cent, Rot0, Rotmat );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
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
    for( vector<double>& point : side ) point = rotmat*point;
    Rotate( cent, rot[0], rotmat );

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );

    t = T;
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Pippo::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;
    cent += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) vec = Rotmat*vec;
    for( vector<double>& point : side ) point = Rotmat*point; 
    Rotate( cent, Rot0, Rotmat );


    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
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

    // Rotate
    Rotate( l1, rot[0], rotmat);
    Rotate( l2, rot[0], rotmat);

    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( delta, rot[0], rotmat );
}

// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodie
void Capsule::BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) {
    // Translate
    for( vector<double>& point : rot ) point += Delta;
    l1 += Delta;
    l2 += Delta;

    // Rotate
    for( vector<double>& point : rot ) Rotate( point, Rot0, Rotmat );
    for( vector<double>& vec : trans ) vec = Rotmat*vec;
    Rotate( l1, Rot0, Rotmat );
    Rotate( l2, Rot0, Rotmat );
    // Move sub-bodies
    for( Body* body : SubBodies ) body->BeMoved( Delta, Rot0, Rotmat );
}




// Complete ManyBody constructor 
ManyBody::ManyBody( const vector<Sphere>& Spheres, const vector<Pippo>& Pippos, const vector<Capsule>& Capsules ): Body() {
    for( Sphere sphere : Spheres ) AddBody(sphere);
    for( Pippo pippo : Pippos ) AddBody(pippo);
    for( Capsule capsule : Capsules ) AddBody(capsule);
}

// Constructor from file
ManyBody::ManyBody( string filename ): Body() {
    ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    string line;
    // Skip header
    for( size_t i = 0; i < 6; i++ ) getline( file, line );

    // Get bodies
    while( getline( file, line ) ) {
        istringstream iss(line);
        string type, name, superbody;
        iss >> type;
    
        if( type == "Sphere" ) {
            vector<vector<double>> rot, trans;
            vector<double> cent(3);
            double rad, w;
            size_t trans_size;
            for( size_t i = 0; i < 3; i++ ) iss >> cent[i];
            iss >> rad;

            iss >> name;
            iss >> superbody;
            iss >> w;
            if( w != 0 ) {
                rot = {{0,0,0}, {0,0,0}};
                for( size_t i = 0; i < 3; i++ ) iss >> rot[0][i];
                for( size_t i = 0; i < 3; i++ ) iss >> rot[1][i];
            }
            iss >> trans_size;
            for( size_t i = 0; i < trans_size; i++ ){
                vector<double> temp(3);
                for( size_t j = 0; j < 3; j++ ) {
                    iss >> temp[j];
                }
                trans.push_back(temp);
            }

            AddBody( Sphere( cent, rad, name, rot, w*M_PI, trans ));
            if( superbody != "None" ) Attach( name, superbody );
        }

        if( type == "Pippo" ) {
            vector<vector<double>> rot, trans, sides;
            vector<double> cent(3);
            double w;
            size_t trans_size;
            for( size_t i = 0; i < 3; i++ ) iss >> cent[i];
            for( size_t i = 0; i < 3; i++ ){
                vector<double> temp(3);
                for( size_t j = 0; j < 3; j++ ) {
                    iss >> temp[j];
                }
                sides.push_back(temp);
            }

            iss >> name;
            iss >> superbody;
            iss >> w;
            if( w != 0 ) {
                rot = {{0,0,0}, {0,0,0}};
                for( size_t i = 0; i < 3; i++ ) iss >> rot[0][i];
                for( size_t i = 0; i < 3; i++ ) iss >> rot[1][i];
            }
            iss >> trans_size;
            for( size_t i = 0; i < trans_size; i++ ){
                vector<double> temp(3);
                for( size_t j = 0; j < 3; j++ ) {
                    iss >> temp[j];
                }
                trans.push_back(temp);
            }

            AddBody( Pippo( cent, sides, name, rot, w*M_PI, trans ));
            if( superbody != "None" ) Attach( name, superbody );
        }

        if( type == "Capsule" ) {
            vector<vector<double>> rot, trans;
            vector<double> cent(3), l1(3), l2(3);
            double rad, w;
            size_t trans_size;
            for( size_t i = 0; i < 3; i++ ) iss >> l1[i];
            for( size_t i = 0; i < 3; i++ ) iss >> l2[i];
            iss >> rad;

            iss >> name;
            iss >> superbody;
            iss >> w;
            if( w != 0 ) {
                rot = {{0,0,0}, {0,0,0}};
                for( size_t i = 0; i < 3; i++ ) iss >> rot[0][i];
                for( size_t i = 0; i < 3; i++ ) iss >> rot[1][i];
            }
            iss >> trans_size;
            for( size_t i = 0; i < trans_size; i++ ){
                vector<double> temp(3);
                for( size_t j = 0; j < 3; j++ ) {
                    iss >> temp[j];
                }
                trans.push_back(temp);
            }

            AddBody( Capsule( l1, l2, rad, name, rot, w*M_PI, trans ));
            if( superbody != "None" ) Attach( name, superbody );
        }
    }

    file.close();
}

// Destructor
ManyBody::~ManyBody() {
    // cout << "Destroying ManyBody" << endl;
    // Delete all dynamically allocated Sphere objects
    for (Sphere* sphere : spheres) {
        delete sphere;
    }

    // Delete all dynamically allocated Pippo objects
    for (Pippo* pippo : pippos) {
        delete pippo;
    }

    // Delete all dynamically allocated Capsule objects
    for (Capsule* capsule : capsules) {
        delete capsule;
    }
}

// Primes the body to be checked. Primes each body
void ManyBody::Prime( vector<double> p, vector<double> v  ) {
    for( size_t i = 0; i < spheres.size(); i++ ) spheres[i]->Prime( p, v );
    for( size_t i = 0; i < pippos.size(); i++ ) pippos[i]->Prime( p, v );
    for( size_t i = 0; i < capsules.size(); i++ ) capsules[i]->Prime( p, v );
}

// Checks if the ManyBody is making contact with a ray
bool ManyBody::Check( Ray& ray ) {
    for( size_t i = 0; i < spheres.size(); i++ ) if(spheres[i]->Check( ray )) return true;
    for( size_t i = 0; i < pippos.size(); i++ ) if(pippos[i]->Check( ray )) return true;
    for( size_t i = 0; i < capsules.size(); i++ ) if(capsules[i]->Check( ray )) return true;
    return false;
}

// Time evolution of the body
void ManyBody::Move( double T ) {
    if( T == t ) return;
    // Move parts
    for( Sphere* sphere : spheres ) sphere->Move(T);
    for( Pippo*  pippo : pippos ) pippo->Move(T);
    for( Capsule* capsule : capsules ) capsule->Move(T);

    t = T;
}

// Pointer to the body with that name
Body* ManyBody::Find( string name ) {
    for( size_t i = 0; i < spheres.size(); i++ ) {
        if( spheres[i]->GetName() == name ) return spheres[i];
    }
    for( size_t i = 0; i < pippos.size(); i++ ) {
        if( pippos[i]->GetName() == name ) return pippos[i];
    }
    for( size_t i = 0; i < capsules.size(); i++ ) {
        if( capsules[i]->GetName() == name ) return capsules[i];
    }
    cout << name << " not found!" << endl;
    return nullptr;
}

// Attaches the sub-body to the super-body
void ManyBody::Attach( Body SubBody, string SuperName ){
    Body* Super = Find(SuperName);
    if(Super) Super->AddSubBody( SubBody );
}

void ManyBody::Attach( string SubName, string SuperName ) { 
    Body* Super = Find(SuperName);
    Body* Sub = Find(SubName);
    if( Super and Sub ) Super->AddSubBody(*Sub);
}

// Prints to file the state (all the bodies and their parameters)
void ManyBody::PrintState( ofstream &fout ) {
    for( Sphere* sphere : spheres ) {
        vector<double> cent = sphere->GetCent();
        double rad = sphere->GetRad();
        fout << "S," << cent[0] << "," << cent[1] << "," << cent[2] << "," << rad << endl;
    }

    for( Pippo* pippo : pippos ) {
        vector<vector<double>> Side = pippo->GetSide();
        vector<double> cent = pippo->GetCent();
        fout << "P," << cent[0] << "," << cent[1] << "," << cent[2] << ",";
        for( size_t i = 0; i < Side.size(); i++ ) {
            fout << Side[i][0] << "," << Side[i][1] << "," << Side[i][2];
            if( i+1 != Side.size() ) fout << ",";
        }
        fout << endl;
    }

    for( Capsule* capsule : capsules ) {
        vector<double> l1 = capsule->GetL1();
        vector<double> l2 = capsule->GetL2();
        double rad = capsule->GetRad();
        fout << "C," << l1[0] << "," << l1[1] << "," << l1[2] << "," << l2[0] << "," << l2[1] << "," << l2[2] << "," << rad << endl;
    }
    
}

void ManyBody::PrintState( string outfile ) {
    ofstream fout(outfile);
    PrintState( fout );
    fout.close();
} 