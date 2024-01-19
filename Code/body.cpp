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
    // Add sin terms
    for( size_t i = 0; 2*i < trans.size(); i++ ){
        delta += trans[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Add cos terms
    for( size_t i = 0; 2*i+1 < trans.size(); i++ ){
        delta += trans[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double>& point : rot ){
        point += delta;
    }

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w.size() > 0 and rot.size() == 2 ) {
        double theta = 0;
        // Add sin terms
        for( size_t i = 0; 2*i < w.size(); i++ ){
        theta += w[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
        }
        // Add cos terms
        for( size_t i = 0; 2*i+1 < w.size(); i++ ){
            theta += w[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
        }
        rotmat =  RotMat( rot[1]-rot[0], theta );
    } else {
        rotmat = IdMat(3);
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

// Finds smallest box around body
void Body::FindBox( vector<double> &min, vector<double> &max ) {
    if( min.size() != 3 or max.size() != 3 ) {
        cout << "Error in FindBox: min and max must be of size 3" << endl;
    }
}

// Prints to file the state of the body
void Body::PrintState( ofstream &fout ) {}

void Body::PrintState( string outfile ) {
    ofstream fout(outfile);
    PrintState( fout );
    fout.close();
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
    // Add sin terms
    for( size_t i = 0; 2*i < trans.size(); i++ ){
        delta += trans[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Add cos terms
    for( size_t i = 0; 2*i+1 < trans.size(); i++ ){
        delta += trans[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double>& point : rot ) point += delta;
    cent += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w.size() > 0 and rot.size() == 2 ) {
        double theta = 0;
        // Add sin terms
        for( size_t i = 0; 2*i < w.size(); i++ ){
        theta += w[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
        }
        // Add cos terms
        for( size_t i = 0; 2*i+1 < w.size(); i++ ){
            theta += w[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
        }
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

// Finds smallest box around body
void Sphere::FindBox( vector<double> &min, vector<double> &max ) {
    if( min.size() != 3 or max.size() != 3 ) {
        cout << "Error in FindBox: min and max must be of size 3" << endl;
    }
    for( size_t i = 0; i < 3; i++ ) {
        if( cent[i] - rad < min[i] ) min[i] = cent[i] - rad;
        if( cent[i] + rad > max[i] ) max[i] = cent[i] + rad;
    }
}

// Prints to file the state of the body
void Sphere::PrintState( ofstream &fout ) {
    fout << setprecision(4);
    fout << "S,"  << cent[0] << "," << cent[1] << "," << cent[2] << "," << rad << endl;
}

void Sphere::PrintState( string outfile ) {
    ofstream fout(outfile);
    PrintState( fout );
    fout.close();
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
    // Add sin terms
    for( size_t i = 0; 2*i < trans.size(); i++ ){
        delta += trans[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Add cos terms
    for( size_t i = 0; 2*i+1 < trans.size(); i++ ){
        delta += trans[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
    }
    // Translate
    for( vector<double>& point : rot ) point += delta;
    cent += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w.size() > 0 and rot.size() == 2 ) {
        double theta = 0;
        // Add sin terms
        for( size_t i = 0; 2*i < w.size(); i++ ){
        theta += w[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
        }
        // Add cos terms
        for( size_t i = 0; 2*i+1 < w.size(); i++ ){
            theta += w[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
        }
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
                vertices.push_back( cent + i*side[0]/2. + j*side[1]/2. + k*side[2]/2. );
            }
        }
    }
    return vertices;
}

// Finds smallest box around body
void Pippo::FindBox( vector<double> &min, vector<double> &max ) {
    if( min.size() != 3 or max.size() != 3 ) {
        cout << "Error in FindBox: min and max must be of size 3" << endl;
    }
    vector<vector<double>> vertices = GetVertices();
    for( vector<double> vertex : vertices ){
        for( size_t i = 0; i < 3; i++ ) {
            if( vertex[i] < min[i] ) min[i] = vertex[i];
            if( vertex[i] > max[i] ) max[i] = vertex[i];
        }
    }
    
}


// Prints to file the state of the body
void Pippo::PrintState( ofstream &fout ) {
    fout << setprecision(4);
    fout << "P," << cent[0] << "," << cent[1] << "," << cent[2] << ",";
        for( size_t i = 0; i < side.size(); i++ ) {
            fout << side[i][0] << "," << side[i][1] << "," << side[i][2];
            if( i+1 != side.size() ) fout << ",";
        }
        fout << endl;
}

void Pippo::PrintState( string outfile ) {
    ofstream fout(outfile);
    PrintState( fout );
    fout.close();
} 




// Primes the body to be checked (projects l1 and l2 onto the surface)
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
    // Add sin terms
    for( size_t i = 0; 2*i < trans.size(); i++ ){
        delta += trans[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
    }
    // Add cos terms
    for( size_t i = 0; 2*i+1 < trans.size(); i++ ){
        delta += trans[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
    }
    // Translate 
    for( vector<double>& point : rot )  point += delta;
    l1 += delta;
    l2 += delta;

    // Generate rotation matrix
    vector<vector<double>> rotmat;
    if( w.size() > 0 and rot.size() == 2 ) {
        double theta = 0;
        // Add sin terms
        for( size_t i = 0; 2*i < w.size(); i++ ){
        theta += w[2*i]*( sin(T*2*M_PI/(i+1)) - sin(t*2*M_PI/(i+1) ));
        }
        // Add cos terms
        for( size_t i = 0; 2*i+1 < w.size(); i++ ){
            theta += w[2*i+1]*( cos(T*2*M_PI/(i+1)) - cos(t*2*M_PI/(i+1) ));
        }
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

    t = T;
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

// Finds smallest box around body
void Capsule::FindBox( vector<double> &min, vector<double> &max ) {
    if( min.size() != 3 or max.size() != 3 ) {
        cout << "Error in FindBox: min and max must be of size 3" << endl;
    }
    for( size_t i = 0; i < 3; i++ ) {
        if( l1[i] - rad < min[i] ) min[i] = l1[i] - rad;
        if( l1[i] + rad > max[i] ) max[i] = l1[i] + rad;
    }
    for( size_t i = 0; i < 3; i++ ) {
        if( l2[i] - rad < min[i] ) min[i] = l2[i] - rad;
        if( l2[i] + rad > max[i] ) max[i] = l2[i] + rad;
    }
}

// Prints to file the state of the body
void Capsule::PrintState( ofstream &fout ) {
    fout << setprecision(4);
    fout << "C," << l1[0] << "," << l1[1] << "," << l1[2] << "," << l2[0] << "," << l2[1] << "," << l2[2] << "," << rad << endl;
}

void Capsule::PrintState( string outfile ) {
    ofstream fout(outfile);
    PrintState( fout );
    fout.close();
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

    // Get bodies
    while( getline( file, line ) ) {
        istringstream iss(line);
        string type, name, superbody;
        iss >> type;
    
        if( type == "Sphere" ) {
            vector<vector<double>> rot, trans;
            vector<double> cent(3);
            double rad;
            size_t trans_size, w_size;
            for( size_t i = 0; i < 3; i++ ) iss >> cent[i];
            iss >> rad;

            iss >> name;
            iss >> superbody;
            iss >> w_size;
            vector<double> w(w_size);
            for( size_t i = 0; i < w_size; i++ ) iss >> w[i];
            w = w*(M_PI/180);
            if( w_size != 0 ) {
                rot = {{0,0,0}, {0,0,0}};
                for( size_t i = 0; i < 3; i++ ) iss >> rot[0][i];
                vector<double> axis(3);
                for( size_t i = 0; i < 3; i++ ) iss >> axis[i];
                rot[1] = rot[0] + axis;
            }
            iss >> trans_size;
            for( size_t i = 0; i < trans_size; i++ ){
                vector<double> temp(3);
                for( size_t j = 0; j < 3; j++ ) {
                    iss >> temp[j];
                }
                trans.push_back(temp);
            }

            AddBody( Sphere( cent, rad, name, rot, w, trans ));
            if( superbody != "None" ) Attach( name, superbody );
        }

        if( type == "Pippo" ) {
            vector<vector<double>> rot, trans, sides;
            vector<double> cent(3);
            size_t trans_size, w_size;
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
            iss >> w_size;
            vector<double> w(w_size);
            for( size_t i = 0; i < w_size; i++ ) iss >> w[i];
            w = w*(M_PI/180);
            if( w_size != 0 ) {
                rot = {{0,0,0}, {0,0,0}};
                for( size_t i = 0; i < 3; i++ ) iss >> rot[0][i];
                vector<double> axis(3);
                for( size_t i = 0; i < 3; i++ ) iss >> axis[i];
                rot[1] = rot[0] + axis;
            }
            iss >> trans_size;
            for( size_t i = 0; i < trans_size; i++ ){
                vector<double> temp(3);
                for( size_t j = 0; j < 3; j++ ) {
                    iss >> temp[j];
                }
                trans.push_back(temp);
            }

            AddBody( Pippo( cent, sides, name, rot, w, trans ));
            if( superbody != "None" ) Attach( name, superbody );
        }

        if( type == "Capsule" ) {
            vector<vector<double>> rot, trans;
            vector<double> cent(3), l1(3), l2(3);
            double rad;
            size_t trans_size, w_size;
            for( size_t i = 0; i < 3; i++ ) iss >> l1[i];
            for( size_t i = 0; i < 3; i++ ) iss >> l2[i];
            iss >> rad;

            iss >> name;
            iss >> superbody;
            iss >> w_size;
            vector<double> w(w_size);
            for( size_t i = 0; i < w_size; i++ ) iss >> w[i];
            w = w*(M_PI/180);
            if( w_size != 0 ) {
                rot = {{0,0,0}, {0,0,0}};
                for( size_t i = 0; i < 3; i++ ) iss >> rot[0][i];
                vector<double> axis(3);
                for( size_t i = 0; i < 3; i++ ) iss >> axis[i];
                rot[1] = rot[0] + axis;
            }
            iss >> trans_size;
            for( size_t i = 0; i < trans_size; i++ ){
                vector<double> temp(3);
                for( size_t j = 0; j < 3; j++ ) {
                    iss >> temp[j];
                }
                trans.push_back(temp);
            }

            AddBody( Capsule( l1, l2, rad, name, rot, w, trans ));
            if( superbody != "None" ) Attach( name, superbody );
        }
    }

    file.close();
}

// Destructor
ManyBody::~ManyBody() {
    // cout << "Destroying ManyBody" << endl;
    // Delete all dynamically allocated Body objects
    for (Body* body : bodies) {
        delete body;
    }
}

// Primes the body to be checked. Primes each body
void ManyBody::Prime( vector<double> p, vector<double> v  ) {
    for(Body* body : bodies) body->Prime( p, v );
}

// Checks if the ManyBody is making contact with a ray
bool ManyBody::Check( Ray& ray ) {
    for(Body* body : bodies) if(body->Check( ray )) return true;
    return false;
}

// Time evolution of the body
void ManyBody::Move( double T ) {
    if( T == t ) return;
    // Move parts
    for(Body* body : bodies) body->Move(T);

    t = T;
}

// Returns a pointer to the body with that name in the ManyBody
Body* ManyBody::Find( string name ) {
    for(Body* body : bodies) if( body->GetName() == name ) return body;

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

// Finds smallest box around body
void ManyBody::FindBox( vector<double> &min, vector<double> &max ) {
    if( min.size() != 3 or max.size() != 3 ) {
        cout << "Error in FindBox: min and max must be of size 3" << endl;
    }
    for( Body* body : bodies ) {
        body->FindBox( min, max );
    }
}

// Prints to file the state (all the bodies and their parameters)
void ManyBody::PrintState( ofstream &fout ) {
    for( Body* body : bodies ) {
        body->PrintState(fout);
    }
}

void ManyBody::PrintState( string outfile ) {
    ofstream fout(outfile);
    PrintState( fout );
    fout.close();
} 