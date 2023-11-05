#include "RainFunctions.h"
#include "body.h"
#include "ray.h"

using namespace std;


// Projects the Point on a plane perpendicular to v and passing through p
vector<double> Project( vector<double> Point, vector<double> p, vector<double> v ){
    vector<double> diff = p - Point;
    return ( Point + v*(diff*v)/(v*v) );
}

// Finds the vertex in the middle of the three seen faces of a parallelepiped defined by a pont and three sides 
vector<double> FindMiddle( vector<double> p, vector<vector<double>> sides, vector<double> v ){
    for( unsigned int i = 0; i < sides.size(); i++ ){
        if( sides[i]*v < 0 ) p += sides[i];
    }
    return p;
}

// Finds the hexagonal projection H of a parallelepiped defined by p nad sides on a plane perpendicular to v and passing through P0
vector<vector<double>> FindHexProj(  vector<double> p, vector<vector<double>> Side, vector<double> v, vector<double> P0){
    vector<vector<double>> H = {FindMiddle( p, Side, v )};
    vector<vector<double>> delta(3, vector<double>(3, 0.0));        // Used to calculate the position of the vertices to project
    for( int i = 0; i < 3; i++ ){
        delta[i] = Side[i]*v < 0 ? ((double)-1)*Side[i] : Side[i];
    }
    H.push_back( H[0] + delta[0] );
    H.push_back( H[0] + delta[0] + delta[1]) ;
    H.push_back( H[0] + delta[1] );
    H.push_back( H[0] + delta[1] + delta[2] );
    H.push_back( H[0] + delta[2] );
    H.push_back( H[0] + delta[2] + delta[0] );
    for( int i = 0; i < 7; i++ ){
        H[i] = Project( H[i], P0, v );         // We project them
    }

    return H;
}

// Returns the highest absolute value of the projections of the vertices of H on a line in direction u1 passing through H[0]
double MaxU(vector<vector<double>> H, vector<double> u ) {
    double result = 0;
    for( int i = 1; i < 7; i++ ) {
        H[i] -= H[0];
        double proj = abs(H[i]*u/Norm(u));
        if(proj > result) result = proj;
    }
    // cout << result << endl;
    return result;
}



// Returns wether the Point is inside the hexagon H using triangles
bool PointIsInsideT( vector<double> Point, vector<vector<double>> H ){
    // Centers on p
    Point -= H[0];
    for(int i = 1; i < 7; i++ ){
        H[i] -= H[0];
    }

    // Checks if Point is inside the rectangle with vertices H[0], H[i], H[i+1]
    for( int i = 1; i < 7; i++ ){
        double epsilon = 1e-10;
        double A = Norm( CrossProduct( H[i], H[PBCH(i+1)]) );
        double alpha = Norm( CrossProduct( Point, H[PBCH(i+1)]) )/A;
        double beta = Norm( CrossProduct( Point, H[i]) )/A;
        double gamma = Norm( CrossProduct( Point-H[i], Point-H[PBCH(i+1)]) )/A;
        if( 0 <= alpha and alpha <= 1 and 0 <= beta and beta <= 1 and 0 <= gamma and gamma <= 1 and abs(alpha + beta + gamma - 1 ) <  epsilon ){
            return true;
        }
    }
    return false;
}


// Periodic Boundary conditions for the index of H, keeps it between 1 and 6
int PBCH( int i ) {
    while( i > 6 ) i -= 6;
    while( i < 1 ) i += 6;
    return i;
}


// Checks rays generation 
void RayGenCheck( string outfile, vector<double> box, vector<double> rel_vel ){
    ofstream Pout("outfile");
    for( int i = 0; i < 10000; i+=50 ){
        ProjSurface temp( box, rel_vel, i+1 );
        Pout << (double)temp.GetNRays()/(i+1) << endl;
    }
    Pout.close();
}


// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<double> body_v(N);
    vector<double> wetness(N);
    for( unsigned int i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        vector<double> relvel = rain_v;
        relvel[0] -= body_v[i];
        wetness[i] = Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body)/body_v[i];
    }
    return {body_v, wetness};
}


// Estimates wetness for N velocities of the body between vmin and vmax (measured as fractions of vertical rain speed), and returns a matrix with the velocities as the first colunmn and the respective theorical wetness as the second column and the estimated wetness as the third
vector<vector<double>> CompareAN( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vmin*=-rain_v[2];
    vmax*=-rain_v[2];
    
    vector<double> body_v(N);
    vector<double> analytical(N);
    vector<double> wetness(N);
    for( unsigned int i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        vector<double> relvel = rain_v;
        relvel[0] -= body_v[i];
        analytical[i] = body.Anal( relvel, body_v[i]);
        wetness[i] = Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body)/body_v[i];
    }
    
    return {body_v, analytical,  wetness};
}

// Estimates wetness for N velocities of two body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the wetness of the first body as the second column and of the second body as the third column
vector<vector<double>> CompareBB( vector<double> box, Body& body1, Body& body2, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx){
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<double> body_v(N);
    vector<double> wetness1(N);
    vector<double> wetness2(N);
    for( unsigned int i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        vector<double> relvel = rain_v;
        relvel[0] -= body_v[i];
        wetness1[i] = Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body1)/body_v[i];
        wetness2[i] = Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body2)/body_v[i];
    }
    return {body_v, wetness1, wetness2};
}



// Returns the minimum distance between the point p and the segment line with extremes l1 and l2
double PointSegDist( vector<double> p, vector<double> l1, vector<double> l2 ) {
    // Changes frame of reference to l1 = 0
    p -= l1;
    l2 -= l1;
    // Calculates projection of p on line passing through l1 and l2
    double proj = p*l2/Norm(l2);
    // Returns distance between p and the closest point belonging to the segment
    if( proj <= 0 )  return Norm(p);
    if( proj >= Norm(l2) ) return Norm(p-l2);
    return Norm( p - l2*proj);
}



// Returns NxN identity matrix
vector<vector<double>> IdMat( unsigned int N ) {
    vector<vector<double>> idmat(N, vector<double>(N));
    for( unsigned int i = 0; i < N; i++ ){
        for( unsigned int j = 0; j < N; j++ ){
            idmat[i][j] = (i==j) ? 1 : 0; 
        }
    }
    return idmat;
}



// Returns the rotation matrix
vector<vector<double>> RotMat( vector<double> axis, double theta ) {
    if ( Norm(axis) == 0 or axis.size() != 3 ) return {}; // Handle zero-length vector to avoid division by zero
    axis = axis/Norm(axis);
    double s = sin(theta);
    double c = cos(theta);
    double G = 1-c;

    return { { axis[0]*axis[0]*G + c,           axis[0]*axis[1]*G - axis[2]*s,  axis[0]*axis[2]*G + axis[1]*s },
             { axis[1]*axis[0]*G + axis[2]*s,   axis[1]*axis[1]*G + c,          axis[1]*axis[2]*G - axis[0]*s },
             { axis[2]*axis[0]*G - axis[1]*s,   axis[2]*axis[1]*G + axis[0]*s,  axis[2]*axis[2]*G + c         } };
}



void Rotate( vector<double>& Point, vector<double> Rot0, const vector<vector<double>>& Rotmat ){
    Point -= Rot0;
    Point = Rotmat*Point;
    Point += Rot0; 
}