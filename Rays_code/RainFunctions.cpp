#include "RainFunctions.h"
#include "body.h"
#include "ray.h"

using namespace std;


// Projects the Point on a plane perpendicular to v and passing through p
vector<long double> Project( vector<long double> Point, vector<long double> p, vector<long double> v ){
    vector<long double> diff = p - Point;
    return ( Point + v*(diff*v)/(v*v) );
}

// Finds the vertex in the middle of the three seen faces of a parallelepiped defined by a pont and three sides 
vector<long double> FindMiddle( vector<long double> p, vector<vector<long double>> sides, vector<long double> v ){
    for( unsigned int i = 0; i < sides.size(); i++ ){
        if( sides[i]*v < 0 ) p += sides[i];
    }
    return p;
}

// Finds the hexagonal projection H of a parallelepiped defined by p nad sides on a plane perpendicular to v and passing through P0
vector<vector<long double>> FindHexProj(  vector<long double> p, vector<vector<long double>> Side, vector<long double> v, vector<long double> P0){
    vector<vector<long double>> H = {FindMiddle( p, Side, v )};
    vector<vector<long double>> delta(3, vector<long double>(3, 0.0));        // Used to calculate the position of the vertices to project
    for( int i = 0; i < 3; i++ ){
        delta[i] = Side[i]*v < 0 ? ((long double)-1)*Side[i] : Side[i];
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
long double MaxU(vector<vector<long double>> H, vector<long double> u ) {
    long double result = 0;
    for( int i = 1; i < 7; i++ ) {
        H[i] -= H[0];
        long double proj = abs(H[i]*u/Norm(u));
        if(proj > result) result = proj;
    }
    // cout << result << endl;
    return result;
}



// Returns wether the Point is inside the hexagon H using triangles
bool PointIsInsideT( vector<long double> Point, vector<vector<long double>> H ){
    // Centers on p
    Point -= H[0];
    for(int i = 1; i < 7; i++ ){
        H[i] -= H[0];
    }

    // Checks if Point is inside the reiangle with vertices H[0], H[i], H[i+1]
    for( int i = 1; i < 7; i++ ){
        long double epsilon = 1e-10;
        long double A = Norm( CrossProduct( H[i], H[PBCH(i+1)]) );
        long double alpha = Norm( CrossProduct( Point, H[PBCH(i+1)]) )/A;
        long double beta = Norm( CrossProduct( Point, H[i]) )/A;
        long double gamma = Norm( CrossProduct( Point-H[i], Point-H[PBCH(i+1)]) )/A;
        if( 0 <= alpha and alpha <= 1 and 0 <= beta and beta <= 1 and 0 <= gamma and gamma <= 1 and abs(alpha + beta + gamma - 1 ) <  epsilon ){
            return true;
        }
    }
    return false;
}

// NOT WORKING PROPERLY: Returns wether the Point is inside the hexagon H using parallelograms 
// bool PointIsInsideP( vector<long double> Point, vector<vector<long double>> H ){
//     // centers the points in p
//     Point -= H[0];
//     for(int i = 0; i < 7; i++ ) H[i] -=H[0];
//     // Calculate stuff
//     long double a, b;
//     long double epsilon = 1e-16;
//     long double H1p = H[1]*Point;
//     long double H1H1 = H[1]*H[1];
//     long double H1H3 = H[1]*H[3];
//     long double H3p = H[3]*Point;
//     long double H3H3 = H[3]*H[3];
//     long double H3H5 = H[3]*H[4];
//     long double H5p = H[5]*Point;
//     long double H5H5 = H[5]*H[5];
//     long double H5H1 = H[5]*H[1];

//     if( H1H1*H3H3 - H1H3*H1H3 > epsilon and H1H1 > epsilon ){
//         b = (H3p*H1H1 - H1p*H1H3)/(H1H1*H3H3 - H1H3*H1H3);
//         a = (H1p*H3H3 - H3p*H1H3)/(H1H1*H3H3 - H1H3*H1H3);
//         if( 0 <= a and a <= 1 and 0 <= b and b <= 1 ) return true;
//     }

//     if( H3H3*H5H5 - H3H5*H3H5 > epsilon and H3H3 > epsilon ){
//         b = (H5p*H3H3 - H3p*H3H5)/(H3H3*H5H5 - H3H5*H3H5);
//         a = (H3p*H5H5 - H5p*H3H5)/(H3H3*H5H5 - H3H5*H3H5);
//         if( 0 <= a and a <= 1 and 0 <= b and b <= 1 ) return true;
//     } 

//     if( H5H5*H1H1 - H5H1*H5H1 > epsilon and H5H5 > epsilon ){
//         b = (H1p*H5H5 - H5p*H5H1)/(H5H5*H5H5 - H5H1*H5H1);
//         a = (H5p*H1H1 - H1p*H5H1)/(H5H5*H5H5 - H5H1*H5H1);
//         if( 0 <= a and a <= 1 and 0 <= b and b <= 1 ) return true;
//     }

//     return false;
// }



// Periodic Boundary conditions for the index of H, keeps it between 1 and 6
int PBCH( int i ) {
    while( i > 6 ) i -= 6;
    while( i < 1 ) i += 6;
    return i;
}


// Checks rays generation 
void RayGenCheck( string outfile, vector<long double> box, vector<long double> rel_vel ){
    ofstream Pout("outfile");
    for( int i = 0; i < 10000; i+=50 ){
        ProjSurface temp( box, rel_vel, i+1 );
        Pout << (double)temp.GetNRays()/(i+1) << endl;
    }
    Pout.close();
}


// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<long double>> Simulate( vector<long double> box, Body& body, vector<long double> rain_v, long double vmin, long double vmax, unsigned int N, unsigned int nrays ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<long double> body_v(N);
    vector<long double> wetness(N);
    for( unsigned int i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(long double)i/((long double)N-1) );
        vector<long double> relvel = rain_v;
        relvel[0] -= body_v[i];
        wetness[i] = Norm(relvel)*ProjSurface( box, relvel, nrays ).BodyProj(body)/body_v[i];
    }
    return {body_v, wetness};
}

// Estimates wetness for N velocities of the body between vmin and vmax (measured as fractions of vertical rain speed), and returns a matrix with the velocities as the first colunmn and the respective theorical wetness as the second column and the estimated wetness as the third
vector<vector<long double>> CompareAN( vector<long double> box, Body& body, vector<long double> rain_v, long double vmin, long double vmax, unsigned int N, unsigned int nrays ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vmin*=-rain_v[2];
    vmax*=-rain_v[2];
    vector<long double> body_v(N);
    vector<long double> analytical(N);
    vector<long double> wetness(N);
    for( unsigned int i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(long double)i/((long double)N-1) );
        vector<long double> relvel = rain_v;
        relvel[0] -= body_v[i];
        analytical[i] = body.Anal( relvel, body_v[i]);
        wetness[i] = Norm(relvel)*ProjSurface( box, relvel, nrays ).BodyProj(body)/body_v[i];
    }
    return {body_v, analytical,  wetness};
}
