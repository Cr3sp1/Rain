#include "RainFunctions.h"
#include "body.h"
#include "ray.h"

using namespace std;


// Projects the Point on a plane perpendicular to v and passing through p
vector<long double> Project( vector<long double> Point, vector<long double> p, vector<long double> v ){
    vector<long double> diff = p - Point;
    return ( Point + v*(diff*v)/(v*v) );
}

// Auxiliary function used only in the constructor that checks wether a point is inside the hexagon
bool PointIsInsideT( vector<long double> Point, vector<long double> p, vector<vector<long double>> H ){
    // Centers on p
    Point -= p;
    for(int i = 0; i < 6; i++ ){
        H[i] -= p;
    }

    // Checks if Point is inside the reiangle with vertices p, H[i], H[i+1]
    for( int i = 0; i < 6; i++ ){
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

// Auxiliary function used only in the constructor that checks wether a point is inside the hexagon
bool PointIsInsideP( vector<long double> Point, vector<long double> p, vector<vector<long double>> H ){
    // centers the points in p
    Point -= p;
    for(int i = 0; i < 6; i++ ) H[i] -=p;
    // Calculate stuff
    long double a, b;
    long double epsilon = 1e-16;
    long double H0p = H[0]*Point;
    long double H0H0 = H[0]*H[0];
    long double H0H2 = H[0]*H[2];
    long double H2p = H[2]*Point;
    long double H2H2 = H[2]*H[2];
    long double H2H4 = H[2]*H[4];
    long double H4p = H[4]*Point;
    long double H4H4 = H[4]*H[4];
    long double H4H0 = H[4]*H[0];

    if( H0H0*H2H2 - H0H2*H0H2 > epsilon and H0H0 > epsilon ){
        b = (H2p*H0H0 - H0p*H0H2)/(H0H0*H2H2 - H0H2*H0H2);
        a = (H0p*H2H2 - H2p*H0H2)/(H0H0*H2H2 - H0H2*H0H2);
        if( 0 <= a and a <= 1 and 0 <= b and b <= 1 ) return true;
    }

    if( H2H2*H4H4 - H2H4*H2H4 > epsilon and H2H2 > epsilon ){
        b = (H4p*H2H2 - H2p*H2H4)/(H2H2*H4H4 - H2H4*H2H4);
        a = (H2p*H4H4 - H4p*H2H4)/(H2H2*H4H4 - H2H4*H2H4);
        if( 0 <= a and a <= 1 and 0 <= b and b <= 1 ) return true;
    } 

    if( H4H4*H0H0 - H4H0*H4H0 > epsilon and H4H4 > epsilon ){
        b = (H0p*H4H4 - H4p*H4H0)/(H4H4*H4H4 - H4H0*H4H0);
        a = (H4p*H0H0 - H0p*H4H0)/(H4H4*H4H4 - H4H0*H4H0);
        if( 0 <= a and a <= 1 and 0 <= b and b <= 1 ) return true;
    }

    return false;
}



// Periodic Boundary conditions for the index of H
int PBCH( int i ) {
    while( i > 5 ) i -= 6;
    while( i < 0 ) i += 6;
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
