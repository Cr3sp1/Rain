#include "RainFunctions.h"

using namespace std;


// Auxiliary function used only in the constructor that checks wether a point is inside the hexagon
bool PointIsInside( vector<long double> Point, vector<long double> p, vector<vector<long double>> H ) {
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



// Projects the Point on a plane perpendicular to v and passing through p
vector<long double> Project( vector<long double> Point, vector<long double> p, vector<long double> v ){
    vector<long double> diff = p - Point;
    return ( Point + v*(diff*v)/(v*v) );
}