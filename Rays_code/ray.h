#ifndef __Ray_h__
#define __Ray_h__


#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "VectorOperations.h"
#include "RainFunctions.h"

using namespace std;


// Forward declaratio
class Body;

// Ray class
class Ray {

  protected:
    // Starting point
	vector<long double> R0;    

    // Direction vector 
    vector<long double> V;      

    // The ray starts active and deactivates when it hits an object
    bool Active;

	
  public:
    // Complete constructor
	  Ray( vector<long double> position, vector<long double> direction);

    // Turns the ray on and off
    void On(){Active = true;}
    void Off(){Active = false;}

    // Gets stuff
    vector<long double> GetR0(){return R0;}    
    vector<long double> GetV(){return V;};     
    bool IsOn(){return Active;}

};


// Class of the projection of the box on the surface perpendicular to the relative velocity of the rain v
class ProjSurface{

  protected:
    
    // Vertices of the hexagon
    vector<vector<long double>> H;
    // Total surface
    long double surf;
    // Rays
    vector<Ray> rays;

  public:

    // Complete constructor
    ProjSurface(vector<long double> box, vector<long double> vel, unsigned int nrays);

    // Resets the surface between steps
    void reset();

    // Get stuff
    int GetNRays(){return rays.size();}

    // Prints all the origins of the rays to file
    void PrintR( string outfile );

    // Prints H to file
    void PrintH( string outfile );

    // Returns an estimate of the projection of the body on the plane
    long double BodyProj( Body& body );

};




#endif // __Ray_h__
