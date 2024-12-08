#ifndef __Ray_h__
#define __Ray_h__


#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "VectorOperations.h"
#include "VectorStat.h"
#include "RainFunctions.h"

using namespace std;


// Forward declaratio
class Body;

// Ray class
class Ray {

  protected:
    // Starting point
	  vector<double> R0;         

    // The ray starts active and deactivates when it hits an object
    bool Active;

    // Weight of the ray, desctibes its proximity to the object
    double Weight;

	
  public:
    // Static direction vector 
    static vector<double> V; 

    // Direction setter
    void SetV( vector<double> v ) { v = V; }

    // Complete static constructor
	  Ray( vector<double> position, vector<double> direction): R0(position), Active(true), Weight(0.) { V = direction; };

    // Efficient constructor
    Ray( vector<double> position ): R0(position), Active(true), Weight(0.) {};

    // Turns the ray on and off
    void On() { Active = true; }
    void Off() { Active = false; }

    // Sets the weight to a value
    void SetWeight( double w ) { Weight = w; }

    // Gets stuff
    vector<double> GetR0() { return R0; }    
    vector<double> GetV() { return V; };     
    bool IsOn() { return Active; }
    double GetWeight() { return Weight; }

};


// Class of the projection of the box on the surface perpendicular to the relative velocity of the rain v
class ProjSurface{

  protected:
    
    // Vertices of the hexagon
    vector<vector<double>> H;
    // Precision (distance between ray origins)
    double dx;
    // Rays
    vector<Ray> rays;

  public:

    // Complete constructor
    ProjSurface(vector<double> box, vector<double> vel, double Dx);

    // Resets the surface between steps (activates each ray)
    void reset();

    // Get stuff
    int GetNRays() { return rays.size(); }
    vector<vector<double>> GetH() { return H; }
    vector<double> GetV() { return rays[0].GetV(); }

    // Prints all the origins of the active rays to file, last column is 1 if active, 0 if not
    void PrintR( ofstream &fout );
    void PrintR( string outfile );

    // Prints all the origins of the rays projected on the x-y plane to file, last column is 1 if active, 0 if not
    void PrintRaysFlat( ofstream &fout );
    void PrintRaysFlat( string outfile );

    // Prints all the origins of the rays projected on the x-y plane and their weight to file
    void PrintRaysFlatSmooth( ofstream &fout );
    void PrintRaysFlatSmooth( string outfile );

    // Prints H to file
    void PrintH( ofstream &fout );
    void PrintH( string outfile );

    // Returns an estimate of the projection of the body on the plane
    double BodyProj( Body& body );

    // Returns an estimate of the projection of the dynamic body on the plane
    double BodyProj( Body& body, double tmin, double tmax, unsigned int nstep );

    // Returns a smooth estimate of the projection of the body on the plane
    double BodyProjSmooth( Body& body );

    // Returns a smooth estimate of the projection of the dynamic body on the plane
    double BodyProjSmooth( Body& body, double tmin, double tmax, unsigned int nstep );

};




#endif // __Ray_h__
