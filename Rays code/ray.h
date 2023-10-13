#ifndef __Ray_h__
#define __Ray_h__


#include <fstream>
#include <cmath>
#include <vector>
#include "VectorOperations.h"
#include "body.h"

using namespace std;


// Ray class
class Ray {

  protected:
    // Starting point
	vector<double> R0;    

    // Direction vector 
    vector<double> V;      

    // The ray starts active and deactivates when it hits an object
    bool Active;

	
  public:
    // Complete constructor
	Ray( vector<double> position, vector<double> direction);

    // Turns the ray on and off
    void On(){Active = true;}
    void Off(){Active = false;}

    // Gets stuff
    vector<double> GetR0(){return R0;}    
    vector<double> GetV(){return V;};     
    bool IsOn(){return Active;}

};


// Class of the projection of the box on the surface perpendicular to the relative velocity of the rain v
class ProjSurface{

  protected:
    
    // Vertices of the hexagon
    vector<double> H;
    // Total surface
    double surf;
    // Number of rays
    int nrays;
    // Number of rays intersecting the object
    int nhit;
    // Rays
    vector<Ray> rays;

  public:

    // Complete constructor
    ProjSurface(vector<double> box, vector<double> vel, int nrays);

    // Resets the surface between steps
    void reset();

    // Returns an estimate of the projection of the body on the plane
    double proj();

};




#endif // __Ray_h__
