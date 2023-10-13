#ifndef __Body_h__
#define __Body_h__


#include <fstream>
#include <cmath>
#include <vector>
#include "ray.h"


using namespace std;

// Body class
class Body {

  public:
    // Default constructor
	Body();
	// Time evolution ( [dt] = [s] )
	virtual void Move( double dt ) {}
	// Checks if the body is making contact with a ray
	// virtual bool Check( Ray ray ) {}

};


// Sphere class
class Sphere: public Body {

  private:
	
	vector<double> cent;	// Position of the center of the sphere (mm)
	double rad;	// Radius of the sphere (mm)


  public:

	// Complete constructor ( [center] = [mm], [radius] = [mm] )
	Sphere( vector<double> center, double radius );
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	// bool Check( Ray ray ) override;

};


// Parallelepiped class
class Pippo: public Body {

  private:
	
	vector<double> cent;	// Position of the center of the sphere (mm)
	vector<double> dim;  // Dimensions along the axes (mm)


  public:

	// Complete constructor ( [center]=[mm], [dimensions]=[mm] )
	Pippo( vector<double> center, vector<double> dim );
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	// bool Check( Ray ray ) override;

};



#endif // __Body_h__
