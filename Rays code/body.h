#ifndef __Body_h__
#define __Body_h__


#include <fstream>
#include <cmath>
#include "body.h"
#include <vector>

using namespace std;

// Body class
class Body {

  protected:

	double wetness;     // Volume of water absorbed (mm^3)

	
  public:
    // Default constructor
	Body();
	// Returns the wetness
	double Wetness() const;
	// Time evolution ( [dt] = [s] )
	virtual void Move( double dt ) {}
	// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
	virtual void Check() {}
	// Checks an array of raindrops	WIP
	// void Check( vector<Raindrop> rain ):

};


// Sphere class
class Sphere: public Body {

  private:
	
	vector<double> cent;	// Position of the center of the sphere (mm)
	double rad;	// Radius of the sphere (mm)


  public:

	// Complete constructor ( [center] = [mm], [radius] = [mm] )
	Sphere( vector<double> center, double radius );
	// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
	void Check( ) override;

};


// Parallelepiped class
class Pippo: public Body {

  private:
	
	vector<double> cent;	// Position of the center of the sphere (mm)
	vector<double> dim;  // Dimensions along the axes (mm)


  public:

	// Complete constructor ( [center]=[mm], [dimensions]=[mm] )
	Pippo( vector<double> center, vector<double> dim );
	// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
	void Check() override;

};



#endif // __Body_h__
