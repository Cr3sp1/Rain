#ifndef __Body_h__
#define __Body_h__

#include <armadillo>
#include "rain.h"

using namespace arma;


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
	virtual void Check( Raindrop& drop ) {}
	// Checks an array of raindrops	WIP
	// void Check( vector<Raindrop> rain ):

};


// Sphere class
class Sphere: public Body {

  private:
	
	vec3 cent;	// Position of the center of the sphere (mm)
	double rad;	// Radius of the sphere (mm)


  public:

	// Complete constructor ( [center] = [mm], [radius] = [mm] )
	Sphere( vec3 center, double radius );
	// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
	void Check( Raindrop& drop ) override;

};


// Parallelepiped class
class Pippo: public Body {

  private:
	
	vec3 cent;	// Position of the center of the sphere (mm)
	vec3 dim;  // Dimensions along the axes (mm)


  public:

	// Complete constructor ( [center]=[mm], [dimensions]=[mm] )
	Pippo( vec3 center, vec3 dim );
	// Checks if the body is making contact with a raindrop and if so adds its the volume to the wetness
	void Check( Raindrop& drop ) override;

};



#endif // __Body_h__
