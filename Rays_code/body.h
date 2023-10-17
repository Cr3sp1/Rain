#ifndef __Body_h__
#define __Body_h__


#include <fstream>
#include <cmath>
#include <vector>



using namespace std;

// Forward declaration
class Ray;
class ProjSurface;

// Body class
class Body {

  public:
    // Default constructor
	Body();
	// Time evolution ( [dt] = [s] )
	virtual void Move( long double dt );
	// Checks if the body is making contact with a ray
	virtual bool Check( Ray& rayy );

};


// Sphere class
class Sphere: public Body {

  private:
	
	vector<long double> cent;	// Position of the center of the sphere (mm)
	long double rad;	// Radius of the sphere (mm)


  public:

	// Complete constructor ( [center] = [mm], [radius] = [mm] )
	Sphere( vector<long double> center, long double radius );
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	bool Check( Ray& ray ) override;
	// Gets stuff
	vector<long double> GetCent(){return cent;}
	long double GetRad(){return rad;}

};



#endif // __Body_h__
