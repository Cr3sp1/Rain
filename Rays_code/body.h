#ifndef __Body_h__
#define __Body_h__


#include <fstream>
#include <cmath>
#include <vector>
#include "VectorOperations.h"



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
	// Primes the body to be checked
	virtual void Prime( vector<long double> p, vector<long double> v );
	// Checks if the body is making contact with a ray
	virtual bool Check( Ray& rayy );
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	virtual long double Anal( vector<long double> RelVel, long double bodyvel  );


};




// Sphere class
class Sphere: public Body {

  private:
	
	vector<long double> cent;	// Position of the center of the sphere 
	long double rad;	// Radius of the sphere
	vector<long double> Hcent;		// Projection of the center of the sphere


  public:

	// Complete constructor ( [center] = [mm], [radius] = [mm] )
	Sphere( vector<long double> center, long double radius );
	// Primes the body to be checked
	void Prime( vector<long double> p, vector<long double> v  ) override;
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	bool Check( Ray& ray ) override;
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	long double Anal( vector<long double> v , long double bodyvel ) override;
	// Gets stuff
	vector<long double> GetCent(){return cent;}
	long double GetRad(){return rad;}

};



// Parallelepiped class
class Pippo: public Body {

  private:
	
	vector<long double> p;	// Position of a vertex 
	vector<vector<long double>> side;	// Sides of the parallelepiped 
	vector<vector<long double>> H;		// Hexagonal projection on surface


  public:

	// Complete constructor 
	Pippo( vector<long double> P, vector<vector<long double>> Side );
	// Primes the body to be checked
	void Prime( vector<long double> p, vector<long double> v  ) override;
	// Checks if the body is making contact with a ray
	bool Check( Ray& ray ) override;
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	long double Anal( vector<long double> v, long double bodyvel  ) override;
	// Gets stuff
	vector<long double> GetP(){return p;}
	vector<vector<long double>>  GetSide(){return side;}

};



// Capsule class
class Capsule: public Body {

  private:
	
	vector<long double> l1,l2;	// Position of the two extremes of the axis
	long double rad;	// Radius of the sphere
	vector<long double> H1, H2;		// Projections of the two extremes of the axise


  public:

	// Complete constructor 
	Capsule( vector<long double> l1, vector<long double> l2, long double radius );
	// Primes the body to be checked
	void Prime( vector<long double> p, vector<long double> v  ) override;
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	bool Check( Ray& ray ) override;
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	long double Anal( vector<long double> v , long double bodyvel ) override;
	// Gets stuff
	vector<long double> GetL1(){return l1;}
	vector<long double> GetL2(){return l2;}
	long double GetRad(){return rad;}

};


#endif // __Body_h__
