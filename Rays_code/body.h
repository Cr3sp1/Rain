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
	
  private:
  	// Time (should go from 0 to 1)
	long double t;
	// Points used for rotation in time evolution: rot[0] is center, rot[1]-rot[0] is the axis of rotation
	vector<vector<long double>> rot;
	// Amplitude of oscillation in radiants (should go from 0 to pi)
	long double w;
	// Vectors used for periodic translation in time evolution (coefficient in sin expasion of periodic motion)
	vector<vector<long double>> trans;
	// Pointers to all the bodies that move relative to this one
	vector<Body*> SubBodies;
	// Pointer to the body this moves relative to 
    Body* SuperBody; 



  public:
    // Default constructor
	Body(): t(0), SuperBody(nullptr) {};
	// Copy constructor
    Body(const Body& other);
	// Copy assignment operator
	Body& operator=(const Body& other);
	// Destructor
    virtual ~Body();
	// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
	virtual void Move( long double T );
	// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodies
	virtual void BeMoved( long double T, vector<vector<long double>> Rot, long double W, vector<vector<long double>> Trans );
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	virtual void Prime( vector<long double> p, vector<long double> v ) {}
	// Checks if the body is making contact with a ray
	virtual bool Check( Ray& rayy ) { return false; }
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	virtual long double Anal( vector<long double> RelVel, long double bodyvel  ) { return -1; }
	// Adds a Body to SubBodies
	virtual void AddSubBody(Body* body) { SubBodies.push_back(body); }
	// Sets SuperBody
	virtual void SetSuperBody(Body* body) { SuperBody = body; }


};




// Sphere class
class Sphere: public Body {

  private:
	
	vector<long double> cent;	// Position of the center of the sphere 
	long double rad;	// Radius of the sphere
	vector<long double> Hcent;		// Projection of the center of the sphere


  public:

	// Complete static constructor
	Sphere( vector<long double> center, long double radius ): Body(), cent(center), rad(radius) {}
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
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
	
	vector<long double> cent;	// Position of the center
	vector<vector<long double>> side;	// Sides of the parallelepiped 
	vector<vector<long double>> H;		// Hexagonal projection on surface


  public:

	// Complete static constructor 
	Pippo( vector<long double> Center, vector<vector<long double>> Side ): Body(), cent(Center), side(Side) {}
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	void Prime( vector<long double> p, vector<long double> v  ) override;
	// Checks if the body is making contact with a ray
	bool Check( Ray& ray ) override;
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	long double Anal( vector<long double> v, long double bodyvel  ) override;
	// Gets stuff
	vector<long double> GetCent(){return cent;}
	vector<vector<long double>>  GetSide(){return side;}

};




// Capsule class
class Capsule: public Body {

  private:
	
	vector<long double> l1,l2;	// Position of the two extremes of the axis
	long double rad;	// Radius of the sphere
	vector<long double> H1, H2;		// Projections of the two extremes of the axise


  public:

	// Complete static constructor 
	Capsule( vector<long double> L1, vector<long double> L2, long double Radius ): Body(), l1(L1), l2(L2), rad(Radius) {} 
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
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




// ManyBody class
class ManyBody: public Body {

  private:
	
	vector<Sphere> spheres;			// Position of the center of the sphere 
	vector<Pippo> pippos;			// Position of the center of the sphere 
	vector<Capsule> capsules;		// Position of the center of the sphere


  public:

	// Complete static constructor 
	ManyBody( vector<Sphere> Spheres, vector<Pippo> Pippos, vector<Capsule> Capsules );
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	void Prime( vector<long double> p, vector<long double> v  ) override;
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	bool Check( Ray& ray ) override;
	// Gets stuff
	vector<long double> GetSphCent( unsigned int index );
	long double GetSphRad( unsigned int index );
	vector<long double> GetPipCent( unsigned int index );
	vector<vector<long double>>  GetPipSide( unsigned int index );
	vector<long double> GetCapL1( unsigned int index );
	vector<long double> GetCapL2( unsigned int index );
	long double GetCapRad( unsigned int index );

};



#endif // __Body_h__
