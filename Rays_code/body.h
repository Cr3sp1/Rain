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
	
  protected:
	// Name (optional)
	string name;
  	// Time (should go from 0 to 1)
	double t;
	// Points used for rotation in time evolution: rot[0] is center, rot[1]-rot[0] is the axis of rotation
	vector<vector<double>> rot;
	// Amplitude of oscillation in radiants (should go from 0 to pi)
	double w;
	// Vectors used for periodic translation in time evolution (coefficient in sin expasion of periodic motion)
	vector<vector<double>> trans;
	// Pointers to all the bodies that move relative to this one
	vector<Body*> SubBodies;
	// Pointer to the body this moves relative to 
    Body* SuperBody; 



  public:
    // Default constructor
	Body(): t(0), w(0), SuperBody(nullptr) {};
	// Dynamic constructor
	Body( string Name, vector<vector<double>> Rot, double W, vector<vector<double>> Trans ): name(Name), t(0), rot(Rot), w(W), trans(Trans), SuperBody(nullptr) {} 
	// Copy constructor
    Body(const Body& other);
	// Copy assignment operator
	Body& operator=(const Body& other);
	// Destructor
    virtual ~Body();
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	virtual void Prime( vector<double> p, vector<double> v ) {}
	// Checks if the body is making contact with a ray
	virtual bool Check( Ray& rayy ) { return false; }
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	virtual double Anal( vector<double> RelVel, double bodyvel  ) { return -1; }
	// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
	virtual void Move( double T );
	// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodies
	virtual void BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat );
	// Sets SuperBody
	virtual void SetSuperBody(Body* body) { SuperBody = body; }
	// Adds a Body to SubBodies
	virtual void AddSubBody(Body* body);


};




// Sphere class
class Sphere: public Body {

  protected:
	
	vector<double> cent;	// Position of the center of the sphere 
	double rad;	// Radius of the sphere
	vector<double> Hcent;		// Projection of the center of the sphere


  public:

	// Complete static constructor
	Sphere( vector<double> center, double radius ): Body(), cent(center), rad(radius) {}
	// Complete dynamic constructor
	Sphere( vector<double> center, double radius, string Name, vector<vector<double>> Rot, double W, vector<vector<double>> Trans  ): Body( Name, Rot, W, Trans ), cent(center), rad(radius) {}
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	void Prime( vector<double> p, vector<double> v  ) override;
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	bool Check( Ray& ray ) override;
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	double Anal( vector<double> v , double bodyvel ) override;
	// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
	void Move( double T ) override;
	// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodies
	void BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) override;
	// Gets stuff
	vector<double> GetCent(){return cent;}
	double GetRad(){return rad;}

};




// Parallelepiped class
class Pippo: public Body {

  protected:
	
	vector<double> cent;	// Position of the center
	vector<vector<double>> side;	// Sides of the parallelepiped 
	vector<vector<double>> H;		// Hexagonal projection on surface


  public:

	// Complete static constructor 
	Pippo( vector<double> Center, vector<vector<double>> Side ): Body(), cent(Center), side(Side) {}
	// Complete dynamic constructor
	Pippo( vector<double> Center, vector<vector<double>> Side, string Name, vector<vector<double>> Rot, double W, vector<vector<double>> Trans  ): Body( Name, Rot, W, Trans ), cent(Center), side(Side) {}
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	void Prime( vector<double> p, vector<double> v  ) override;
	// Checks if the body is making contact with a ray
	bool Check( Ray& ray ) override;
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	double Anal( vector<double> v, double bodyvel  ) override;
	// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
	void Move( double T ) override;
	// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodies
	void BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) override;
	// Gets stuff
	vector<double> GetCent() {return cent;}
	vector<vector<double>>  GetSide() {return side;}
	vector<vector<double>>  GetVertices();

};




// Capsule class
class Capsule: public Body {

  protected:
	
	vector<double> l1,l2;	// Position of the two extremes of the axis
	double rad;	// Radius of the sphere
	vector<double> H1, H2;		// Projections of the two extremes of the axise


  public:

	// Complete static constructor 
	Capsule( vector<double> L1, vector<double> L2, double Radius ): Body(), l1(L1), l2(L2), rad(Radius) {} 
	// Complete dynamic constructor
	Capsule( vector<double> L1, vector<double> L2, double Radius, string Name, vector<vector<double>> Rot, double W, vector<vector<double>> Trans  ): Body( Name, Rot, W, Trans ), l1(L1), l2(L2), rad(Radius) {}
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	void Prime( vector<double> p, vector<double> v  ) override;
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	bool Check( Ray& ray ) override;
	// Analytical solution of rain intercepted. v is relative velocity, bodyvel is body velocity
	double Anal( vector<double> v , double bodyvel ) override;
	// Time evolution of the body in its own frame of reference, also propagates to the sub-bodies
	void Move( double T ) override;
	// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodies
	void BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) override;
	// Gets stuff
	vector<double> GetL1(){return l1;}
	vector<double> GetL2(){return l2;}
	double GetRad(){return rad;}

};




// ManyBody class
class ManyBody: public Body {

  public:

	// Bodies contained in ManyBody
	vector<Sphere> spheres;
	vector<Pippo> pippos;
	vector<Capsule> capsules;
	// Complete constructor 
	ManyBody( const vector<Sphere>& Spheres, const vector<Pippo>& Pippos, const vector<Capsule>& Capsules );
	// Primes the body to be checked. p is a point on the surface containing the ray origins and v is the relative velocity
	void Prime( vector<double> p, vector<double> v  ) override;
	// Checks if the body is making contact with a ray and if so adds its the volume to the wetness
	bool Check( Ray& ray ) override;
	// Time evolution of all the bodies
	void Move( double T ) override;
	// Time evolution caused by the super-body, affects the whole frame of reference, also propagates to the sub-bodies
	void BeMoved( vector<double> Delta, vector<double> Rot0, vector<vector<double>> Rotmat ) override {};
	// Add bodies
	void AddSphere( const Sphere& sphere ) { spheres.push_back(sphere); }
	void AddPippo( const Pippo& pippo ) { pippos.push_back(pippo); }
	void AddSphere( const Capsule& capsule ) { capsules.push_back(capsule); }
	// Attaches the sub-body to the super-body
	void Attach( Body& SubBody, Body& SuperBody );
	// Prints to file the state (all the bodies and their parameters)
	void PrintState( ofstream &fout );
	void PrintState( string outfile );

};



#endif // __Body_h__
