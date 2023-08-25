#ifndef __Rain_h__
#define __Rain_h__

#include <armadillo>

using namespace arma;


// Single raindrop class
class Raindrop {

  private:

	double rad;     // Radius (mm)
	double vol;		// Volume (mm^3)
	bool hit;		// Whether it has come in contact with the body
	vec3 pos;		// Cartesian Coordinates (mm)
	vec3 vel;		// Velocity in cartesian coordinates (mm/s)
	

  public:

    // Default constructor
	Raindrop();
	// Complete constructor ([wind_vel] = [body_vel] = [mm/s] )
	Raindrop( double radius, vec3 position, vec3 wind_vel, double body_vel );
	// Sets the hit status
	void SetHit(bool status);
	// Returns the radius
	double Rad() const;
	// Returns the volume
	double Vol() const;
	// Returns the hit status
	bool Hit() const;
	// Returns the position
	vec3 Pos() const;
	// Returns the velocity
	vec3 Vel() const;
	// Time evolution without boundaries ([dt] = [s])
	void Move( double dt );
	// Time evolution with periodic boundary conditions ([dt] = [s], [max] = [mm])
	void Move( double dt, vec3 max );


};

#endif // __Rain_h__
