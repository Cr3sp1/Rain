#ifndef __rainfun_h__
#define __rainfun_h__


#include <armadillo>
#include "rain.h"
#include "body.h"

using namespace arma;
using namespace std;


// Rain builder
vector<Raindrop> Rain( int nrain, vec3 box, vec wind, double bodyvel, Random& rand );

// Enacts the simulation



#endif // __rainfun_h__