#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <cmath>
#include "random.h"
#include "rainfun.h"

using namespace arma;
using namespace std;


// Rain builder
vector<Raindrop> Rain( int nrain, vec3 box, vec wind, double bodyvel, Random& rand ) {
    vector<Raindrop> rain(nrain);
    for( int i = 0; i < (int) rain.size(); i++ ) {
        rain[i] = Raindrop(rand.Rannyu(0, 3), {rand.Rannyu(0, box(0)),rand.Rannyu(0, box(1)),rand.Rannyu(0, box(1))}, wind, bodyvel);
    }
    return rain;
}

// Enacts the simulation