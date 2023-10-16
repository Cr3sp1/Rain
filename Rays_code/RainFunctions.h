#ifndef __RainFunctions_h__
#define __RainFunctions_h__


#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "VectorOperations.h"
#include "body.h"

using namespace std;


// Auxiliary function used only in the surface constructor that checks wether a point is inside the hexagon
bool PointIsInside( vector<long double> Point, vector<long double> p, vector<vector<long double>> H );


// Projects the Point on a plane perpendicular to v and passing through p
vector<long double> Project( vector<long double> Point, vector<long double> p, vector<long double> v );

#endif