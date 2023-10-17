#ifndef __RainFunctions_h__
#define __RainFunctions_h__


#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "VectorOperations.h"


using namespace std;

// Forward declaration
class Body;
class Ray;
class ProjSurface;


// Projects the Point on a plane perpendicular to v and passing through p
vector<long double> Project( vector<long double> Point, vector<long double> p, vector<long double> v );


// (NOTE: DOESN'T ALWAYS WORK)  Auxiliary function used only in the surface constructor that checks wether a point is inside the hexagon checking parallelograms
bool PointIsInsideP( vector<long double> Point, vector<long double> p, vector<vector<long double>> H );

// Auxiliary function used only in the surface constructor that checks wether a point is inside the hexagon checking triangles
bool PointIsInsideT( vector<long double> Point, vector<long double> p, vector<vector<long double>> H );

// Periodic Boundary conditions for the index of H
int PBCH( int i );

// Checks rays generation 
void RayGenCheck( string outfile, vector<long double> box, vector<long double> rel_vel );

#endif