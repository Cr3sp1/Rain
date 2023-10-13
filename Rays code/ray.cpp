#include "ray.h"

using namespace std;


// Complete constructor
Ray::Ray( vector<double> position, vector<double> direction){
    R0 = position;
    V = direction;
    Active = true;
}