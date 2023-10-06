#ifndef __Random__
#define __Random__

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>

using namespace std;

// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // Default constructor
  Random();
  // Constructor that sets the seed from file
  Random(const char* PrimesFile, const char* SeedFile);
  // Destructor
  ~Random();
  // Method to set the seed for the RNG
  void SetRandom(int * , int, int);
  // Method to set the seed for the RNG taking them from files
  void SetRandFromFile( const char* primes, const char* seed );
  // Method to save the seed to a file
  void SaveSeed();
  // Method to generate a random number in the range [0,1)
  double Rannyu(void);
  // Method to generate a random number in the range [min,max)
  double Rannyu(double min, double max);
  // Method to generate a random number with a Gaussian distribution
  double Gauss(double mean, double sigma);
  // Method to generate a random number with an exponential distribution
  double Exponential(double lambda);
  // Method to generate a random number with a Cauchy-Lorentz distribution
  double CauchyLorentz(double mean, double gamma);
  // Method to generate a discrete random number in the range [min,max]
  int Discrete( int min, int max );
  // Method to generate an equiprobable bool
  bool Bool();
  // Method to generate a point distributed uniformly on the surface of an n-dimensional sphere.
  vector<double> Sphere( int dimentions, double radius );
};

#endif // __Random__

