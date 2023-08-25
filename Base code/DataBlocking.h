#ifndef __DataBlocking__
#define __DataBlocking__

#include <armadillo>
using namespace arma;

// Evalues the error on the firs N blocks
double block_error(const vec& BlockAvg, int N);
// Same but for rows
double block_error(const rowvec& BlockAvg, int N);
// Builds the error vector
vec block_error(const vec& BlockAvg);
// Builds the average of first i blocks vector
vec tot_avg(const vec& BlockAvg);
// Builds and prints in the same file the averages and the errors as vectors 
void print_avg_err(const vec& BlockAvg, const char* namefile);



#endif