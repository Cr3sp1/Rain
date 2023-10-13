#include "ray.h"

using namespace std;


// Complete ray constructor
Ray::Ray( vector<double> position, vector<double> direction){
    R0 = position;
    V = direction;
    Active = true;
}

// Complete constructor
    ProjSurface::ProjSurface(vector<double> box, vector<double> vel, unsigned int nrays){

        // Checks validity of arguments
        if( box.size() != 3 or box[0] <= 0 or box[1] <= 0 or box[2] <= 0  ) {
            throw std::invalid_argument("Box must be of size 3 with all components positive!");
        }
        if( vel.size() != 3 or vel[2] <= 0  ) {
            throw std::invalid_argument("Box must be of size 3 with third component negative!");
        }

        // Finds vertex between three "seen" faces
        vector<double> p(3);
        for( int i = 0; i < 3; i++ ){
            p[i] = vel[i] < 0 ? box[i] : 0;
        }
        cout << "Seen vertex is " << p[0] << p[1]  << p[2] << endl;


        // Temporary
        H = box;
        surf = 0;

    }