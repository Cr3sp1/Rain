#include "ray.h"
#include "body.h"

using namespace std;


// Declares V
vector<double> Ray::V(3);

// Complete static ray constructor
Ray::Ray( vector<double> position, vector<double> direction){
    R0 = position;
    V = direction;
    Active = true;
}

// Efficient constructor
Ray::Ray( vector<double> position ) {
    R0 = position;
    Active = true;
}

// Resets the surface between steps (turns on all rays)
void ProjSurface::reset() {
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        rays[i].On();
    }
}


// Complete static constructor
ProjSurface::ProjSurface(vector<double> box, vector<double> vel, double dx){

    // Checks validity of arguments
    if( box.size() != 3 or box[0] <= 0 or box[1] <= 0 or box[2] <= 0  ) {
        throw invalid_argument("Box must be of size 3 with all components positive!");
    }
    if( vel.size() != 3 or vel[2] >= 0  ) {
        throw invalid_argument("Velocity must be of size 3 with third component negative!");
    }

    
    // Finds vertex between three "seen" faces and projects the other vertices on the plane passing by it perpendicular to vel
    vector<vector<double>> sides = { {box[0],0,0}, {0,box[1],0}, {0,0,box[2]} };
    H = FindHexProj( {0,0,0}, sides, vel, FindMiddle({0,0,0}, sides, vel) );


    // Evaluates the surface
    surf = 0;
    for( int i = 1; i < 7; i+=2 ) {
        surf += Norm( CrossProduct( H[0]-H[i], H[i]-H[i+1] ) );
    }
    // cout << "Tot Surface = " << surf << endl;


    // Generates the rays on a square grid along directions u1 and u2 (u1 and ud perpendicular and belonging to surface)
    // cout << "Generating rays" << endl;
    rays = {};
    // cout << "dx = " << dx << endl;
    vector<double> u1 = (Norm(H[5] - H[0]) > Norm(H[3] - H[0])) ? H[5] - H[0] : H[3] - H[0] ;      // makes sure that u1 isn't infinitesimal
    u1 = (u1/Norm(u1))*dx;
     double maxu1 = MaxU( H, u1 );
    // cout << "|u1| = " << Norm(u1) << endl;
    // cout << "steps1 = " << Norm(u1) << endl;
    vector<double> u2 = CrossProduct(u1, vel);
    u2 = (u2/Norm(u2))*dx;
    double maxu2 = MaxU( H, u2 );
    // cout << "|u2| = " << Norm(u2) << endl;

    for( double sign1 = -1; sign1 < 2; sign1+=2 ) {
        vector<double> point1 = (sign1 == -1) ? H[0] : H[0] + u1;
        
        while( Norm( point1 - H[0]) < maxu1 ) {
            for( double sign2 = -1; sign2 < 2; sign2+=2 ){
                vector<double> point2 = (sign2 == -1) ? point1 : point1 + u2;

                while( Norm(point2 - point1) < maxu2 ){
                    if( PointIsInsideT( point2, H )) {
                        Ray temp( point2 );
                        rays.push_back(temp);
                    }
                    point2 += (sign2*u2);
                }
            }
            point1 += (sign1*u1); 
        }
    }


    // Sets the rain speed
    Ray::V = vel;

    // cout << "Number of rays generated: " <<  rays.size()  << endl;
}


// Prints all the origins of the rays to file
void ProjSurface::PrintR( ofstream &fout ){
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        vector<double> r = rays[i].GetR0();
        fout << r[0] << ", " << r[1] << ", " << r[2] << endl;
    }
}

void ProjSurface::PrintR( string outfile ){
    ofstream fout(outfile);
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        vector<double> r = rays[i].GetR0();
        fout << r[0] << ", " << r[1] << ", " << r[2] << endl;
    }
    fout.close();
}


// Prints H to file
void ProjSurface::PrintH( ofstream &fout ){
    for( long unsigned int i = 0; i < H.size(); i++ ){
        fout << H[i][0] << ", " << H[i][1] << ", " << H[i][2] << endl;
    }
}

void ProjSurface::PrintH( string outfile ){
    ofstream fout(outfile);
    for( long unsigned int i = 0; i < H.size(); i++ ){
        fout << H[i][0] << ", " << H[i][1] << ", " << H[i][2] << endl;
    }
    fout.close();
}


// Returns an estimate of the projection of the body on the plane
double ProjSurface::BodyProj( Body& body ) {
    unsigned int nhit = 0;
    // cout << "Projecting on " << rays.size() << " rays" << endl;
    body.Prime( H[0], Ray::V );
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        if( body.Check( rays[i]) ) nhit++;
    }
    // cout << "nhit = " << nhit << endl;
    return surf*nhit/rays.size();
}