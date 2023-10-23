#include "ray.h"
#include "body.h"

using namespace std;


// Declares V
vector<long double> Ray::V(3);

// Complete ray constructor
Ray::Ray( vector<long double> position, vector<long double> direction){
    R0 = position;
    V = direction;
    Active = true;
}

// Efficient constructor
Ray::Ray( vector<long double> position ) {
    R0 = position;
    Active = true;
}

// Resets the surface between steps (turns on all rays)
void ProjSurface::reset() {
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        rays[i].On();
    }
}


// Complete constructor
ProjSurface::ProjSurface(vector<long double> box, vector<long double> vel, unsigned int n_rays){

    // Checks validity of arguments
    if( box.size() != 3 or box[0] <= 0 or box[1] <= 0 or box[2] <= 0  ) {
        throw invalid_argument("Box must be of size 3 with all components positive!");
    }
    if( vel.size() != 3 or vel[2] >= 0  ) {
        throw invalid_argument("Velocity must be of size 3 with third component negative!");
    }

    
    // Finds vertex between three "seen" faces and projects the other vertices on the plane passing by it perpendicular to vel
    vector<vector<long double>> sides = { {box[0],0,0}, {0,box[1],0}, {0,0,box[2]} };
    H = FindHexProj( {0,0,0}, sides, vel, FindMiddle({0,0,0}, sides, vel) );


    // Evaluates the surface
    surf = 0;
    for( int i = 1; i < 7; i+=2 ) {
        surf += Norm( CrossProduct( H[0]-H[i], H[i]-H[i+1] ) );
    }
    // cout << "Tot Surface = " << surf << endl;


    // Generates the rays on a square grid along directions u1 and u2
    // cout << "Generating rays" << endl;
    rays = {};
    long double d = sqrt(surf/n_rays);
    // cout << "d = " << d << endl;
    vector<long double> u1 = (Norm(H[5] - H[0]) > Norm(H[3] - H[0])) ? H[5] - H[0] : H[3] - H[0] ;      // makes sure that u1 isn't infinitesimal
    u1 = (u1/Norm(u1))*d;
    // cout << "|u1| = " << Norm(u1) << endl;
    vector<long double> u2 = CrossProduct(u1, vel);
    u2 = (u2/Norm(u2))*d;
    // cout << "|u2| = " << Norm(u2) << endl;
    vector<long double> point1 = H[0];
    bool still_inside1, still_inside2;

    do{
        vector<long double> point2 = point1;

        do{
            Ray temp( point2);
            rays.push_back(temp);
            point2 = Project( point2 + u2, H[0], vel);
            still_inside2 = PointIsInsideT(point2, H );
        } 
        while(still_inside2);

        point2 = Project(point1 - u2, H[0], vel);
        still_inside2 = PointIsInsideT(point2, H );
        while(still_inside2){
            Ray temp( point2);
            rays.push_back(temp);
            point2 = Project( point2 - u2, H[0], vel);
            still_inside2 = PointIsInsideT(point2, H );
        } 
        

        point1 = Project( point1 + u1, H[0], vel);
        still_inside1 = PointIsInsideT(point1, H );
    }
    while(still_inside1);

    point1 = Project( H[0] - u1, H[0], vel);
    still_inside1 = PointIsInsideT(point1, H );
    while(still_inside1){
        vector<long double> point2 = point1;
        do{
            Ray temp( point2);
            rays.push_back(temp);
            point2 = Project( point2 + u2, H[0], vel);
            still_inside2 = PointIsInsideT(point2,  H );
        } 
        while(still_inside2);

        point2 = Project(point1 - u2, H[0], vel);
        still_inside2 = PointIsInsideT(point2, H );;
        while(still_inside2){
            Ray temp( point2 );
            rays.push_back(temp);
            point2 = Project( point2 - u2, H[0], vel);
            still_inside2 = PointIsInsideT(point2, H );
        } 

        point1 = Project( point1 - u1, H[0], vel);
        still_inside1 = PointIsInsideT(point1, H );
    }
    
    // Sets the rain speed
    Ray::V = vel;

    // cout << "Number of rays generated: " <<  rays.size() << "/" << n_rays << endl;
}


// Prints all the origins of the rays to file
void ProjSurface::PrintR( string outfile ){
    ofstream fout(outfile);
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        vector<long double> r = rays[i].GetR0();
        fout << r[0] << ", " << r[1] << ", " << r[2] << endl;
    }
    fout.close();
}


// Prints H to file
void ProjSurface::PrintH( string outfile ){
    ofstream fout(outfile);
    for( long unsigned int i = 0; i < H.size(); i++ ){
        fout << H[i][0] << ", " << H[i][1] << ", " << H[i][2] << endl;
    }
    fout.close();
}


// Returns an estimate of the projection of the body on the plane
long double ProjSurface::BodyProj( Body& body ) {
    unsigned int nhit = 0;
    // cout << "Projecting on " << rays.size() << " rays" << endl;
    body.Prime( H[0], Ray::V );
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        if( body.Check( rays[i]) ) nhit++;
    }
    // cout << "nhit = " << nhit << endl;
    return surf*nhit/rays.size();
}