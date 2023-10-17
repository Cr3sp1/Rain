#include "ray.h"
#include "body.h"

using namespace std;


// Complete ray constructor
Ray::Ray( vector<long double> position, vector<long double> direction){
    R0 = position;
    V = direction;
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
        throw std::invalid_argument("Box must be of size 3 with all components positive!");
    }
    if( vel.size() != 3 or vel[2] >= 0  ) {
        throw std::invalid_argument("Velocity must be of size 3 with third component negative!");
    }


    // Finds vertex between three "seen" faces
    vector<long double> p(3);
    for( int i = 0; i < 3; i++ ){
        p[i] = vel[i] < 0 ? box[i] : 0;
    }
    // cout << "p = (" << p[0] << ", " << p[1] << ", " << p[2] << ")" << endl;


    // Finds the projections H of the six vertices adjacent to P
    vector<vector<long double>> delta(3, vector<long double>(3, 0.0));        // Used to calculate the position of the vertices
    for( int i = 0; i < 3; i++ ){
        delta[i][i] = vel[i] < 0 ? (-box[i]) : box[i];
    }
    vector<vector<long double>> h(6, vector<long double>(3));       // h[i] is the position of the vertex i
    h[0] = p + delta[0];
    h[1] = p + delta[0] + delta[1];
    h[2] = p + delta[1];
    h[3] = p + delta[1] + delta[2];
    h[4] = p + delta[2];
    h[5] = p + delta[2] + delta[0];
    for( int i = 0; i < 6; i++ ){
        h[i] = Project( h[i], p, vel );         // We project them
    }
    H = h;
    H.push_back(p);
    // for( int i = 0; i < 6; i++) cout << "H["<<i<<"] = (" << H[i][0] << ", " << H[i][1] << ", " << H[i][2] << ")" << endl;
    // for( int i = 0; i < 6; i++) cout << (H[i]-p)*vel << endl;


    // Evaluates the surface
    surf = 0;
    for( int i = 0; i < 6; i+=2 ) {
        surf += Norm( CrossProduct( p-H[i], H[i]-H[i+1] ) );
    }
    // cout << "Tot Surface = " << surf << endl;


    // Generates the rays on a square grid along directions u1 and u2
    cout << "Generating rays" << endl;
    rays = {};
    long double d = sqrt(surf/n_rays);
    // cout << "d = " << d << endl;
    vector<long double> u1 = (Norm(H[4] - p) > Norm(H[2] - p)) ? H[4] - p : H[2] - p ;      // makes sure that u1 isn't infinitesimal
    u1 = (u1/Norm(u1))*d;
    // cout << "|u1| = " << Norm(u1) << endl;
    vector<long double> u2 = CrossProduct(u1, vel);
    u2 = (u2/Norm(u2))*d;
    // cout << "|u2| = " << Norm(u2) << endl;
    vector<long double> point1 = p;
    bool still_inside1, still_inside2;

    do{
        vector<long double> point2 = point1;

        do{
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2 = Project( point2 + u2, p, vel);
            still_inside2 = PointIsInsideT(point2, p, H );
        } 
        while(still_inside2);

        point2 = Project(point1 - u2, p, vel);
        still_inside2 = PointIsInsideT(point2, p, H );
        while(still_inside2){
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2 = Project( point2 - u2, p, vel);
            still_inside2 = PointIsInsideT(point2, p, H );
        } 
        

        point1 = Project( point1 + u1, p, vel);
        still_inside1 = PointIsInsideT(point1, p, H );
    }
    while(still_inside1);

    point1 = Project( p - u1, p, vel);
    still_inside1 = PointIsInsideT(point1, p, H );
    while(still_inside1){
        vector<long double> point2 = point1;
        do{
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2 = Project( point2 + u2, p, vel);
            still_inside2 = PointIsInsideT(point2, p, H );
        } 
        while(still_inside2);

        point2 = Project(point1 - u2, p, vel);
        still_inside2 = PointIsInsideT(point2, p, H );;
        while(still_inside2){
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2 = Project( point2 - u2, p, vel);
            still_inside2 = PointIsInsideT(point2, p, H );
        } 

        point1 = Project( point1 - u1, p, vel);
        still_inside1 = PointIsInsideT(point1, p, H );
    }
    

    cout << "Number of rays generated: " <<  rays.size() << "/" << n_rays << endl;
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
    cout << "Projecting on " << rays.size() << " rays" << endl;
    for( long unsigned int i = 0; i < rays.size(); i++ ){
        if( body.Check( rays[i]) ) {
            nhit++;
            rays[i].Off();
        }
    }
    cout << "nhit = " << nhit << endl;
    return surf*nhit/rays.size();
}