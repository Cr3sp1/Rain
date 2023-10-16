#include "ray.h"

using namespace std;


// Complete ray constructor
Ray::Ray( vector<long double> position, vector<long double> direction){
    R0 = position;
    V = direction;
    Active = true;
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
    for( int i = 0; i < 6; i++) cout << "H["<<i<<"] = (" << H[i][0] << ", " << H[i][1] << ", " << H[i][2] << ")" << endl;
    for( int i = 0; i < 6; i++) cout << (H[i]-p)*vel << endl;


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
    vector<long double> u1 = (Norm(H[0] - p) > Norm(H[2] - p)) ? H[0] - p : H[2] - p ;
    u1 = (u1/Norm(u1))*d;
    // cout << "|u1| = " << Norm(u1) << endl;
    vector<long double> u2 = CrossProduct(u1, vel);
    u2 = (u2/Norm(u2))*d;
    // cout << "|u2| = " << Norm(u2) << endl;
    vector<long double> point1 = p;
    bool end_reached1, end_reached2;

    do{
        vector<long double> point2 = point1;

        do{
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2+=u2;
            end_reached2 = PointIsInside(point2, p, H );
        } 
        while(end_reached2);

        point2 = point1;
        do{
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2-=u2;
            end_reached2 = PointIsInside(point2, p, H );
        } 
        while(end_reached2);

        point1 += u1;
        end_reached1 = PointIsInside(point1, p, H );
    }
    while(end_reached1);

    point1 = p;

    do{
        vector<long double> point2 = point1;
        do{
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2+=u2;
            end_reached2 = PointIsInside(point2, p, H );
        } 
        while(end_reached2);

        point2 = point1;
        do{
            Ray temp( point2, vel );
            rays.push_back(temp);
            point2-=u2;
            end_reached2 = PointIsInside(point2, p, H );
        } 
        while(end_reached2);

        point1 -= u1;
        end_reached1 = PointIsInside(point1, p, H );
    }
    while(end_reached1);

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