#include "ray.h"

using namespace std;


// Complete ray constructor
Ray::Ray( vector<double> position, vector<double> direction){
    R0 = position;
    V = direction;
    Active = true;
}


// Auxiliary function used only in the constructor that checks wether a point is inside the hexagon
bool PointIsInside( vector<double> Point, vector<double> p, vector<vector<double>> H ) {
    // centers the points in p
    Point -= p;
    for(int i = 0; i < 6; i++ ) H[i] -=p;
    // Calculate stuff
    double a, b;
    double epsilon = 1e-14;
    double H0p = H[0]*p;
    double H0H0 = H[0]*H[0];
    double H0H2 = H[0]*H[2];
    double H2p = H[2]*p;
    double H2H2 = H[2]*H[2];
    double H2H4 = H[2]*H[4];
    double H4p = H[4]*p;
    double H4H4 = H[4]*H[4];
    double H4H0 = H[4]*H[0];

    if( H0H0*H2H2 - H0H2*H0H0 > epsilon and H0H0 > epsilon ){
        b = (H2p*H0H0 - H0p*H0H2)/(H0H0*H2H2 - H0H2*H0H0);
        a = (H0p - b*H0H2)/H0H0;
        if( 0 < a and a <= 1 and 0 < b and b <= 1 ) return true;
    }

    if( H2H2*H4H4 - H2H4*H2H2 > epsilon and H2H2 > epsilon ){
        b = (H4p*H2H2 - H2p*H2H4)/(H2H2*H4H4 - H2H4*H2H2);
        a = (H2p - b*H2H4)/H2H2;
        if( 0 < a and a <= 1 and 0 < b and b <= 1 ) return true;
    }

    if( H4H4*H4H4 - H4H0*H4H4 > epsilon and H4H4 > epsilon ){
        b = (H0p*H4H4 - H4p*H4H0)/(H4H4*H4H4 - H4H0*H4H4);
        a = (H4p - b*H4H0)/H4H4;
        if( 0 < a and a <= 1 and 0 < b and b <= 1 ) return true;
    }

    return false;
}



// Complete constructor
    ProjSurface::ProjSurface(vector<double> box, vector<double> vel, unsigned int n_rays){

        // Checks validity of arguments
        if( box.size() != 3 or box[0] <= 0 or box[1] <= 0 or box[2] <= 0  ) {
            throw std::invalid_argument("Box must be of size 3 with all components positive!");
        }
        if( vel.size() != 3 or vel[2] >= 0  ) {
            throw std::invalid_argument("Velocity must be of size 3 with third component negative!");
        }


        // Finds vertex between three "seen" faces
        vector<double> p(3);
        for( int i = 0; i < 3; i++ ){
            p[i] = vel[i] < 0 ? box[i] : 0;
        }
        cout << "p = (" << p[0] << ", " << p[1] << ", " << p[2] << ")" << endl;


        // Finds the projections H of the six vertices adjacent to P
        vector<vector<double>> delta(3, vector<double>(3, 0.0));        // Used to calculate the Delta
        for( int i = 0; i < 3; i++ ){
            delta[i][i] = vel[i] < 0 ? (-box[i]) : box[i];
        }
        vector<vector<double>> Delta(6);        // Delta[i] is the displacement from p of the vertex V[i] that we want to project
        Delta[0] = delta[0];
        Delta[1] = delta[0] + delta[1];
        Delta[2] = delta[1];
        Delta[3] = delta[1] + delta[2];
        Delta[4] = delta[2];
        Delta[5] = delta[2] + delta[0];
        // for( int i = 0; i < 6; i++) cout << "D["<<i<<"] = (" << Delta[i][0] << ", " << Delta[i][1] << ", " << Delta[i][2] << ")" << endl;
        vector<vector<double>> h(6, vector<double>(3));
        for( int i = 0; i < 6; i++ ){

            h[i] = p + Delta[i] - vel*(Delta[i]*vel)/(vel*vel);
        }
        H = h;
        for( int i = 0; i < 6; i++) cout << "H["<<i<<"] = (" << H[i][0] << ", " << H[i][1] << ", " << H[i][2] << ")" << endl;
        // for( int i = 0; i < 6; i++) cout << (H[i]-p)*vel << endl;


        // Evaluates the surface
        surf = 0;
        for( int i = 0; i < 6; i+=2 ) {
            surf += Norm( CrossProduct( p-H[i], H[i]-H[i+1] ) );
        }
        cout << "Tot Surface = " << surf << endl;


        // Generates the rays on a square grid along directions u1 and u2
        double d = sqrt(surf/n_rays);
        vector<double> u1 = p-H[0];
        u1 = u1*d/Norm(u1);
        vector<double> u2 = CrossProduct(u1, vel);
        u2 = u2*d/Norm(u2);
        vector<double> point1 = p;
        rays = {};
        bool end_reached1, end_reached2;
        do{
            vector<double> point2 = point1;
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
            vector<double> point2 = point1;
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


        

        cout << "Number of rays: " <<  rays.size() << endl;
    }



