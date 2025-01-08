#include "RainFunctions.h"
#include "body.h"
#include "ray.h"

using namespace std;


// Projects the Point on a plane perpendicular to v and passing through p
vector<double> Project( vector<double> Point, vector<double> p, vector<double> v ){
    vector<double> diff = p - Point;
    return ( Point + v*(diff*v)/(v*v) );
}



// Finds the vertex in the middle of the three seen faces of a parallelepiped defined by a point p and three sides 
vector<double> FindMiddle( vector<double> p, vector<vector<double>> sides, vector<double> v ){
    for( size_t i = 0; i < sides.size(); i++ ){
        if( sides[i]*v < 0 ) p += sides[i];
    }
    return p;
}



// Finds the hexagonal projection H of a parallelepiped defined by a vertex p and sides on a plane perpendicular to v and passing through P0
vector<vector<double>> FindHexProj(  vector<double> p, vector<vector<double>> Side, vector<double> v, vector<double> P0){
    vector<vector<double>> H = {FindMiddle( p, Side, v )};
    vector<vector<double>> delta(3, vector<double>(3, 0.0));        // Used to calculate the position of the vertices to project
    for( int i = 0; i < 3; i++ ){
        delta[i] = Side[i]*v < 0 ? ((double)-1)*Side[i] : Side[i];
    }
    H.push_back( H[0] + delta[0] );
    H.push_back( H[0] + delta[0] + delta[1]) ;
    H.push_back( H[0] + delta[1] );
    H.push_back( H[0] + delta[1] + delta[2] );
    H.push_back( H[0] + delta[2] );
    H.push_back( H[0] + delta[2] + delta[0] );
    for( int i = 0; i < 7; i++ ){
        H[i] = Project( H[i], P0, v );         // We project them
    }

    return H;
}



// Returns the highest absolute value of the projections of the vertices of H on a line in direction u1 passing through H[0]
double MaxU(vector<vector<double>> H, vector<double> u ) {
    double result = 0;
    for( int i = 1; i < 7; i++ ) {
        H[i] -= H[0];
        double proj = abs(H[i]*u/Norm(u));
        if(proj > result) result = proj;
    }
    // cout << result << endl;
    return result;
}



// Returns wether the Point is inside the hexagon H using triangles and barycentric coordinates
bool PointIsInsideT( vector<double> Point, vector<vector<double>> H ){
    // Centers all other points on H[1]
    Point -= H[1];
    for(int i = 2; i < 7; i++ ){
        H[i] -= H[1];
    }

    // Checks if Point is inside the triangle with vertices H[1], H[i], H[i+1]
    for( int i = 2; i < 6; i++ ){
        int i_next = i + 1;
        double epsilon = 1e-10;
        double A = Norm( CrossProduct( H[i], H[i_next]) );
        double alpha = Norm( CrossProduct( Point, H[i_next]) )/A;
        double beta = Norm( CrossProduct( Point, H[i]) )/A;
        double gamma = Norm( CrossProduct( Point-H[i], Point-H[i_next]) )/A;
        if( 0 <= alpha and alpha <= 1 and
            0 <= beta and beta <= 1 and
            0 <= gamma and gamma <= 1 and
            abs(alpha + beta + gamma - 1 ) <  epsilon ){
            return true;
        }
    }
    return false;
}


// Periodic Boundary conditions for the index of H, keeps it between 1 and 6
int PBCH( int i ) {
    while( i > 6 ) i -= 6;
    while( i < 1 ) i += 6;
    return i;
}


// Checks rays generation 
void RayGenCheck( string outfile, vector<double> box, vector<double> rel_vel ){
    ofstream Pout("outfile");
    for( int i = 0; i < 10000; i+=50 ){
        ProjSurface temp( box, rel_vel, i+1 );
        Pout << (double)temp.GetNRays()/(i+1) << endl;
    }
    Pout.close();
}



// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
double Wetness( vector<double> box, Body& body, vector<double> rain_v, double vb, double dx ) {
    vector<double> relvel = rain_v;
    relvel[0] -= vb;
    return Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body)/vb;
}



// Estimates wetness for N velocities of the dynamic body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
double Wetness( vector<double> box, Body& body, vector<double> rain_v, double vb, double dx, double tmin, double tmax, unsigned int nstep ) {
    vector<double> relvel = rain_v;
    relvel[0] -= vb;
    return Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body, tmin, tmax, nstep )/vb;
}



// Estimates wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<double> body_v(N);
    vector<double> wetness(N);
    for( size_t i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        wetness[i] = Wetness( box, body, rain_v, body_v[i], dx );
    }
    return Transpose(vector<vector<double>>{ body_v, wetness});
}



// Estimates wetness for N velocities of the dynamic body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> Simulate( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx, double tmin, double tmax, unsigned int nstep ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<double> body_v(N);
    vector<double> wetness(N);
    for( size_t i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        wetness[i] = Wetness( box, body, rain_v, body_v[i], dx, tmin, tmax, nstep );
    }
    return Transpose(vector<vector<double>>{ body_v, wetness});
}



// Estimates wetness for N values of dx between dxmin and dxmax, and returns a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr( vector<double> box, Body& body, vector<double> rain_vel, double body_vel, unsigned int N, double dxmin, double dxmax ) {
    vector<double> d = {dxmax};
    vector<double> S = {Wetness( box, body, rain_vel, body_vel, dxmax)};
    double k =  N == 0 ? 0  :  pow(dxmin/dxmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        d.push_back(d[i-1]*k);
        cout << "dx = " << d[i] << endl;
        S.push_back(Wetness( box, body, rain_vel, body_vel, d[i]));
    }

    return Transpose(vector<vector<double>>{ d, S});
}



// Estimates wetness for N values of dx between dxmin and dxmax for dynamic body, and returns a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr( vector<double> box, Body& body, vector<double> rain_vel, double body_vel, unsigned int N, double dxmin, double dxmax, double tmin, double tmax, unsigned int nstep) {
    vector<double> dx = {dxmax};
    vector<double> S = {Wetness( box, body, rain_vel, body_vel, dxmax, tmin, tmax, nstep )};
    double k =  N == 0 ? 0  :  pow(dxmin/dxmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        dx.push_back(dx[i-1]*k);
        cout << "dx = " << dx[i] << endl;
        S.push_back(Wetness( box, body, rain_vel, body_vel, dx[i], tmin, tmax, nstep ));
    }

    return Transpose(vector<vector<double>>{ dx, S});
}



// Estimates wetness for N values of nstep between nstepmin and nstepmax, and returns a matrix with nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErrT( vector<double> box, Body& body, vector<double> rain_vel, double body_vel, double dx, double tmin, double tmax, unsigned int N, unsigned int nstepmin, unsigned int nstepmax) {
    double dtmin = 1./nstepmax;
    double dtmax = 1./nstepmin;
    vector<double> dt = {dtmax};
    vector<double> nstep = { (double) nstepmin};
    vector<double> S = {Wetness( box, body, rain_vel, body_vel, dx, tmin, tmax, nstepmin )};
    double k =  N == 0 ? 0  :  pow(dtmin/dtmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        dt.push_back(dt[i-1]*k);
        double nstepnew = round(1./dt[i]);
        if( nstepnew > nstep.back()) {
            cout << "Nstep = " << nstepnew << endl;
            nstep.push_back(nstepnew);
            S.push_back(Wetness( box, body, rain_vel, body_vel, dx, tmin, tmax, nstepnew ));
        }
    }

    return Transpose(vector<vector<double>>{ nstep, S});
}



// Estimate wetness for N velocities of the body between vmin and vmax (measured as fractions of vertical rain speed), and returns a matrix with the velocities as the first colunmn and the respective theorical wetness as the second column and the estimated wetness as the third
vector<vector<double>> CompareAN( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    
    vector<double> body_v(N);
    vector<double> analytical(N);
    vector<double> wetness(N);
    for( size_t i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        vector<double> relvel = rain_v;
        relvel[0] -= body_v[i];
        cout << "relvel = (" << relvel[0] << ", " << relvel[1] << ", "<< relvel[2] << ")" << endl;
        analytical[i] = body.Anal( relvel, body_v[i]);
        cout << "anal = " << analytical[i];
        wetness[i] = Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body)/body_v[i];
    }
    
    vector<vector<double>> mat{body_v, analytical, wetness};
    return Transpose(mat);
}



// Returns the minimum distance between the point p and the segment line with extremes l1 and l2
double PointSegDist( vector<double> p, vector<double> l1, vector<double> l2 ) {
    // Changes frame of reference to l1 = 0
    p -= l1;
    l2 -= l1;
    double l2_norm2 = l2*l2;
    // Makes sure not to divide by zero
    if( l2_norm2 == 0 ) return Norm(p);
    // Calculates projection of p on line passing through l1 (0) and l2, value of projection is multiplied by l2_norm
    double proj = p*l2;
    // Returns distance between p and the closest point belonging to the segment
    if( proj <= 0 ) return Norm(p);                    // p closest to l1
    if( proj >= l2_norm2 ) return Norm(p-l2);           // p closest to l2
    return sqrt( p*p - proj*proj/l2_norm2 );               // p closest to its projection on the segment, calculates distance with pythagoras
}



// Returns NxN identity matrix
vector<vector<double>> IdMat( unsigned int N ) {
    vector<vector<double>> idmat(N, vector<double>(N));
    for( size_t i = 0; i < N; i++ ){
        for( size_t j = 0; j < N; j++ ){
            idmat[i][j] = (i==j) ? 1 : 0; 
        }
    }
    return idmat;
}



// Returns the rotation matrix
vector<vector<double>> RotMat( vector<double> axis, double theta ) {
    if ( Norm(axis) == 0 or axis.size() != 3 ) return IdMat(axis.size()); // Handle zero-length vector to avoid division by zero
    axis = axis/Norm(axis);
    double s = sin(theta);
    double c = cos(theta);
    double G = 1-c;

    return { { axis[0]*axis[0]*G + c,           axis[0]*axis[1]*G - axis[2]*s,  axis[0]*axis[2]*G + axis[1]*s },
             { axis[1]*axis[0]*G + axis[2]*s,   axis[1]*axis[1]*G + c,          axis[1]*axis[2]*G - axis[0]*s },
             { axis[2]*axis[0]*G - axis[1]*s,   axis[2]*axis[1]*G + axis[0]*s,  axis[2]*axis[2]*G + c         } };
}



// Rotates a Point relative to the point Rot0
void Rotate( vector<double>& Point, const vector<double>& Rot0, const vector<vector<double>>& Rotmat ){
    if( Rotmat == IdMat(3) ) return;
    Point -= Rot0;
    Point = Rotmat*Point;
    Point += Rot0; 
}



// Prints the shadow of a body at nstep different time steps in [tmin, tmax)
void PrintDynShadow( vector<double> box, Body& body, vector<double> relvel, double dx, double tmin, double tmax, unsigned int nstep, string outfile){
    ProjSurface canvas(box, relvel, dx);
    double dt = nstep < 2 ? 0 : ( tmax - tmin )/nstep;
    double t = tmin;

    for( unsigned int i = 0; i < nstep; i++ ){
        body.Move(t);
        canvas.reset();
        canvas.BodyProj(body);
        string out = outfile + to_string(t) + ".dat";
        cout << "Printing to " << out << endl;
        canvas.PrintRaysFlat(out);
        t += dt;
    }
}



// Prints the smooth shadow of a body at nstep different time steps in [tmin, tmax)
void PrintDynShadowSmooth( vector<double> box, Body& body, vector<double> relvel, double dx, double tmin, double tmax, unsigned int nstep, string outfile){
    ProjSurface canvas(box, relvel, dx);
    double dt = nstep < 2 ? 0 : ( tmax - tmin )/nstep;
    double t = tmin;

    for( unsigned int i = 0; i < nstep; i++ ){
        body.Move(t);
        canvas.reset();
        canvas.BodyProjSmooth(body);
        string out = outfile + to_string(t) + ".dat";
        cout << "Printing to " << out << endl;
        canvas.PrintRaysFlatSmooth(out);
        t += dt;
    }
}



// Prints the state of a body at nstep different time steps in [tmin, tmax)
void PrintDynState( Body& body, double tmin, double tmax, unsigned int nstep, string outfile ) {
    double dt = nstep < 2 ? 0 : ( tmax - tmin )/nstep;
    double t = tmin;

    for( unsigned int i = 0; i < nstep; i++ ){
        body.Move(t);
        string out = outfile + to_string(t) + ".dat";
        cout << "Printing to " << out << endl;
        body.PrintState(out);
        t += dt;
    }
}



// Looks for the minimum of wetness between vmin and vmax with N steps and returns its value if it finds it, else returns -1.
double FindMin( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep ) {
    // Inintialize with limit for vb to inf
    vector<double> vlim = { 1., 0., 0., };
    double wbest = ProjSurface( box, vlim, dx ).BodyProj(body, 0, 1, nstep);
    double vbest = -1;
    
    // Finds minimum in [vmin, vmax] 
    for( size_t i = 0; i < N; i++ ){
        double vb = N < 2 ? vmin : vmin + i*(vmax-vmin)/(N-1);
        vector<double> vrel = rain_v;
        vrel[0] -= vb;
        double wetness = Norm(vrel)*ProjSurface( box, vrel, dx ).BodyProj(body, 0, 1, nstep)/vb;
        if ( wetness < wbest ) {
            wbest = wetness;
            vbest = vb;
        }
    }

    return vbest;
}



// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max], and calculates wetness for n_fit values around it, returns all these values
vector<vector<double>> FindMinFit(vector<double> box, Body& body, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep, unsigned int n_fit, double vcross, double vtail_min, double vtail_max, unsigned int n_tail ) {
    vector<double> vtail, vb, wetness;

    for( size_t i = 0; i < n_tail; i++ ) {
        double vtail_i = n_tail > 1 ?  vtail_min + i*(vtail_max-vtail_min)/(n_tail-1) : vtail_min;
        vector<double> vb_i, wetness_i;

        // Calculates all wetnesses
        for( size_t j = 0; j < N; j++ ) {
            double vb_ij = N > 1 ? vmin + j*(vmax-vmin)/(N-1) : vmin;
            vb_i.push_back(vb_ij);

            vector<double> vrel = { vtail_i, vcross, -1};
            vrel[0] -= vb_ij;
            wetness_i.push_back(Norm(vrel)*ProjSurface( box, vrel, dx ).BodyProj(body, 0, 1, nstep)/vb_ij);
        }

        // Finds the index of the minimum
        auto min_it = min_element(wetness_i.begin(), wetness_i.end());
        int min_index = distance(wetness_i.begin(), min_it);

        // Saves up to n_fit points around min
        int jmin = min_index - n_fit/2;
        for( int j = jmin; j < jmin + (int) n_fit; j++ ) {
            // recycles wetness from those found
            if( 0 <= j and j < (int) vb_i.size() ) {
                vtail.push_back(vtail_i);
                vb.push_back(vb_i[j]);
                wetness.push_back(wetness_i[j]);

            } else if( j >= (int) vb_i.size() ) {       // finds wetness if not in range
                double vbnew = vmax + ( j - (int) vb_i.size() + 1)*(vmax-vmin)/(N-1);
                vector<double> vrel = { vtail_i, vcross, -1};
                vrel[0] -= vbnew;
                double wetnew = Norm(vrel)*ProjSurface( box, vrel, dx ).BodyProj(body, 0, 1, nstep)/vbnew;

                vtail.push_back(vtail_i);
                vb.push_back(vbnew);
                wetness.push_back(wetnew);
            }else if( j < 0 ) {
                double vbnew = vmin + j*(vmax-vmin)/(N-1);
                if ( vbnew > 0 ) {
                    vector<double> vrel = { vtail_i, vcross, -1};
                    vrel[0] -= vbnew;
                    double wetnew = Norm(vrel)*ProjSurface( box, vrel, dx ).BodyProj(body, 0, 1, nstep)/vbnew;

                    vtail.push_back(vtail_i);
                    vb.push_back(vbnew);
                    wetness.push_back(wetnew);
                }
            }
        }
    }

    vector<vector<double>> results{vtail, vb, wetness};
    return Transpose(results);
}



// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max]x[vcross_min, vcross_max], and calculates wetness for n_fit values around it, returns all these values
vector<vector<double>> OptMapFit(vector<double> box, Body& body, double vmin, double vmax, unsigned int N, double dx, unsigned int nstep, unsigned int n_fit, double vtail_min, double vtail_max, unsigned int n_tail, double vcross_min, double vcross_max, unsigned int n_cross ) {
    vector<double> vtail, vcross, vb, wetness;


    for( size_t i = 0; i < n_tail; i++ ) {
        double vtail_i = n_tail > 1 ?  vtail_min + i*(vtail_max-vtail_min)/(n_tail-1) : vtail_min;


        for( size_t j = 0; j < n_cross; j++ ) {
            double vcross_j = n_cross > 1 ?  vcross_min + j*(vcross_max-vcross_min)/(n_cross-1) : vcross_min;
            vector<double> vb_ij, wetness_ij;

            // Calculates all wetnesses
            for( size_t k = 0; k < N; k++ ) {
                double vb_ijk = N > 1 ? vmin + k*(vmax-vmin)/(N-1) : vmin;
                vb_ij.push_back(vb_ijk);

                vector<double> vrel = { vtail_i, vcross_j, -1};
                vrel[0] -= vb_ijk;
                wetness_ij.push_back(Norm(vrel)*ProjSurface( box, vrel, dx ).BodyProj(body, 0, 1, nstep)/vb_ijk);
            }

            // Finds the index of the minimum
            auto min_it = min_element(wetness_ij.begin(), wetness_ij.end());
            int min_index = distance(wetness_ij.begin(), min_it);

            // Saves up to n_fit points around min
            int kmin = min_index - n_fit/2;
            for( int k = kmin; k < kmin + (int) n_fit; k++ ) {
                // recycles wetness from those found
                if( 0 <= k and k < (int) vb_ij.size() ) {
                    vtail.push_back(vtail_i);
                    vcross.push_back(vcross_j);
                    vb.push_back(vb_ij[k]);
                    wetness.push_back(wetness_ij[k]);

                } else if( k >= (int) vb_ij.size() ) {       // finds wetness if not in range
                    double vbnew = vmax + ( k - (int) vb_ij.size() + 1)*(vmax-vmin)/(N-1);
                    vector<double> vrel = { vtail_i, vcross_j, -1};
                    vrel[0] -= vbnew;
                    double wetnew = Norm(vrel)*ProjSurface( box, vrel, dx ).BodyProj(body, 0, 1, nstep)/vbnew;

                    vtail.push_back(vtail_i);
                    vcross.push_back(vcross_j);
                    vb.push_back(vbnew);
                    wetness.push_back(wetnew);
                } else if( k < 0 ) {
                    double vbnew = vmin + k*(vmax-vmin)/(N-1);
                    if ( vbnew > 0 ) {
                        vector<double> vrel = { vtail_i, vcross_j, -1};
                        vrel[0] -= vbnew;
                        double wetnew = Norm(vrel)*ProjSurface( box, vrel, dx ).BodyProj(body, 0, 1, nstep)/vbnew;

                        vtail.push_back(vtail_i);
                        vcross.push_back(vcross_j);
                        vb.push_back(vbnew);
                        wetness.push_back(wetnew);
                    }
                }
            }
        }
    }

    vector<vector<double>> results{vtail, vcross, vb, wetness};
    return Transpose(results);
}



// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max] using Brent algorithm
vector<vector<double>> FindMinBrent(vector<double> box, Body& body, double vmin, double vmax, double dx, unsigned int nstep, double tol, unsigned int n_tries, double vcross, double vtail_min, double vtail_max, unsigned int n_tail ) {
    vector<double> vtail, vb, wetness;
    Brent mins(tol);

    for( size_t i = 0; i < n_tail; i++ ) {
        double vtail_i = n_tail > 1 ?  vtail_min + i*(vtail_max-vtail_min)/(n_tail-1) : vtail_min;
        vtail.push_back(vtail_i);
        vector<double> rain_vel = {vtail_i, vcross, -1};
        auto wetfunc = [&box, &body, &rain_vel, dx, nstep] (double x) {return Wetness( box, body, rain_vel, x, dx, 0., 1., nstep );};

        if( mins.bracket( vmin, vmax, 3, wetfunc ) ) {
            mins.minimize( wetfunc);
            vb.push_back( mins.xmin );
            wetness.push_back( mins.fmin );
        } else {
            vb.push_back( mins.cx );
            wetness.push_back( mins.fc );
        }
        cout << "Step " << i+1 << "/" << n_tail << " completed" << endl;
    }

    return Transpose(vector<vector<double>>{vtail, vb, wetness});
}


// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max]x[vcross_min, vcross_max] using Brent algorithm
vector<vector<double>> OptMapBrent(vector<double> box, Body& body, double vmin, double vmax, double dx, unsigned int nstep, double tol, unsigned int n_tries, double vtail_min, double vtail_max, unsigned int n_tail, double vcross_min, double vcross_max, unsigned int n_cross ) {
    vector<double> vtail, vcross, vb, wetness;
    Brent mins(tol);


    for( size_t i = 0; i < n_tail; i++ ) {
        double vtail_i = n_tail > 1 ?  vtail_min + i*(vtail_max-vtail_min)/(n_tail-1) : vtail_min;

        for( size_t j = 0; j < n_cross; j++ ) {
            double vcross_j = n_cross > 1 ?  vcross_min + j*(vcross_max-vcross_min)/(n_cross-1) : vcross_min;
            vtail.push_back(vtail_i);
            vcross.push_back(vcross_j);
            vector<double> rain_vel = {vtail_i, vcross_j, -1};
            auto wetfunc = [&box, &body, &rain_vel, dx, nstep] (double x) {return Wetness( box, body, rain_vel, x, dx, 0., 1., nstep );};

            if( mins.bracket( vmin, vmax, n_tries, wetfunc ) ) {
                mins.minimize( wetfunc);
                vb.push_back( mins.xmin );
                wetness.push_back( mins.fmin );
            } else {
                vb.push_back( mins.cx );
                wetness.push_back( mins.fc );
            }

            cout << "Step " << i*n_cross+j+1 << "/" << n_tail*n_cross << " completed" << endl;
        }
    }

    return Transpose(vector<vector<double>>{vtail, vcross, vb, wetness});
}


// Transforms distance into a value in [0, 1]
double smooth_w( double delta_r, double dx ) {
    delta_r /= dx;
    if( delta_r >= 1 ) return 0.;
    if( delta_r <= -1 ) return 1.;
    return ( 1 - sin(delta_r*M_PI/2) )/2;
}


// Estimates smooth wetness
double WetnessSmooth( vector<double> box, Body& body, vector<double> rain_v, double vb, double dx ) {
    vector<double> relvel = rain_v;
    relvel[0] -= vb;
    return Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProjSmooth(body)/vb;
}


// Estimates smooth wetness of the dynamic body
double WetnessSmooth( vector<double> box, Body& body, vector<double> rain_v, double vb, double dx, double tmin, double tmax, unsigned int nstep ) {
    vector<double> relvel = rain_v;
    relvel[0] -= vb;
    return Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProjSmooth(body, tmin, tmax, nstep )/vb;
}


// Estimates smooth wetness for N velocities of the dynamic body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimulateSmooth( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx, double tmin, double tmax, unsigned int nstep ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<double> body_v(N);
    vector<double> wetness(N);
    for( size_t i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        wetness[i] = WetnessSmooth( box, body, rain_v, body_v[i], dx, tmin, tmax, nstep );
    }
    return Transpose(vector<vector<double>>{ body_v, wetness});
}


// Estimates smooth wetness for N velocities of the body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimulateSmooth( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx ) {
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<double> body_v(N);
    vector<double> wetness(N);
    for( size_t i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        wetness[i] = WetnessSmooth( box, body, rain_v, body_v[i], dx);
    }
    return Transpose(vector<vector<double>>{ body_v, wetness});
}


// Estimates smooth wetness for N values of dx between dxmin and dxmax, and returns a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimSmoothErr( vector<double> box, Body& body, vector<double> rain_vel, double body_vel, unsigned int N, double dxmin, double dxmax ) {
    vector<double> d = {dxmax};
    vector<double> S = {WetnessSmooth( box, body, rain_vel, body_vel, dxmax)};
    double k =  N == 0 ? 0  :  pow(dxmin/dxmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        d.push_back(d[i-1]*k);
        cout << "dx = " << d[i] << endl;
        S.push_back(WetnessSmooth( box, body, rain_vel, body_vel, d[i]));
    }

    return Transpose(vector<vector<double>>{ d, S});
}

// Estimates smooth wetness for N values of dx between dxmin and dxmax for dynamic body, and returns a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimSmoothErr( vector<double> box, Body& body, vector<double> rain_vel, double body_vel, unsigned int N, double dxmin, double dxmax, double tmin, double tmax, unsigned int nstep) {
    vector<double> dx = {dxmax};
    vector<double> S = {WetnessSmooth( box, body, rain_vel, body_vel, dxmax, tmin, tmax, nstep )};
    double k =  N == 0 ? 0  :  pow(dxmin/dxmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        dx.push_back(dx[i-1]*k);
        cout << "dx = " << dx[i] << endl;
        S.push_back(WetnessSmooth( box, body, rain_vel, body_vel, dx[i], tmin, tmax, nstep ));
    }

    return Transpose(vector<vector<double>>{ dx, S});
}

// Estimates smooth wetness for N values of nstep between nstepmin and nstepmax, and returns a matrix with nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimSmoothErrT( vector<double> box, Body& body, vector<double> rain_vel, double body_vel, double dx, double tmin, double tmax, unsigned int N, unsigned int nstepmin, unsigned int nstepmax) {
    double dtmin = 1./nstepmax;
    double dtmax = 1./nstepmin;
    vector<double> dt = {dtmax};
    vector<double> nstep = { (double) nstepmin};
    vector<double> S = {WetnessSmooth( box, body, rain_vel, body_vel, dx, tmin, tmax, nstepmin )};
    double k =  N == 0 ? 0  :  pow(dtmin/dtmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        dt.push_back(dt[i-1]*k);
        double nstepnew = round(1./dt[i]);
        if( nstepnew > nstep.back()) {
            cout << "Nstep = " << nstepnew << endl;
            nstep.push_back(nstepnew);
            S.push_back(WetnessSmooth( box, body, rain_vel, body_vel, dx, tmin, tmax, nstepnew ));
        }
    }

    return Transpose(vector<vector<double>>{ nstep, S});
}


// Finds minimums of wetness for a fixed vcross and [vtail_min, vtail_max] using Brent algorithm
vector<vector<double>> FindMinBrentSmooth(vector<double> box, Body& body, double vmin, double vmax, double dx, unsigned int nstep, double tol, unsigned int n_tries, double vcross, double vtail_min, double vtail_max, unsigned int n_tail ) {
    vector<double> vtail, vb, wetness;
    Brent mins(tol);

    for( size_t i = 0; i < n_tail; i++ ) {
        double vtail_i = n_tail > 1 ?  vtail_min + i*(vtail_max-vtail_min)/(n_tail-1) : vtail_min;
        vtail.push_back(vtail_i);
        vector<double> rain_vel = {vtail_i, vcross, -1};
        auto wetfunc = [&box, &body, &rain_vel, dx, nstep] (double x) {return WetnessSmooth( box, body, rain_vel, x, dx, 0., 1., nstep );};

        if( mins.bracket( vmin, vmax, 3, wetfunc ) ) {
            mins.minimize( wetfunc);
            vb.push_back( mins.xmin );
            wetness.push_back( mins.fmin );
        } else {
            vb.push_back( mins.cx );
            wetness.push_back( mins.fc );
        }
        cout << "Step " << i+1 << "/" << n_tail << " completed" << endl;
    }

    return Transpose(vector<vector<double>>{vtail, vb, wetness});
}


// Estimates smooth wetness for N_v velocities of the dynamic body between vmin and vmax, and N_nstep values of nstep between nstep_min and nstep_max and returns a matrix with the nsteps velocities as the first colunmn, the velocities as the second column and the respective wetness as the third column
vector<vector<double>> SimulateNstepSmooth( vector<double> box, Body& body, vector<double> rain_v, double vmin, double vmax, unsigned int N_v, double dx, unsigned int nstep_min, unsigned int nstep_max, unsigned int N_nstep ) {
    double dtmin = 1./nstep_max;
    double dtmax = 1./nstep_min;
    double k =  N_nstep < 2 ? 1  :  pow(dtmin/dtmax, (double)1/(N_nstep-1));
    vector<double> nstep, vb, wetness;
    
    for( size_t i = 0; i < N_nstep; i++ ) {
        double dt_i = dtmax * pow( k, i );
        unsigned int nstep_i = round(1./dt_i);
        for( size_t j = 0; j < N_v; j++ ) {
            double vb_j = ( N_v < 2 ? vmin : vmin + (vmax - vmin)*(double)j/((double)N_v-1) );
            nstep.push_back(nstep_i);
            vb.push_back(vb_j);
            wetness.push_back( WetnessSmooth( box, body, rain_v, vb_j, dx, 0., 1., nstep_i) );
        }
    }

    return Transpose(vector<vector<double>>{ nstep, vb, wetness });
}


// Fits points with a parabola y = ax^2 + bx + c, using least squares minimization, and returns a tuple containing (a, b, c)
tuple<double, double, double> QuadraticRegression( vector<double> x_vals, vector<double> y_vals ){
    size_t n = x_vals.size();
    assert( y_vals.size() == n && "Number of x and y values must be equal!" );

    // Build design matrix
    vector<vector<double>> X( n, vector<double>(3) );
    for( size_t i = 0; i < n; i++ ) {
        for( size_t j = 0; j < 3; j++ ){
           X[i][j] = pow( x_vals[i], j );
        }
    }

    // Evaluate coefficients
    vector<vector<double>> X_t = Transpose(X);
    vector<double> beta = Inverse( X_t*X ) * X_t * y_vals;

    return make_tuple( beta[2], beta[1], beta[0] );
}


// Finds minimums of smooth wetness for a fixed vcross and [vtail_min, vtail_max] using Brent algorithm, and calculates wetness for n_fit values spaced dv around it, returns all these values
vector<vector<double>> FindMinFitSmooth(vector<double> box, Body& body, double vmin, double vmax, double dx, unsigned int nstep, double vcross, double vtail_min, double vtail_max, unsigned int n_tail, int n_fit, double dv ) {
    vector<double> vtail, vb, wetness;
    Brent mins(dv);

    for( size_t i = 0; i < n_tail; i++ ) {
        double vtail_i = n_tail > 1 ?  vtail_min + i*(vtail_max-vtail_min)/(n_tail-1) : vtail_min;
        vector<double> rain_vel = {vtail_i, vcross, -1};
        auto wetfunc = [&box, &body, &rain_vel, dx, nstep] (double x) {return WetnessSmooth( box, body, rain_vel, x, dx, 0., 1., nstep );};

        if( mins.bracket( vmin, vmax, 3, wetfunc ) ) {
            vector<double> vtail_vec( n_fit, vtail_i );
            vector<double> x, fx;
            mins.minimize( wetfunc);

            // Evaluate n_fit points around minimum
            for( int j = -(n_fit-1)/2; j <= n_fit/2; j++ ) {  
                double x_j = mins.xmin + j*dv;
                double fx_j = ( j == 0 ) ? mins.fmin : wetfunc(x_j);
                x.push_back( x_j );
                fx.push_back( fx_j );
            }

            // Try quadratic fit, if resulting minimum isn't in the middle of evaluated points move towards it and repeat until it is
            int max_iter = 100;
            int n_half = n_fit/2;
            for( int iter = 0; iter < max_iter; iter++ ) {
                double a_fit, b_fit, c_fit, min_fit;
                tie( a_fit, b_fit, c_fit ) = QuadraticRegression( x, fx );

                min_fit = ( a_fit != 0 )?  -b_fit/(2*a_fit) : -sgn(b_fit)*2*vmax;        // moving towards -sgn(b_fit)*2*vmax means moving following first derivative

                // Avoid following minimum outside of range
                if( a_fit > 0 && ( min_fit <= vmin || min_fit >= vmax ) ) break;

                // Find number of values of x that are lower (and higher) than min_fit
                int n_lower = 0;
                while( n_lower < n_fit && x[n_lower] < min_fit ) n_lower ++;
                int n_higher = n_fit - n_lower;

                // Account for failed fit, if a_fit < 0 move in  opposite direction
                if( a_fit < 0 ) swap( n_lower, n_higher );

                cout << "Iteration " << iter + 1 << ", n_lower = " << n_lower << ", n_higher = " << n_higher << endl;

                // Move towards minimum
                if( n_lower < n_half ) {
                    double new_x = x[0] - dv;
                    x.insert( x.begin(), new_x );
                    fx.insert( fx.begin(), wetfunc( new_x ) );
                    x.pop_back();
                    fx.pop_back();
                } else if( n_higher < n_half ) {
                    double new_x = x.back() + dv;
                    x.push_back( new_x );
                    fx.push_back( wetfunc( new_x ) );
                    x.erase( x.begin() );
                    fx.erase( fx.begin() );
                } else break;
            }

            // Add resulting points to output
            vb.insert( vb.end(), x.begin(), x.end() );
            wetness.insert( wetness.end(), fx.begin(), fx.end() );
            vtail.insert( vtail.end(), vtail_vec.begin(), vtail_vec.end() );

        } else {    // If bracketing fails the bound minimum is in vb = vmax
            vb.push_back( mins.cx );
            wetness.push_back( mins.fc );
            vtail.push_back( vtail_i );
        }
        cout << "Step " << i+1 << "/" << n_tail << " completed" << endl;
    }

    return Transpose(vector<vector<double>>{vtail, vb, wetness});
}


// Finds minimums of smooth wetness for a fixed vcross and vtail_min using Brent algorithm with nstep in [nstep_min, nstep_max], and calculates wetness for n_fit values spaced dv around it, returns all these values
vector<vector<double>> FindMinFitSmooth(vector<double> box, Body& body, double vmin, double vmax, double dx, unsigned int nstep_min, unsigned int nstep_max, unsigned int N_nstep, double vcross, double vtail, int n_fit, double dv ) {
    double dtmin = 1./nstep_max;
    double dtmax = 1./nstep_min;
    double k =  N_nstep < 2 ? 1  :  pow(dtmin/dtmax, (double)1/(N_nstep-1));

    vector<double> nstep, vb, wetness;
    Brent mins(dv);
    vector<double> rain_vel = {vtail, vcross, -1};

    for( size_t i = 0; i < N_nstep; i++ ) {
        double dt_i = dtmax * pow( k, i );
        unsigned int nstep_i = round(1./dt_i);
        auto wetfunc = [&box, &body, &rain_vel, dx, nstep_i] (double x) {return WetnessSmooth( box, body, rain_vel, x, dx, 0., 1., nstep_i );};

        if( mins.bracket( vmin, vmax, 3, wetfunc ) ) {
            vector<double> nstep_vec( n_fit, nstep_i );
            vector<double> x, fx;
            mins.minimize( wetfunc);
            
            // Evaluate n_fit points around minimum
            for( int j = -(n_fit-1)/2; j <= n_fit/2; j++ ) {  
                double x_j = mins.xmin + j*dv;
                double fx_j = ( j == 0 ) ? mins.fmin : wetfunc(x_j);
                x.push_back( x_j );
                fx.push_back( fx_j );
            }

            // Try quadratic fit, if resulting minimum isn't in the middle of evaluated points move towards it and repeat until it is
            int max_iter = 100;
            int n_half = n_fit/2;
            for( int iter = 0; iter < max_iter; iter++ ) {
                double a_fit, b_fit, c_fit, min_fit;
                tie( a_fit, b_fit, c_fit ) = QuadraticRegression( x, fx );

                min_fit = ( a_fit != 0 )?  -b_fit/(2*a_fit) : -sgn(b_fit)*2*vmax;        // moving towards -sgn(b_fit)*2*vmax means moving following first derivative

                // Avoid following minimum outside of range
                if( a_fit > 0 && ( min_fit <= vmin || min_fit >= vmax ) ) break;

                // Find number of values of x that are lower (and higher) than min_fit
                int n_lower = 0;
                while( n_lower < n_fit && x[n_lower] < min_fit ) n_lower ++;
                int n_higher = n_fit - n_lower;

                // Account for failed fit, if a_fit < 0 move in  opposite direction
                if( a_fit < 0 ) swap( n_lower, n_higher );

                cout << "Iteration " << iter + 1 << ", n_lower = " << n_lower << ", n_higher = " << n_higher << endl;

                // Move towards minimum
                if( n_lower < n_half ) {
                    double new_x = x[0] - dv;
                    x.insert( x.begin(), new_x );
                    fx.insert( fx.begin(), wetfunc( new_x ) );
                    x.pop_back();
                    fx.pop_back();
                } else if( n_higher < n_half ) {
                    double new_x = x.back() + dv;
                    x.push_back( new_x );
                    fx.push_back( wetfunc( new_x ) );
                    x.erase( x.begin() );
                    fx.erase( fx.begin() );
                } else break;
            }

            // Add remaining points to output
            vb.insert( vb.end(), x.begin(), x.end() );
            wetness.insert( wetness.end(), fx.begin(), fx.end() );
            nstep.insert( nstep.end(), nstep_vec.begin(), nstep_vec.end() );

        } else {
            vb.push_back( mins.cx );
            wetness.push_back( mins.fc );
            nstep.push_back( nstep_i );
        }
        cout << "Step " << i+1 << "/" << N_nstep << " completed" << endl;
    }

    return Transpose(vector<vector<double>>{nstep, vb, wetness});
}