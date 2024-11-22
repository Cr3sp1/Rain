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
    // Centers on p
    Point -= H[0];
    for(int i = 1; i < 7; i++ ){
        H[i] -= H[0];
    }

    // Checks if Point is inside the rectangle with vertices H[0], H[i], H[i+1]
    for( int i = 1; i < 7; i++ ){
        double epsilon = 1e-10;
        double A = Norm( CrossProduct( H[i], H[PBCH(i+1)]) );
        double alpha = Norm( CrossProduct( Point, H[PBCH(i+1)]) )/A;
        double beta = Norm( CrossProduct( Point, H[i]) )/A;
        double gamma = Norm( CrossProduct( Point-H[i], Point-H[PBCH(i+1)]) )/A;
        if( 0 <= alpha and alpha <= 1 and 0 <= beta and beta <= 1 and 0 <= gamma and gamma <= 1 and abs(alpha + beta + gamma - 1 ) <  epsilon ){
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
vector<vector<double>> SimErr( vector<double> box, Body& body, vector<double> relvel, double bodyvel, unsigned int N, double dxmin, double dxmax) {
    vector<double> d = {dxmax};
    vector<double> S = {ProjSurface( box, relvel, d[0] ).BodyProj(body)*Norm(relvel)/bodyvel};
    double k =  N == 0 ? 0  :  pow(dxmin/dxmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        d.push_back(d[i-1]*k);
        cout << "dx = " << d[i] << endl;
        S.push_back(ProjSurface( box, relvel, d[i] ).BodyProj(body)*Norm(relvel)/bodyvel);
    }

    return Transpose(vector<vector<double>>{ d, S});
}

// Estimates wetness for N values of dx between dxmin and dxmax for dynamic body, and returns a matrix with dx as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErr( vector<double> box, Body& body, vector<double> relvel, double bodyvel, unsigned int N, double dxmin, double dxmax, double tmin, double tmax, unsigned int nstep) {
    vector<double> dx = {dxmax};
    vector<double> S = {ProjSurface( box, relvel, dx[0] ).BodyProj(body)*Norm(relvel)/bodyvel};
    double k =  N == 0 ? 0  :  pow(dxmin/dxmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        dx.push_back(dx[i-1]*k);
        cout << "dx = " << dx[i] << endl;
        S.push_back(ProjSurface( box, relvel, dx[i] ).BodyProj(body, tmin, tmax, nstep )*Norm(relvel)/bodyvel);
    }

    return Transpose(vector<vector<double>>{ dx, S});
}

// Estimates wetness for N values of nstep between nstepmin and nstepmax, and returns a matrix with nstep as the first colunmn and the respective wetness as the second column
vector<vector<double>> SimErrT( vector<double> box, Body& body, vector<double> relvel, double bodyvel, double dx, double tmin, double tmax, unsigned int N, unsigned int nstepmin, unsigned int nstepmax) {
    double dtmin = 1./nstepmax;
    double dtmax = 1./nstepmin;
    vector<double> dt = {dtmax};
    vector<double> nstep = { (double) nstepmin};
    vector<double> S = {ProjSurface( box, relvel, dx ).BodyProj(body, tmin, tmax, nstep[0] )*Norm(relvel)/bodyvel};
    double k =  N == 0 ? 0  :  pow(dtmin/dtmax, (double)1/(N-1));
    
    for( size_t i = 1; i < N; i++ ) {
        dt.push_back(dt[i-1]*k);
        double nstepnew = round(1./dt[i]);
        if( nstepnew > nstep.back()) {
            cout << "Nstep = " << nstepnew << endl;
            nstep.push_back(nstepnew);
            S.push_back(ProjSurface( box, relvel, dx ).BodyProj(body, tmin, tmax, nstepnew )*Norm(relvel)/bodyvel);
        }
    }

    return Transpose(vector<vector<double>>{ nstep, S});
}

// Estimate wetness for N_dx values of dx between dxmin and dxmax, and N_t values of nstep between nstepmin and nstepmax, and returns a matrix where first column is the nstep and first row is the dx
vector<vector<double>> SimErrTdx( vector<double> box, Body& body, vector<double> relvel, double bodyvel, unsigned int N_dx, double dxmin, double dxmax, unsigned int N_t, unsigned int nstepmin, unsigned int nstepmax) {
    double k_dx =  N_dx == 0 ? 0  :  pow(dxmin/dxmax, (double)1/(N_dx-1));
    double k_t =  N_t == 0 ? 0  :  (double)(nstepmax-nstepmin)/(N_t-1);

    // Build indices
    vector<double> index_dx = { 0, dxmax };
    for( size_t i = 1; i < N_dx; i++ ) index_dx.push_back(index_dx[i]*k_dx);

    vector<double> index_t = { (double)nstepmin };
    for( size_t i = 1; i < N_t; i++ ) index_t.push_back(index_t[i-1]+k_t);

    vector<vector<double>> results = {index_dx};
    for( double& nstep : index_t ) results.push_back({floor(nstep)});

    // Fill the matrix
    for( size_t i = 1; i <= N_t; i++ ) {
        cout << "Nstep = " << results[i][0] << endl;
        for( size_t j = 1; j <= N_dx; j++) {
            results[i].push_back(ProjSurface( box, relvel, results[0][j] ).BodyProj(body, 0, 1, results[i][0] )*Norm(relvel)/bodyvel);
        }
    }

    return results;
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

// Estimates wetness for N velocities of two body between vmin and vmax, and returns a matrix with the velocities as the first colunmn and the wetness of the first body as the second column and of the second body as the third column
vector<vector<double>> CompareBB( vector<double> box, Body& body1, Body& body2, vector<double> rain_v, double vmin, double vmax, unsigned int N, double dx){
    if( vmin > vmax or vmin < 0 ) cout << "Error: Vmin and Vmax have to be positive and Vmax > Vmin!" << endl;
    vector<double> body_v(N);
    vector<double> wetness1(N);
    vector<double> wetness2(N);
    for( size_t i = 0; i < N; i++ ){
        body_v[i] = ( N == 1 ? vmin : vmin + (vmax - vmin)*(double)i/((double)N-1) );
        vector<double> relvel = rain_v;
        relvel[0] -= body_v[i];
        wetness1[i] = Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body1)/body_v[i];
        wetness2[i] = Norm(relvel)*ProjSurface( box, relvel, dx ).BodyProj(body2)/body_v[i];
    }
    vector<vector<double>> mat{body_v, wetness1, wetness2};
    return Transpose(mat);
}



// Returns the minimum distance between the point p and the segment line with extremes l1 and l2
double PointSegDist( vector<double> p, vector<double> l1, vector<double> l2 ) {
    // Changes frame of reference to l1 = 0
    p -= l1;
    l2 -= l1;
    if( Norm(l2) == 0 ) return Norm(p);
    // Calculates projection of p on line passing through l1 (0) and l2
    double proj = p*l2/Norm(l2);
    // Returns distance between p and the closest point belonging to the segment
    if( proj <= 0 )  return Norm(p);
    if( proj >= Norm(l2) ) return Norm(p-l2);
    // cout << "ok" << endl;
    return Norm( p - proj*l2/Norm(l2));
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
        canvas.BodyProj(body);
        string out = outfile + to_string(t) + ".dat";
        cout << "Printing to " << out << endl;
        canvas.PrintRaysFlat(out);
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
vector<vector<double>> FindMinBrent(vector<double> box, Body& body, double vmin, double vmax, double dx, unsigned int nstep, double tol, double vcross, double vtail_min, double vtail_max, unsigned int n_tail ) {
    vector<double> vtail, vb, wetness;
    Brent mins(tol);

    for( size_t i = 0; i < n_tail; i++ ) {
        double vtail_i = n_tail > 1 ?  vtail_min + i*(vtail_max-vtail_min)/(n_tail-1) : vtail_min;
        vtail.push_back(vtail_i);
        vector<double> rain_vel = {vtail_i, vcross, -1};
        auto wetfunc = [&box, &body, &rain_vel, dx, nstep] (double x) {return Wetness( box, body, rain_vel, x, dx, 0., 1., nstep );};

        if( mins.bracket( vmin + tol, vmax - 2*tol, vmax, wetfunc ) ) {
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
vector<vector<double>> OptMapBrent(vector<double> box, Body& body, double vmin, double vmax, double dx, unsigned int nstep, double tol, double vtail_min, double vtail_max, unsigned int n_tail, double vcross_min, double vcross_max, unsigned int n_cross ) {
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

            if( mins.bracket( vmin + tol, vmax - 2*tol, vmax, wetfunc ) ) {
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



