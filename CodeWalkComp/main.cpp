#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "body.h"
#include "ray.h"
# include "VectorStat.h"


using namespace std;


int main (int argc, char *argv[]){

    // Declare stuff
    vector<double> rain_vel(3);
    vector<double> rel_vel(3);
    vector<double> box(3); 
    double body_vel, dx;
    unsigned int nstep_v, nstep_t;
    int N_vb, N_fit;
    

    // Reads the input file
    ifstream ReadInput;
    ReadInput.open("input.in"); 
    ReadInput >> box[0] >> box[1] >> box[2];
    ReadInput >> dx;
    ReadInput >> nstep_t;
    ReadInput >> body_vel;
    ReadInput >> nstep_v;
    ReadInput >> N_vb;
    ReadInput >> N_fit;
    ReadInput >> rain_vel[0] >> rain_vel[1] >> rain_vel[2];
    rel_vel = rain_vel;
    rel_vel[0] -= body_vel;

    ReadInput.close();


    // Outputs settings and stuff
    cout << "relative velocity = [ " << rel_vel[0] << ", " << rel_vel[1] << ", " << rel_vel[2] << " ]" << endl;
    cout << "Object velocity = " << body_vel  << ", nstep_v = " << nstep_v << ", dx = " << dx << endl;


    // Walking man
    ManyBody Walk("../Bodies/WalkingMan.in");
    vector<double> boxW = {0.90, 0.56, 1.74 };

    // Running man
    ManyBody Run("../Bodies/RunningMan.in");
    vector<double> boxR = {1.13, 0.58, 1.82};


    // Check boxes
    // PrintDynShadow(boxW, Walk, {0, 0, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_xy" );
    // PrintDynShadow(boxW, Walk, {0, -100, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_xz" );
    // PrintDynShadow(boxW, Walk, {-100, 0, -1}, dx, 0, 1, 60, "../data/Walk/Proj/Walk_yz" );
    // PrintDynShadow(boxR, Run, {0, 0, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_xy" );
    // PrintDynShadow(boxR, Run, {0, -100, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_xz" );
    // PrintDynShadow(boxR, Run, {-100, 0, -1}, dx, 0, 1, 60, "../data/Run/Proj/Run_yz" );

    // // Print dynamic state
    // PrintDynState( Walk, 0, 1, 60, "../data/Walk/Status/Walk");
    // PrintDynState( Run, 0, 1, 60, "../data/Run/Status/Run");

    // Walking & running error analysis
    // vector<vector<double>> WalkResDxT = SimErrTdx(boxW, Walk, rel_vel, body_vel, nstep_v, dx, 0.1, nstep_t, 1, nstep_t );
    // Print( "../data/Walk/ErrDxT.dat", WalkResDxT, 12 );
    // vector<vector<double>> WalkResT = SimErrT(boxW, Walk, rel_vel, body_vel, dx, 0, 1, 100, 1, nstep_t);
    // Print( "../data/Walk/ErrT.dat", WalkResT, 12 );
    // vector<vector<double>> WalkResDx = SimErr(boxW, Walk, rel_vel, body_vel, nstep_v, dx, 1, 0, 1, nstep_t);
    // Print( "../data/Walk/ErrDx.dat", WalkResDx, 12 );

    // vector<vector<double>> RunResDxT = SimErrTdx(boxR, Run, rel_vel, body_vel, nstep_v, dx, 0.1, nstep_t, 1, nstep_t );
    // Print( "../data/Run/ErrDxT.dat", RunResDxT, 12 );
    // vector<vector<double>> RunResT = SimErrT(boxR, Run, rel_vel, body_vel, dx, 0, 1, 100, 1, nstep_t);
    // Print( "../data/Run/ErrT.dat", RunResT, 12 );
    // vector<vector<double>> RunResDx = SimErr(boxW, Run, rel_vel, body_vel, nstep_v, dx, 1, 0, 1, nstep_t);
    // Print( "../data/Run/ErrDx.dat", RunResDx, 12 );
    


    // Error analysis sphere
    // ManyBody TrialS("../Bodies/Sphere.in");
    // vector<vector<double>> resultsS;
    // resultsS = SimErr( box, TrialS, rel_vel, body_vel, 200, 0.0001, 1 );
    // Print( "../data/Sphere/ErrorS.dat", resultsS, 15);
        
    // Error analysis pippo
    // ManyBody TrialP("../Bodies/Pippo.in");
    // vector<vector<double>> resultsP;
    // resultsP = SimErr( box, TrialP, rel_vel, body_vel, 200, 0.0001, 1 );
    // Print( "../data/Pippo/ErrorP.dat", resultsP, 15);

    // Error analysis capsule
    // ManyBody TrialC("../Bodies/Capsule.in");
    // vector<vector<double>> resultsC;
    // resultsC = SimErr( box, TrialC, rel_vel, body_vel, 200, 0.0001, 1 );
    // Print( "../data/Capsule/ErrorC.dat", resultsC, 15);
    


    // // Draw shadow
    // ProjSurface Plz( box, rel_vel, dx );
    // Plz.BodyProj(trialM1);
    // Plz.PrintR("../data/RayOriginsCap.dat");
    // Plz.PrintRaysFlat("../data/RayOriginsCapF.dat");



    // Simulate different body velocities
    // resultsS = CompareAN( box, *TrialS.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Sphere/CompareS.dat", resultsS, 15);
    // resultsP = CompareAN( box, *TrialP.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Pippo/CompareP.dat", resultsP, 15);
    // resultsC = CompareAN( box, *TrialC.Find("Name"), rain_vel, 1, 10, nstep_v, dx );
    // Print("../data/Capsule/CompareC.dat", resultsC, 15);



    // Simulation of two pippos compenetrating
    // ManyBody Trial2P("../Bodies/DoublePippo.in");
    // vector<double> dist;
    // vector<double> wet2P;
    // ProjSurface Surf2P( box, rel_vel, dx );
    // for( size_t i = 0; i < nstep_t; i++ ) {
    //     Trial2P.Move( asin( (double)i/(nstep_t-1))/(2*M_PI));   // it just works ;)
    //     vector<double> cent1 = dynamic_cast<Pippo*>(Trial2P.Find("Still"))->GetCent();
    //     vector<double> cent2 = dynamic_cast<Pippo*>(Trial2P.Find("Moving"))->GetCent();
    //     dist.push_back( Norm( cent1 - cent2 ));
    //     Surf2P.reset();
    //     wet2P.push_back(Surf2P.BodyProj(Trial2P)*Norm(rel_vel)/body_vel);
    //     // Surf2P.PrintRaysFlat("../data/Pippo/Proj2P/dist" + to_string(dist[i]) + ".dat");
    // }

    // vector<vector<double>> results2P = { dist, wet2P};
    // results2P = Transpose(results2P);
    // Print("../data/Pippo/DoubleP.dat", results2P, 12 );
    


    // Find v_opt fit map
    
    // vector<vector<double>> WalkOptMapFit = OptMapFit( boxW, Walk, 0.7.*30./49., 0.7, N_vb-30, dx, nstep_t, N_fit,  0, 0.7, nstep_v, 0, 0.35, nstep_v );
    // Print( "../data/Walk/OptMapFitW.dat", WalkOptMapFit, 12 );

    // vector<vector<double>> RunOptMapFit = OptMapFit( boxR, Run, 2.*17./49., 2, N_vb-17, dx, nstep_t, N_fit, 0, 2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Run/OptMapFitR.dat", RunOptMapFit, 12 ); 



    // Compare walking and running
    vector<vector<double>> WalkOptMapComp = OptMapFit( boxW, Walk, 0.7*30./49., 0.7, N_vb-30, dx, nstep_t, N_fit,  -0.5, 2, nstep_v, 0, 1.2, nstep_v );
    Print( "../data/Walk/OptMapCompareW.dat", WalkOptMapComp, 12 );

    // vector<vector<double>> RunOptMapComp = OptMapFit( boxR, Run, 2.*17./49., 2, N_vb-17, dx, nstep_t, N_fit,  -0.5, 2, nstep_v, 0, 1.2, nstep_v );
    // Print( "../data/Run/OptMapCompareR.dat", RunOptMapComp, 12 );



    // V_opt fit graphs
    // vector<vector<double>> WalkMins1 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit, 0., 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW0.dat", WalkMins1, 12 );

    // vector<vector<double>> WalkMins2 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit,  0.15, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW015.dat", WalkMins2, 12 );

    // vector<vector<double>> WalkMins3 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit, 0.25, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW025.dat", WalkMins3, 12 );

    // vector<vector<double>> WalkMins4 = FindMinFit( boxW, Walk, 0, 0.7, N_vb, dx, nstep_t, N_fit, 0.30, 0., 0.7, nstep_v );
    // Print( "../data/Walk/OptFitW030.dat", WalkMins4, 12 );



    // vector<vector<double>> RunMins1 = FindMinFit( boxR, Run, 2.*17./49., 2, N_vb-17, dx, nstep_t, N_fit,  0, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR0.dat", RunMins1, 12 );

    // vector<vector<double>> RunMins2 = FindMinFit( boxR, Run, 2.*17./49., 2, N_vb-17, dx, nstep_t, N_fit,  0.5, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR05.dat", RunMins2, 12 );

    // vector<vector<double>> RunMins3 = FindMinFit( boxR, Run, 2.*17./49., 2, N_vb-17, dx, nstep_t, N_fit,  0.75, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR075.dat", RunMins3, 12 );

    // vector<vector<double>> RunMins4 = FindMinFit( boxR, Run, 2.*17./49., 2, N_vb-17, dx, nstep_t, N_fit,  1, 0, 2, nstep_v );
    // Print( "../data/Run/OptFitR1.dat", RunMins4, 12 );





    // For fitting graph 
    // vector<vector<double>> RunFitGraph = WetFit(boxR, Run, 0, 2, N_vb, dx, nstep_t, rain_vel[0], rain_vel[1] );
    // Print( "../data/Run/GraphFit.dat", RunFitGraph, 12 );


    
    return 0;
}

