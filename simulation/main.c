/*
Magneto-Optical Trapping Simulation - Header file
Code release version: v1.2


Monte Carlo trajectory simulation for narrow-line MOTs.
Generates 2D x-z histogram files and trapping time files.



Ramon Gabriel Teixeira Rosa, PhD
Optics and Photonics Research Center (CePOF)
University of Sao Paulo - Sao Carlos Institute of Physics
ramongabriel.tr@usp.br
+55(16)3373.9810 (Ext: 225)

July 13, 2018
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "MOTsimHeader_v1.h"


///-- TITLE FOR OUTPUT FILE HEADER --///
char TITLE[] = "TESTE_GEOM_OPTIM_0123456789012345678901234567890123456789012345678901234567890123456789";
//without_spaces
//will be used on output folder name



///------------------------------------------------ CONSTANTS ----------------------------------------------------///
const double h = 6.62607004e-34;
const double mub = 9.274009994e-24;
const double e = 1.60217662e-19;
const double kb = 1.38064852e-23;
const double c = 2.99792458e8;
// SI units
///---------------------------------------------------------------------------------------------------------------///
///------------------------------------------------ ENVIRONMENT---------------------------------------------------///
double gx=+0.0;
double gy=+0.0;
double gz=-9.8;
// gravity [m/s2]
double GBx = 0.046/2;
double GBy = 0.046/2;
double GBz = -0.046;
// Magnetic field gradient [T/m]   (1 G/cm = 0.01 T/m)
double B0[] = {0,0,0}; // DO NOT CHANGE UNTIL MODIFYING MAGNETIC FORCE ON DIPOLE!!!
// Magnetic field bias [T]   (1 G = 1e-4 T)
double detuning = -70;
// units of Gamma
#define NUMBEAMS 3
//double BEAMS[][3] = {
//        {-1, 0, 0},
//        {+1, 0, 0},
//        { 0,-1, 0},
//        { 0,+1, 0},
//        { 0, 0,-1},
//        { 0, 0,+1},
//};
//double BEAMS[][3] = {
//        {-1, 0,+1},
//        {+1, 0,+1},
//        { 0,-1, 0},
//        { 0,+1, 0},
//};
//double BEAMS[][3] = {
//        {-1,0,sqrt(2)/2},
//        {+0.5,-sqrt(3)/2,sqrt(2)/2},
//        {+0.5,+sqrt(3)/2,sqrt(2)/2},
//        {+1,0,-sqrt(2)/2},
//        {-0.5,+sqrt(3)/2,-sqrt(2)/2},
//        {-0.5,-sqrt(3)/2,-sqrt(2)/2}
//};
double BEAMS[][3] = {
        {-1.0,0,0.4*sqrt(2)/2.0},
        {+0.5,-sqrt(3)/2.0,0.4*sqrt(2)/2.0},
        {+0.5,+sqrt(3)/2.0,0.4*sqrt(2)/2.0}
};
// All beams cross the point (0,0,0);
const char BEAMS_POL[NUMBEAMS] = "lll";
//const char BEAMS_POL[NUMBEAMS] = "rrrrrr";
//const char BEAMS_POL[NUMBEAMS] = "llllrr";
//const char BEAMS_POL[NUMBEAMS] = "llrr";
//const char BEAMS_POL[NUMBEAMS] = "llrrrr";  // 'r': RHCP
                                            // 'l': LHCP
                                            // 'x': polarization vector on plane defined by X and ek (beam propagation vector)
                                            // 'y': [...]
                                            // 'z': [...]
//const double so[] = {160,160,160};
//const double so[] = {0.6,0.6,0.6};
//const double so[] = {50,50,50,50,50,50};
//const double so[] = {160,160,160,160,160,160};
const double so[] = {160,160,160};
//const double so[] = {0.6,0.6,0.6,0.6,0,0.6};
//const double so[] = {0.6,0.6,0.6,0.6};
//const double so[] = {0.17,0.17,0.17,0.17,0.17,0.17};
// peak intensity/saturation_intensity
//const double BEAMS_WAIST[] = {0.036,0.036,0.036};
//const double BEAMS_WAIST[] = {0.036,0.036,0.036,0.036,0.036,0.036};
const double BEAMS_WAIST[] = {0.036,0.036,0.036};
//const double BEAMS_WAIST[] = {0.036,0.036,0.036,0.036};
// 1/e2 gaussian beam waist [m]

const int NUM_SIDEBANDS = 30; //30
// Each side of central band
const double FreqStep_SIDEBANDS = 105e3;
// [Hz]


double complex BEAMS_POL_VECTORS[NUMBEAMS][3];
///---------------------------------------------------------------------------------------------------------------///
///------------------------------------------------ ATOM ---------------------------------------------------------///
// Dy
const double m = (1.660539040e-27)*163.9291748;
const double Gamma = 136e3;
const double lambda = 626e-9;
const double Jgnd = 8;
const double Jexc = 9;
const double gjg = 1.24;
const double gje = 1.29;

//Er
//const double m = (1.660539040e-27)*166.0;
//const double Gamma = 190e3;
//const double lambda = 583e-9;
//const double Jgnd = 6;
//const double Jexc = 7;
//const double gjg = 1.167;
//const double gje = 1.195;

// Cs
//const double m = (1.660539040e-27)*132.90545;
//const double Gamma = 5.2e6;
//const double lambda = 852e-9;
//const double Jgnd = 4;
//const double Jexc = 5;
//const double gjg = 0.25;
//const double gje = 0.40;
///---------------------------------------------------------------------------------------------------------------///
///------------------------------------------------ INITIAL CONDITIONS -------------------------------------------///
double xx0c = 0;
double yy0c = 0;
double zz0c = 0;
// initial position
double vx0c = 0;
double vy0c = 0;
double vz0c = 0;
// inital velocity
double STDxx0 = 2e-3;
double STDyy0 = 2e-3;
double STDzz0 = 2e-3;
// standard deviation of initial position gaussian distribution (around [xx0c,yy0c,zz0c])
double T = 50e-6;
// Temperature [K]
// Velocity distribution following Maxwell-Boltzmann distribution with offset [vx0c,vy0c,vz0c];
double xx0,yy0,zz0,vx0,vy0,vz0;
///---------------------------------------------------------------------------------------------------------------///
///------------------------------------------------ SIMULATION ---------------------------------------------------///
const int NUMITER=1e5;
//maximum number of iterations
const double BOUNDARY = 30e-3;
const double xc = 0;
const double yc = 0;
const double zc = 0;
// if( sqrt((x-xc)^2 + (y-yc)^2 + (z-zc)^2) > BOUNDARY ), assume atom escaped trap
// [m]
const int NumVar = 31;
// number of different configurations (TRAP_DEPTH_ROUTINE() and TRAPPING_TIME_ROUTINE())
const int NumVel = 30;
// numer of initial velocities (TRAP_DEPTH_ROUTINE())
const double MaxVel_TrapDepth = 7;
// maximum velocity for trap depth routine [m/s]
const int NumAvg = 50;
// averages per configuration
const double DetuningRange[] = {-0,-200};
// detuning range in units of gamma;

const int PRINT_STEPBYSTEP = 0;
// print simulation results step by step (test and debug) [0,1]
///---------------------------------------------------------------------------------------------------------------///
///------------------------------------------------ RESULTS ------------------------------------------------------///
#define NumVoxels 200
double PositionHistogram[NumVoxels][NumVoxels][NumVoxels];
double xbins[NumVoxels];
double ybins[NumVoxels];
double zbins[NumVoxels];

#define NumVelocityBins 5000
#define NumTimeWindows_VelHist 50
const double TimeStepLength_VelHist = 10e-3;
// Length of time windows on velocity evolution histogram
const double MaxVelocityBin = 5.0;
// Last value for VelocityBins [m/s]
double VelocityHistogram[NumTimeWindows_VelHist][NumVelocityBins];
double VelocityBins[NumVelocityBins];
///---------------------------------------------------------------------------------------------------------------///














/// END OF USER DEFINED PARAMETERS -------------------------------------------------------------------------------///
///---------------------------------------------------------------------------------------------------------------///
///---------------------------------------------------------------------------------------------------------------///
///---------------------------------------------------------------------------------------------------------------///
///---------------------------------------------------------------------------------------------------------------///
///---------------------------------------------------------------------------------------------------------------///


// Definition of structs for return of functions
struct BeamSelection{
    int BEAMindex;
    double CycleTime;
    int dmj;
};
struct SimulationResults{
    int FLAG;
    double TrapTime;
    int Iterations;
};
struct CollisionMomentum{
    double dvxc;
    double dvyc;
    double dvzc;
};





//Calculate magnetic fields given x,y,z
//Returns field strength B and modifies Bn[] (unit vector of field direction)
double MagneticField (double x,double y, double z,double Bn[]){

    double Bx = GBx*x + B0[0];
    double By = GBy*y + B0[1];
    double Bz = GBz*z + B0[2];
    double B;
    double r0,r1,r2;

    B = sqrt(Bx*Bx + By*By + Bz*Bz);

    if(B==0){
        r0 = (randr()-0.5)*2;
        r1 = (randr()-0.5)*2;
        r2 = (randr()-0.5)*2;
        Bn[0] = r0 / sqrt(r0*r0 + r1*r1 + r2*r2);
        Bn[1] = r1 / sqrt(r0*r0 + r1*r1 + r2*r2);
        Bn[2] = r2 / sqrt(r0*r0 + r1*r1 + r2*r2);
        return B;
    }

    Bn[0]=Bx/B;
    Bn[1]=By/B;
    Bn[2]=Bz/B;

    return B;
}

//Calculate polarization vectors on the lab frame
void CalculatePolarizationVectors(void){
    int i;
    int ERROR = 0;
    double k[3],ex[3],ey[3];
    double r[3];
    double exek;

    char p;
        // 'r': RHCP
        // 'l': LHCP
        // 'x': polarization vector on plane defined by X and ek (beam propagation vector)
        // 'y': [...]
        // 'z': [...]


    for (i=0;i<NUMBEAMS;i++){
        k[0] = BEAMS[i][0];
        k[1] = BEAMS[i][1];
        k[2] = BEAMS[i][2];
        p = BEAMS_POL[i];

        if(p == 'r'){
            // random ex not parallel to k
            do{
                r[0] = (randr()-0.5)*2;
                r[1] = (randr()-0.5)*2;
                r[2] = (randr()-0.5)*2;
                ex[0] = r[0] / sqrt(dotproduct(r,r));
                ex[1] = r[1] / sqrt(dotproduct(r,r));
                ex[2] = r[2] / sqrt(dotproduct(r,r));
                exek = dotproduct(ex,k);
            }while(abs(exek) == 1);
            // orthogonalization
            r[0] = ex[0] - exek*k[0];
            r[1] = ex[1] - exek*k[1];
            r[2] = ex[2] - exek*k[2];
            // renormalization
            ex[0] = r[0] / sqrt(dotproduct(r,r));
            ex[1] = r[1] / sqrt(dotproduct(r,r));
            ex[2] = r[2] / sqrt(dotproduct(r,r));
            // ey = k(*)ex
            ey[0] = k[1]*ex[2] - k[2]*ex[1];
            ey[1] = k[2]*ex[0] - k[0]*ex[2];
            ey[2] = k[0]*ex[1] - k[1]*ex[0];
            // RHCP
            BEAMS_POL_VECTORS[i][0] = sqrt(0.5)*(ex[0] - I*ey[0]);
            BEAMS_POL_VECTORS[i][1] = sqrt(0.5)*(ex[1] - I*ey[1]);
            BEAMS_POL_VECTORS[i][2] = sqrt(0.5)*(ex[2] - I*ey[2]);
        }
        else if(p == 'l'){
            // random ex not parallel to k
            do{
                r[0] = (randr()-0.5)*2;
                r[1] = (randr()-0.5)*2;
                r[2] = (randr()-0.5)*2;
                ex[0] = r[0] / sqrt(dotproduct(r,r));
                ex[1] = r[1] / sqrt(dotproduct(r,r));
                ex[2] = r[2] / sqrt(dotproduct(r,r));
                exek = ex[0]*k[0] + ex[1]*k[1] + ex[2]*k[2];
            }while(exek == 1);
            // orthogonalization
            r[0] = ex[0] - exek*k[0];
            r[1] = ex[1] - exek*k[1];
            r[2] = ex[2] - exek*k[2];
            // renormalization
            ex[0] = r[0] / sqrt(dotproduct(r,r));
            ex[1] = r[1] / sqrt(dotproduct(r,r));
            ex[2] = r[2] / sqrt(dotproduct(r,r));
            // ey = k(*)ex
            ey[0] = k[1]*ex[2] - k[2]*ex[1];
            ey[1] = k[2]*ex[0] - k[0]*ex[2];
            ey[2] = k[0]*ex[1] - k[1]*ex[0];
            // LHCP
            BEAMS_POL_VECTORS[i][0] = sqrt(0.5)*(ex[0] + I*ey[0]);
            BEAMS_POL_VECTORS[i][1] = sqrt(0.5)*(ex[1] + I*ey[1]);
            BEAMS_POL_VECTORS[i][2] = sqrt(0.5)*(ex[2] + I*ey[2]);
        }
        else if(p == 'x' || p == 'y' || p == 'z'){
            r[0]=0;
            r[1]=0;
            r[2]=0;
            if(p=='x') r[0] = 1;
            if(p=='y') r[1] = 1;
            if(p=='z') r[2] = 1;

            if(abs(dotproduct(r,k))==1){
                ERROR = 1;
                printf("\n\n\nERROR: UNCONSISTENT POLARIZATION [BEAM #%d]\n\n",i);
                printf("Polarization: %c\n",p);
                printf("Beam propagation: [%+.2e %+.2e %+.2e]\n\n",k[0],k[1],k[2]);
                printf("Polarization vector cannot be parallel to beam propagation direction!\n\n\n\n");
                printf("Press any key to exit");
                getchar();
                printf("\n\n");
            }
            r[0] = r[0] - dotproduct(r,k)*k[0];
            r[1] = r[1] - dotproduct(r,k)*k[1];
            r[2] = r[2] - dotproduct(r,k)*k[2];
            BEAMS_POL_VECTORS[i][0] = r[0]/sqrt(dotproduct(r,r));
            BEAMS_POL_VECTORS[i][1] = r[1]/sqrt(dotproduct(r,r));
            BEAMS_POL_VECTORS[i][2] = r[2]/sqrt(dotproduct(r,r));

        }
        else{
            ERROR = 1;
            printf("\n\n\nERROR: POLARIZATION DEFINITION NOT RECOGNIZED [BEAM #%d]\n\n",i);
            printf("Polarization: %c\n\n",p);
            printf("Polarization must be defined as one of the following options:\n['r','l','x','y','z']\n\n\n\n");
            printf("Press any key to exit");
            getchar();
            printf("\n\n");
        }
    }
    if(ERROR) exit(-1);
}

// Calculate transition vectors given magnetic field direction
void CalculateTransitionVectors (double Bn[],double complex polsm[],double complex polsp[],double complex polpi[]){
    int i;
    double r[3],ex[3],ey[3];
    double exek;

    // polpi = Bn
    for (i=0;i<3;i++)
        polpi[i]=Bn[i];

    // random ex not parallel to k
    do{
        for(i=0;i<3;i++)
            r[i] = (randr()-0.5)*2;
        for(i=0;i<3;i++)
            ex[i] = r[i] / sqrt(dotproduct(r,r));
        exek = dotproduct(ex,Bn);
    } while(abs(exek)==1);
    // orthogonalization
    for(i=0;i<3;i++)
        r[i] = ex[i] - exek*Bn[i];
    // renormalization
    for(i=0;i<3;i++)
        ex[i] = r[i] / sqrt(dotproduct(r,r));
    // ey = k(*)ex
    ey[0] = Bn[1]*ex[2] - Bn[2]*ex[1];
    ey[1] = Bn[2]*ex[0] - Bn[0]*ex[2];
    ey[2] = Bn[0]*ex[1] - Bn[1]*ex[0];
    // SIGMA-
    polsm[0] = sqrt(0.5)*(ex[0] - I*ey[0]);
    polsm[1] = sqrt(0.5)*(ex[1] - I*ey[1]);
    polsm[2] = sqrt(0.5)*(ex[2] - I*ey[2]);
    // SIGMA+
    polsp[0] = sqrt(0.5)*(ex[0] + I*ey[0]);
    polsp[1] = sqrt(0.5)*(ex[1] + I*ey[1]);
    polsp[2] = sqrt(0.5)*(ex[2] + I*ey[2]);
}

// Selects beam to be absorbed given x,y,z,vx,vy,vz;
// Returns beam index, time step and transition type (delta mj)
struct BeamSelection ChooseBeam (double x,double y,double z,double vx, double vy, double vz){
    int BEAMindex=0,i=0,j=0,k=0,iT=0;
    double B=0.,Bn[]={0.,0.,0.},BEAM[3];
    double dZeeman=0.,mje=0.,mjg=0.,sop=0.,delta=0.,Rscatt=0.,doppler;
    double complex polsm[3];
    double complex polsp[3];
    double complex polpi[3];
    double complex TransitionVectors[3][3]; //TransitionVectors[transition index][x,y,z]
    double complex trpolvec[3];             //trpolvec[x,y,z]
    double complex bPOL[3];
    double cycletime[NUMBEAMS][3];          //TransitionVectors[beam index][transition index]
    double soeff,soeffx,deltax,cycletimex,trpolvec_beampol;

//    double DEBUG1[NUMBEAMS][3];
//    double DEBUG2[NUMBEAMS][3];
//    double DEBUG3[NUMBEAMS][3];


    B = MagneticField(x,y,z,Bn);
    CalculateTransitionVectors(Bn,polsm,polsp,polpi);


    for (i=0;i<3;i++){
        TransitionVectors[0][i] = polsm[i];
        TransitionVectors[1][i] = polpi[i];
        TransitionVectors[2][i] = polsp[i];
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // SET GROUND STATE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    mjg = -Jgnd; // <---------------------- [#REVIEW]

    for (i=0;i<NUMBEAMS;i++){
        //BEAM
        for (j=0;j<3;j++){
            BEAM[j]=BEAMS[i][j];
            bPOL[j]=BEAMS_POL_VECTORS[i][j];
        }

        //Gaussian beam;
        sop = so[i]*exp(-2*(pow(BEAM[1]*z - BEAM[2]*y,2) + pow(BEAM[2]*x - BEAM[0]*z,2) + pow(BEAM[0]*y - BEAM[1]*x,2))/pow(BEAMS_WAIST[i],2));

        // Loop over possible transitions
        mje=mjg-1;
        for (iT=0;iT<3;iT++){
            if(fabs(mjg)>Jexc)
                cycletime[i][iT] = 1./0.;   // State does not exist
            else{
                // State does exist
                for(j=0;j<3;j++)
                    trpolvec[j] = TransitionVectors[iT][j];
                trpolvec_beampol = pow(AbsDotProductComplex(bPOL,trpolvec),2);//////// ESCREVER AQUI
                dZeeman = mub*B*(gjg*mjg - gje*mje )/h;
                doppler = -(c/lambda)*(vx*BEAM[0] + vy*BEAM[1] + vz*BEAM[2])/c;
                delta = detuning*Gamma + dZeeman + doppler;
                soeff = sop*trpolvec_beampol;

                cycletime[i][iT] = 1./0.;
                for (k=-NUM_SIDEBANDS;k<=NUM_SIDEBANDS;k++){
                    soeffx = soeff/(2*NUM_SIDEBANDS+1);
                    deltax = delta + k*FreqStep_SIDEBANDS;
                    Rscatt = (Gamma/2)*soeffx/( 1+soeffx+(4*(deltax*deltax)/(Gamma*Gamma)) );
                    cycletimex = (1/Rscatt)*RandomExpDist();
                    if(cycletimex < cycletime[i][iT])
                        cycletime[i][iT] = cycletimex;
                }


                //cycletime[i][iT] = (1/Rscatt)*RandomExpDist();

//                DEBUG1[i][iT] = Rscatt;
//                DEBUG2[i][iT] = (Gamma/2)*soeff/( 1+soeff+(4*((delta-doppler)*(delta-doppler))/(Gamma*Gamma)) );;
//                DEBUG3[i][iT] = delta/Gamma;

            }
            mje++;
        }

    }


    int dmj[]={-1,0,1};
    int xdmj = 0;
    double dt = 1/0.;

    for (i=0;i<NUMBEAMS;i++){
        for (iT=0;iT<3;iT++){
            if(cycletime[i][iT]<dt){
                dt = cycletime[i][iT];
                xdmj = dmj[iT];
                BEAMindex = i;
            }
        }
    }



    /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ///
    /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ///
    ///  #DEBUG
//    if(PRINT_STEPBYSTEP){
//                double AUX;
//                static int numA=0,numB=0;
//                printf("\n\n\n");
//                printf("\n\n\nr = %+.2e  %+.2e  %+.2e\n",x,y,z);
//                printf("v = %+.2e  %+.2e  %+.2e\n\n",vx,vy,vz);
//                for (i=0;i<NUMBEAMS;i++){
//                    for (iT=0;iT<3;iT++){
//                        //printf("%+.2e  ",cycletime[i][iT]);
//                        printf("%+.8e  ",DEBUG1[i][iT]);
//                    }
//                    if(i==BEAMindex){
//                        printf("(dmj = %+d)  ",xdmj);
//                        AUX = BEAMS[i][0]*x + BEAMS[i][1]*y + BEAMS[i][2]*z;
//                        if(AUX>0){
//                            printf("[+]");
//                            numA++;
//                        }
//                        if(AUX<0){
//                            printf("[-]");
//                            numB++;
//                        }
//                        if(AUX==0)
//                            printf("[0]");
//                    }
//                    printf("\n");
//                    for (iT=0;iT<3;iT++){
//                        //printf("%+.2e  ",cycletime[i][iT]);
//                        printf("%+.8e  ",DEBUG2[i][iT]);
//                    }
//                    printf("\n");
//                    for (iT=0;iT<3;iT++){
//                        //printf("%+.2e  ",cycletime[i][iT]);
//                        printf("%+.8e  ",DEBUG3[i][iT]);
//                    }
//                    printf("\n\n");
//                }
//                printf("\n\n[%d\t\t%d]\n",numA,numB);
//                getchar();
//    }
    /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ///
    /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ///


    struct BeamSelection result;
    result.BEAMindex = BEAMindex;
    result.CycleTime = dt;
    result.dmj = xdmj;
    return result;
}

// Calculate collision momentum transfer assuming collision of the atom with another atom with a temperature T
struct CollisionMomentum CalculateCollisions(double x, double y, double z, double vx, double vy, double vz, double dt){
    // Canonical ensemble
    struct CollisionMomentum pc;
//    double dv[3];
//    double v, tau;
//
//    v = sqrt(vx*vx + vy*vy + vz*vz);
//    //tau = (0.40/v)*RandomExpDist();
//    tau = (5e-6/v)*RandomExpDist();
//
//    if(tau<=dt){
//        dv[0] = gaussian()*sqrt(kb*T/m);
//        dv[1] = gaussian()*sqrt(kb*T/m);
//        dv[2] = gaussian()*sqrt(kb*T/m);
//    }
//    else{
//        dv[0] = 0;
//        dv[1] = 0;
//        dv[2] = 0;
//    }
//
//    pc.dvxc=dv[0];
//    pc.dvyc=dv[1];
//    pc.dvzc=dv[2];
    pc.dvxc=0.0;
    pc.dvyc=0.0;
    pc.dvzc=0.0;

    return pc;
};

//Simulates path of one atom with defined configuration
struct SimulationResults RunSimulation (){
    struct SimulationResults SimRes;
    if(abs(Jexc - Jgnd)!=1){
        printf("\n\n\t\tError!\n\n\t\tJexc != Jgnd +- 1 \n\n");
        SimRes.FLAG = -1;
        SimRes.Iterations=0;
        SimRes.TrapTime=0;
        exit(-1);
    }

    int i,BEAMindex;
    double x,y,z,vx,vy,vz,t;
    struct BeamSelection BEAMresult;
    struct CollisionMomentum pCollision_Struct;
    double r0,r1,r2,pxe,pye,pze,pxa,pya,pza;
    int ix,iy,iz,iv,itwvh; //Histogram indexes
    double TRAPTIME=0;
    double dvxc,dvyc,dvzc;
    int dmj=0;

    dmj=dmj; // TEMP SÓ PRA SUMIR O WARNING ENQUANTO NÃO USO DE VERDADE


    double dvx,dvy,dvz,dt;
    double Amagx,Amagy,Amagz;

    // Initial conditions
    t=0;
    vx = vx0;
    vy = vy0;
    vz = vz0;
    x = xx0;
    y = yy0;
    z = zz0;

    int FLAG = 0;

    for (i=0;i<NUMITER && FLAG==0;i++){
        BEAMresult = ChooseBeam (x,y,z,vx,vy,vz);
        BEAMindex = BEAMresult.BEAMindex;
        dt = BEAMresult.CycleTime;
        dmj = BEAMresult.dmj;
        //B = MagneticField(x,y,z,Bn);

        // Emission in random direction
        r0 = (randr()-0.5)*2;
        r1 = (randr()-0.5)*2;
        r2 = (randr()-0.5)*2;
        pxe = (h/lambda)*r0/sqrt(r0*r0 + r1*r1 + r2*r2);
        pye = (h/lambda)*r1/sqrt(r0*r0 + r1*r1 + r2*r2);
        pze = (h/lambda)*r2/sqrt(r0*r0 + r1*r1 + r2*r2);
        // Absorption
        pxa = (h/lambda)*BEAMS[BEAMindex][0];
        pya = (h/lambda)*BEAMS[BEAMindex][1];
        pza = (h/lambda)*BEAMS[BEAMindex][2];
        // Collisions
        pCollision_Struct = CalculateCollisions(x,y,z,vx,vy,vz,dt);
        dvxc=pCollision_Struct.dvxc;
        dvyc=pCollision_Struct.dvyc;
        dvzc=pCollision_Struct.dvzc;

        dvx = (pxe + pxa)/m + dvxc;
        dvy = (pye + pya)/m + dvyc;
        dvz = (pze + pza)/m + dvzc;

        Amagx = -gjg*Jgnd*mub*(GBx*GBx*x)/sqrt(pow(GBx*x,2) + pow(GBy*y,2) + pow(GBz*z,2));
        Amagy = -gjg*Jgnd*mub*(GBy*GBy*y)/sqrt(pow(GBx*x,2) + pow(GBy*y,2) + pow(GBz*z,2));
        Amagz = -gjg*Jgnd*mub*(GBz*GBz*z)/sqrt(pow(GBx*x,2) + pow(GBy*y,2) + pow(GBz*z,2));

         x += vx*dt + (gx/2)*dt*dt + (Amagx/2)*dt*dt;
         y += vy*dt + (gy/2)*dt*dt + (Amagy/2)*dt*dt;
         z += vz*dt + (gz/2)*dt*dt + (Amagz/2)*dt*dt;
        vx += gx*dt + dvx + Amagx*dt;
        vy += gy*dt + dvy + Amagy*dt;
        vz += gz*dt + dvz + Amagz*dt;
         t += dt;



         if((vx*x + vy*y + vz*z)<0)
            TRAPTIME = t;



        if(sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc))>BOUNDARY)
            FLAG = 1;
        else{
            ix = round((x+BOUNDARY)*(NumVoxels-1)/(2*BOUNDARY));
            iy = round((y+BOUNDARY)*(NumVoxels-1)/(2*BOUNDARY));
            iz = round((z+BOUNDARY)*(NumVoxels-1)/(2*BOUNDARY));
            if(ix<=NumVoxels-1 || iz<=NumVoxels-1 || iy<=NumVoxels-1)
                PositionHistogram[ix][iy][iz] += dt;
            iv = round((NumVelocityBins-1)*sqrt(vx*vx + vy*vy + vz*vz)/VelocityBins[NumVelocityBins-1]);
            itwvh = floor(t/TimeStepLength_VelHist);
            if(itwvh > NumTimeWindows_VelHist-1)
                itwvh = NumTimeWindows_VelHist-1;
            if(iv<=NumVelocityBins-1)
                VelocityHistogram[itwvh][iv] += dt;

//            printf("\n\n### %f %f %d ###\n\n",sqrt(vx*vx + vy*vy + vz*vz),VelocityBins[NumVelocityBins-1],iv);
//            getchar();
        }
    }

    SimRes.FLAG = FLAG;
    SimRes.TrapTime = TRAPTIME;
    SimRes.Iterations = i;
    return SimRes;
}

// Write output file header
void WriteResultsHeader(FILE *fid){
    int i;
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);

    fprintf(fid,"%s;\n",TITLE);
    fprintf(fid,"Data starts after exclamation mark;\n");
    fprintf(fid,"%d_%02d_%02d %02d:%02d:%02d;\n",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    fprintf(fid,"Gravity = (%f,%f,%f) m/s2;\n",gx,gy,gz);
    fprintf(fid,"Magnetic field gradient = (%e,%e,%e) T/m;\n",GBx,GBy,GBz);
    fprintf(fid,"Magnetic field bias = (%e %e %e) T;\n",B0[0],B0[1],B0[2]);
    fprintf(fid,"Detuning = %f (units of gamma);\n",detuning);
    fprintf(fid,"Number of beams = %d;\n",NUMBEAMS);

    fprintf(fid,"BEAMS: propagation direction = \n");
    fprintf(fid,"((%+f,%+f,%+f)",BEAMS[0][0],BEAMS[0][1],BEAMS[0][2]);
    for (i=1;i<NUMBEAMS;i++)
        fprintf(fid,",\n(%+f,%+f,%+f)",BEAMS[i][0],BEAMS[i][1],BEAMS[i][2]);
    fprintf(fid,");\n");

    fprintf(fid,"BEAMS: polarization definition = (%c",BEAMS_POL[0]);
    for (i=1;i<NUMBEAMS;i++)
        fprintf(fid,",%c",BEAMS_POL[i]);
    fprintf(fid,");\n");

    fprintf(fid,"BEAMS: polarization vectors = \n");
    fprintf(fid,"((%+f%+fi,%+f%+fi,%+f%+fi)",creal(BEAMS_POL_VECTORS[0][0]),cimag(BEAMS_POL_VECTORS[0][0]),creal(BEAMS_POL_VECTORS[0][1]),cimag(BEAMS_POL_VECTORS[0][1]),creal(BEAMS_POL_VECTORS[0][2]),cimag(BEAMS_POL_VECTORS[0][2]));;
    for (i=1;i<NUMBEAMS;i++)
        fprintf(fid,",\n(%+f%+fi,%+f%+fi,%+f%+fi)",creal(BEAMS_POL_VECTORS[i][0]),cimag(BEAMS_POL_VECTORS[i][0]),creal(BEAMS_POL_VECTORS[i][1]),cimag(BEAMS_POL_VECTORS[i][1]),creal(BEAMS_POL_VECTORS[i][2]),cimag(BEAMS_POL_VECTORS[i][2]));;
    fprintf(fid,");\n");

    fprintf(fid,"BEAMS: peak intensity/saturation intensity = (%f",so[0]);
    for (i=1;i<NUMBEAMS;i++)
        fprintf(fid,",%f",so[i]);
    fprintf(fid,");\n");
    fprintf(fid,"BEAMS: 1/e2 waist = (%e",BEAMS_WAIST[0]);
    for (i=1;i<NUMBEAMS;i++)
        fprintf(fid,",%e",BEAMS_WAIST[i]);
    fprintf(fid,") m;\n");
    fprintf(fid,"ATOM: mass = %e kg;\n",m);
    fprintf(fid,"ATOM: transition gamma = %e Hz;\n",Gamma);
    fprintf(fid,"ATOM: transition wavelength = %e m;\n",lambda);
    fprintf(fid,"ATOM: J (ground) = %d;\n",(int)Jgnd);
    fprintf(fid,"ATOM: J (excited) = %d;\n",(int)Jexc);
    fprintf(fid,"ATOM: g_lande (ground) = %f;\n",gjg);
    fprintf(fid,"ATOM: g_lande (excited) = %f;\n",gje);
    fprintf(fid,"SIMULATION: maximum number of iterations = %d;\n",NUMITER);
    fprintf(fid,"SIMULATION: boundary radius = %e m;\n",BOUNDARY);
    fprintf(fid,"SIMULATION: boundary center = (%e,%e,%e) m;\n",xc,yc,zc);
    fprintf(fid,"SIMULATION: initial position = (%e,%e,%e) m;\n",xx0,yy0,zz0);
    fprintf(fid,"SIMULATION: initial velocity = (%e,%e,%e) m/s;\n",vx0,vy0,vz0);
    fprintf(fid,"<<<DATA>>>!\n");
}


int TRAPPING_TIME_ROUTINE (void){
    srand(time(NULL));
    struct SimulationResults SimRes;
    int i,j,k;
    int ix,iy,iz;
    double PosHist2Dxy[NumVoxels][NumVoxels],PosHist2Dxz[NumVoxels][NumVoxels],PosHist2Dyz[NumVoxels][NumVoxels];

    // Magnetic field check
    if(GBx+GBy+GBz!=0){
        printf("Magnetic field divergence != 0\n\n GBx + GBy + GBz !=0\n\n");
        getchar();
        exit(-1);
    }

    //------------------------------------------------------- Files --------------------------------------------------------//
    time_t datetime;
    struct tm tm;
    double clock1,clock2;
    datetime = time(NULL);
    tm = *localtime(&datetime);
    clock1 = clock();

    char OUTPUTFOLDER[200];
    char TRAPTIMEFILE[400];
    char INFOFILE[400];

    sprintf(OUTPUTFOLDER,"./Results/%s_%04d%02d%02d_%02d%02d%02d",TITLE,tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    sprintf(TRAPTIMEFILE,"%s/Trapping_Time.dat",OUTPUTFOLDER);
    sprintf(INFOFILE,"%s/info.txt",OUTPUTFOLDER);

    mkdir(OUTPUTFOLDER);

    FILE *fLOG;
    fLOG = fopen("outputlog.txt","a");
    // log file
    FILE *fTT;
    fTT = fopen(TRAPTIMEFILE,"w");
    // traptime file

    if(fTT == NULL){
        printf("\n\n\n\n FAILED TO CREATE:\n%s \n\n\n\n\n",TRAPTIMEFILE);
        printf("Press <enter> to continue");
        getchar();
        printf("\n\n\n\n\n\n\n\n");
    }





    fprintf(fLOG,"\n\n\n\n\n\n\n\n");
    for (k=0;k<100;k++)
        fprintf(fLOG,"#");
    fprintf(fLOG,"\n<<<Simulation started (%d_%02d_%02d %02d:%02d:%02d)>>>\n",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    //----------------------------------------------------------------------------------------------------------------------//


    // Generate xbins,ybins,zbins
    for (k=0;k<NumVoxels;k++){
        xbins[k] = -BOUNDARY + (2*BOUNDARY)*k/NumVoxels;
        ybins[k] = -BOUNDARY + (2*BOUNDARY)*k/NumVoxels;
        zbins[k] = -BOUNDARY + (2*BOUNDARY)*k/NumVoxels;
    }
    // Normalize BEAMS
    double S;
    for (i=0;i<NUMBEAMS;i++){
        S=0;
        for (j=0;j<3;j++)
            S += BEAMS[i][j]*BEAMS[i][j];
        S = sqrt(S);
        for (j=0;j<3;j++)
            BEAMS[i][j] = BEAMS[i][j]/S;
    }


    CalculatePolarizationVectors();


    FILE *fINFO;
    fINFO = fopen(INFOFILE,"w");
    WriteResultsHeader(fINFO);
    fclose(fINFO);

    //*
    double detuningVar[NumVar];
    linspace(DetuningRange[0],DetuningRange[1],NumVar,detuningVar);

    // Execute simulation multiple times for averaging and varying detuning
    for (j=0;j<NumVar;j++){
        detuning = detuningVar[j];

        fprintf(fTT,"%e\t",detuning);

        // Zero position histograms (every
        for (ix=0;ix<NumVoxels;ix++){
            for (iy=0;iy<NumVoxels;iy++){
                PosHist2Dxy[ix][iy]=0;
                PosHist2Dxz[ix][iy]=0;
                PosHist2Dyz[ix][iy]=0;
                for (iz=0;iz<NumVoxels;iz++){
                    PositionHistogram[ix][iy][iz]=0;
                }
            }
        }

        printf("\n Detuning = %f\n\n",detuning);
        fprintf(fLOG,"\n Detuning = %f (units of Gamma)\n\n",detuning);


        for (i=0;i<NumAvg;i++){
            //Random starting position and velocity
            xx0 = xx0c + gaussian()*STDxx0;
            yy0 = yy0c + gaussian()*STDyy0;
            zz0 = zz0c + gaussian()*STDzz0;
            vx0 = vx0c + gaussian()*sqrt(kb*T/m);
            vy0 = vy0c + gaussian()*sqrt(kb*T/m);
            vz0 = vz0c + gaussian()*sqrt(kb*T/m);


            SimRes = RunSimulation();
                  printf("[%-5d/%d var; %5d/%d avg]  (%d)  %-8d  %.3e s\n",j+1,NumVar,i+1,NumAvg,SimRes.FLAG,SimRes.Iterations,SimRes.TrapTime);
            fprintf(fLOG,"[%-5d/%d var; %5d/%d avg]  (%d)  %-8d  %.3e s\n",j+1,NumVar,i+1,NumAvg,SimRes.FLAG,SimRes.Iterations,SimRes.TrapTime);
            fprintf(fTT,"%e\t",SimRes.TrapTime);
        }
        fprintf(fTT,"\n");

        for (ix=0;ix<NumVoxels;ix++){
            for (iy=0;iy<NumVoxels;iy++){
                for (iz=0;iz<NumVoxels;iz++){
                    PosHist2Dxy[ix][iy] += PositionHistogram[ix][iy][iz];
                    PosHist2Dxz[ix][iz] += PositionHistogram[ix][iy][iz];
                    PosHist2Dyz[iy][iz] += PositionHistogram[ix][iy][iz];
                }
            }
        }


        ///-----------------------------------------------------------------------------------------------------------------------///
        /// ---------- 2D HISTOGRAMS FILES ---------------------------------------------------------------------------------------///
            FILE *fid;
            char filename[200];
            sprintf(filename,"%s/det%f.dat",OUTPUTFOLDER,detuning);
            fid = fopen(filename,"w");
            // 2D histograms files

            if(fid == NULL){
                printf("\n\n\n\n FAILED TO CREATE:\n%s \n\n\n\n\n",filename);
                printf("Press <enter> to continue");
                getchar();
                printf("\n\n\n\n\n\n\n\n");
            }

            WriteResultsHeader(fid);
            fprintf(fid,"%e\t",.0);

            for (ix=0;ix<NumVoxels;ix++){
                fprintf(fid,"%e\t",xbins[ix]);
            }
            fprintf(fid,"\n");
            for (iy=0;iy<NumVoxels;iy++){
                fprintf(fid,"%e\t",ybins[iy]);
                for (ix=0;ix<NumVoxels;ix++){
                        fprintf(fid,"%e\t",PosHist2Dxz[ix][iy]);
                }
                fprintf(fid,"\n");
            }
            fflush(fid);
            fclose(fid);
        ///-----------------------------------------------------------------------------------------------------------------------///
        ///-----------------------------------------------------------------------------------------------------------------------///



    } //var

    datetime = time(NULL);
    tm = *localtime(&datetime);
    clock2 = clock();
    fprintf(fLOG,"\n\nElapsed time: %f s\n",(clock2 - clock1)/CLOCKS_PER_SEC);
    fprintf(fLOG,"<<<Simulation ended (%d_%02d_%02d %02d:%02d:%02d)>>>\n",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    for (k=0;k<100;k++)
        fprintf(fLOG,"#");
    fprintf(fLOG,"\n\n\n\n\n\n\n\n");
    fflush(fLOG);
    fclose(fLOG);

    //*/

    return 0;
}






int TRAP_DEPTH_ROUTINE (void){
    srand(time(NULL));
    struct SimulationResults SimRes;
    int i,j,k;

    // Magnetic field check
    if(GBx+GBy+GBz!=0){
        printf("Magnetic field divergence != 0\n\n GBx + GBy + GBz !=0\n\n");
        printf("\n%f\n%f\n%f",GBx,GBy,GBz);
        getchar();
        exit(-1);
    }

    //------------------------------------------------------- Files --------------------------------------------------------//
    time_t datetime;
    struct tm tm;
    datetime = time(NULL);
    tm = *localtime(&datetime);

    char OUTPUTFOLDER[200];
    char TRAPDEPTHFILE[400];
    char COPYMAINCOM[400];
    char INFOFILE[400];

    sprintf(OUTPUTFOLDER,"./Results/%s_%04d%02d%02d_%02d%02d%02d",TITLE,tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    sprintf(TRAPDEPTHFILE,"%s/TrapDepth.dat",OUTPUTFOLDER);
    sprintf(COPYMAINCOM,"copy .\\main.c .\\Results\\%s_%04d%02d%02d_%02d%02d%02d\\ ",TITLE,tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    sprintf(INFOFILE,"%s/info.txt",OUTPUTFOLDER);
    mkdir(OUTPUTFOLDER);
    system(COPYMAINCOM);
    printf("\n\n>>> %s\n\n",COPYMAINCOM);



    FILE *fLOG;
    fLOG = fopen("outputlog.txt","a");
    // log file
    FILE *fTD;
    fTD = fopen(TRAPDEPTHFILE,"w");
    // traptime file

    if(fTD == NULL){
        printf("\n\n\n\n FAILED TO CREATE:\n%s \n\n\n\n\n",TRAPDEPTHFILE);
        printf("Press <enter> to continue");
        getchar();
        printf("\n\n\n\n\n\n\n\n");
    }





    fprintf(fLOG,"\n\n\n\n\n\n\n\n");
    for (k=0;k<100;k++)
        fprintf(fLOG,"#");
    fprintf(fLOG,"\n<<<Simulation started (%d_%02d_%02d %02d:%02d:%02d)>>>\n",tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    //----------------------------------------------------------------------------------------------------------------------//


    // Generate xbins,ybins,zbins
    for (k=0;k<NumVoxels;k++){
        xbins[k] = -BOUNDARY + (2*BOUNDARY)*k/NumVoxels;
        ybins[k] = -BOUNDARY + (2*BOUNDARY)*k/NumVoxels;
        zbins[k] = -BOUNDARY + (2*BOUNDARY)*k/NumVoxels;
    }
    // Normalize BEAMS
    double S;
    for (i=0;i<NUMBEAMS;i++){
        S=0;
        for (j=0;j<3;j++)
            S += BEAMS[i][j]*BEAMS[i][j];
        S = sqrt(S);
        for (j=0;j<3;j++)
            BEAMS[i][j] = BEAMS[i][j]/S;
    }


    CalculatePolarizationVectors();

    FILE *fINFO;
    fINFO = fopen(INFOFILE,"w");
    WriteResultsHeader(fINFO);
    fclose(fINFO);

    //*/
    int iVar,iVel,iAvg;
    double r[3];
    double velocity[NumVel];

    linspace(0,MaxVel_TrapDepth,NumVel,velocity);

    double detuningVar[NumVar];
    linspace(DetuningRange[0],DetuningRange[1],NumVar,detuningVar);

    int Ntrp;


    fprintf(fTD,"%.8e\t",0.0);
    for (iVel=0;iVel<NumVel;iVel++)
        fprintf(fTD,"%.8e\t",velocity[iVel]);
    fprintf(fTD,"\n");


    for (iVar=0;iVar<NumVar;iVar++){
        detuning = detuningVar[iVar];
        printf("Detuning = %10.2f\n",detuning);
        fprintf(fTD,"%.8e\t",detuning);
        for (iVel=0;iVel<NumVel;iVel++){
            printf("    vel = %.2e m/s\n",velocity[iVel]);
            Ntrp = 0;
            for (iAvg=0;iAvg<NumAvg;iAvg++){
                RandomUnitVector(r);
                vx0 = velocity[iVel]*r[0];
                vy0 = velocity[iVel]*r[1];
                vz0 = velocity[iVel]*r[2];
                xx0 = xx0c + gaussian()*STDxx0;
                yy0 = yy0c + gaussian()*STDyy0;
                zz0 = zz0c + gaussian()*STDzz0;
                SimRes = RunSimulation();
                if(SimRes.FLAG == 0)
                    Ntrp++;
                printf("[%5d/%d var; %5d/%d vel; %5d/%d avg]  (%d)  %-8d  %.3e s\n",iVar+1,NumVar,iVel+1,NumVel,iAvg+1,NumAvg,SimRes.FLAG,SimRes.Iterations,SimRes.TrapTime);
                //fprintf(fLOG,"[%5d/%d var; %5d/%d vel; %5d/%d avg]  (%d)  %-8d  %.3e s\n",iVar+1,NumVar,iVel+1,NumVel,iAvg+1,NumAvg,SimRes.FLAG,SimRes.Iterations,SimRes.TrapTime);
            }
            //getchar();
            printf("\n\n%e\n\n",((double)Ntrp)/((double)NumAvg));
            fprintf(fTD,"%.8e\t",((double)Ntrp)/((double)NumAvg));
        }
        fprintf(fTD,"\n");
    }
    //*/
    return 0;
}

int VELOCITY_HISTOGRAM_EVOLUTION_ROUTINE (void){
    ///------------------------------///
    ///------------------------------///
    ///------------------------------///
    // Check magnetic field gradients
    if(GBx+GBy+GBz!=0){
        printf("Magnetic field divergence != 0\n\n GBx + GBy + GBz !=0\n\n");
        getchar();
        exit(-1);
    }
    //Random generator seed
    srand(time(NULL));
    // Normalize BEAMS
    double S;
    int i,j;
    for (i=0;i<NUMBEAMS;i++){
        S=0;
        for (j=0;j<3;j++)
            S += BEAMS[i][j]*BEAMS[i][j];
        S = sqrt(S);
        for (j=0;j<3;j++)
            BEAMS[i][j] = BEAMS[i][j]/S;
    }
    //Polarization Vectors on lab frame
    CalculatePolarizationVectors();
    //Velocity histogram = zeros
    for(j=0;j<NumTimeWindows_VelHist;j++)
        for (i=0;i<NumVelocityBins;i++)
            VelocityHistogram[j][i]=0;
    //Create velocity histogram bins
    linspace(0,MaxVelocityBin,NumVelocityBins,VelocityBins);
    // Files
    time_t datetime;
    struct tm tm;
    datetime = time(NULL);
    tm = *localtime(&datetime);
    char OUTPUTFOLDER[200];
    char VELHISTFILE[400];
    sprintf(OUTPUTFOLDER,"./Results/%s_%04d%02d%02d_%02d%02d%02d",TITLE,tm.tm_year+1900,tm.tm_mon+1,tm.tm_mday,tm.tm_hour,tm.tm_min,tm.tm_sec);
    sprintf(VELHISTFILE,"%s/VelocityHistogram.dat",OUTPUTFOLDER);
    mkdir(OUTPUTFOLDER);
    //
    ///------------------------------///
    ///------------------------------///
    ///------------------------------///

    struct SimulationResults SimRes;

    for (i=0;i<NumAvg;i++){
        xx0 = xx0c + gaussian()*STDxx0;
        yy0 = yy0c + gaussian()*STDyy0;
        zz0 = zz0c + gaussian()*STDzz0;
        vx0 = vx0c + gaussian()*sqrt(kb*T/m);
        vy0 = vy0c + gaussian()*sqrt(kb*T/m);
        vz0 = vz0c + gaussian()*sqrt(kb*T/m);

        SimRes = RunSimulation();
        printf("[%3d/%d]  (%d)  %-8d  %.3e s\n",i+1,NumAvg,SimRes.FLAG,SimRes.Iterations,SimRes.TrapTime);
    }

    FILE *fout;
    fout = fopen(VELHISTFILE,"w");
    fprintf(fout,"%e",0.0);
    for (j=0;j<NumTimeWindows_VelHist;j++)
        fprintf(fout,"\t%e",TimeStepLength_VelHist*j);
    for (i=0;i<NumVelocityBins;i++){
        fprintf(fout,"\n%e",VelocityBins[i]);
        for (j=0;j<NumTimeWindows_VelHist;j++){
            fprintf(fout,"\t%e",VelocityHistogram[j][i]);
        }
    }
    fflush(fout);
    fclose(fout);
//
    return 0;
}

int GEOMETRY_OPTIMIZATION (void){

    double Beam_zComponent;
    double beam_angle;
    for (beam_angle=2.5;beam_angle<=25.0;beam_angle+=2.5)
    {
        sprintf(TITLE,"GEOMOPTIM_FINE_3B_ANGLE_%f",beam_angle);
        Beam_zComponent = tan(beam_angle*pi/180);

        BEAMS[0][0] = -1.0;
        BEAMS[0][1] = +0.0;
        BEAMS[0][2] = Beam_zComponent;

        BEAMS[1][0] = +0.5;
        BEAMS[1][1] = -0.5*sqrt(3);
        BEAMS[1][2] = Beam_zComponent;

        BEAMS[2][0] = +0.5;
        BEAMS[2][1] = +0.5*sqrt(3);
        BEAMS[2][2] = Beam_zComponent;

        printf("\n\n\n BEAM ANGLE = %f degrees\n\n\n",beam_angle);
        TRAP_DEPTH_ROUTINE();
    }


    return 0;
}

int main (void){
    //TRAPPING_TIME_ROUTINE();
    //TRAP_DEPTH_ROUTINE();
    //VELOCITY_HISTOGRAM_EVOLUTION_ROUTINE();
    GEOMETRY_OPTIMIZATION();

    return 0;
}
