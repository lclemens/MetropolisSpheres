/*** Allard Lab jun.allard@uci.edu                    ***/

#define TWISTER genrand_real3()
#define NTMAX           10
#define NMAX            30
#define NTADAPT         20000
#define NTCHECK         200000
#define DCHIMIN         1e-4
#define PI              3.14159265359
#define INF             1e14
#define DCHIINIT        0.1
#define KSCRITICAL      0.01
#define MEMBRANE        1
#define MULTIPLE        1
#define CPMAX           1e8
#define TALKATIVE       1
#define VISUALIZE       1

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "twister.c"


/*******************************************************************************/
//  GLOBAL VARIABLES
/*******************************************************************************/

/* General Global Variables */
char listName[100];
FILE *fList;
//
char paramsFilename[100], filamentFilename[100], iSiteFilename[100], bSiteFilename[100], basicSiteFilename[100];
FILE *paramsFile, *filList, *iSiteList, *bSiteList, *basicSiteList;

long NumberiSites;
long Ncurrent;
double c0, c1, irLigand;

int i1,i2;

long iseed;

double contourLength;
double rSphere[NMAX][3],rSpherePropose[NMAX][3];
double rAnchor[NMAX][3];
double sRadius;
int NSphere;
double norm;

long st;
long proposals[2], accepts[2], nt, iChi, i, iPropose, iChiPropose, ix, iParam, ntNextStationarityCheck, iStart;
double E, ENew, rate[2], dChi[2], dChiHere, Force;
long constraintProposalsTotal;

int filamentInputMethod, iSiteInputMethod;
long commandiSites;
char *iSiteLocations;
char input[4*NMAX];
long j,m;
int verboseTF;

/* Convergence Global Variables */
int convergedTF, constraintSatisfiedTF;

/* MULTIPLE Global Variables*/
double brLigand;

long ib, ib2;

/* MULTIPLE FILAMENT Variables*/
double baseSepDistance;

double distCurrent;

/* Filament tail dimerization force */
double dimerDistCurrent, dimerDist0;
double kdimer;

/* MULTIPLE ligands energy variables */
double kBound;
double boundCentertoJointDistance, boundCentertoBaseDistance, boundCentertoBaseLigandDistance, boundCentertoBoundDistance;

int endConstraintPassedTF;

double kSphere;
/*******************************************************************************/
//  INCLUDES
/*******************************************************************************/

#include "outputControl.c"
#include "getParameters.c"
#include "metropolisJoint.c"


/*******************************************************************************/
//  MAIN
/*******************************************************************************/

// arguments: listName,
int main( int argc, char *argv[] )
{

    if(argv[1]) //get name of parameters file
        strcpy(paramsFilename,argv[1]);
    if (TALKATIVE) printf("This is the parameter filename: %s\n", paramsFilename);

    // use text file to set parameters
    getParameters();
    
    if(argv[2]) //get name of parameters file
        strcpy(listName,argv[2]);
    if (TALKATIVE) printf("This is the output filename: %s\n", listName);
    
    if(argv[3]) //get size of spheres
        sRadius = atof(argv[3]);
    if (TALKATIVE) printf("This is the sphere radius: %f\n", sRadius);


    
    
    // initialize random seed
	iseed = RanInitReturnIseed(0);
    if (TALKATIVE) printf("iseed: %ld\n",iseed);
	
    // run metropolis algorithm
	metropolisJoint();

	return 0;
	
} // finished main
