#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#pragma warning(disable:4996)  // get rid of some annoying warnings

//----- PROGRAM CONSTANTS ---------------------------------------------------------------
#define X         0                     // x-direction
#define Y         1                     // y-direction
#define MAXLINE   1024                  // string buffer size
#define PI        (4.0*atan(1.0))
#define DEG2RAD   (PI/180.0)
#define TOL       (1.0e-8)              // determine if two points are same location

// headers for data file parsing
#define JOINT_COORDINATE_HEADER      "[JOINT COORDINATES]"
#define MEMBER_CONNECTIVITY_HEADER   "[MEMBER JOINT CONNECTIVITY]"
#define REACTIONS_HEADER             "[REACTIONS AT NODES]"
#define EXTERNAL_FORCES_HEADER       "[EXTERNAL FORCES]"
#define FORCE_UNITS_HEADER           "[FORCE UNITS]"

//----- STRUCTURE DEFINITIONS -----------------------------------------------------------

typedef struct JOINT
{
	double p[2];  // joint x,y position.  p[0]=x, p[1]=y
}
JOINT;

typedef struct MEMBER
{
   int j[2];      // joint indexes
   double F, L;   // internal Force (+ve is tension), member Length
}
MEMBER;



typedef struct FORCE
{
   double F;      // Force

   int j;         // joint index (reactions at nodes)
   double theta;  // if force is given as angle from positive x axis (reaction at node angles)

   int idir;      // if true y else x


   FORCE()  // initialize with bad values
   {
      F=1.0e30;
      j=-1, idir=-1;
      theta=1.0e30;
   };
}
FORCE;

typedef struct TRUSS_SIZES
{
   int NJ, NM, NR, NEF;  // number of Joints, Members, Reactions, External Forces
}
TRUSS_SIZES;

//----- FUNCTION PROTOTYPES -------------------------------------------------------------

// parse the file to get number of joints,members,reactions,forces
TRUSS_SIZES getTrussSizes(FILE *f);

// gets the truss data from a file
bool getTrussData(FILE *f, TRUSS_SIZES, JOINT*, MEMBER*, FORCE*, FORCE*, char*);

// free up memory and close file
void cleanup(FILE *f, JOINT*, MEMBER*, FORCE*, FORCE*, double**, double**, double*, double*, int);

bool isEmpty(char *);  // checks for empty line in file

// function to invert an NxN matrix
bool InverseNxN(double **A, int N, double **Ainv);
// helper function for InverseNxN
bool ELGS(double **A, int N, int *indx);

