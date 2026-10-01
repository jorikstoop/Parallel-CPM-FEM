// def.h
#ifndef _DEF
#define _DEF

//#define NULL 0
#define FALSE 0
#define TRUE 1
typedef int BOOL;
#define SEED 3

#define NVX 100        // domain size (number of voxels)
#define NVY NVX
#define NV  (NVX*NVY)
#define NNX (NVX+1)     // number of nodes
#define NNY (NVY+1)
#define NN  (NNX*NNY)
#define NDOF (2*NN)      // number of degrees of freedom
#define VOXSIZE .0000025 // [m]
#define NRINC 100      // 100 for profiling, 3000 for defualt sim

#define MAXNRITER 1000
#define ACCURACY .00001

// material properties
#define YOUNGS 32000 // [Pa]
#define POISSON .45//

// loading
#define LOAD 0 //3E6
#define FORCE (LOAD*VOXSIZE)

// cells
#define IMMOTILITY 1.0//50

#define CELLDIAM .000020 // cell diameter [m]
#define CELLRVOX (CELLDIAM/2/VOXSIZE) // cell radius [pixels]
#define TARGETVOLUME  (3.1415*CELLRVOX*CELLRVOX) // targetvolume [pixels]
#define INELASTICITY 500.0// [-] 1.0E20 // [/m4]

#define NOSTICKJ 500000.0 //10000// [/m] contact penalty for none-adhesive surface
#define JCM (NOSTICKJ*VOXSIZE)  // cell-medium
#define JCC (2.0*JCM) // cell-cell

#define MAXDHSTR 10.0 // unscaled at the moment,
#define THRESHOLDSTIFF 15E3    // threshold stiffness for durotaxis
#define STIFFSENSITIVITY .0005 // steepness of durotaxis sigmoid
#define STIFFENINGSTIFF .1 // steepness of strain-stiffening

#define THICKNESS 10E-6 // effective thickness of substrate
#define CELLFORCE (1.0E-5/THICKNESS) // [N]

#define SQ05 .707107 //sqrt(.5), used often enough to make this convenient


#endif
