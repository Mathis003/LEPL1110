#ifndef _FEM_H_
#define _FEM_H_

/* Import */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>

/* Preprocessor directive */

#define ErrorScan(a) femErrorScan(a, __LINE__, __FILE__)
#define Error(a)     femError(a, __LINE__, __FILE__)
#define Warning(a)   femWarning(a, __LINE__, __FILE__)
#define ErrorGmsh(a) femErrorGmsh(a, __LINE__, __FILE__)

#define TRUE 1
#define FALSE 0
#define MAXNAME 256

/* Enumerations */

typedef enum { DIRICHLET_X, DIRICHLET_Y, DIRICHLET_XY, DIRICHLET_N, DIRICHLET_T, DIRICHLET_NT, NEUMANN_X, NEUMANN_Y, NEUMANN_N, NEUMANN_T, UNDEFINED=-1} femBoundaryType;
typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE } femElementType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } femElasticCase;
typedef enum {FEM_FULL,FEM_BAND} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

/* Structures */

typedef struct {
    int nNodes;
    double *X;
    double *Y;
} femNodes;

typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;

typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

typedef struct {
    femDomain *domain;
    femBoundaryType type;
    double value1;
    double value2;
} femBoundaryCondition;

typedef struct {
    femBoundaryType type;
    double nx, ny;
    double value1, value2;
} femConstrainedNode;

typedef struct {
    double widthSpanBridge, heightBridge;
    double widthWindow, heightWindow;
    double widthSubRoadWay, heightSubRoadWay;
    double rxArc, ryArc;
    double rxLongArc, ryLongArc;
    double widthPiles;
    double widthPillars, heightPillars;
    double widthPylons, heightPylons;
    double widthStayCables, heightStayCables, distStayCables, angleStayCables;

    double defaultSize;
    
    double (*geoSize)(double x, double y);
    double * (*getMaterialProperties)(char *material);
    char * (*getMaterials)(double x, double y);

    // For the example
    double LxPlate, LyPlate;

    femElementType elementType;
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeometry;

typedef struct {
    int n;
    femElementType type;
    void (*x2) (double *xsi, double *eta);
    void (*phi2) (double xsi, double eta, double *phi);
    void (*dphi2dx) (double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x) (double *xsi);
    void (*phi) (double xsi, double *phi);
    void (*dphidx) (double xsi, double *dphidxsi);
} femDiscrete;

typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femSystemCopy;

typedef struct {
    double *B;
    double **A;
} femFullSystem;

typedef struct {
    double *B;
    double **A;        
    int band;        
} femBandSystem;

typedef struct {
    femSolverType type;
    void *solver;
    int size;
} femSolver;

typedef struct {
    double E, nu, rho, gx, gy;
    double A, B, C;
    int planarStrainStress;
    int nBoundaryConditions;
    femBoundaryCondition **conditions;
    double *soluce;
    double *residuals;
    femGeometry *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femDiscrete *spaceEdge;
    femIntegration *ruleEdge;
    femSolver *solver;
    femSystemCopy *copySystem;
    femConstrainedNode *constrainedNodes;
} femProblem;

/* Variables */

extern double *positionMeshNodes; // To renumber the mesh nodes
extern femGeometry theGeometry;

/* Functions declaration */

/**********/
/* Solver */
/**********/

femFullSystem *femFullSystemCreate(int size);
void femFullSystemFree(femFullSystem *mySystem);
void femFullSystemAlloc(femFullSystem *mySystem, int size);
void femFullSystemInit(femFullSystem *mySystem, int size);
void femFullSystemAssemble(femFullSystem *system, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, int xLoc, int nLoc);
void femFullSystemConstrainXY(femFullSystem *mySystem, int size, int myNode, double value);
void femFullSystemConstrainNT(femFullSystem *system, int size, int node1, int node2, double a, double b);
double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol);
double *femFullSystemEliminate(femFullSystem *mySystem, int size);
void femFullSystemPrint(femFullSystem *mySystem, int size);
void femFullSystemPrintInfos(femFullSystem *mySystem, int size);

femBandSystem *femBandSystemCreate(int size, int band);
void femBandSystemFree(femBandSystem *myBandSystem);
void femBandSystemInit(femBandSystem *myBandSystem, int size);
void femBandSystemAlloc(femBandSystem *system, int size, int band);
int comparPosNode(const void *a, const void *b);
void femMeshRenumber(femMesh *theMesh, femRenumType renumType);
int femMeshComputeBand(femMesh *theMesh);
void femBandSystemAssemble(femBandSystem *system, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, int xLoc, int nLoc);
double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol);
void femBandSystemConstrainXY(femBandSystem *system, int size, int node, double value); 
void femBandSystemConstrainNT(femBandSystem *system, int size, int node1, int node2, double a, double b);
double *femBandSystemEliminate(femBandSystem *myBand, int size);
void femBandSystemPrint(femBandSystem *myBand, int size);
void femBandSystemPrintInfos(femBandSystem *myBand, int size);
int comparPositionNode(const void *a, const void *b);
void femMeshRenumber(femMesh *theMesh, femRenumType renumType);
int femMeshComputeBand(femMesh *theMesh);

femSolver *femSolverCreate(int sizeLoc);
femSolver *femSolverFullCreate(int size, int sizeLoc);
femSolver *femSolverBandCreate(int size, int sizeLoc, int band);
femSolver *femSolverIterativeCreate(int size, int sizeLoc);
void femSolverFree(femSolver *mySolver);
void femSolverInit(femSolver *mySolver);
double femSolverGet(femSolver *mySolver, int i, int j);
void femSolverPrint(femSolver *mySolver);
void femSolverPrintInfos(femSolver *mySolver);
void femSolverAssemble(femSolver *mySolver, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, int xLoc, int nLoc);
void femSolverSystemConstrainXY(femSolver *mySolver, int node, double value);
void femSolverSystemConstrainNT(femSolver *mySolver, int node1, int node2, double a, double b);
double *femSolverEliminate(femSolver *mySolver);
int femSolverConverged(femSolver *mySolver);
double **getMatrixA(femSolver *mySolver);
double *getVectorB(femSolver *mySolver);
double getSizeMatrix(femSolver *mySolver);

/**************************************/
/* Discrete space + Integration space */
/**************************************/

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = {0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = {0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = {0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = {0.166666666666667, 0.166666666666667, 0.166666666666667};
static const double _gaussEdge2Xsi[2]    = { 0.577350269189626,-0.577350269189626};
static const double _gaussEdge2Weight[2] = { 1.000000000000000, 1.000000000000000};

femIntegration *femIntegrationCreate(int n, femElementType type);
void femIntegrationFree(femIntegration *theRule);

femDiscrete *femDiscreteCreate(int n, femElementType type);
void femDiscreteFree(femDiscrete *mySpace);
void femDiscretePrint(femDiscrete *mySpace);
void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta);
void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi);
void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void femDiscreteXsi(femDiscrete *mySpace, double *xsi);
void femDiscretePhi(femDiscrete *mySpace, double xsi, double *phi);
void femDiscreteDphi(femDiscrete *mySpace, double xsi, double *dphidxsi);

/************/
/* Geometry */
/************/

femGeometry *geoGetGeometry(void);
double geoSize(double x, double y);
void geoFree(void);
void geoSetSizeCallback(double (*geoSize)(double x, double y));
void geoMeshPrint(void);
void geoMeshWrite(const char *filename);
void geoMeshRead(const char *filename);
void geoSetDomainName(int iDomain, char *name);
int geoGetDomain(char *name);

/************************/
/* Elasticity functions */
/************************/

femProblem *femElasticityCreate(femGeometry *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase);
void femElasticityFree(femProblem *theProblem);
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2);
double *femElasticitySolve(femProblem *theProblem);
double *femElasticityForces(femProblem *theProblem);

void femElasticityPrint(femProblem *theProblem);
femProblem *femElasticityRead(femGeometry *theGeometry, femSolverType typeSolver, const char *filename);
void femElasticityWrite(femProblem *theProbconst, const char *filename);
void femSolutionWrite(int nNodes, int nfields, double *data, const char *filename);
int femSolutiondRead(int allocated_size, double *value, const char *filename);

/***********/
/* General */
/***********/

double femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y));
double femMin(double *x, int n);
double femMax(double *x, int n);
void femError(char *text, int line, char *file);
void femErrorScan(int test, int line, char *file);
void femWarning(char *text, int line, char *file);

#endif // _FEM_H_