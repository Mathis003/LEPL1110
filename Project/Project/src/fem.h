#ifndef _FEM_H_
#define _FEM_H_

/*
* Note:
* This header file is used to define the functions and structures used in the project.
* The functions are defined in the fem.c file.
* These functions are used in the Pronect/ folder but also in the ProjectPreProcessor and the ProjectPostProcessor/ folder.
*/

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
typedef enum {FEM_DISCRETE_TYPE_LINEAR, FEM_DISCRETE_TYPE_QUADRATIC} femDiscreteType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } femElasticCase;
typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE } femElementType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM,FEM_RCMK} femRenumType;
typedef enum {FEM_FULL,FEM_BAND} femSolverType;

/* Structures */

typedef struct {
    int nNodes;
    int *number;
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
    double * (*getMaterialProperties)(char *material); // Not used yet
    char * (*getMaterials)(double x, double y); // Not used yet

    // For the example only (beam or U form)
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
    femConstrainedNode *constrainedNodes;
} femProblem;

typedef struct Queue
{
	int capacity;
	int size;
	int front;
	int rear;
	int *elements;
} Queue;

/* Variables */

extern double *positionMeshNodes;
extern femGeometry theGeometry;

/* Functions declaration */

/**********/
/* Solver */
/**********/

femFullSystem *femFullSystemCreate(int size);
void femFullSystemFree(femFullSystem *mySystem);
void femFullSystemAlloc(femFullSystem *mySystem, int size);
void femFullSystemInit(femFullSystem *mySystem, int size);
double femFullSystemGetA_Entry(femFullSystem *mySystem, int myRow, int myCol);
double femFullSystemGetB_Entry(femFullSystem *mySystem, int myRow);
void femFullSystemAssemble(femFullSystem *system, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, double xLoc, int nLoc);
void femFullSystemConstrainXY(femFullSystem *mySystem, int myNode, double value, int size);
void femFullSystemConstrainNT(femFullSystem *theSystem, double vect1_x, double vect1_y, double vect2_x, double vect2_y, double myValue, int Ux, int Uy, int size);
double *femFullSystemEliminate(femFullSystem *mySystem, int size);
void femFullSystemGetResidual(femFullSystem *mySystem, int size, double *residuals, double *theSoluce);
void femFullSystemPrint(femFullSystem *mySystem, int size);
void femFullSystemPrintInfos(femFullSystem *mySystem, int size);

femBandSystem *femBandSystemCreate(int size, int band);
void femBandSystemFree(femBandSystem *myBandSystem);
void femBandSystemInit(femBandSystem *myBandSystem, int size);
void femBandSystemAlloc(femBandSystem *system, int size, int band);
int isInBand(int band, int myRow, int myCol);
double femBandSystemGetA_Entry(femBandSystem *mySystem, int myRow, int myCol);
double femBandSystemGetB_Entry(femBandSystem *mySystem, int myRow);
void femBandSystemAssemble(femBandSystem *system, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, double xLoc, int nLoc);
void femBandSystemConstrainXY(femBandSystem *system, int node, double value, int size);
void femBandSystemConstrainNT(femBandSystem *theSystem, double vect1_x, double vect1_y, double vect2_x, double vect2_y, double myValue, int Ux, int Uy, int size);
double *femBandSystemEliminate(femBandSystem *myBand, int size);
void femBandSystemGetResidual(femBandSystem *mySystem, int size, double *residuals, double *theSoluce);
void femBandSystemPrint(femBandSystem *myBand, int size);
void femBandSystemPrintInfos(femBandSystem *myBand, int size);
void femMeshRenumber(femMesh *theMesh, femRenumType renumType);
int femMeshComputeBand(femMesh *theMesh);

femSolver *femSolverCreate(int size);
femSolver *femSolverFullCreate(int size);
femSolver *femSolverBandCreate(int size, int band);
void femSolverFree(femSolver *mySolver);
void femSolverInit(femSolver *mySolver);
double femSolverGetA_Entry(femSolver *mySolver, int myRow, int myCol);
double femSolverGetB_Entry(femSolver *mySolver, int myRow);
double **femSolverGetA(femSolver *mySolver);
double *femSolverGetB(femSolver *mySolver);
void femSolverSet(femSolver *mySolver, double **newA, double *newB);
void femSolverGetResidual(femSolver *mySolver, double *residuals, double *theSoluce);
void femSolverPrint(femSolver *mySolver);
void femSolverPrintInfos(femSolver *mySolver);
void femSolverAssemble(femSolver *mySolver, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, double xLoc, int nLoc);
void femSolverSystemConstrainXY(femSolver *mySolver, int node, double value);
void femSolverSystemConstrainNT(femSolver *mySolver, double vect1_x, double vect1_y, double vect2_x, double vect2_y, double myValue, int Ux, int Uy);
double *femSolverEliminate(femSolver *mySolver);

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

void _q1c0_x_linear(double *xsi, double *eta);
void _q1c0_phi_linear(double xsi, double eta, double *phi);
void _q1c0_dphidx_linear(double xsi, double eta, double *dphidxsi, double *dphideta);
void _p1c0_x_linear(double *xsi, double *eta);
void _p1c0_phi_linear(double xsi, double eta, double *phi);
void _p1c0_dphidx_linear(double xsi, double eta, double *dphidxsi, double *dphideta);
void _e1c0_x_linear(double *xsi);
void _e1c0_phi_linear(double xsi,  double *phi);
void _e1c0_dphidx_linear(double xsi, double *dphidxsi);

void _q1c0_x_quadratic(double *xsi, double *eta);
void _q1c0_phi_quadratic(double xsi, double eta, double *phi);
void _q1c0_dphidx_quadratic(double xsi, double eta, double *dphidxsi, double *dphideta);
void _p1c0_x_quadratic(double *xsi, double *eta);
void _p1c0_phi_quadratic(double xsi, double eta, double *phi);
void _p1c0_dphidx_quadratic(double xsi, double eta, double *dphidxsi, double *dphideta);
void _e1c0_x_quadratic(double *xsi);
void _e1c0_phi_quadratic(double xsi,  double *phi);
void _e1c0_dphidx_quadratic(double xsi, double *dphidxsi);

femDiscrete *femDiscreteCreate(femElementType type, femDiscreteType dType);
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
void geoMeshWrite(const char *filename, femDiscreteType dType);
void geoMeshRead(const char *filename, femDiscreteType dType);
void geoSetDomainName(int iDomain, char *name);
int geoGetDomain(char *name);

/************************/
/* Elasticity functions */
/************************/

femProblem *femElasticityCreate(femGeometry *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase, femDiscreteType dType);
void femElasticityFree(femProblem *theProblem);
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2);
void femElasticityAssembleElements(femProblem *theProblem);
void femElasticityApplyDirichlet(femProblem *theProblem, double FACTOR);
double *femElasticitySolve(femProblem *theProblem, femRenumType renumType, double FACTOR, int currAnim);
double *femElasticityForces(femProblem *theProblem);

void femElasticityPrint(femProblem *theProblem);
femProblem *femElasticityRead(femGeometry *theGeometry, femSolverType typeSolver, const char *filename, femRenumType renumType, femDiscreteType dType);
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

/***************/
/* Renumbering */
/***************/

int comparPositionNode(const void *a, const void *b);
void femMeshRenumber(femMesh *theMesh, femRenumType renumType);

void swap(int *a, int *b);
void reverse_array(int *array, int n);

Queue *createQueue(int maxCapacity);
void enqueue(Queue *Q, int element);
void dequeue(Queue *Q);
int peek(Queue *Q);
int isEmpty(Queue *Q);
int isFull(Queue *Q);

int partition(int arr1[], int arr2[], int low, int high);
void quickSort(int arr1[], int arr2[], int low, int high);

int *createAdjacencyMatrix(femMesh *mesh);
void add_neighbors_to_queue(int *adj, int n, int *degrees, int *inserted, Queue *Q, int element_idx);
Queue *rcm(femMesh *theMesh, int nNodes);

/***********************/
/* Animation functions */
/***********************/

double adaptForceForMotionCar(double force, double node_position, int currAnim, int nTotalAnim);
double adaptForceForMotionCarReversed(double force, double node_position, int currAnim, int nTotalAnim);

#endif // _FEM_H_