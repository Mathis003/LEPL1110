#ifndef _FEM_OTHER_H_
#define _FEM_OTHER_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define ErrorScan(a) femErrorScan(a, __LINE__, __FILE__)
#define Error(a)     femError(a, __LINE__, __FILE__)
#define Warning(a)   femWarning(a, __LINE__, __FILE__)

#define TRUE 1
#define FALSE 0

#define MAXNAME 256

typedef enum { DIRICHLET_X, DIRICHLET_Y, DIRICHLET_XY, DIRICHLET_N, DIRICHLET_T, DIRICHLET_NT, NEUMANN_X, NEUMANN_Y, NEUMANN_N, NEUMANN_T, UNDEFINED=-1} femBoundaryType;
typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE } femElementType;

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

double femMin(double *x, int n);
double femMax(double *x, int n);
void femError(char *text, int line, char *file);
void femErrorScan(int test, int line, char *file);
void femWarning(char *text, int line, char *file);

#endif // _FEM_OTHER_H_