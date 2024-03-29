/*
*  fem.c
*  Library for LEPL1110 : Finite Elements for dummies
*
*  Copyright (C) 2023 UCL-IMMC : Vincent Legat
*  All rights reserved.
*
*/

#ifndef _FEM_H_
#define _FEM_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ErrorScan(a) femErrorScan(a, __LINE__, __FILE__)
#define ErrorGmsh(a) femErrorGmsh(a, __LINE__, __FILE__)
#define Error(a)     femError(a, __LINE__, __FILE__)
#define Warning(a)   femWarning(a, __LINE__, __FILE__)
#define FALSE 0
#define TRUE 1
#define MAXNAME 256

typedef enum { FEM_TRIANGLE, FEM_QUAD, FEM_EDGE } femElementType;
typedef enum { DIRICHLET_X, DIRICHLET_Y, DIRICHLET_XY, DIRICHLET_N, DIRICHLET_T, DIRICHLET_NT, NEUMANN_X, NEUMANN_Y, NEUMANN_N, NEUMANN_T, UNDEFINED=-1} femBoundaryType;
typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } femElasticCase;

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
    double LxPlate, LyPlate;
    double h;
    femElementType elementType;
    double (*geoSize)(double x, double y);
    femNodes *theNodes;
    femMesh *theElements;
    femMesh *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;

typedef struct {
    int n;
    femElementType type;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x)(double *xsi);
    void (*phi)(double xsi, double *phi);
    void (*dphidx)(double xsi, double *dphidxsi);
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
    double E, nu, rho, gx, gy;
    double A, B, C;
    int planarStrainStress;
    int nBoundaryConditions;
    femBoundaryCondition **conditions;
    double *soluce;
    double *residuals;
    femGeo *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femDiscrete *spaceEdge;
    femIntegration *ruleEdge;
    femFullSystem *system;
    femConstrainedNode *constrainedNodes;
} femProblem;

void geoInitialize(void);
femGeo *geoGetGeometry(void);
double geoSize(double x, double y);
double geoSizeDefault(double x, double y);
void geoSetSizeCallback(double (*geoSize)(double x, double y));
void geoMeshGenerate(void);
void geoMeshGenerateGeo(void);
void geoMeshGenerateGeoFile(const char *filename);
void geoMeshGenerateMshFile(const char *filename);
void geoMeshImport(void);
void geoMeshPrint(void);
void geoMeshWrite(const char *filename);
void geoMeshRead(const char *filename);
void geoSetDomainName(int iDomain, char *name);
int geoGetDomain(char *name);
void geoFinalize(void);
void geoFree(void);

femProblem *femElasticityCreate(femGeo *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase);
void femElasticityFree(femProblem *theProblem);
void femElasticityPrint(femProblem *theProblem);
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2);
double *femElasticitySolve(femProblem *theProblem);
void femElasticityWrite(femProblem *theProbconst, const char *filename);
femProblem *femElasticityRead(femGeo *theGeometry, const char *filename);

void femSolutionWrite(int nNodes, int nfields, double *data, const char *filename);
int femSolutiondRead(int allocated_size, double *value, const char *filename);

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

femFullSystem *femFullSystemCreate(int size);
void femFullSystemFree(femFullSystem *mySystem);
void femFullSystemPrint(femFullSystem *mySystem);
void femFullSystemInit(femFullSystem *mySystem);
void femFullSystemAlloc(femFullSystem *mySystem, int size);
double *femFullSystemEliminate(femFullSystem *mySystem);
void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double value);

double femMin(double *x, int n);
double femMax(double *x, int n);
void femError(char *text, int line, char *file);
void femErrorScan(int test, int line, char *file);
void femErrorGmsh(int test, int line, char *file);
void femWarning(char *text, int line, char *file);

#endif