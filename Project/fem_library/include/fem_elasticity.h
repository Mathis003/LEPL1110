#ifndef _FEM_ELASTICITY_H_
#define _FEM_ELASTICITY_H_

#include <math.h>
#include <string.h>

#include "fem_solver.h"
#include "fem_integrate.h"

typedef enum { PLANAR_STRESS, PLANAR_STRAIN, AXISYM } femElasticCase;

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
    femFullSystem *system;
    femConstrainedNode *constrainedNodes;
} femProblem;

femProblem *femElasticityCreate(femGeometry *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase);
void femElasticityFree(femProblem *theProblem);
void femElasticityPrint(femProblem *theProblem);
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2);
double *femElasticitySolve(femProblem *theProblem);

femProblem *femElasticityRead(femGeometry *theGeometry, const char *filename);
void femElasticityWrite(femProblem *theProbconst, const char *filename);

void femSolutionWrite(int nNodes, int nfields, double *data, const char *filename);
int femSolutiondRead(int allocated_size, double *value, const char *filename);

void geoSetDomainName(int iDomain, char *name);
int geoGetDomain(char *name);

#endif // _FEM_ELASTICITY_H_