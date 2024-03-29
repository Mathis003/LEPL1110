/*
*  fem.c
*  Library for LEPL1110 : Finite Elements for dummies
*
*  Copyright (C) 2021 UCL-IMMC : Vincent Legat
*  All rights reserved.
*
*/

#include "fem.h"

femGeo theGeometry;

femGeo *geoGetGeometry(void) { return &theGeometry; }

double geoSizeDefault(double x, double y) { return theGeometry.h; }

void geoFree(void)
{
    if (theGeometry.theNodes)
    {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes);
    }
    if (theGeometry.theElements)
    {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements);
    }
    if (theGeometry.theEdges)
    {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges);
    }
    for (int i = 0; i < theGeometry.nDomains; i++)
    {
        free(theGeometry.theDomains[i]->elem);
        free(theGeometry.theDomains[i]);
    }
    free(theGeometry.theDomains);
}

void geoSetSizeCallback(double (*geoSize)(double x, double y)) { theGeometry.geoSize = geoSize; }

void geoMeshPrint(void)
{
    femNodes *theNodes = theGeometry.theNodes;
    if (theNodes != NULL)
    {
        printf("Number of nodes %d \n", theNodes->nNodes);
        for (int i = 0; i < theNodes->nNodes; i++) { printf("%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]); }
    }
    femMesh *theEdges = theGeometry.theEdges;
    if (theEdges != NULL)
    {
        printf("Number of edges %d \n", theEdges->nElem);
        int *elem = theEdges->elem;
        for (int i = 0; i < theEdges->nElem; i++) { printf("%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]); }
    }
    femMesh *theElements = theGeometry.theElements;
    if (theElements != NULL)
    {
        if (theElements->nLocalNode == 3)
        {
            printf("Number of triangles %d \n", theElements->nElem);
            int *elem = theElements->elem;
            for (int i = 0; i < theElements->nElem; i++) { printf("%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]); }
        }
        else if (theElements->nLocalNode == 4)
        {
            printf("Number of quads %d \n", theElements->nElem);
            int *elem = theElements->elem;
            for (int i = 0; i < theElements->nElem; i++) { printf("%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]); }
        }
    }
    int nDomains = theGeometry.nDomains;
    printf("Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        printf("  Domain : %6d \n", iDomain);
        printf("  Name : %s\n", theDomain->name);
        printf("  Number of elements : %6d\n", theDomain->nElem);
        for (int i = 0; i < theDomain->nElem; i++)
        {
            // if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
            printf("%6d", theDomain->elem[i]);
            if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0) { printf("\n"); }
        }
        printf("\n");
    }
}

void geoMeshWrite(const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }

    femNodes *theNodes = theGeometry.theNodes;
    fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
    for (int i = 0; i < theNodes->nNodes; i++) { fprintf(file, "%6d : %14.7e %14.7e \n", i, theNodes->X[i], theNodes->Y[i]); }

    femMesh *theEdges = theGeometry.theEdges;
    fprintf(file, "Number of edges %d \n", theEdges->nElem);
    int *elem = theEdges->elem;
    for (int i = 0; i < theEdges->nElem; i++) { fprintf(file, "%6d : %6d %6d \n", i, elem[2 * i], elem[2 * i + 1]); }

    femMesh *theElements = theGeometry.theElements;
    if (theElements->nLocalNode == 3)
    {
        fprintf(file, "Number of triangles %d \n", theElements->nElem);
        elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) { fprintf(file, "%6d : %6d %6d %6d\n", i, elem[3 * i], elem[3 * i + 1], elem[3 * i + 2]); }
    }
    else if (theElements->nLocalNode == 4)
    {
        fprintf(file, "Number of quads %d \n", theElements->nElem);
        elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) { fprintf(file, "%6d : %6d %6d %6d %6d\n", i, elem[4 * i], elem[4 * i + 1], elem[4 * i + 2], elem[4 * i + 3]); }
    }

    int nDomains = theGeometry.nDomains;
    fprintf(file, "Number of domains %d\n", nDomains);
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        fprintf(file, "  Domain : %6d \n", iDomain);
        fprintf(file, "  Name : %s\n", theDomain->name);
        fprintf(file, "  Number of elements : %6d\n", theDomain->nElem);
        for (int i = 0; i < theDomain->nElem; i++)
        {
            fprintf(file, "%6d", theDomain->elem[i]);
            if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0)  { fprintf(file, "\n"); }
            fprintf(file, "\n");
        }
    }

    fclose(file);
}

void geoMeshRead(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }

    int trash, *elem;

    femNodes *theNodes = malloc(sizeof(femNodes));
    if (theNodes == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theGeometry.theNodes = theNodes;
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
    theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->Y == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < theNodes->nNodes; i++) { ErrorScan(fscanf(file, "%d : %le %le \n", &trash, &theNodes->X[i], &theNodes->Y[i])); }

    femMesh *theEdges = malloc(sizeof(femMesh));
    if (theEdges == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theGeometry.theEdges = theEdges;
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
    theEdges->elem = malloc(sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
    if (theEdges->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < theEdges->nElem; ++i)
    {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash, &elem[2 * i], &elem[2 * i + 1]));
    }

    femMesh *theElements = malloc(sizeof(femMesh));
    if (theElements == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theGeometry.theElements = theElements;
    theElements->nLocalNode = 0;
    theElements->nodes = theNodes;
    char elementType[MAXNAME];
    ErrorScan(fscanf(file, "Number of %s %d \n", elementType, &theElements->nElem));
    if (strncasecmp(elementType, "triangles", MAXNAME) == 0)
    {
        theElements->nLocalNode = 3;
        theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
        if (theElements->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theElements->nElem; ++i)
        {
            elem = theElements->elem;
            ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", &trash, &elem[3 * i], &elem[3 * i + 1], &elem[3 * i + 2]));
        }
    }
    if (strncasecmp(elementType, "quads", MAXNAME) == 0)
    {
        theElements->nLocalNode = 4;
        theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
        if (theElements->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theElements->nElem; ++i)
        {
            elem = theElements->elem;
            ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", &trash, &elem[4 * i], &elem[4 * i + 1], &elem[4 * i + 2], &elem[4 * i + 3]));
        }
    }

    ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
    int nDomains = theGeometry.nDomains;
    theGeometry.theDomains = malloc(sizeof(femDomain *) * nDomains);
    if (theGeometry.theDomains == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = malloc(sizeof(femDomain));
        if (theDomain == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theGeometry.theDomains[iDomain] = theDomain;
        theDomain->mesh = theEdges;
        ErrorScan(fscanf(file, "  Domain : %6d \n", &trash));
        ErrorScan(fscanf(file, "  Name : %[^\n]s \n", (char *)&theDomain->name));
        ErrorScan(fscanf(file, "  Number of elements : %6d\n", &theDomain->nElem));
        theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
        if (theDomain->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theDomain->nElem; i++)
        {
            ErrorScan(fscanf(file, "%6d", &theDomain->elem[i]));
            if ((i + 1) != theDomain->nElem && (i + 1) % 10 == 0) { ErrorScan(fscanf(file, "\n")); }
        }
    }

    fclose(file);
}

void geoSetDomainName(int iDomain, char *name)
{
    if (iDomain >= theGeometry.nDomains) { Error("Illegal domain number"); }
    if (geoGetDomain(name) != -1) { Error("Cannot use the same name for two domains"); }
    sprintf(theGeometry.theDomains[iDomain]->name, "%s", name);
}

int geoGetDomain(char *name)
{
    int theIndex = -1;
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name, theDomain->name, MAXNAME) == 0) { theIndex = iDomain; }
    }
    return theIndex;
}

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626, -0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = {0.577350269189626, -0.577350269189626, -0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = {1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = {0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = {0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = {0.166666666666667, 0.166666666666667, 0.166666666666667};

femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (theRule == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    if (type == FEM_QUAD && n == 4)
    {
        theRule->n = 4;
        theRule->xsi = _gaussQuad4Xsi;
        theRule->eta = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theRule->n = 3;
        theRule->xsi = _gaussTri3Xsi;
        theRule->eta = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight;
    }
    else { Error("Cannot create such an integration rule !"); }
    return theRule;
}

void femIntegrationFree(femIntegration *theRule) { free(theRule); }

void _q1c0_x(double *xsi, double *eta)
{
    xsi[0] = 1.0;
    eta[0] = 1.0;
    xsi[1] = -1.0;
    eta[1] = 1.0;
    xsi[2] = -1.0;
    eta[2] = -1.0;
    xsi[3] = 1.0;
    eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =  (1.0 + eta) / 4.0;
    dphidxsi[1] = -(1.0 + eta) / 4.0;
    dphidxsi[2] = -(1.0 - eta) / 4.0;
    dphidxsi[3] =  (1.0 - eta) / 4.0;
    dphideta[0] =  (1.0 + xsi) / 4.0;
    dphideta[1] =  (1.0 - xsi) / 4.0;
    dphideta[2] = -(1.0 - xsi) / 4.0;
    dphideta[3] = -(1.0 + xsi) / 4.0;
}

void _p1c0_x(double *xsi, double *eta)
{
    xsi[0] = 0.0;
    eta[0] = 0.0;
    xsi[1] = 1.0;
    eta[1] = 0.0;
    xsi[2] = 0.0;
    eta[2] = 1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;
}

femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (theSpace == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    if (type == FEM_QUAD && n == 4)
    {
        theSpace->n = 4;
        theSpace->x2 = _q1c0_x;
        theSpace->phi2 = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theSpace->n = 3;
        theSpace->x2 = _p1c0_x;
        theSpace->phi2 = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
    }
    else { Error("Cannot create such a discrete space !"); }
    
    return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace) { free(theSpace); }

void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta) { mySpace->x2(xsi, eta); }

void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi) { mySpace->phi2(xsi, eta, phi); }

void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta) { mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta); }

void femDiscretePrint(femDiscrete *mySpace)
{
    int i, j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];

    femDiscreteXsi2(mySpace, xsi, eta);
    for (i = 0; i < n; i++)
    {
        femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
        femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);

        for (j = 0; j < n; j++)
        {
            printf("(xsi=%+.1f,eta=%+.1f) : ", xsi[i], eta[i]);
            printf(" phi(%d)=%+.1f", j, phi[j]);
            printf("   dphidxsi(%d)=%+.1f", j, dphidxsi[j]);
            printf("   dphideta(%d)=%+.1f \n", j, dphideta[j]);
        }
        printf(" \n");
    }
}

femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    if (theSystem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    femFullSystemAlloc(theSystem, size);
    femFullSystemInit(theSystem);
    return theSystem;
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemAlloc(femFullSystem *mySystem, int size)
{
    int i;
    double *elem = malloc(sizeof(double) * size * (size + 1));
    if (elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    mySystem->A = malloc(sizeof(double *) * size);
    if (mySystem->A == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    mySystem->B = elem;
    mySystem->A[0] = elem + size;
    mySystem->size = size;
    for (i = 1; i < size; i++) { mySystem->A[i] = mySystem->A[i - 1] + size; }
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int i, size = mySystem->size;
    for (i = 0; i < size * (size + 1); i++) { mySystem->B[i] = 0; }
}

void femFullSystemPrint(femFullSystem *mySystem)
{
    double **A, *B;
    int i, j, size;

    A = mySystem->A;
    B = mySystem->B;
    size = mySystem->size;

    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            if (A[i][j] == 0) { printf("         "); }
            else              { printf(" %+.1e", A[i][j]); }
        }
        printf(" :  %+.1e \n", B[i]);
    }
}

double *femFullSystemEliminate(femFullSystem *mySystem)
{
    double **A, *B, factor;
    int i, j, k, size;

    A = mySystem->A;
    B = mySystem->B;
    size = mySystem->size;

    /* Gauss elimination */
    for (k = 0; k < size; k++)
    {
        if (fabs(A[k][k]) <= 1e-16) { printf("Pivot index %d  ", k); printf("Pivot value %e  ", A[k][k]); Error("Cannot eliminate with such a pivot"); }
        for (i = k + 1; i < size; i++)
        {
            factor = A[i][k] / A[k][k];
            for (j = k + 1; j < size; j++) { A[i][j] = A[i][j] - A[k][j] * factor; }
            B[i] = B[i] - B[k] * factor;
        }
    }

    /* Back-substitution */
    for (i = size - 1; i >= 0; i--)
    {
        factor = 0;
        for (j = i + 1; j < size; j++) { factor += A[i][j] * B[j]; }
        B[i] = (B[i] - factor) / A[i][i];
    }

    return mySystem->B;
}

void femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue)
{
    double **A, *B;
    int i, size;

    A = mySystem->A;
    B = mySystem->B;
    size = mySystem->size;

    for (i = 0; i < size; i++) { B[i] -= myValue * A[i][myNode]; A[i][myNode] = 0; }
    for (i = 0; i < size; i++) { A[myNode][i] = 0; }

    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

femProblem *femElasticityCreate(femGeo *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    if (theProblem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    theProblem->E = E;
    theProblem->nu = nu;
    theProblem->gx = gx;
    theProblem->gy = gy;
    theProblem->rho = rho;

    if (iCase == PLANAR_STRESS)
    {
        theProblem->A = E / (1 - nu * nu);
        theProblem->B = E * nu / (1 - nu * nu);
        theProblem->C = E / (2 * (1 + nu));
    }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM)
    {
        theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
        theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
        theProblem->C = E / (2 * (1 + nu));
    }

    theProblem->planarStrainStress = iCase;
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;

    int nNodes = theGeometry->theNodes->nNodes;
    int size = 2 * nNodes;
    theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
    if (theProblem->constrainedNodes == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    for (int i = 0; i < nNodes; i++)
    {
        theProblem->constrainedNodes[i].type = UNDEFINED;
        theProblem->constrainedNodes[i].nx = NAN;
        theProblem->constrainedNodes[i].ny = NAN;
        theProblem->constrainedNodes[i].value2 = NAN;
        theProblem->constrainedNodes[i].value2 = NAN;
    }

    theProblem->geometry = theGeometry;
    if (theGeometry->theElements->nLocalNode == 3)
    {
        theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
    }
    else if (theGeometry->theElements->nLocalNode == 4)
    {
        theProblem->space = femDiscreteCreate(4, FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
    }
    theProblem->system = femFullSystemCreate(size);
    return theProblem;
}

void femElasticityFree(femProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    for (int i = 0; i < theProblem->nBoundaryConditions; i++) { free(theProblem->conditions[i]); }
    free(theProblem->conditions);
    free(theProblem->soluce);
    free(theProblem->residuals);
    free(theProblem->constrainedNodes);
    free(theProblem);
}

/*
`value2` is only used for `DIRICHLET_XY` and `DIRICHLET_NT` boundary conditions. Otherwise it is ignored and set to NAN.
*/
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2)
{
    int iDomain = geoGetDomain(nameDomain);
    if (iDomain == -1) { Error("Undefined domain :-("); }
    value2 = ((type != DIRICHLET_XY) && (type != DIRICHLET_NT)) ? NAN : value2;

    femBoundaryCondition *theBoundary = malloc(sizeof(femBoundaryCondition));
    if (theBoundary == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theBoundary->domain = theProblem->geometry->theDomains[iDomain];
    theBoundary->value1 = value1;
    theBoundary->value2 = value2;
    theBoundary->type = type;
    theProblem->nBoundaryConditions++;
    int nBoundaryConditions = theProblem->nBoundaryConditions;

    if (theProblem->conditions == NULL)
    {
        theProblem->conditions = malloc(nBoundaryConditions * sizeof(femBoundaryCondition *));
        if (theProblem->conditions == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    }

    femNodes *theNodes = theProblem->geometry->theNodes;
    if (type == DIRICHLET_X || type == DIRICHLET_Y || type == DIRICHLET_XY || type == DIRICHLET_N || type == DIRICHLET_T || type == DIRICHLET_NT)
    {
        // Ensure that there is only one Dirichlet boundary condition per domain
        for (int i = 0; i < nBoundaryConditions - 1; i++)
        {
            if (theProblem->conditions[i]->domain != theBoundary->domain) { continue; }
            femBoundaryType type_i = theProblem->conditions[i]->type;
            if (type_i == DIRICHLET_X || type_i == DIRICHLET_Y || type_i == DIRICHLET_XY || type_i == DIRICHLET_N || type_i == DIRICHLET_T || type_i == DIRICHLET_NT)
            {
                printf("\nTrying to set a second Dirichlet boundary condition on domain \"%s\"", nameDomain);
                Error("Only one Dirichlet boundary condition is allowed per domain");
            }
        }

        femDomain *theDomain = theProblem->geometry->theDomains[iDomain];
        int *elem = theDomain->elem;
        int nElem = theDomain->nElem;
        femConstrainedNode constrainedNode;
        constrainedNode.type = type;
        constrainedNode.value1 = value1;
        constrainedNode.value2 = value2;
        constrainedNode.nx = NAN;
        constrainedNode.ny = NAN;
        if (type == DIRICHLET_X || type == DIRICHLET_Y || type == DIRICHLET_XY)
        {
            for (int iElem = 0; iElem < nElem; iElem++)
            {
                for (int i = 0; i < 2; i++)
                {
                    int node = theDomain->mesh->elem[2 * elem[iElem] + i];
                    theProblem->constrainedNodes[node] = constrainedNode;
                }
            }
        }
        else // Need to compute normals
        {
            int nNodes = theNodes->nNodes;
            double *NX = malloc(nNodes * sizeof(double));
            if (NX == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
            double *NY = malloc(nNodes * sizeof(double));
            if (NY == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
            for (int iElem = 0; iElem < nElem; iElem++)
            {
                int node0 = theDomain->mesh->elem[2 * elem[iElem] + 0];
                int node1 = theDomain->mesh->elem[2 * elem[iElem] + 1];
                NX[node0] = 0;
                NY[node0] = 0;
                NX[node1] = 0;
                NY[node1] = 0;
            }
            for (int iElem = 0; iElem < nElem; iElem++)
            {
                int node0 = theDomain->mesh->elem[2 * elem[iElem] + 0];
                int node1 = theDomain->mesh->elem[2 * elem[iElem] + 1];
                double tx = theNodes->X[node1] - theNodes->X[node0];
                double ty = theNodes->Y[node1] - theNodes->Y[node0];
                double nx = ty;
                double ny = -tx;
                NX[node0] += nx;
                NY[node0] += ny;
                NX[node1] += nx;
                NY[node1] += ny;
            }

            for (int iElem = 0; iElem < nElem; iElem++)
            {
                for (int i = 0; i < 2; i++)
                {
                    int node = theDomain->mesh->elem[2 * elem[iElem] + i];
                    double nx = NX[node];
                    double ny = NY[node];
                    double norm = hypot(nx, ny);
                    theProblem->constrainedNodes[node] = constrainedNode;
                    theProblem->constrainedNodes[node].nx = nx / norm;
                    theProblem->constrainedNodes[node].ny = ny / norm;
                }
            }
            free(NX);
            free(NY);
        }
    }

    theProblem->conditions = realloc(theProblem->conditions, nBoundaryConditions * sizeof(femBoundaryCondition *));
    theProblem->conditions[nBoundaryConditions - 1] = theBoundary;
}

void femElasticityPrint(femProblem *theProblem)
{
    printf("\n\n ======================================================================================= \n\n");
    printf(" Linear elasticity problem \n");
    printf("   Young modulus   E   = %14.7e [N/m2]\n", theProblem->E);
    printf("   Poisson's ratio nu  = %14.7e [-]\n", theProblem->nu);
    printf("   Density         rho = %14.7e [kg/m3]\n", theProblem->rho);
    printf("   Gravity-X       gx  = %14.7e [m/s2]\n", theProblem->gx);
    printf("   Gravity-Y       gy  = %14.7e [m/s2]\n", theProblem->gy);

    if (theProblem->planarStrainStress == PLANAR_STRAIN)      { printf("   Planar strains formulation \n"); }
    else if (theProblem->planarStrainStress == PLANAR_STRESS) { printf("   Planar stresses formulation \n"); }
    else if (theProblem->planarStrainStress == AXISYM)        { printf("   Axisymmetric formulation \n"); }

    printf("   Boundary conditions : \n");
    for (int i = 0; i < theProblem->nBoundaryConditions; i++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[i];
        double value1 = theCondition->value1;
        double value2 = theCondition->value2;
        printf("  %20s :", theCondition->domain->name);
        if (theCondition->type == DIRICHLET_X)       { printf(" imposing %9.2e as the horizontal displacement  \n", value1); }
        else if (theCondition->type == DIRICHLET_Y)  { printf(" imposing %9.2e as the vertical displacement  \n", value1); }
        else if (theCondition->type == DIRICHLET_XY) { printf(" imposing %9.2e, %9.2e as the displacement  \n", value1, value2); }
        else if (theCondition->type == DIRICHLET_N)  { printf(" imposing %9.2e as the normal displacement  \n", value1); }
        else if (theCondition->type == DIRICHLET_T)  { printf(" imposing %9.2e as the tangential displacement  \n", value1); }
        else if (theCondition->type == DIRICHLET_NT) { printf(" imposing (%9.2e, %9.2e) as the normal and tangential displacement  \n", value1, value2); }
        else if (theCondition->type == NEUMANN_X)    { printf(" imposing %9.2e as the horizontal force  \n", value1); }
        else if (theCondition->type == NEUMANN_Y)    { printf(" imposing %9.2e as the vertical force  \n", value1); }
        else if (theCondition->type == NEUMANN_N)    { printf(" imposing %9.2e as the normal force  \n", value1); }
        else if (theCondition->type == NEUMANN_T)    { printf(" imposing %9.2e as the tangential force  \n", value1); }
    }
    printf(" ======================================================================================= \n\n");
}

void femElasticityWrite(femProblem *theProblem, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }

    if (theProblem->planarStrainStress == PLANAR_STRESS)      { fprintf(file, "Type of problem    :  Planar stresses  \n"); }
    else if (theProblem->planarStrainStress == PLANAR_STRAIN) { fprintf(file, "Type of problem    :  Planar strains \n"); }
    else if (theProblem->planarStrainStress == AXISYM)        { fprintf(file, "Type of problem    :  Axi-symetric problem \n"); }
    else                                                      { fprintf(file, "Type of problem    :  Undefined  \n"); }

    fprintf(file, "Young modulus      : %14.7e  \n", theProblem->E);
    fprintf(file, "Poisson ratio      : %14.7e  \n", theProblem->nu);
    fprintf(file, "Mass density       : %14.7e  \n", theProblem->rho);
    fprintf(file, "Gravity-X          : %14.7e  \n", theProblem->gx);
    fprintf(file, "Gravity-Y          : %14.7e  \n", theProblem->gy);

    for (int i = 0; i < theProblem->nBoundaryConditions; i++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[i];
        double value1 = theCondition->value1;
        double value2 = theCondition->value2;
        fprintf(file, "Boundary condition : ");

        if (theCondition->type == DIRICHLET_X)       { fprintf(file, " Dirichlet-X        = %14.7e, %14.7e ", value1, NAN); }
        else if (theCondition->type == DIRICHLET_Y)  { fprintf(file, " Dirichlet-Y        = %14.7e, %14.7e ", value1, NAN); }
        else if (theCondition->type == DIRICHLET_XY) { fprintf(file, " Dirichlet-XY       = %14.7e, %14.7e ", value1, value2); }
        else if (theCondition->type == DIRICHLET_N)  { fprintf(file, " Dirichlet-N        = %14.7e, %14.7e ", value1, NAN); }
        else if (theCondition->type == DIRICHLET_T)  { fprintf(file, " Dirichlet-T        = %14.7e, %14.7e ", value1, NAN); }
        else if (theCondition->type == DIRICHLET_NT) { fprintf(file, " Dirichlet-NT       = %14.7e, %14.7e ", value1, value2); }
        else if (theCondition->type == NEUMANN_X)    { fprintf(file, " Neumann-X          = %14.7e, %14.7e ", value1, NAN); }
        else if (theCondition->type == NEUMANN_Y)    { fprintf(file, " Neumann-Y          = %14.7e, %14.7e ", value1, NAN); }
        else if (theCondition->type == NEUMANN_N)    { fprintf(file, " Neumann-N          = %14.7e, %14.7e ", value1, NAN); }
        else if (theCondition->type == NEUMANN_T)    { fprintf(file, " Neumann-T          = %14.7e, %14.7e ", value1, NAN); }
        else                                         { fprintf(file, " Undefined          = %14.7e, %14.7e ", NAN, NAN); }
    
        fprintf(file, ": %s\n", theCondition->domain->name);
    }
    fclose(file);
}

femProblem *femElasticityRead(femGeo *theGeometry, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }
    femProblem *theProblem = malloc(sizeof(femProblem));
    if (theProblem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;

    int nNodes = theGeometry->theNodes->nNodes;
    int size = 2 * nNodes;
    theProblem->soluce = malloc(size * sizeof(double));
    if (theProblem->soluce == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    theProblem->residuals = malloc(size * sizeof(double));
    if (theProblem->residuals == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    for (int i = 0; i < size; i++) { theProblem->soluce[i] = 0.0; theProblem->residuals[i] = 0.0; }

    theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
    if (theProblem->constrainedNodes == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    for (int i = 0; i < nNodes; i++)
    {
        theProblem->constrainedNodes[i].type   = UNDEFINED;
        theProblem->constrainedNodes[i].nx     = NAN;
        theProblem->constrainedNodes[i].ny     = NAN;
        theProblem->constrainedNodes[i].value2 = NAN;
        theProblem->constrainedNodes[i].value2 = NAN;
    }

    theProblem->geometry = theGeometry;
    if (theGeometry->theElements->nLocalNode == 3)
    {
        theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
    }
    else if (theGeometry->theElements->nLocalNode == 4)
    {
        theProblem->space = femDiscreteCreate(4, FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4, FEM_QUAD);
    }
    theProblem->system = femFullSystemCreate(size);

    char theLine[MAXNAME];
    char theDomain[MAXNAME];
    char theArgument[MAXNAME];
    double value1, value2;
    femBoundaryType typeCondition;

    while (!feof(file))
    {
        ErrorScan(fscanf(file, "%19[^\n]s \n", (char *)&theLine));
        if (strncasecmp(theLine, "Type of problem     ", 19) == 0)
        {
            ErrorScan(fscanf(file, ":  %[^\n]s \n", (char *)&theArgument));
            if (strncasecmp(theArgument, "Planar stresses", 13) == 0)           { theProblem->planarStrainStress = PLANAR_STRESS; }
            if (strncasecmp(theArgument, "Planar strains", 13) == 0)       { theProblem->planarStrainStress = PLANAR_STRAIN; }
            if (strncasecmp(theArgument, "Axi-symetric problem", 13) == 0) { theProblem->planarStrainStress = AXISYM; }
        }
        if (strncasecmp(theLine, "Young modulus       ", 19) == 0)     { ErrorScan(fscanf(file, ":  %le\n", &theProblem->E)); }
        if (strncasecmp(theLine, "Poisson ratio       ", 19) == 0) { ErrorScan(fscanf(file, ":  %le\n", &theProblem->nu)); }
        if (strncasecmp(theLine, "Mass density        ", 19) == 0) { ErrorScan(fscanf(file, ":  %le\n", &theProblem->rho)); }
        if (strncasecmp(theLine, "Gravity-X           ", 19) == 0) { ErrorScan(fscanf(file, ":  %le\n", &theProblem->gx)); }
        if (strncasecmp(theLine, "Gravity-Y           ", 19) == 0) { ErrorScan(fscanf(file, ":  %le\n", &theProblem->gy)); }
        if (strncasecmp(theLine, "Boundary condition  ", 19) == 0)
        {
            ErrorScan(fscanf(file, ":  %19s = %le, %le : %[^\n]s\n", (char *)&theArgument, &value1, &value2, (char *)&theDomain));
            if (strncasecmp(theArgument, "Dirichlet-X", 19) == 0)  { typeCondition = DIRICHLET_X; }
            if (strncasecmp(theArgument, "Dirichlet-Y", 19) == 0)  { typeCondition = DIRICHLET_Y; }
            if (strncasecmp(theArgument, "Dirichlet-XY", 19) == 0) { typeCondition = DIRICHLET_XY; }
            if (strncasecmp(theArgument, "Dirichlet-N", 19) == 0)  { typeCondition = DIRICHLET_N; }
            if (strncasecmp(theArgument, "Dirichlet-T", 19) == 0)  { typeCondition = DIRICHLET_T; }
            if (strncasecmp(theArgument, "Dirichlet-NT", 19) == 0) { typeCondition = DIRICHLET_NT; }
            if (strncasecmp(theArgument, "Neumann-X", 19) == 0)    { typeCondition = NEUMANN_X; }
            if (strncasecmp(theArgument, "Neumann-Y", 19) == 0)    { typeCondition = NEUMANN_Y; }
            if (strncasecmp(theArgument, "Neumann-N", 19) == 0)    { typeCondition = NEUMANN_N; }
            if (strncasecmp(theArgument, "Neumann-T", 19) == 0)    { typeCondition = NEUMANN_T; }
            femElasticityAddBoundaryCondition(theProblem, theDomain, typeCondition, value1, value2);
        }
        ErrorScan(fscanf(file, "\n"));
    }

    int iCase = theProblem->planarStrainStress;
    double E = theProblem->E;
    double nu = theProblem->nu;

    if (iCase == PLANAR_STRESS)
    {
        theProblem->A = E / (1 - nu * nu);
        theProblem->B = E * nu / (1 - nu * nu);
        theProblem->C = E / (2 * (1 + nu));
    }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM)
    {
        theProblem->A = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
        theProblem->B = E * nu / ((1 + nu) * (1 - 2 * nu));
        theProblem->C = E / (2 * (1 + nu));
    }

    fclose(file);
    return theProblem;
}

void femSolutionWrite(int nNodes, int nfields, double *data, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }
    fprintf(file, "Size %d,%d\n", nNodes, nfields);
    for (int i = 0; i < nNodes; i++)
    {
        for (int j = 0; j < nfields-1; j++) { fprintf(file, "%.18le,", data[i * nfields + j]); }
        fprintf(file, "%.18le", data[i * nfields + nfields-1]);
        fprintf(file, "\n");
    }
    fclose(file);
}

int femSolutiondRead(int allocated_size, double *value, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }
    int nNodes,nFields;
    ErrorScan(fscanf(file, "Size %d,%d\n", &nNodes, &nFields));
    if (nNodes * nFields > allocated_size)
    {
        printf("Error: allocated size is %d, but the solution file has %d nodes and %d fields", allocated_size, nNodes, nFields);
        Error("The allocated size is too small for femSolutiondRead");
    }
    for (int i = 0; i < nNodes; i++)
    {
        for (int j = 0; j < nFields; j++) { ErrorScan(fscanf(file, "%le,", &value[i * nFields + j])); }
        ErrorScan(fscanf(file, "\n"));
    }
    printf("Reading solution of shape (%d,%d)\n", nNodes, nFields);
    fclose(file);
    return nNodes * nFields;
}

double femMin(double *x, int n)
{
    double myMin = x[0];
    int i;
    for (i = 1; i < n; i++) { myMin = fmin(myMin, x[i]); }
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i = 1; i < n; i++) { myMax = fmax(myMax, x[i]); }
    return myMax;
}

void femError(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s:%d at line %d : \n  %s\n", file, line, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femErrorScan(int test, int line, char *file)
{
    if (test >= 0) { return; }

    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s:%d at line %d : \n", file, line, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femWarning(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s:%d at line %d : \n  %s\n", file, line, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
}