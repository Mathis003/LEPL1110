#include "../include/fem.h"


/* Variables (needed for linux compilation)*/

femGeometry theGeometry;

/**********************************/
/******* Solvers functions *******/
/**********************************/

/* Full System */
femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *system = malloc(sizeof(femFullSystem));
    if (system == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    femFullSystemAlloc(system, size);
    femFullSystemInit(system, size);
    return system;
}

void femFullSystemFree(femFullSystem *system)
{   
    if (system != NULL)
    {
        if (system->A != NULL) { free(system->A); system->A = NULL; }
        if (system->B != NULL) { free(system->B); system->B = NULL; }
        free(system); system = NULL;
    }
}

void femFullSystemAlloc(femFullSystem *system, int size)
{
    double *elem = malloc(sizeof(double) * size * (size + 1));
    if (elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    system->A = malloc(sizeof(double *) * size);
    if (system->A == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    system->B = elem;
    system->A[0] = elem + size;
    for (int i = 1; i < size; i++) { system->A[i] = system->A[i - 1] + size; }
}

void femFullSystemInit(femFullSystem *system, int size)
{
    for (int i = 0; i < size * (size + 1); i++) { system->B[i] = 0; }
}

void femFullSystemAssemble(femFullSystem *system, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, double xLoc, int nLoc, const double FACTOR)
{
    double **A = system->A;
    double *B  = system->B;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double gx  = theProblem->gx;
    double gy  = theProblem->gy;

    if (theProblem->planarStrainStress == PLANAR_STRAIN || theProblem->planarStrainStress == PLANAR_STRESS)
    {
        for (int i = 0; i < nLoc; i++)
        {
            for (int j = 0; j < nLoc; j++)
            {
                A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac;
                A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac;
                A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac;
                A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac;
            }
            B[mapX[i]] += phi[i] * gx * rho * weightedJac * FACTOR;
            B[mapY[i]] += phi[i] * gy * rho * weightedJac * FACTOR;
        }
    }
    else if (theProblem->planarStrainStress == AXISYM)
    {
        for (int i = 0; i < nLoc; i++)
        {
            for (int j = 0; j < nLoc; j++)
            {
                A[mapX[i]][mapX[j]] += (dphidx[i] * a * xLoc * dphidx[j] + dphidy[i] * c * xLoc * dphidy[j] + dphidx[i] * b * phi[j] + phi[i] * (b * dphidx[j] + a * phi[j] / xLoc)) * weightedJac;
                A[mapX[i]][mapY[j]] += (dphidx[i] * b * xLoc * dphidy[j] + dphidy[i] * c * xLoc * dphidx[j] + phi[i] * b * dphidy[j]) * weightedJac;
                A[mapY[i]][mapX[j]] += (dphidy[i] * b * xLoc * dphidx[j] + dphidx[i] * c * xLoc * dphidy[j] + dphidy[i] * b * phi[j]) * weightedJac;
                A[mapY[i]][mapY[j]] += (dphidy[i] * a * xLoc * dphidy[j] + dphidx[i] * c * xLoc * dphidx[j]) * weightedJac;
            }
            B[mapX[i]] += phi[i] * xLoc * gx * rho * weightedJac * FACTOR;
            B[mapY[i]] += phi[i] * xLoc * gy * rho * weightedJac * FACTOR;
        }
    }
    else { Error("Unexpected problem type"); }
}

void femFullSystemConstrainXY(femFullSystem *system, int node, double value, int size)
{
    double **A, *B;
    int i;

    A = system->A;
    B = system->B;

    for (i = 0; i < size; i++) { B[i] -= value * A[i][node]; A[i][node] = 0; }
    for (i = 0; i < size; i++) { A[node][i] = 0; }

    A[node][node] = 1;
    B[node] = value;
}

void femFullSystemConstrainNT(femFullSystem *system, int size, int node1, int node2, double a, double b)
{
    double **A, *B;
    int i;

    A = system->A;
    B = system->B;

    // Force this constraint to be applied : node2 = b + a * node1
    for (int i = 0; i < size; i++)
    {
        if      (i == node2) { A[node2][i] = 1.0; }
        else if (i == node1) { A[node2][i] = - a; }
        else                 { A[node2][i] = 0.0; }
    }
    B[node2] = b;
}

void femFullSystemSetSystem(femFullSystem *system, double **A, double *B)
{
    system->A = A;
    system->B = B;
}

double femFullSystemGet(femFullSystem *system, int row, int col) { return system->A[row][col]; }

double *femFullSystemEliminate(femFullSystem *system, int size)
{
    double **A, *B, factor;
    int i, j, k;

    A = system->A;
    B = system->B;

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

    return system->B;
}

void femFullSystemPrint(femFullSystem *system, int size)
{
    double **A, *B;
    int i, j;

    A = system->A;
    B = system->B;

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

void femFullSystemPrintInfos(femFullSystem *system, int size)
{
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int) sizeof(double) * size * (size + 1));     
}


/* Band System */
femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *system = malloc(sizeof(femBandSystem));
    if (system == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    femBandSystemAlloc(system, size, band);
    femBandSystemInit(system, size);
    return system;
}

void femBandSystemAlloc(femBandSystem *system, int size, int band)
{   
    system->B = malloc(sizeof(double) * size * (band + 1));
    if (system->B == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    system->A = malloc(sizeof(double *) * size);
    if (system->A == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }     
    system->band = band;
    system->A[0] = system->B + size;
    for (int i = 1 ; i < size ; i++) { system->A[i] = system->A[i-1] + band - 1; }
}

void femBandSystemInit(femBandSystem *system, int size)
{
    double **A = system->A;
    double *B  = system->B;
    int band   = system->band;
    for (int i = 0 ; i < size * (band + 1) ; i++) { B[i] = 0; }
}

void femBandSystemFree(femBandSystem *system)
{
    free(system->B);
    free(system->A); 
    free(system);
}

void femBandSystemSetSystem(femBandSystem *system, double **A, double *B)
{
    system->A = A;
    system->B = B;
}

int isInBand(int band, int myRow, int myCol) { return myCol >= myRow && myCol < myRow + band; }

void femBandSystemPrint(femBandSystem *system, int size)
{
    double  **A, *B;
    int band;
    A    = system->A;
    B    = system->B;
    band = system->band;

    for (int i = 0; i < size; i++)
    {
        for (int j = i; j < i + band; j++)
        {
            if (A[i][j] == 0) { printf("         "); }
            else              { printf(" %+.1e",A[i][j]); }
        }
        printf(" :  %+.1e \n",B[i]);
    }
}

void femBandSystemPrintInfos(femBandSystem *system, int size)
{
    int band = system->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int) sizeof(double) * size * (band + 1));
}

// TODO
void femBandSystemAssemble(femBandSystem *system, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, double xLoc, int nLoc, const double FACTOR)
{
    double **A = system->A;
    double *B  = system->B;
    int band   = system->band;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;
    double rho = theProblem->rho;
    double gx  = theProblem->gx;
    double gy  = theProblem->gy;

    if (theProblem->planarStrainStress == PLANAR_STRAIN || theProblem->planarStrainStress == PLANAR_STRESS)
    {
        for (int i = 0; i < nLoc; i++)
        {
            for (int j = 0; j < nLoc; j++)
            {
                A[mapX[i]][mapX[j]] += isInBand(band, mapX[i], mapX[j]) ? (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac : 0.0;
                A[mapX[i]][mapY[j]] += isInBand(band, mapX[i], mapY[j]) ? (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac : 0.0;
                A[mapY[i]][mapX[j]] += isInBand(band, mapY[i], mapX[j]) ? (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac : 0.0;
                A[mapY[i]][mapY[j]] += isInBand(band, mapY[i], mapY[j]) ? (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac : 0.0;
            }
            B[mapX[i]] += phi[i] * gx * rho * weightedJac * FACTOR;
            B[mapY[i]] += phi[i] * gy * rho * weightedJac * FACTOR;
        }
    }
    else if (theProblem->planarStrainStress == AXISYM)
    {
        for (int i = 0; i < nLoc; i++)
        {
            for (int j = 0; j < nLoc; j++)
            {
                A[mapX[i]][mapX[j]] += isInBand(band, mapX[i], mapX[j]) ? (dphidx[i] * a * xLoc * dphidx[j] + dphidy[i] * c * xLoc * dphidy[j] + phi[i] * ((b * dphidx[j]) + (a * phi[j] / xLoc)) + dphidx[i] * b * phi[j]) * weightedJac : 0.0;
                A[mapX[i]][mapY[j]] += isInBand(band, mapX[i], mapY[j]) ? (dphidx[i] * b * xLoc * dphidy[j] + dphidy[i] * c * xLoc * dphidx[j] + phi[i] * b * dphidy[j]) * weightedJac : 0.0;
                A[mapY[i]][mapX[j]] += isInBand(band, mapY[i], mapX[j]) ? (dphidy[i] * b * xLoc * dphidx[j] + dphidx[i] * c * xLoc * dphidy[j] + dphidy[i] * b * phi[j]) * weightedJac : 0.0;
                A[mapY[i]][mapY[j]] += isInBand(band, mapY[i], mapY[j]) ? (dphidy[i] * a * xLoc * dphidy[j] + dphidx[i] * c * xLoc * dphidx[j]) * weightedJac : 0.0;
            }
            B[mapX[i]] += phi[i] * xLoc * gx * rho * weightedJac * FACTOR;
            B[mapY[i]] += phi[i] * xLoc * gy * rho * weightedJac * FACTOR;
        }
    }
    else { Error("Unexpected problem type"); }
}

// TODO
double femBandSystemGet(femBandSystem *system, int myRow, int myCol)
{
    double **A = system->A;
    int band   = system->band;
    // return A[myRow][myCol];
    return isInBand(band, myRow, myCol) ? A[myRow][myCol] : 0.0;
}

// TODO : Is this good ?
void femBandSystemConstrainXY(femBandSystem *system, int node, double value, int size)
{
    double **A, *B;
    int i, band;

    A = system->A;
    B = system->B;
    band = system->band;

    for (i = 0; i < size; i++)
    {
        if (isInBand(band, i, node))
        {
            B[i] -= value * A[i][node];
            A[i][node] = 0;
        }
    }
    for (i = 0; i < size; i++) { if (isInBand(band, i, node)) { A[node][i] = 0; }}

    A[node][node] = 1.0;
    B[node] = value;


    // int lowerBound = (node >= band) ? node - band + 1 : 0;
    // for (i = lowerBound ; i < node; i++)
    // {
    //     B[i] -= value * A[i][node];
    //     A[i][node] = 0;
    // }

    // int upperBound = (node + band >= size) ? size : node + band;
    // for (i = node + 1; i < upperBound ; i++)
    // {
    //     B[i] -= value * A[node][i];
    //     A[node][i] = 0;
    // }

    // A[node][node] = 1;
    // B[node] = value;
}

// TODO : Implement this function (Not good for now) 
void femBandSystemConstrainNT(femBandSystem *system, int size, int node1, int node2, double a, double b)
{
    double **A, *B;
    int i, band;

    A = system->A;
    B = system->B;
    band = system->band;

    // Force this constraint to be applied : node2 = b + a * node1
    for (int i = 0; i < size; i++)
    {
        if (!(isInBand(band, node2, i))) { continue; }

        if      (i == node2) { A[node2][i] = 1.0; }
        else if (i == node1) { A[node2][i] = - a; }
        else                 { A[node2][i] = 0.0; }
    }
    B[node2] = b;
}

double *femBandSystemEliminate(femBandSystem *system, int size)
{
    double  **A, *B, factor;
    int     i, j, k, jend, band;
    A    = system->A;
    B    = system->B;
    band = system->band;

    /* Gauss elimination */
    for (k = 0; k < size; k++)
    {
        if (fabs(A[k][k]) <= 1e-16) { Error("Cannot eliminate with such a pivot.\n"); }
        jend = (k + band < size) ? k + band : size;
        for (i = k + 1; i < jend; i++)
        {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) { A[i][j] -= factor * A[k][j]; }
            B[i] -= factor * B[k];
        }    
    }
    
    /* Back-substitution */
    for (i = size - 1; i >= 0 ; i--)
    {
        factor = 0;
        jend = (i + band < size) ? i + band : size;
        for (j = i + 1 ; j < jend; j++) { factor += A[i][j] * B[j]; }
        B[i] = ( B[i] - factor) / A[i][i];
    }
    return system->B;
}

double *positionMeshNodes;

int comparPositionNode(const void *a, const void *b)
{
    const int *nodePos_a = (const int *) a;
    const int *nodePos_b = (const int *) b;
    
    double diff = positionMeshNodes[*nodePos_a] - positionMeshNodes[*nodePos_b];
    return (diff < 0) - (diff > 0);
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int nNodes = theMesh->nodes->nNodes;
    int *mapper = (int *) malloc(nNodes * sizeof(int));
    if (mapper == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < nNodes; i++) { mapper[i] = i; }

    switch (renumType)
    {        
        case FEM_NO :
            break;

        case FEM_XNUM :
            positionMeshNodes = theMesh->nodes->X;
            qsort(mapper, nNodes, sizeof(int), comparPositionNode);
            break;

        case FEM_YNUM :
            positionMeshNodes = theMesh->nodes->Y;
            qsort(mapper, nNodes, sizeof(int), comparPositionNode);
            break;    

        default : Error("Unexpected renumbering option");
    }

    for (int i = 0; i < nNodes; i++) { theMesh->nodes->number[mapper[i]] = i; }
    free(mapper); 
}

int femMeshComputeBand(femMesh *theMesh)
{
    int maxNum, minNum, nodeNum, elemNum;
    int band = 0;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++)
    {   
        maxNum = INT_MIN;
        minNum = INT_MAX;

        for (int j = 0; j < theMesh->nLocalNode; j++)
        {
            elemNum = theMesh->elem[iElem * theMesh->nLocalNode + j];
            nodeNum = theMesh->nodes->number[elemNum];

            maxNum = (nodeNum > maxNum) ? nodeNum : maxNum;
            minNum = (nodeNum < minNum) ? nodeNum : minNum;
        }
        if (band < maxNum - minNum) { band = maxNum - minNum; }
    }
    return ++band;
}


/* Solver Abstraction */

femSolver *femSolverCreate(int size)
{
    femSolver *solver = malloc(sizeof(femSolver));
    if (solver == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    solver->size = size;
    if (solver == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    return solver;
}

femSolver *femSolverFullCreate(int size)
{
    femSolver *mySolver = femSolverCreate(size);
    mySolver->type = FEM_FULL;
    mySolver->solver = (femSolver *) femFullSystemCreate(size);
    return mySolver;
}

femSolver *femSolverBandCreate(int size, int band)
{
    femSolver *mySolver = femSolverCreate(size);
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *) femBandSystemCreate(size, band);
    return mySolver;
}

void femSolverFree(femSolver *mySolver)
{
    if(mySolver == NULL) { return; }
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemFree((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemFree((femBandSystem *)mySolver->solver); break;
        default :       Error("Unexpected solver type");
    }
    free(mySolver);
}

void femSolverSetSystem(femSolver *mySolver, double **A, double *B)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemSetSystem((femFullSystem *) mySolver->solver, A, B); break;
        case FEM_BAND : femBandSystemSetSystem((femBandSystem *) mySolver->solver, A, B); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverInit(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemInit((femFullSystem *) mySolver->solver, mySolver->size); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *) mySolver->solver, mySolver->size); break;
        default :       Error("Unexpected solver type");
    }
}

double femSolverGet(femSolver *mySolver, int i, int j)
{
    double value = 0.0;
    switch (mySolver->type)
    {  
        case FEM_FULL : value = femFullSystemGet((femFullSystem *) mySolver->solver, i, j); break;
        case FEM_BAND : value = femBandSystemGet((femBandSystem *) mySolver->solver, i, j); break;
        default :       Error("Unexpected solver type");
    }
    return value;
}

void femSolverPrint(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemPrint((femFullSystem *) mySolver->solver, mySolver->size); break;
        case FEM_BAND : femBandSystemPrint((femBandSystem *) mySolver->solver, mySolver->size); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverPrintInfos(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemPrintInfos((femFullSystem *) mySolver->solver, mySolver->size); break;
        case FEM_BAND : femBandSystemPrintInfos((femBandSystem *) mySolver->solver, mySolver->size); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverAssemble(femSolver *mySolver, femProblem *theProblem, int *mapX, int *mapY, double *phi, double *dphidx, double *dphidy, double weightedJac, double xLoc, int nLoc, const double FACTOR)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *) mySolver->solver, theProblem, mapX, mapY, phi, dphidx, dphidy, weightedJac, xLoc, nLoc, FACTOR); break;
        case FEM_BAND : femBandSystemAssemble((femBandSystem *) mySolver->solver, theProblem, mapX, mapY, phi, dphidx, dphidy, weightedJac, xLoc, nLoc, FACTOR); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverSystemConstrainXY(femSolver *mySolver, int node, double value)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemConstrainXY((femFullSystem *) mySolver->solver, node, value, mySolver->size); break;
        case FEM_BAND : femBandSystemConstrainXY((femBandSystem *) mySolver->solver, node, value, mySolver->size); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverSystemConstrainNT(femSolver *mySolver, int node1, int node2, double a, double b)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemConstrainNT((femFullSystem *) mySolver->solver, mySolver->size, node1 ,node2, a, b); break;
        case FEM_BAND : femBandSystemConstrainNT((femBandSystem *) mySolver->solver, mySolver->size, node1, node2, a, b); break;
        default :       Error("Unexpected solver type");
    }
}

double *femSolverEliminate(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : return femFullSystemEliminate((femFullSystem *) mySolver->solver, mySolver->size); break;
        case FEM_BAND : return femBandSystemEliminate((femBandSystem *) mySolver->solver, mySolver->size); break;
        default :       Error("Unexpected solver type");
    }
}

int femSolverConverged(femSolver *mySolver)
{
    int testConvergence;
    switch (mySolver->type)
    {
        case FEM_FULL : testConvergence = 1; break;
        case FEM_BAND : testConvergence = 1; break;
        default :       Error("Unexpected solver type");
    }
    return testConvergence;
}

double **getMatrixA(femSolver *mySolver)
{
    switch(mySolver->type)
    {
        case FEM_FULL : return ((femFullSystem *) mySolver->solver)->A;
        case FEM_BAND : return ((femBandSystem *) mySolver->solver)->A;
        default :       Error("Unexpected solver type");
    }
}

double *getVectorB(femSolver *mySolver)
{
    switch(mySolver->type)
    {
        case FEM_FULL : return ((femFullSystem *) mySolver->solver)->B;
        case FEM_BAND : return ((femBandSystem *) mySolver->solver)->B;
        default :       Error("Unexpected solver type");
    }
}

double getSizeMatrix(femSolver *mySolver)
{
    return mySolver->size;
}


/**********************************************************/
/******* Discrete space + Integrate space functions *******/
/**********************************************************/

femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (theRule == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    if (type == FEM_QUAD && n == 4)
    {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight;
    }
    else if (type == FEM_EDGE && n == 2) {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight;
        }
    else { Error("Cannot create such an integration rule !"); }
    return theRule;
}

void femIntegrationFree(femIntegration *theRule) { free(theRule); }

void _q1c0_x(double *xsi, double *eta)
{
    xsi[0] = 1.0;  eta[0] = 1.0;
    xsi[1] = -1.0; eta[1] = 1.0;
    xsi[2] = -1.0; eta[2] = -1.0;
    xsi[3] = 1.0;  eta[3] = -1.0;
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
    dphidxsi[0] =  (1.0 + eta) / 4.0; dphidxsi[1] = -(1.0 + eta) / 4.0;
    dphidxsi[2] = -(1.0 - eta) / 4.0; dphidxsi[3] =  (1.0 - eta) / 4.0;
    dphideta[0] =  (1.0 + xsi) / 4.0; dphideta[1] =  (1.0 - xsi) / 4.0;
    dphideta[2] = -(1.0 - xsi) / 4.0; dphideta[3] = -(1.0 + xsi) / 4.0;
}

void _p1c0_x(double *xsi, double *eta)
{
    xsi[0] = 0.0; eta[0] = 0.0;
    xsi[1] = 1.0; eta[1] = 0.0;
    xsi[2] = 0.0; eta[2] = 1.0;
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

void _e1c0_x(double *xsi) 
{
    xsi[0] = -1.0;  
    xsi[1] =  1.0;  
}

void _e1c0_phi(double xsi,  double *phi)
{
    phi[0] = (1 - xsi) / 2.0;  
    phi[1] = (1 + xsi) / 2.0;
}

void _e1c0_dphidx(double xsi, double *dphidxsi)
{
    dphidxsi[0] = -0.5;  
    dphidxsi[1] =  0.5;
}

femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (theSpace == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    if (type == FEM_QUAD && n == 4)
    {
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
    }
    else if (type == FEM_EDGE && n == 2)
    {
        theSpace->n       = 2;
        theSpace->x       = _e1c0_x;
        theSpace->phi     = _e1c0_phi;
        theSpace->dphidx  = _e1c0_dphidx;
        }
    else { Error("Cannot create such a discrete space !"); }
    return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace) { free(theSpace); }

void femDiscreteXsi2(femDiscrete *mySpace, double *xsi, double *eta) { mySpace->x2(xsi, eta); }

void femDiscretePhi2(femDiscrete *mySpace, double xsi, double eta, double *phi) { mySpace->phi2(xsi, eta, phi); }

void femDiscreteDphi2(femDiscrete *mySpace, double xsi, double eta, double *dphidxsi, double *dphideta) { mySpace->dphi2dx(xsi, eta, dphidxsi, dphideta); }

void femDiscreteXsi(femDiscrete *mySpace, double *xsi) { mySpace->x(xsi); }

void femDiscretePhi(femDiscrete *mySpace, double xsi, double *phi) { mySpace->phi(xsi, phi); }

void femDiscreteDphi(femDiscrete *mySpace, double xsi, double *dphidxsi) { mySpace->dphidx(xsi, dphidxsi); }

void femDiscretePrint(femDiscrete *mySpace)
{
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    int n = mySpace->n;

    femDiscreteXsi2(mySpace, xsi, eta);
    for (int i = 0; i < n; i++)
    {
        femDiscretePhi2(mySpace, xsi[i], eta[i], phi);
        femDiscreteDphi2(mySpace, xsi[i], eta[i], dphidxsi, dphideta);

        for (int j = 0; j < n; j++)
        {
            printf("(xsi = %+.1f, eta = %+.1f) : ", xsi[i], eta[i]);
            printf(" phi(%d) = %+.1f", j, phi[j]);
            printf("   dphidxsi(%d) = %+.1f", j, dphidxsi[j]);
            printf("   dphideta(%d) = %+.1f \n", j, dphideta[j]);
        }
        printf(" \n");
    }
}

/************************************/
/******* Elasticity functions *******/
/************************************/

femProblem *femElasticityCreate(femGeometry *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    if (theProblem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
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

    theProblem->planarStrainStress  = iCase;
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;

    int nNodes = theGeometry->theNodes->nNodes;
    int size = 2 * nNodes;
    theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
    if (theProblem->constrainedNodes == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
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
        theProblem->rule  = femIntegrationCreate(3, FEM_TRIANGLE);
    }
    else if (theGeometry->theElements->nLocalNode == 4)
    {
        theProblem->space = femDiscreteCreate(4, FEM_QUAD);
        theProblem->rule  = femIntegrationCreate(4, FEM_QUAD);
    }
    theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
    theProblem->ruleEdge  = femIntegrationCreate(2, FEM_EDGE); 
    theProblem->solver    = femSolverCreate(size);

    return theProblem;
}

void femElasticityFree(femProblem *theProblem)
{
    femSolverFree(theProblem->solver);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    for (int i = 0; i < theProblem->nBoundaryConditions; i++) { free(theProblem->conditions[i]); theProblem->conditions[i] = NULL; }
    free(theProblem->conditions);
    theProblem->conditions = NULL;
    free(theProblem->soluce);
    theProblem->soluce = NULL;
    free(theProblem->residuals);
    theProblem->residuals = NULL;
    free(theProblem->constrainedNodes);
    theProblem->constrainedNodes = NULL;
    free(theProblem);
    theProblem = NULL;
}

void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2)
{
    int iDomain = geoGetDomain(nameDomain);
    if (iDomain == -1) { Error("Undefined domain :-("); }

    // This variable is only used for 'DIRICHLET_XY' and 'DIRICHLET_NT' boundary conditions. Otherwise it is ignored (set to NAN).
    value2 = ((type != DIRICHLET_XY) && (type != DIRICHLET_NT)) ? NAN : value2;

    femBoundaryCondition *theBoundary = malloc(sizeof(femBoundaryCondition));
    if (theBoundary == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theBoundary->domain = theProblem->geometry->theDomains[iDomain];
    theBoundary->value1 = value1;
    theBoundary->value2 = value2;
    theBoundary->type = type;
    theProblem->nBoundaryConditions++;
    int nBoundaryConditions = theProblem->nBoundaryConditions;

    if (theProblem->conditions == NULL)
    {
        theProblem->conditions = malloc(nBoundaryConditions * sizeof(femBoundaryCondition *));
        if (theProblem->conditions == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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
        else
        {
            // Need to compute normals
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

femProblem *femElasticityRead(femGeometry *theGeometry, femSolverType typeSolver, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }
    femProblem *theProblem = malloc(sizeof(femProblem));
    if (theProblem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;

    int nNodes = theGeometry->theNodes->nNodes;
    int size = 2 * nNodes;
    theProblem->soluce = malloc(size * sizeof(double));
    if (theProblem->soluce == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    theProblem->residuals = malloc(size * sizeof(double));
    if (theProblem->residuals == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    for (int i = 0; i < size; i++) { theProblem->soluce[i] = 0.0; theProblem->residuals[i] = 0.0; }

    theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
    if (theProblem->constrainedNodes == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
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
    theProblem->spaceEdge = femDiscreteCreate(2, FEM_EDGE);
    theProblem->ruleEdge  = femIntegrationCreate(2, FEM_EDGE); 

    if (typeSolver == FEM_FULL) { theProblem->solver = femSolverFullCreate(size); }
    else if (typeSolver == FEM_BAND)
    {
        int band = femMeshComputeBand(theGeometry->theElements);
        printf("Band = %d, size = %d\n", band, size);
        theProblem->solver = femSolverBandCreate(size, band);
    }
    else { Error("Unknown solver type"); }

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

void femSystemWrite(double **A, double *B, int size, const char *filename)
{
    FILE *file = fopen(filename, "w");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }
    fprintf(file, "Size %d\n", size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { fprintf(file, "%.18le\t", A[i][j]); }
        fprintf(file, "\n");
    }
    for (int i = 0; i < size; i++) { fprintf(file, "%.18le\t", B[i]); }
}

int femSystemRead(double ***A, double **B, int *size, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }

    ErrorScan(fscanf(file, "Size %d\n", size));

    *A = (double **) malloc((*size) * sizeof(double *));
    *B = (double *) malloc((*size) * sizeof(double));
    if (*A == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    if (*B == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    for (int i = 0; i < *size; i++)
    {
        (*A)[i] = (double *) malloc((*size) * sizeof(double));
        if ((*A)[i] == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    }

    for (int i = 0; i < *size; i++)
    {
        for (int j = 0; j < *size; j++) { fscanf(file, "%le\t", &(*A)[i][j]); }
        fscanf(file, "\n");
    }
    for (int i = 0; i < *size; i++) { fscanf(file, "%le\t", &(*B)[i]); }

    fclose(file);
    return EXIT_SUCCESS;
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


/**********************************/
/******* Geometry functions *******/
/**********************************/

femGeometry *geoGetGeometry(void) { return &theGeometry; }

void geoSetSizeCallback(double (*geoSize)(double x, double y)) { theGeometry.geoSize = geoSize; }

void geoFree(void)
{
    if (theGeometry.theNodes)
    {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes->number);
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

void geoMeshPrint(void)
{
    femNodes *theNodes = theGeometry.theNodes;
    if (theNodes != NULL)
    {
        printf("Number of nodes %d \n", theNodes->nNodes);
        for (int i = 0; i < theNodes->nNodes; i++) { printf("%6d : %6d : %14.7e %14.7e \n", i, theNodes->number[i], theNodes->X[i], theNodes->Y[i]); }
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
    for (int i = 0; i < theNodes->nNodes; i++) { fprintf(file, "%6d : %6d %14.7e %14.7e \n", i, theNodes->number[i], theNodes->X[i], theNodes->Y[i]); }

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
    if (theNodes == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theGeometry.theNodes = theNodes;
    ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
    theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->Y == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->number = malloc(sizeof(int) * theNodes->nNodes);
    if (theNodes->number == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < theNodes->nNodes; i++) { ErrorScan(fscanf(file, "%d : %d %le %le \n", &trash, &theNodes->number[i], &theNodes->X[i], &theNodes->Y[i])); }

    femMesh *theEdges = malloc(sizeof(femMesh));
    if (theEdges == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theGeometry.theEdges = theEdges;
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
    theEdges->elem = malloc(sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
    if (theEdges->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < theEdges->nElem; ++i)
    {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash, &elem[2 * i], &elem[2 * i + 1]));
    }

    femMesh *theElements = malloc(sizeof(femMesh));
    if (theElements == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theGeometry.theElements = theElements;
    theElements->nLocalNode = 0;
    theElements->nodes = theNodes;
    char elementType[MAXNAME];
    ErrorScan(fscanf(file, "Number of %s %d \n", elementType, &theElements->nElem));
    if (strncasecmp(elementType, "triangles", MAXNAME) == 0)
    {
        theElements->nLocalNode = 3;
        theElements->elem = malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
        if (theElements->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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
        if (theElements->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theElements->nElem; ++i)
        {
            elem = theElements->elem;
            ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", &trash, &elem[4 * i], &elem[4 * i + 1], &elem[4 * i + 2], &elem[4 * i + 3]));
        }
    }

    ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
    int nDomains = theGeometry.nDomains;
    theGeometry.theDomains = malloc(sizeof(femDomain *) * nDomains);
    if (theGeometry.theDomains == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = malloc(sizeof(femDomain));
        if (theDomain == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theGeometry.theDomains[iDomain] = theDomain;
        theDomain->mesh = theEdges;
        ErrorScan(fscanf(file, "  Domain : %6d \n", &trash));
        ErrorScan(fscanf(file, "  Name : %[^\n]s \n", (char *)&theDomain->name));
        ErrorScan(fscanf(file, "  Number of elements : %6d\n", &theDomain->nElem));
        theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
        if (theDomain->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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
    if (iDomain >= theGeometry.nDomains || iDomain < 0) { Error("Illegal domain number"); }
    if (geoGetDomain(name) != -1) { Error("Cannot use the same name for two domains"); }
    sprintf(theGeometry.theDomains[iDomain]->name, "%s", name);
}

int geoGetDomain(char *name)
{
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++)
    {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name, theDomain->name, MAXNAME) == 0) { return iDomain; }
    }
    return -1;
}

/*********************************/
/******* General functions *******/
/*********************************/


double femElasticityIntegrate(femProblem *theProblem, double (*f)(double x, double y))
{
    femIntegration *theRule     = theProblem->rule;
    femGeometry    *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theMesh     = theGeometry->theElements;
    femDiscrete    *theSpace    = theProblem->space;

    double x[4] ,y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, i, map[4];
    int nLocal = theMesh->nLocalNode;
    double value = 0.0;

    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (i = 0; i < nLocal; i++)
        {
            map[i] = theMesh->elem[iElem * nLocal + i];
            x[i]   = theNodes->X[map[i]];
            y[i]   = theNodes->Y[map[i]];
        }
        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        { 
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theProblem->space,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);

            double dxdxsi = 0.0; double dxdeta = 0.0;
            double dydxsi = 0.0; double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {  
                dxdxsi += x[i] * dphidxsi[i];       
                dxdeta += x[i] * dphideta[i];   
                dydxsi += y[i] * dphidxsi[i];   
                dydeta += y[i] * dphideta[i];
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            double weightedJac = jac * weight;

            for (i = 0; i < theProblem->space->n; i++) { value += phi[i] * f(x[i], y[i]) * weightedJac; }
        }
    }
    return value;
}

double femMin(double *x, int n)
{
    double myMin = x[0];
    for (int i = 1; i < n; i++) { myMin = fmin(myMin, x[i]); }
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    for (int i = 1; i < n; i++) { myMax = fmax(myMax, x[i]); }
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