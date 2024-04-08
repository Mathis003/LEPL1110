#include "../include/fem_solver.h"


/*********************************************************************************************************************/
/******* FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER *******/
/*********************************************************************************************************************/
double *positionMeshNodes; // To renumber the mesh nodes

femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *system = malloc(sizeof(femFullSystem));
    if (system == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    femFullSystemAlloc(system, size);
    femFullSystemInit(system);
    return system;
}

void femFullSystemFree(femFullSystem *system)
{
    free(system->A);
    free(system->B);
    free(system);
}

void femFullSystemAlloc(femFullSystem *system, int size)
{
    double *elem = malloc(sizeof(double) * size * (size + 1));
    if (elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    system->A = malloc(sizeof(double *) * size);
    if (system->A == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    system->B = elem;
    system->A[0] = elem + size;
    system->size = size;
    for (int i = 1; i < size; i++) { system->A[i] = system->A[i - 1] + size; }
}

void femFullSystemInit(femFullSystem *system)
{
    int size = system->size;
    for (int i = 0; i < size * (size + 1); i++) { system->B[i] = 0; }
}

void femFullSystemAssemble(femFullSystem *system, double *Aloc, double *Bloc, int *map, int nLoc)
{
    for (int i = 0; i < nLoc; i++)
    {
        for (int j = 0; j < nLoc; j++) { system->A[map[i]][map[j]] += Aloc[i * nLoc+ j]; }
        system->B[map[i]] += Bloc[i];
    }
}

void femFullSystemConstrain(femFullSystem *system, int node, double value)
{
    double **A, *B;
    int i, size;

    A = system->A;
    B = system->B;
    size = system->size;

    for (i = 0; i < size; i++) { B[i] -= value * A[i][node]; A[i][node] = 0; }
    for (i = 0; i < size; i++) { A[node][i] = 0; }

    A[node][node] = 1;
    B[node] = value;
}

double femFullSystemGet(femFullSystem *system, int row, int col) { return system->A[row][col]; }

double *femFullSystemEliminate(femFullSystem *system)
{
    double **A, *B, factor;
    int i, j, k, size;

    A = system->A;
    B = system->B;
    size = system->size;

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

void femFullSystemPrint(femFullSystem *system)
{
    double **A, *B;
    int i, j, size;

    A = system->A;
    B = system->B;
    size = system->size;

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

void femFullSystemPrintInfos(femFullSystem *system)
{
    int size = system->size;
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int) sizeof(double) * size * (size + 1));     
}


/*********************************************************************************************************************/
/******* BAND SOLVER ***** BAND SOLVER ***** BAND SOLVER ***** BAND SOLVER ***** BAND SOLVER ***** BAND SOLVER *******/
/*********************************************************************************************************************/

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *system = malloc(sizeof(femBandSystem));
    if (system == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    femBandSystemAlloc(system, size, band);
    femBandSystemInit(system);
    return system;
}

void femBandSystemAlloc(femBandSystem *system, int size, int band)
{   
    system->B = malloc(sizeof(double) * size * (band + 1));
    if (system->B == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    system->A = malloc(sizeof(double *) * size);
    if (system->A == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }     
    system->size = size;
    system->band = band;
    system->A[0] = system->B + size;
    for (int i = 1 ; i < size ; i++) { system->A[i] = system->A[i-1] + band - 1; }
}

void femBandSystemInit(femBandSystem *system)
{
    int size = system->size;
    int band = system->band;
    for (int i = 0 ; i < size * (band + 1) ; i++) { system->B[i] = 0; }
}

void femBandSystemFree(femBandSystem *system)
{
    free(system->B);
    free(system->A); 
    free(system);
}

void femBandSystemPrint(femBandSystem *system)
{
    double  **A, *B;
    int    band, size;
    A    = system->A;
    B    = system->B;
    size = system->size;
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
  
void femBandSystemPrintInfos(femBandSystem *system)
{
    int size = system->size;
    int band = system->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int) sizeof(double) * size * (band + 1));     
}

void femBandSystemAssemble(femBandSystem *system, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int currRow, currCol;
    double **A = system->A;
    double *B  = system->B;
    
    for (int i = 0; i < nLoc; i++)
    {
        currRow = map[i];
        for (int j = 0; j < nLoc; j++)
        {
            currCol = map[j];
            if (currRow <= currCol) { A[currRow][currCol] += Aloc[i * nLoc + j]; }
        }
        B[currRow] += Bloc[i];
    }
}

double femBandSystemGet(femBandSystem *system, int myRow, int myCol)
{
    return (myCol >= myRow && myCol < myRow + system->band) ? system->A[myRow][myCol] : 0.0;
}

void femBandSystemConstrain(femBandSystem *system, int node, double value)
{
    double **A, *B;
    int i, size, band;

    A = system->A;
    B = system->B;
    size = system->size;
    band = system->band;

    for (i = 0; i < size; i++)
    {
        if (i >= node && i < node + band) { A[node][i] = (i == node) ? 1 : 0; }
        else { B[i] -= value * A[i][node]; A[i][node] = 0; }
    }
    A[node][node] = 1;
    B[node] = value;
}

double *femBandSystemEliminate(femBandSystem *system)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = system->A;
    B    = system->B;
    size = system->size;
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
    if (mapper == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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

    // TODO : Check if this is correct
    for (int i = 0; i < nNodes; i++) { theMesh->elem[mapper[i]] = i; }

    free(mapper); 
}

int femMeshComputeBand(femMesh *theMesh)
{

    /*
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
    */
    int maxNum, minNum, nodeNum, elemNum;
    int band = 0;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++)
    {   
        maxNum = INT_MIN;
        minNum = INT_MAX;

        for (int j = 0; j < theMesh->nLocalNode; j++)
        {
            elemNum = theMesh->elem[iElem * theMesh->nLocalNode + j];
            // TODO : Check if this is correct
            nodeNum = theMesh->elem[elemNum];

            maxNum = (nodeNum > maxNum) ? nodeNum : maxNum;
            minNum = (nodeNum < minNum) ? nodeNum : minNum;
        }
        if (band < maxNum - minNum) { band = maxNum - minNum; }
    }
    return ++band;
}


/*********************************************************************************************************************/
/*************** ITERATIVE SOLVER ***** ITERATIVE SOLVER ***** ITERATIVE SOLVER ***** ITERATIVE SOLVER ***************/
/*********************************************************************************************************************/

femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *iterativeSolver = malloc(sizeof(femIterativeSolver));
    if (iterativeSolver == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    femIterativeSolverAlloc(iterativeSolver, size);
    femIterativeSolverInit(iterativeSolver);
    return iterativeSolver;
}

void femIterativeSolverFree(femIterativeSolver *iterativeSolver)
{
    free(iterativeSolver->R);
    free(iterativeSolver);
}

void femIterativeSolverAlloc(femIterativeSolver *iterativeSolver, int size)
{
    iterativeSolver->R = malloc(sizeof(double) * size * 4);
    if (iterativeSolver->R == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    iterativeSolver->D = iterativeSolver->R + size;       
    iterativeSolver->S = iterativeSolver->R + size * 2;       
    iterativeSolver->X = iterativeSolver->R + size * 3;       
    iterativeSolver->size = size;
}

void femIterativeSolverInit(femIterativeSolver *iterativeSolver)
{
    iterativeSolver->iter  = 0;
    iterativeSolver->error = 10.0e+12;
    for (int i = 0 ; i < iterativeSolver->size * 4 ; i++) { iterativeSolver->R[i] = 0; }        
}
 
void femIterativeSolverPrint(femIterativeSolver *iterativeSolver)
{
    double *R;
    int size;

    R    = iterativeSolver->R;
    size = iterativeSolver->size;

    for (int i = 0; i < size; i++) { printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *iterativeSolver)
{
    if (iterativeSolver->iter == 1) { printf("\n    Iterative solver \n"); }
    printf("    Iteration %4d : %14.7e\n",iterativeSolver->iter, iterativeSolver->error);
}

int femIterativeSolverConverged(femIterativeSolver *iterativeSolver)
{
    int  testConvergence = 0;
    if (iterativeSolver->iter  > 3000)    { testConvergence = -1; }
    if (iterativeSolver->error < 10.0e-6) { testConvergence = 1; }
    return testConvergence;
}

void femIterativeSolverAssemble(femIterativeSolver *iterativeSolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    for (int i = 0; i < nLoc; i++)
    {
        int row = map[i];
        iterativeSolver->R[row] -= Bloc[i];
        for(int j = 0; j < nLoc; j++) { iterativeSolver->R[row] += Aloc[i * nLoc + j] * Uloc[j]; }
    }
}

double *femIterativeSolverEliminate(femIterativeSolver *iterativeSolver)
{
    iterativeSolver->iter++;
    double error = 0.0;
    for (int i = 0; i < iterativeSolver->size; i++)
    {
        error += iterativeSolver->R[i] * iterativeSolver->R[i];
        iterativeSolver->X[i] = - iterativeSolver->R[i] / 5.0; 
        iterativeSolver->R[i] = 0.0;
    }
    iterativeSolver->error = sqrt(error);
    return iterativeSolver->X;
}


/*********************************************************************************************************************/
/******* ABSTRACTION : FEMSOLVER ********** ABSTRACTION : FEMSOLVER ********** ABSTRACTION : FEMSOLVER ***************/
/*********************************************************************************************************************/

femSolver *femSolverCreate(int sizeLoc)
{
    femSolver *solver = malloc(sizeof(femSolver));
    if (solver == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    solver->local = femFullSystemCreate(sizeLoc);
    return solver;
}

femSolver *femSolverFullCreate(int size, int sizeLoc)
{
    femSolver *mySolver = femSolverCreate(sizeLoc);
    mySolver->type = FEM_FULL;
    mySolver->solver = (femSolver *) femFullSystemCreate(size);
    return mySolver;
}

femSolver *femSolverBandCreate(int size, int sizeLoc, int band)
{
    femSolver *mySolver = femSolverCreate(sizeLoc);
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *) femBandSystemCreate(size,band);
    return mySolver;
}

femSolver *femSolverIterativeCreate(int size, int sizeLoc)
{
    femSolver *mySolver = femSolverCreate(sizeLoc);
    mySolver->type = FEM_ITER;
    mySolver->solver = (femSolver *) femIterativeSolverCreate(size);
    return mySolver;
}

void femSolverFree(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemFree((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemFree((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverFree((femIterativeSolver *) mySolver->solver); break;
        default :       Error("Unexpected solver type");
    }
    free(mySolver->local);
    free(mySolver);
}

void femSolverInit(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemInit((femFullSystem *) mySolver->solver); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *) mySolver->solver); break;
        case FEM_ITER : femIterativeSolverInit((femIterativeSolver *) mySolver->solver); break;
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
        case FEM_ITER : value = (i == j); break;
        default :       Error("Unexpected solver type");
    }
    return value;
}

void femSolverPrint(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemPrint((femFullSystem *) mySolver->solver); break;
        case FEM_BAND : femBandSystemPrint((femBandSystem *) mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrint((femIterativeSolver *) mySolver->solver); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverPrintInfos(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemPrintInfos((femFullSystem *) mySolver->solver); break;
        case FEM_BAND : femBandSystemPrintInfos((femBandSystem *) mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrintInfos((femIterativeSolver *) mySolver->solver); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverAssemble(femSolver *mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *) mySolver->solver, Aloc, Bloc, map, nLoc); break;
        case FEM_BAND : femBandSystemAssemble((femBandSystem *) mySolver->solver, Aloc, Bloc, map, nLoc); break;
        case FEM_ITER : femIterativeSolverAssemble((femIterativeSolver *) mySolver->solver, Aloc, Bloc, Uloc, map, nLoc); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverSystemConstrain(femSolver *mySolver, double node, double value)
{
    switch (mySolver->type)
    {
        case FEM_FULL : femFullSystemConstrain((femFullSystem *) mySolver->solver, node, value); break;
        case FEM_BAND : femBandSystemConstrain((femBandSystem *) mySolver->solver, node, value); break;
        case FEM_ITER : femFullSystemConstrain((femFullSystem *) mySolver->solver, node, value); break;
        default :       Error("Unexpected solver type");
    }
}

void femSolverConstrain(femSolver *solver, int node, double value)
{
    switch (solver->type)
    {
        case FEM_FULL : femFullSystemConstrain((femFullSystem *) solver->solver, node, value); break;
        case FEM_BAND : femBandSystemConstrain((femBandSystem *) solver->solver, node, value); break;
        case FEM_ITER : break;
        default :       Error("Unexpected solver type");
    }
}

double *femSolverEliminate(femSolver *mySolver)
{
    switch (mySolver->type)
    {
        case FEM_FULL : return femFullSystemEliminate((femFullSystem *) mySolver->solver); break;
        case FEM_BAND : return femBandSystemEliminate((femBandSystem *) mySolver->solver); break;
        case FEM_ITER : return femIterativeSolverEliminate((femIterativeSolver *) mySolver->solver); break;
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
        case FEM_ITER : testConvergence = femIterativeSolverConverged((femIterativeSolver *) mySolver->solver); break;
        default :       Error("Unexpected solver type");
    }
    return testConvergence;
}