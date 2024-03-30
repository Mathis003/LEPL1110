#include "../include/fem_solver.h"


/*********************************************************************************************************************/
/******* FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER ***** FULL SOLVER *******/
/*********************************************************************************************************************/


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
    double *elem = malloc(sizeof(double) * size * (size + 1));
    if (elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    mySystem->A = malloc(sizeof(double *) * size);
    if (mySystem->A == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    mySystem->B = elem;
    mySystem->A[0] = elem + size;
    mySystem->size = size;
    for (int i = 1; i < size; i++) { mySystem->A[i] = mySystem->A[i - 1] + size; }
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int size = mySystem->size;
    for (int i = 0; i < size * (size + 1); i++) { mySystem->B[i] = 0; }
}

void femFullSystemAssemble(femFullSystem *mySystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    for (int i = 0; i < nLoc; i++)
    {
        for (int j = 0; j < nLoc; j++) { mySystem->A[map[i]][map[j]] += Aloc[i * nLoc+ j]; }
        mySystem->B[map[i]] += Bloc[i];
    }
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

double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol)
{
    return(myFullSystem->A[myRow][myCol]); 
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

void femFullSystemPrintInfos(femFullSystem *mySystem)
{
    int size = mySystem->size;
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
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    if (myBandSystem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    myBandSystem->B = malloc(sizeof(double) * size * (band + 1));
    if (myBandSystem->B == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    myBandSystem->A = malloc(sizeof(double *) * size);
    if (myBandSystem->A == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }     
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    for (int i = 1 ; i < size ; i++) { myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1; }
    femBandSystemInit(myBandSystem);
    return myBandSystem;
}
 
void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}
 
void femBandSystemInit(femBandSystem *myBandSystem)
{
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (int i = 0 ; i < size * (band + 1) ; i++) { myBandSystem->B[i] = 0; }        
}

void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int    band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

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
  
void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int) sizeof(double) * size * (band + 1));     
}

double *posMeshNodes;
int comparPosNode(const void *a, const void *b)
{
    const int *nodePos_a = (const int *) a;
    const int *nodePos_b = (const int *) b;
    
    double diff = posMeshNodes[*nodePos_a] - posMeshNodes[*nodePos_b];
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
        // Sorting the nodes along the x-direction
        case FEM_XNUM :
            posMeshNodes = theMesh->nodes->X;
            qsort(mapper, nNodes, sizeof(int), comparPosNode);
            break;
        // Sorting the nodes along the y-direction
        case FEM_YNUM :
            posMeshNodes = theMesh->nodes->Y;
            qsort(mapper, nNodes, sizeof(int), comparPosNode);
            break;    

        default : Error("Unexpected renumbering option");
    }

    // Renumbering the nodes
    for (int i = 0; i < nNodes; i++)
    {
        theMesh->elem[mapper[i]] = i; // TODO : Check if this is correct
        // theMesh->nodes->number[mapper[i]] = i; // Version in Homework 'BandSolver'
    }

    // Free the memory
    free(mapper); 
}

int femMeshComputeBand(femMesh *theMesh)
{
    int maxNum, minNum, nodeNum, elemNum;
    int myBand = 0;

    // Looping through all elements
    for (int iElem = 0; iElem < theMesh->nElem; iElem++)
    {   
        // Initializing the max and min numbers
        maxNum = INT_MIN;
        minNum = INT_MAX;

        // Looping through all nodes of the element
        for (int j = 0; j < theMesh->nLocalNode; j++)
        {
            // Getting the node number
            elemNum = theMesh->elem[iElem * theMesh->nLocalNode + j];
            // nodeNum = theMesh->nodes->number[elemNum]; // Version in Homework 'BandSolver'
            nodeNum = theMesh->elem[elemNum]; // TODO : Check if this is correct

            // Updating the max and min node numbers
            maxNum = (nodeNum > maxNum) ? nodeNum : maxNum;
            minNum = (nodeNum < minNum) ? nodeNum : minNum;
        }

        // Updating the band
        if (myBand < maxNum - minNum) { myBand = maxNum - minNum; }
    }
    return ++myBand;
}

void femBandSystemAssemble(femBandSystem *myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int currRow, currCol;
    // Looping through the local matrix elements
    for (int i = 0; i < nLoc; i++)
    {
        currRow = map[i];
        for (int j = 0; j < nLoc; j++)
        {
            currCol = map[j];
            // Assemble the local matrix A_e with the global matrix A
            // The condition ensures that we only assemble the upper part of the matrix
            if (currRow <= currCol) { myBandSystem->A[currRow][currCol] += Aloc[i * nLoc + j]; }
        }
        // Assemble the local vector B_e with the global vector B
        myBandSystem->B[currRow] += Bloc[i];
    }
}

double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol]; 
    return(value);
}

double *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Gauss elimination */
    for (k = 0; k < size; k++)
    {
        if (fabs(A[k][k]) <= 1e-8) { Error("Cannot eliminate with such a pivot.\n"); }
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
    return myBand->B;
}



/*********************************************************************************************************************/
/*************** ITERATIVE SOLVER ***** ITERATIVE SOLVER ***** ITERATIVE SOLVER ***** ITERATIVE SOLVER ***************/
/*********************************************************************************************************************/

femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *mySolver = malloc(sizeof(femIterativeSolver));
    if (mySolver == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    mySolver->R = malloc(sizeof(double) * size * 4);
    if (mySolver->R == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    mySolver->D = mySolver->R + size;       
    mySolver->S = mySolver->R + size * 2;       
    mySolver->X = mySolver->R + size * 3;       
    mySolver->size = size;
    femIterativeSolverInit(mySolver);
    return mySolver;
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    mySolver->iter  = 0;
    mySolver->error = 10.0e+12;
    for (int i = 0 ; i < mySolver->size * 4 ; i++) { mySolver->R[i] = 0; }        
}
 
void femIterativeSolverPrint(femIterativeSolver *mySolver)
{
    double *R;
    int size;

    R    = mySolver->R;
    size = mySolver->size;

    for (int i = 0; i < size; i++) { printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *mySolver)
{
    if (mySolver->iter == 1) { printf("\n    Iterative solver \n"); }
    printf("    Iteration %4d : %14.7e\n",mySolver->iter,mySolver->error);
}

int femIterativeSolverConverged(femIterativeSolver *mySolver)
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)    { testConvergence = -1; }
    if (mySolver->error < 10.0e-6) { testConvergence = 1; }
    return testConvergence;
}

void femIterativeSolverAssemble(femIterativeSolver *mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    for (int i = 0; i < nLoc; i++)
    {
        int myRow = map[i];
        mySolver->R[myRow] -= Bloc[i];
        for(int j = 0; j < nLoc; j++) { mySolver->R[myRow] += Aloc[i*nLoc+j]*Uloc[j]; }
    }
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0;
    for (int i = 0; i < mySolver->size; i++)
    {
        error += mySolver->R[i] * mySolver->R[i];
        mySolver->X[i] = - mySolver->R[i] / 5.0; 
        mySolver->R[i] = 0.0;
    }
    mySolver->error = sqrt(error);
    return mySolver->X;
}


/*********************************************************************************************************************/
/******* ABSTRACTION : FEMSOLVER ********** ABSTRACTION : FEMSOLVER ********** ABSTRACTION : FEMSOLVER ***************/
/*********************************************************************************************************************/

femSolver *femSolverCreate(int sizeLoc)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    if (mySolver == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    mySolver->local = femFullSystemCreate(sizeLoc);
    return mySolver;
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
    femFullSystemFree(mySolver->local);
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

double *femSolverEliminate(femSolver *mySolver)
{
    double *soluce;
    switch (mySolver->type)
    {
        case FEM_FULL : soluce = femFullSystemEliminate((femFullSystem *) mySolver->solver); break;
        case FEM_BAND : soluce = femBandSystemEliminate((femBandSystem *) mySolver->solver); break;
        case FEM_ITER : soluce = femIterativeSolverEliminate((femIterativeSolver *) mySolver->solver); break;
        default :       Error("Unexpected solver type");
    }
    return soluce;
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