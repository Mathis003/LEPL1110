#include "fem.h"
#include <limits.h>

// Strip : BEGIN

// Usefull for the renumbering (along x and y directions) method
double *positionMeshNodes;

int comparPositionNode(const void *a, const void *b)
{
    const int *nodePos_a = (const int *) a;
    const int *nodePos_b = (const int *) b;
    
    double diff = positionMeshNodes[*nodePos_a] - positionMeshNodes[*nodePos_b];
    return (diff < 0) - (diff > 0);
}
// Strip : END

#ifndef NORENUMBER 

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;

    // Strip : BEGIN
    int nNodes = theMesh->nodes->nNodes;
    int *mapper = (int *) malloc(nNodes * sizeof(int));
    if (mapper == NULL) { Error("Memory allocation failed !"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < nNodes; i++) { mapper[i] = i; }

    switch (renumType)
    {
        case FEM_NO :
            break;

        // Sorting the nodes along the x-direction
        case FEM_XNUM :
            positionMeshNodes = theMesh->nodes->X;
            qsort(mapper, nNodes, sizeof(int), comparPositionNode);
            break;

        // Sorting the nodes along the y-direction
        case FEM_YNUM :
            positionMeshNodes = theMesh->nodes->Y;
            qsort(mapper, nNodes, sizeof(int), comparPositionNode);
            break;    

        default : Error("Unexpected renumbering option"); }

    // Renumbering the nodes
    for (i = 0; i < nNodes; i++) { theMesh->nodes->number[mapper[i]] = i; }

    // Free the memory
    free(mapper);

    // Strip : END   
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    // Strip : BEGIN
    int myBand = 0;
    int maxNum, minNum, nodeNum, elemNum;

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
            nodeNum = theMesh->nodes->number[elemNum];

            // Updating the max and min node numbers
            maxNum = (nodeNum > maxNum) ? nodeNum : maxNum;
            minNum = (nodeNum < minNum) ? nodeNum : minNum;
        }

        // Updating the band
        if (myBand < maxNum - minNum) { myBand = maxNum - minNum; }
    }
    
    return ++myBand;
    // Strip : END
}

#endif
#ifndef NOBANDASSEMBLE

void femBandSystemAssemble(femBandSystem *myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    // Strip : BEGIN
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
    // Strip : END
}

#endif
#ifndef NOBANDELIMINATE

double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    // Strip : BEGIN

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
    // Strip : END

    return myBand->B;
}

#endif