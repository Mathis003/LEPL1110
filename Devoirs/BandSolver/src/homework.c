#include"fem.h"

// Strip : BEGIN

// Usefull for the renumbering (along x and y directions) method
double *posMeshNodes;

int comparPosNode(const void *a, const void *b)
{
    const int *nodePos_a = (const int *) a;
    const int *nodePos_b = (const int *) b;
    
    return posMeshNodes[*nodePos_a] - posMeshNodes[*nodePos_b];
}
// Strip : END

#ifndef NORENUMBER 

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;

    // Strip : BEGIN

    int nNodes = theMesh->nodes->nNodes;
    int *vect_unit = (int *) malloc(nNodes * sizeof(int));
    for (int i = 0; i < nNodes; i++) { vect_unit[i] = i; }

    switch (renumType)
    {        
        case FEM_NO :
            break;

        // Sorting the nodes along the x-direction
        case FEM_XNUM :
            posMeshNodes = theMesh->nodes->X;
            qsort(vect_unit, nNodes, sizeof(int), &comparPosNode);
            break;

        // Sorting the nodes along the y-direction
        case FEM_YNUM :
            posMeshNodes = theMesh->nodes->Y;
            qsort(vect_unit, nNodes, sizeof(int), &comparPosNode);
            break;    

        default : Error("Unexpected renumbering option"); }

    // Renumbering the nodes
    for (i = 0; i < nNodes; i++) { theMesh->nodes->number[vect_unit[i]] = i; }

    // Free the memory
    free(vect_unit);

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
        maxNum = INT32_MIN;
        minNum = INT32_MAX;

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

    /*
    This is the same function as in the full solver,
    but with a limit on the number of columns to eliminate
    limit = MIN(size, band + idx)
    */

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