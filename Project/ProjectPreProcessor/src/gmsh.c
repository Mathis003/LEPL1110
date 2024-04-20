#include "gmsh.h"

/**********************************/
/******* Gmsh functions *******/
/**********************************/

double geoSizeDefault(double x, double y) { return theGeometry.defaultSize; }

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data) { return theGeometry.geoSize(x, y); }

void geoInitialize(void)
{
    int ierr;
    theGeometry.geoSize = geoSizeDefault;
    gmshInitialize(0, NULL, 1, 0, &ierr);
    ErrorGmsh(ierr);
    gmshModelAdd("MyGeometry", &ierr);
    ErrorGmsh(ierr);
    gmshModelMeshSetSizeCallback(geoGmshSize, NULL, &ierr);
    ErrorGmsh(ierr);
    theGeometry.theNodes    = NULL;
    theGeometry.theElements = NULL;
    theGeometry.theEdges    = NULL;
    theGeometry.theDomains  = NULL;
    theGeometry.nDomains    = 0;
}

void geoFinalize(void)
{
    int ierr;
    geoFree();
    gmshFinalize(&ierr);
    ErrorGmsh(ierr);
}

void geoMeshImport(femDiscreteType discreteType)
{
    int ierr;

    /* Importing nodes */

    size_t nNode, n, m, *node;
    double *xyz, *trash;
    gmshModelMeshGetNodes(&node, &nNode, &xyz, &n, &trash, &m, -1, -1, 0, 0, &ierr); ErrorGmsh(ierr);

    femNodes *theNodes = (femNodes *) malloc(sizeof(femNodes));
    if (theNodes == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->nNodes = nNode;
    theNodes->X = (double *) malloc(sizeof(double) * (theNodes->nNodes));
    theNodes->Y = (double *) malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    if (theNodes->Y == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        theNodes->X[i] = xyz[3 * node[i] - 3];
        theNodes->Y[i] = xyz[3 * node[i] - 2];
    }
    theNodes->number = (int *) malloc(sizeof(int) * (theNodes->nNodes));
    if (theNodes->number == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < theNodes->nNodes; i++) { theNodes->number[i] = i; }

    theGeometry.theNodes = theNodes;
    
    gmshFree(node);
    gmshFree(xyz);
    gmshFree(trash);
    printf("Geo     : Importing %d nodes \n", theGeometry.theNodes->nNodes);

    /* Importing elements */

    // Triangles
    size_t nElem, *elem;
    int elementType = 2;
    int localNode = 3;
    if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { elementType = 9; localNode = 6; }
    gmshModelMeshGetElementsByType(elementType, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    if (nElem != 0)
    {
        femMesh *theElements = (femMesh *) malloc(sizeof(femMesh));
        if (theElements == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theElements->nLocalNode = localNode;
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = (int *) malloc(sizeof(int) * localNode * theElements->nElem);
        if (theElements->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theElements->nElem; i++)
        {
            for (int j = 0; j < localNode; j++) { theElements->elem[localNode * i + j] = node[localNode * i + j] - 1; }
        }
        theGeometry.theElements = theElements;
        gmshFree(node);
        gmshFree(elem);
        printf("Geo     : Importing %d triangles \n", theElements->nElem);
    }

    // Quads
    int nElemTriangles = nElem;
    elementType = 3;
    if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { elementType = 3; }
    gmshModelMeshGetElementsByType(elementType, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    if (nElem != 0 && nElemTriangles != 0) { Error("Cannot consider hybrid geometry with triangles and quads :-("); }

    if (nElem != 0)
    {
        femMesh *theElements = (femMesh *) malloc(sizeof(femMesh));
        if (theElements == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

        theElements->nLocalNode = 4;
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { theElements->nLocalNode = 9; }
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = (int *) malloc(sizeof(int) * theElements->nLocalNode * theElements->nElem);
        if (theElements->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theElements->nElem; i++)
        {
            for (int j = 0; j < theElements->nLocalNode; j++) { theElements->elem[theElements->nLocalNode * i + j] = node[theElements->nLocalNode * i + j] - 1; }
        }        
        theGeometry.theElements = theElements;
        gmshFree(node);
        gmshFree(elem);
        printf("Geo     : Importing %d quads \n", theElements->nElem);
    }

    // Compute node renumbering
    femMesh *theElements = theGeometry.theElements;
    int *connectedNodes = (int *) calloc(theNodes->nNodes, sizeof(int));
    if (connectedNodes == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (int iElem = 0; iElem < theElements->nElem; iElem++)
    {
        for (int i = 0; i < theElements->nLocalNode; i++) { connectedNodes[theElements->elem[iElem * theElements->nLocalNode + i]] = 1; }
    }

    int *nodeRenumber = (int *) malloc(theNodes->nNodes * sizeof(int));
    if (nodeRenumber == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    int countNodes = 0;
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        if (connectedNodes[i]) { nodeRenumber[i] = countNodes; countNodes++; }
        else                   { nodeRenumber[i] = -2147483648; }
    }

    // Condensing nodes
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        if (nodeRenumber[i] < 0) { continue; }
        theNodes->X[nodeRenumber[i]] = theNodes->X[i];
        theNodes->Y[nodeRenumber[i]] = theNodes->Y[i];
    }
    theNodes->nNodes = countNodes;
    theNodes->X = (double *) realloc(theNodes->X, sizeof(double) * (theNodes->nNodes));
    theNodes->Y = (double *) realloc(theNodes->Y, sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    if (theNodes->Y == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    // Renumbering elements
    int nLocalNode = theElements->nLocalNode;
    for (int i = 0; i < theElements->nElem; i++)
    {
        for (int j = 0; j < nLocalNode; j++) { theElements->elem[nLocalNode * i + j] = nodeRenumber[theElements->elem[nLocalNode * i + j]]; }
    }

    elementType = 1;
    localNode = 2;
    if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { elementType = 8; localNode = 3; }
    gmshModelMeshGetElementsByType(elementType, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);

    femMesh *theEdges = (femMesh *) malloc(sizeof(femMesh));
    if (theEdges == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theEdges->nLocalNode = localNode;
    theEdges->nElem = nElem;
    theEdges->nodes = theNodes;
    theEdges->elem  = (int *) malloc(sizeof(int) * localNode * theEdges->nElem);
    if (theEdges->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    int countEdges = 0;
    int *connectedEdges = (int *) calloc(nElem, sizeof(int));
    int *edgeRenumber   = (int *) malloc(nElem * sizeof(int));
    if (connectedEdges == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    if (edgeRenumber == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (int i = 0; i < nElem; i++)
    {
        int map[3];
        for (int j = 0; j < theEdges->nLocalNode; j++)
        {
            map[j] = node[localNode * i + j] - 1;
            if(connectedNodes[map[j]] == 0) { edgeRenumber[j] = -2147483648; continue; }
            theEdges->elem[localNode * countEdges + j] = nodeRenumber[map[j]];
        }
        edgeRenumber[i] = countEdges;
        connectedEdges[i] = 1;
        countEdges++;
    }

    theEdges->nElem = countEdges;
    theEdges->elem = (int *) realloc(theEdges->elem, sizeof(int) * theEdges->nLocalNode * theEdges->nElem);
    if (theEdges->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theGeometry.theEdges = theEdges;
    int shiftEdges = elem[0];
    
    gmshFree(node);
    gmshFree(elem);
    printf("Geo     : Importing %d edges \n", theEdges->nElem);

    /* Importing 1D entities */

    int *dimTags;
    gmshModelGetEntities(&dimTags, &n, 1, &ierr);
    ErrorGmsh(ierr);
    theGeometry.nDomains = n / 2;
    theGeometry.theDomains = (femDomain **) malloc(sizeof(femDomain *) * n / 2);
    if (theGeometry.theDomains == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    printf("Geo     : Importing %d entities \n", theGeometry.nDomains);

    for (int i = 0; i < n / 2; i++)
    {
        int dim = dimTags[2 * i + 0];
        int tag = dimTags[2 * i + 1];
        femDomain *theDomain =(femDomain *) malloc(sizeof(femDomain));
        if (theDomain == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theGeometry.theDomains[i] = theDomain;
        theDomain->mesh = theEdges;
        sprintf(theDomain->name, "Entity %d ", tag - 1);

        int *elementType;
        size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags;
        gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
        theDomain->nElem = nElementTags[0];
        theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
        if (theDomain->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        int nElemCount = 0;
        for (int j = 0; j < theDomain->nElem; j++)
        {
            int mapEdge = edgeRenumber[elementTags[0][j] - shiftEdges];
            if (mapEdge < 0) { continue; }
            theDomain->elem[nElemCount] = mapEdge;
            nElemCount++;
        }
        theDomain->nElem = nElemCount;
        
        printf("Geo     : Entity %d : %d elements \n", i, theDomain->nElem);
        gmshFree(nElementTags);
        gmshFree(nNodesTags);
        gmshFree(elementTags);
        gmshFree(nodesTags);
        gmshFree(elementType);
    }
    gmshFree(dimTags);

    free(connectedNodes); connectedNodes = NULL;
    free(nodeRenumber); nodeRenumber = NULL;
    free(connectedEdges); connectedEdges = NULL;
    free(edgeRenumber); edgeRenumber = NULL;

    // Filter out empty domains
    int countDomains = 0;
    for (int i = 0; i < theGeometry.nDomains; i++)
    {
        if (theGeometry.theDomains[i]->nElem != 0) { theGeometry.theDomains[countDomains] = theGeometry.theDomains[i]; countDomains++; }
        else { free(theGeometry.theDomains[i]->elem); theGeometry.theDomains[i]->elem = NULL; free(theGeometry.theDomains[i]); theGeometry.theDomains[i] = NULL; }
    }
    theGeometry.nDomains = countDomains;

    return;
}

void femErrorGmsh(int ierr, int line, char *file)
{ 
    if (ierr == 0) { return; }
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n", file, line, ierr);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    gmshFinalize(NULL);
    exit(69);
}