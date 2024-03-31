#include "../include/fem_gmsh.h"

double geoSizeDefault(double x, double y) { return theGeometry.h; }

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

void geoMeshImport(void)
{
    int ierr;

    /* Importing nodes */

    size_t nNode, n, m, *node;
    double *xyz, *trash;
    //gmshModelMeshRenumberNodes(&ierr); ErrorGmsh(ierr);
    gmshModelMeshGetNodes(&node, &nNode, &xyz, &n, &trash, &m, -1, -1, 0, 0, &ierr); ErrorGmsh(ierr);

    femNodes *theNodes = malloc(sizeof(femNodes));
    if (theNodes == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->nNodes = nNode;
    theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->Y == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        theNodes->X[i] = xyz[3 * node[i] - 3];
        theNodes->Y[i] = xyz[3 * node[i] - 2];
    }
    theGeometry.theNodes = theNodes;
    
    gmshFree(node);
    gmshFree(xyz);
    gmshFree(trash);
    printf("Geo     : Importing %d nodes \n", theGeometry.theNodes->nNodes);

    /* Importing elements */

    // Triangles
    size_t nElem, *elem;
    gmshModelMeshGetElementsByType(2, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    if (nElem != 0)
    {
        femMesh *theElements = malloc(sizeof(femMesh));
        if (theElements == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theElements->nLocalNode = 3;
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = malloc(sizeof(int) * 3 * theElements->nElem);
        if (theElements->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theElements->nElem; i++)
        {
            for (int j = 0; j < theElements->nLocalNode; j++)
            {
                theElements->elem[3 * i + j] = node[3 * i + j] - 1;
            }
        }
        theGeometry.theElements = theElements;
        gmshFree(node);
        gmshFree(elem);
        printf("Geo     : Importing %d triangles \n", theElements->nElem);
    }

    // Quads
    int nElemTriangles = nElem;
    gmshModelMeshGetElementsByType(3, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    if (nElem != 0 && nElemTriangles != 0) { Error("Cannot consider hybrid geometry with triangles and quads :-("); }

    if (nElem != 0)
    {
        femMesh *theElements = malloc(sizeof(femMesh));
        if (theElements == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theElements->nLocalNode = 4;
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = malloc(sizeof(int) * 4 * theElements->nElem);
        if (theElements->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        for (int i = 0; i < theElements->nElem; i++)
        {
            for (int j = 0; j < theElements->nLocalNode; j++)
            {
                theElements->elem[4 * i + j] = node[4 * i + j] - 1;
            }
        }        
        theGeometry.theElements = theElements;
        gmshFree(node);
        gmshFree(elem);
        printf("Geo     : Importing %d quads \n", theElements->nElem);
    }

    // Compute node renumbering
    femMesh *theElements = theGeometry.theElements;
    int *connectedNodes = calloc(theNodes->nNodes, sizeof(int));
    if (connectedNodes == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (int iElem = 0; iElem < theElements->nElem; iElem++)
    {
        for (int i = 0; i < theElements->nLocalNode; i++)
        {
            connectedNodes[theElements->elem[iElem*theElements->nLocalNode+i]] = 1;
        }
    }

    int *nodeRenumber = malloc(theNodes->nNodes * sizeof(int));
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
    theNodes->X = realloc(theNodes->X, sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->Y = realloc(theNodes->Y, sizeof(double) * (theNodes->nNodes));
    if (theNodes->Y == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    // Renumbering elements
    int nLocalNode = theElements->nLocalNode;
    for (int i = 0; i < theElements->nElem; i++)
    {
        for (int j = 0; j < nLocalNode; j++)
        {
            theElements->elem[nLocalNode * i + j] = nodeRenumber[theElements->elem[nLocalNode * i + j]];
        }
    }

    gmshModelMeshGetElementsByType(1, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);

    femMesh *theEdges = malloc(sizeof(femMesh));
    if (theEdges == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theEdges->nLocalNode = 2;
    theEdges->nElem = nElem;
    theEdges->nodes = theNodes;
    theEdges->elem = malloc(sizeof(int) * 2 * theEdges->nElem);
    if (theEdges->elem == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    int countEdges = 0;
    int *connectedEdges = calloc(nElem, sizeof(int));
    if (connectedEdges == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    int *edgeRenumber = malloc(nElem * sizeof(int));
    if (edgeRenumber == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (int i = 0; i < nElem; i++)
    {
        int map[2] = {node[2*i+0] - 1, node[2*i+1] - 1};
        if (!connectedNodes[map[0]] || !connectedNodes[map[1]]) { edgeRenumber[i] = -2147483648; continue; }
        for (int j = 0; j < theEdges->nLocalNode; j++)
        {
            connectedEdges[i] = 1;
            theEdges->elem[2 * countEdges + j] = nodeRenumber[map[j]];
        }
        edgeRenumber[i] = countEdges;
        countEdges++;
    }

    theEdges->nElem = countEdges;
    theEdges->elem = realloc(theEdges->elem, sizeof(int) * 2 * theEdges->nElem);
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
    theGeometry.theDomains = malloc(sizeof(femDomain *) * n / 2);
    if (theGeometry.theDomains == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    printf("Geo     : Importing %d entities \n", theGeometry.nDomains);

    for (int i = 0; i < n / 2; i++)
    {
        int dim = dimTags[2 * i + 0];
        int tag = dimTags[2 * i + 1];
        femDomain *theDomain = malloc(sizeof(femDomain));
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

    free(connectedNodes);
    free(nodeRenumber);
    free(connectedEdges);
    free(edgeRenumber);

    // Filter out empty domains
    int countDomains = 0;
    for (int i = 0; i < theGeometry.nDomains; i++)
    {
        if (theGeometry.theDomains[i]->nElem != 0) { theGeometry.theDomains[countDomains] = theGeometry.theDomains[i]; countDomains++; }
        else { free(theGeometry.theDomains[i]->elem); free(theGeometry.theDomains[i]); theGeometry.theDomains[i] = NULL; }
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