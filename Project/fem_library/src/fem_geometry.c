#include "../include/fem_geometry.h"

femGeo theGeometry;

femGeo *geoGetGeometry(void) { return &theGeometry; }

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
    theGeometry.nDomains    = 0;
    theGeometry.theDomains  = NULL;
}

void geoFinalize(void)
{
    int ierr;
    geoFree();
    gmshFinalize(&ierr);
    ErrorGmsh(ierr);
}

void geoSetSizeCallback(double (*geoSize)(double x, double y)) { theGeometry.geoSize = geoSize; }

void geoMeshImport(void)
{
    int ierr;

    /* Importing nodes */

    size_t nNode, n, m, *node;
    double *xyz, *trash;
    //gmshModelMeshRenumberNodes(&ierr); ErrorGmsh(ierr);
    gmshModelMeshGetNodes(&node, &nNode, &xyz, &n, &trash, &m, -1, -1, 0, 0, &ierr); ErrorGmsh(ierr);

    femNodes *theNodes = malloc(sizeof(femNodes));
    if (theNodes == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->nNodes = nNode;
    theNodes->X = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->Y = malloc(sizeof(double) * (theNodes->nNodes));
    if (theNodes->Y == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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

    /* Importing elements (pas super joli : a ameliorer pour eviter la triple copie) */

    // Triangles
    size_t nElem, *elem;
    gmshModelMeshGetElementsByType(2, &elem, &nElem, &node, &nNode, -1, 0, 1, &ierr);
    ErrorGmsh(ierr);
    if (nElem != 0)
    {
        femMesh *theElements = malloc(sizeof(femMesh));
        if (theElements == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theElements->nLocalNode = 3;
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = malloc(sizeof(int) * 3 * theElements->nElem);
        if (theElements->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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
        if (theElements == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theElements->nLocalNode = 4;
        theElements->nodes = theNodes;
        theElements->nElem = nElem;
        theElements->elem = malloc(sizeof(int) * 4 * theElements->nElem);
        if (theElements->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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
    if (connectedNodes == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (int iElem = 0; iElem < theElements->nElem; iElem++)
    {
        for (int i = 0; i < theElements->nLocalNode; i++)
        {
            connectedNodes[theElements->elem[iElem*theElements->nLocalNode+i]] = 1;
        }
    }

    int *nodeRenumber = malloc(theNodes->nNodes * sizeof(int));
    if (nodeRenumber == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    int countNodes = 0;
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        if (connectedNodes[i]) { nodeRenumber[i] = countNodes; countNodes++; }
        else                   { nodeRenumber[i] = -2147483648; }
    }

    // condensing nodes
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        if (nodeRenumber[i] < 0) { continue; }
        theNodes->X[nodeRenumber[i]] = theNodes->X[i];
        theNodes->Y[nodeRenumber[i]] = theNodes->Y[i];
    }
    theNodes->nNodes = countNodes;
    theNodes->X = realloc(theNodes->X, sizeof(double) * (theNodes->nNodes));
    if (theNodes->X == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theNodes->Y = realloc(theNodes->Y, sizeof(double) * (theNodes->nNodes));
    if (theNodes->Y == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    // renumbering elements
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
    if (theEdges == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theEdges->nLocalNode = 2;
    theEdges->nElem = nElem;
    theEdges->nodes = theNodes;
    theEdges->elem = malloc(sizeof(int) * 2 * theEdges->nElem);
    if (theEdges->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    int countEdges = 0;
    int *connectedEdges = calloc(nElem, sizeof(int));
    if (connectedEdges == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    int *edgeRenumber = malloc(nElem * sizeof(int));
    if (edgeRenumber == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

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
    if (theEdges->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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
    if (theGeometry.theDomains == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    printf("Geo     : Importing %d entities \n", theGeometry.nDomains);

    for (int i = 0; i < n / 2; i++)
    {
        int dim = dimTags[2 * i + 0];
        int tag = dimTags[2 * i + 1];
        femDomain *theDomain = malloc(sizeof(femDomain));
        if (theDomain == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        theGeometry.theDomains[i] = theDomain;
        theDomain->mesh = theEdges;
        sprintf(theDomain->name, "Entity %d ", tag - 1);

        int *elementType;
        size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags;
        gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
        theDomain->nElem = nElementTags[0];
        theDomain->elem = malloc(sizeof(int) * 2 * theDomain->nElem);
        if (theDomain->elem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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

    // filter out empty domains
    int countDomains = 0;
    for (int i = 0; i < theGeometry.nDomains; i++)
    {
        if (theGeometry.theDomains[i]->nElem != 0) { theGeometry.theDomains[countDomains] = theGeometry.theDomains[i]; countDomains++; }
        else { free(theGeometry.theDomains[i]->elem); free(theGeometry.theDomains[i]); theGeometry.theDomains[i] = NULL; }
    }
    theGeometry.nDomains = countDomains;

    return;
}

double geoSizeDefault(double x, double y) { return theGeometry.h; }

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data) { return theGeometry.geoSize(x, y); }

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

void femErrorGmsh(int ierr, int line, char *file)                                  
{ 
    if (ierr == 0)  return;
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n", file, line, ierr);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    gmshFinalize(NULL);                                        
    exit(69);                                                 
}