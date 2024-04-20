#include "../../Project/src/fem.h"
#include "gmsh.h"

/**********************************/
/********* Gmsh functions *********/
/**********************************/

int createRectangle(double x, double y, double width, double height)
{
    int ierr;
    return gmshModelOccAddRectangle(x, y, 0.0, width, height, -1, 0, &ierr);
}

int createDisk(double xc, double yc, double rx, double ry)
{
    int ierr;
    return gmshModelOccAddDisk(xc, yc, 0.0, rx, ry, -1, NULL, 0, NULL, 0, &ierr);
}

void cutElement(int *mainElement, int *cutElement)
{
    int ierr;
    gmshModelOccCut(mainElement, 2, cutElement, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
}

void fuseElement(int *mainElement, int *fuseElement)
{
    int ierr;
    gmshModelOccFuse(mainElement, 2, fuseElement, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
}

void rotateElement(int *element, double posX, double posY, double angle)
{
    int ierr;
    gmshModelOccRotate(element, 2, posX, posY, 0.0, 0.0, 0.0, 1.0, angle, &ierr);
}


/*************************************/
/********* Position Geometry *********/
/*************************************/

int isRoadWay(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();

    double yBridge = theGeometry->heightPillars + theGeometry->heightBridge;
    double heigthRoadWay = 1.0;

    return ((y >= yBridge - heigthRoadWay) && (y <= yBridge));
}

int isSubRoadWay(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();

    if ((y < theGeometry->heightPillars) || (y > theGeometry->heightSubRoadWay + theGeometry->heightPillars)) { return FALSE; }
    
    if ((x > theGeometry->rxLongArc) && (x < theGeometry->rxLongArc + theGeometry->widthPillars)) { return FALSE; }
    if ((x > theGeometry->rxLongArc + theGeometry->widthPillars + 2 * theGeometry->rxArc) &&
        (x < theGeometry->rxLongArc + 2 * theGeometry->widthPillars + 2 * theGeometry->rxArc))
        { return FALSE; }

    if ((x < -theGeometry->rxLongArc) && (x > -theGeometry->rxLongArc - theGeometry->widthPillars)) { return FALSE; }
    if ((x < -theGeometry->rxLongArc - theGeometry->widthPillars - 2 * theGeometry->rxArc) &&
        (x > -theGeometry->rxLongArc - 2 * theGeometry->widthPillars - 2 * theGeometry->rxArc))
        { return FALSE; }

    return TRUE;
}

int isStayCables(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();

    if ((y <= theGeometry->heightBridge + theGeometry->heightPillars) || (y >= theGeometry->heightBridge + theGeometry->heightPillars + 5 * theGeometry->heightPylons / 3)) { return FALSE; }

    if ((x < theGeometry->rxLongArc + theGeometry->widthPillars / 2 - theGeometry->widthPylons / 2) &&
        (x > - theGeometry->rxLongArc - theGeometry->widthPillars / 2 + theGeometry->widthPylons / 2))
        { return TRUE; }

    if (x > theGeometry->rxLongArc + theGeometry->widthPillars / 2 + theGeometry->widthPylons / 2)  { return TRUE; }
    if (x < -theGeometry->rxLongArc - theGeometry->widthPillars / 2 - theGeometry->widthPylons / 2) { return TRUE; }

    if (y > theGeometry->heightBridge + theGeometry->heightPillars + theGeometry->heightPylons)
    {
        if ((x < theGeometry->rxLongArc + theGeometry->widthPillars / 2 - theGeometry->widthPylons / 4) &&
            (x > -theGeometry->rxLongArc - theGeometry->widthPillars / 2 + theGeometry->widthPylons / 4))
            { return TRUE; }

        if (x > theGeometry->rxLongArc + theGeometry->widthPillars / 2 + theGeometry->widthPylons / 4)  { return TRUE; }
        if (x < -theGeometry->rxLongArc - theGeometry->widthPillars / 2 - theGeometry->widthPylons / 4) { return TRUE; }
    }
    return FALSE;
}

/*****************************/
/********* MATERIALS *********/
/*****************************/

double *getMaterialProperties(char *material)
{
    double *properties = malloc(3 * sizeof(double));
    if (properties == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return NULL; }
    if (strcmp(material, "steel") == 0)
    {
        properties[0] = 200.0e9; // E [Pa]
        properties[1] = 0.3;     // nu [/]
        properties[2] = 7850.0;  // rho [kg/m^3]
    }
    else if (strcmp(material, "reinforced_concrete") == 0)
    {
        properties[0] = 30.0e9; // E [Pa]
        properties[1] = 0.2;    // nu [/]
        properties[2] = 2400.0; // rho [kg/m^3]
    }
    else { Error("Material Unknown."); }
    return properties;
}

char *getMaterials(double x, double y) { return (isSubRoadWay(x, y) == TRUE || isRoadWay(x, y) == TRUE || isStayCables(x, y) == TRUE) ? "steel" : "reinforced_concrete"; }


/********************************************/
/********* Create Bridge components *********/
/********************************************/

void createWindows(femGeometry *theGeometry, int *idWindows)
{
    double width  = theGeometry->widthWindow;
    double height = theGeometry->heightWindow;

    double offsetY = 0.5;
    double x1 = theGeometry->rxLongArc + theGeometry->widthPillars / 4;
    double x2 = theGeometry->rxLongArc + 2 * theGeometry->rxArc + 5 * theGeometry->widthPillars / 4;
    double y = theGeometry->heightPillars + theGeometry->heightSubRoadWay + offsetY;

    idWindows[0] = createRectangle(- x1 - width, y, width, height);
    idWindows[1] = createRectangle(- x2 - width, y, width, height);
    idWindows[2] = createRectangle(x1, y, width, height);
    idWindows[3] = createRectangle(x2, y, width, height);
}

void createPiles(femGeometry *theGeometry, int *idPiles)
{
    double width = theGeometry->widthPiles;
    double height = theGeometry->ryArc;

    double offsetX = theGeometry->widthPiles / 2;
    double y = theGeometry->heightPillars + theGeometry->heightSubRoadWay;

    double x1 = theGeometry->rxLongArc / 4 - offsetX;
    double x2 = 3 * theGeometry->rxLongArc / 4 - offsetX;
    double x3 = theGeometry->rxLongArc + theGeometry->widthPillars + theGeometry->rxArc / 2 - offsetX;
    double x4 = theGeometry->rxLongArc + theGeometry->widthPillars + 3 * theGeometry->rxArc / 2 - offsetX;
    double x5 = theGeometry->widthSpanBridge / 2 - theGeometry->rxArc / 2 - offsetX;

    idPiles[0] = createRectangle(x1, y, width, height);
    idPiles[1] = createRectangle(- x1, y, width, height);
    idPiles[2] = createRectangle(x2, y, width, height);
    idPiles[3] = createRectangle(- x2, y, width, height);
    idPiles[4] = createRectangle(x3, y, width, height);
    idPiles[5] = createRectangle(- x3, y, width, height);
    idPiles[6] = createRectangle(x4, y, width, height);
    idPiles[7] = createRectangle(- x4, y, width, height);
    idPiles[8] = createRectangle(x5, y, width, height);
    idPiles[9] = createRectangle(- x5, y, width, height);
}

void createPillars(femGeometry *theGeometry, int *idPillars)
{
    double width  = theGeometry->widthPillars;
    double height = theGeometry->heightPillars;

    double x1 = theGeometry->rxLongArc;
    double x2 = theGeometry->widthSpanBridge / 2 - theGeometry->rxArc;
    double y = 0.0;
    
    idPillars[0] = createRectangle(- x2, y, width, height);
    idPillars[1] = createRectangle(- x1 - width, y, width, height);
    idPillars[2] = createRectangle(x1, y, width, height);
    idPillars[3] = createRectangle(x2 - width, y, width, height);
}

void createArcs(femGeometry *theGeometry, int *idArcs)
{
    double rxArc = theGeometry->rxArc;
    double ryArc = theGeometry->ryArc;
    double rxLongArc = theGeometry->rxLongArc;
    double ryLongArc = theGeometry->ryLongArc;

    double x1 = theGeometry->rxLongArc + 2 * theGeometry->widthPillars + 3 * theGeometry->rxArc;
    double x2 = theGeometry->rxLongArc + theGeometry->widthPillars;
    double y = theGeometry->heightPillars + theGeometry->heightSubRoadWay;

    idArcs[0] = createDisk(- x1, y, rxArc, ryArc);
    idArcs[1] = createDisk(- x2 - rxArc, y, rxArc, ryArc);
    idArcs[2] = createDisk(0.0, y, rxLongArc, ryLongArc);
    idArcs[3] = createDisk(x2 + rxArc, y, rxArc, ryArc);
    idArcs[4] = createDisk(x1, y, rxArc, ryArc);
}

void createPylons(femGeometry *theGeometry, int *idPylons)
{
    double width  = theGeometry->widthPylons;
    double height = theGeometry->heightPylons;

    double x = theGeometry->rxLongArc + theGeometry->widthPillars / 2;
    double y = theGeometry->heightPillars + theGeometry->heightBridge;

    idPylons[0] = createRectangle(-x - width / 2, y, width, height);
    idPylons[1] = createRectangle(-x - width / 4, y + height, width / 2, 2 * height / 3);
    idPylons[2] = createRectangle(-x - width / 8, y + 5 * height / 3, width / 4, height / 3);
    idPylons[3] = createRectangle(x - width / 2, y, width, height);
    idPylons[4] = createRectangle(x - width / 4, y + height, width / 2, 2 * height / 3);
    idPylons[5] = createRectangle(x - width / 8, y + 5 * height / 3, width / 4, height / 3);
}

void createTopBall(femGeometry *theGeometry, int *idTopBall)
{
    double r = theGeometry->widthPylons / 3;
    double x = theGeometry->rxLongArc + theGeometry->widthPillars / 2;
    double y = theGeometry->heightPillars + theGeometry->heightBridge + 2 * theGeometry->heightPylons + r / 3 + 0.35; 

    idTopBall[0] = createDisk(- x, y, r, r);
    idTopBall[1] = createDisk(x, y, r, r);
}

void createStayCables(femGeometry *theGeometry, int *idStayCables, double *positionX, double *positionY)
{
    double width = theGeometry->widthStayCables;
    double height = theGeometry->heightStayCables;
    double dist = theGeometry->distStayCables;

    double centerPile = theGeometry->rxLongArc + theGeometry->widthPillars / 2;
    const double y  = theGeometry->heightPillars + theGeometry->heightBridge + 5 * theGeometry->heightPylons / 3;
    double x1 = centerPile + theGeometry->widthPylons / 4;
    double x2 = centerPile + theGeometry->widthPylons / 2;

    // Right stay cables on the right pile
    for (int i = 0; i < 4; i++)
    {
        idStayCables[i] = createRectangle(x1, y - i * dist, width, height - i);
        positionX[i] = x1;
        positionY[i] = y - i * dist;
    }
    for (int i = 4; i < 9; i++)
    {
        idStayCables[i] = createRectangle(x2, y - (i + 1) * dist, width, height - (i + 1));
        positionX[i] = x2;
        positionY[i] = y - (i + 1) * dist;
    }

    x1 += width;
    x2 += width;

    // Left stay cables on the left pile
    for (int i = 0; i < 4; i++)
    {
        idStayCables[18 + i] = createRectangle(-x1, y - i * dist, width, height - i);
        positionX[18 + i] = -x1 + width;
        positionY[18 + i] = y - i * dist;
    }
    for (int i = 4; i < 9; i++)
    {
        idStayCables[18 + i] = createRectangle(-x2, y - (i + 1) * dist, width, height - (i + 1));
        positionX[18 + i] = -x2 + width;
        positionY[18 + i] = y - (i + 1) * dist;
    }

    x1 = centerPile - theGeometry->widthPylons / 4 - width;
    x2 = centerPile - theGeometry->widthPylons / 2 - width;

    // Left stay cables on the right pile
    for (int i = 0; i < 4; i++)
    {
        idStayCables[27 + i] = createRectangle(x1, y - i * dist, width, height - i);
        positionX[27 + i] = x1 + width;
        positionY[27 + i] = y - i * dist;
    }
    for (int i = 4; i < 9; i++)
    {
        idStayCables[27 + i] = createRectangle(x2, y - (i + 1) * dist, width, height - (i + 1));
        positionX[27 + i] = x2 + width;
        positionY[27 + i] = y - (i + 1) * dist;
    }

    x1 += width;
    x2 += width;

    // Right stay cables on the left pile
    for (int i = 0; i < 4; i++)
    {
        idStayCables[9 + i] = createRectangle(-x1, y - i * dist, width, height - i);
        positionX[9 + i] = -x1;
        positionY[9 + i] = y - i * dist;
    }
    for (int i = 4; i < 9; i++)
    {
        idStayCables[9 + i] = createRectangle(-x2, y - (i + 1) * dist, width, height - (i + 1));
        positionX[9 + i] = -x2;
        positionY[9 + i] = y - (i + 1) * dist;
    }
}

void createBridgeMainSpan(femGeometry *theGeometry, int *idBridge)
{
    int x = - theGeometry->widthSpanBridge / 2;
    int y = theGeometry->heightPillars;
    int width = theGeometry->widthSpanBridge;
    int height = theGeometry->heightBridge;
    *idBridge = createRectangle(x, y, width, height);
}

void createSubRoadWay(femGeometry *theGeometry, int *idSubRoadWay)
{
    int x = - theGeometry->widthSubRoadWay / 2;
    int y = theGeometry->heightPillars;
    int width = theGeometry->widthSubRoadWay;
    int height = theGeometry->heightSubRoadWay;
    *idSubRoadWay = createRectangle(x, y, width, height);
}

void cutHalfGeometryBySymmetry(femGeometry *theGeometry, int *bridge)
{
    int width = theGeometry->widthSpanBridge / 2;
    int height = theGeometry->heightPillars + theGeometry->heightBridge + 3 * theGeometry->heightPylons;

    int idRectFilterLeft = createRectangle(0.0, 0.0, width, height);
    int *filterLeftRect = malloc(2 * sizeof(int));
    if (filterLeftRect == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    filterLeftRect[0] = 2; filterLeftRect[1] = idRectFilterLeft;

    cutElement(bridge, filterLeftRect);
    
    free(filterLeftRect);
    filterLeftRect = NULL;
}


/*************************************************/
/********* Create Bridge (Generate Mesh) *********/
/*************************************************/

void geoMeshGenerate(femDiscreteType discreteType, int bridgeSimplified)
{
    femGeometry *theGeometry = geoGetGeometry();

    // Define the geometry parameters
    theGeometry->widthSpanBridge  = 62.0;
    theGeometry->heightBridge     = 6.0;
    theGeometry->widthWindow      = 1.5;
    theGeometry->heightWindow     = 0.6;
    theGeometry->widthSubRoadWay  = 62.0;
    theGeometry->heightSubRoadWay = 1.0;
    theGeometry->rxArc            = 5.0;
    theGeometry->ryArc            = 4.0;
    theGeometry->rxLongArc        = 10.0;
    theGeometry->ryLongArc        = 3.0;
    theGeometry->widthPiles       = 0.5;
    theGeometry->widthPillars     = 3.0;
    theGeometry->heightPillars    = 5.0;
    theGeometry->widthPylons      = 2.0;
    theGeometry->heightPylons     = 6.0;
    theGeometry->angleStayCables  = 135 * M_PI / 180;
    theGeometry->widthStayCables  = 0.2;
    theGeometry->heightStayCables = 14.4; // Be carefull (if too long, the stay cables will overlap an arc)
    theGeometry->distStayCables   = 0.8;
    theGeometry->defaultSize      = 1.4;

    theGeometry->geoSize = geoSize;
    theGeometry->getMaterialProperties = getMaterialProperties;
    theGeometry->getMaterials = getMaterials;
    
    int ierr;
    int idBridge, idSubRoadWay;

    int *idWindows    = (int *) malloc(4 * sizeof(int));
    int *idPiles      = (int *) malloc(10 * sizeof(int));
    int *idPillars    = (int *) malloc(4 * sizeof(int));
    int *idArcs       = (int *) malloc(5 * sizeof(int));
    int *idPylons     = (int *) malloc(6 * sizeof(int));
    int *idTopBall    = (int *) malloc(2 * sizeof(int));
    int *idStayCables = (int *) malloc(36 * sizeof(int));
    double *positionStayCablesX = (double *) malloc(36 * sizeof(double));
    double *positionStayCablesY = (double *) malloc(36 * sizeof(double));

    if (idWindows == NULL)           { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idPiles == NULL)             { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idPillars == NULL)           { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idArcs == NULL)              { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idPylons == NULL)            { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idTopBall == NULL)           { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idStayCables == NULL)        { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (positionStayCablesX == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (positionStayCablesY == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }

    int bridge[2], subRoadWay[2], pillars[4][2], pylons[6][2], stayCables[36][2], arcs[5][2], windows[4][2], topBall[2][2], piles[10][2];

    // Create each part of the bridge
    createBridgeMainSpan(theGeometry, &idBridge);
    createPillars(theGeometry, idPillars);

    if (!bridgeSimplified)
    {
        createPylons(theGeometry, idPylons);
        createTopBall(theGeometry, idTopBall);
        createStayCables(theGeometry, idStayCables, positionStayCablesX, positionStayCablesY);
    }
    createArcs(theGeometry, idArcs);
    createWindows(theGeometry, idWindows);
    createSubRoadWay(theGeometry, &idSubRoadWay);
    createPiles(theGeometry, idPiles);

    bridge[0] = 2; bridge[1] = idBridge;
    subRoadWay[0] = 2; subRoadWay[1] = idSubRoadWay;

    for (int i = 0; i < 4; i++)  { pillars[i][0] = 2; pillars[i][1] = idPillars[i]; }
    if (!bridgeSimplified) { for (int i = 0; i < 6; i++)  { pylons[i][0] = 2; pylons[i][1] = idPylons[i]; } }
    for (int i = 0; i < 5; i++)  { arcs[i][0] = 2; arcs[i][1] = idArcs[i]; }
    for (int i = 0; i < 4; i++)  { windows[i][0] = 2; windows[i][1] = idWindows[i]; }
    if (!bridgeSimplified) { for (int i = 0; i < 2; i++)  { topBall[i][0] = 2; topBall[i][1] = idTopBall[i]; } }
    for (int i = 0; i < 10; i++) { piles[i][0] = 2; piles[i][1] = idPiles[i]; }
    if (!bridgeSimplified)
    {
        for (int i = 0; i < 36; i++) { stayCables[i][0] = 2; stayCables[i][1] = idStayCables[i]; }
        for (int i = 0; i < 36; i++) { stayCables[i][0] = 2; stayCables[i][1] = idStayCables[i]; }
    }

    // Cut the arcs and the windows
    for (int i = 0; i < 5; i++) { cutElement(bridge, arcs[i]); }
    for (int i = 0; i < 4; i++) { cutElement(bridge, windows[i]); }

    // Rotate the stay cables
    if (!bridgeSimplified)
    {
        for (int i = 0; i < 18; i++)  { rotateElement(stayCables[i], positionStayCablesX[i], positionStayCablesY[i], -theGeometry->angleStayCables); }
        for (int i = 18; i < 36; i++) { rotateElement(stayCables[i], positionStayCablesX[i], positionStayCablesY[i], theGeometry->angleStayCables); }
    }
    
    // Fuse all the elements together to create the bridge
    fuseElement(bridge, subRoadWay);
    if (!bridgeSimplified)
    {
        fuseElement(pylons[2], topBall[0]);
        fuseElement(pylons[5], topBall[1]);
    }
    for (int i = 0; i < 4; i++)   { fuseElement(bridge, pillars[i]); }
    for (int i = 0; i < 10; i++)  { fuseElement(bridge, piles[i]); }
    if (!bridgeSimplified)
    {
        for (int i = 1; i < 3; i++)   { fuseElement(pylons[0], pylons[i]); }
        for (int i = 4; i < 6; i++)   { fuseElement(pylons[3], pylons[i]); }
        fuseElement(bridge, pylons[0]);
        fuseElement(bridge, pylons[3]);
        for (int i = 0; i < 9; i++)   { fuseElement(bridge, stayCables[i]); }
        for (int i = 27; i < 36; i++) { fuseElement(bridge, stayCables[i]); }
        for (int i = 9; i < 27; i++)  { fuseElement(bridge, stayCables[i]); }
    }

    // Cut the half of the bridge by symmetry
    cutHalfGeometryBySymmetry(theGeometry, bridge);

    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);

    // Generate quads meshing
    if (theGeometry->elementType == FEM_QUAD)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC){ gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    // Generate triangles meshing
    else if (theGeometry->elementType == FEM_TRIANGLE)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    gmshFltkInitialize(&ierr);

    // Free the memory
    free(idWindows);
    free(idPiles);
    free(idPillars);
    free(idArcs);
    free(idPylons);
    free(idTopBall);
    free(idStayCables);
    free(positionStayCablesX);
    free(positionStayCablesY);

    idWindows = NULL;
    idPiles   = NULL;
    idPillars = NULL;
    idArcs    = NULL;
    idPylons  = NULL;
    idTopBall = NULL;
    idStayCables = NULL;
    positionStayCablesX = NULL;
    positionStayCablesY = NULL;
}


/******************************************************/
/********* Name Domains + Boundary Conditions *********/
/******************************************************/

void setDomainsName(int bridgeSimplified)
{
    typedef struct domainMapping {
        int id;
        char *name;
    } domain_Mapping_t;

    if (bridgeSimplified == TRUE)
    {
        domain_Mapping_t domain_mapping[] = {
            {1, "PILAR R 1"},
            {2, "PILAR D 1"},
            {3, "PILAR L 1"},
            {5, "PILAR R 2"},
            {6, "PILAR D 2"},
            {7, "PILAR L 2"},
            {13, "ROADWAY L"},
            {15, "ROADWAY R"},
            {14, "ROADWAY U"},
            {9, "SUB ROADWAY L"},
            {19, "SUB ROADWAY R"},
            {10, "SUB ROADWAY U 1"},
            {28, "SUB ROADWAY U 2"},
            {33, "SUB ROADWAY U 3"},
            {47, "SUB ROADWAY U 4"},
            {35, "SUB ROADWAY U 5"},
            {39, "SUB ROADWAY U 6"},
            {43, "SUB ROADWAY U 7"},
            {18, "SUB ROADWAY U 8"},
            {8, "SUB ROADWAY D 1"},
            {4, "SUB ROADWAY D 2"},
            {0, "SUB ROADWAY D 3"}
        };

        const int NB_DOMAINS = sizeof(domain_mapping) / sizeof(domain_mapping[0]);
        for (int i = 0; i < NB_DOMAINS; ++i)
        {
            domain_Mapping_t domain = domain_mapping[i];
            geoSetDomainName(domain.id, domain.name);
        }
    }
    else
    {
        domain_Mapping_t domain_mapping[] = {
            {5, "PILAR R 1"},
            {6, "PILAR D 1"},
            {7, "PILAR L 1"},
            {1, "PILAR R 2"},
            {2, "PILAR D 2"},
            {3, "PILAR L 2"},
            {10, "SUB ROADWAY U 1"},
            {33, "SUB ROADWAY U 2"},
            {38, "SUB ROADWAY U 3"},
            {45, "SUB ROADWAY U 4"},
            {40, "SUB ROADWAY U 5"},
            {52, "SUB ROADWAY U 6"},
            {56, "SUB ROADWAY U 7"},
            {27, "SUB ROADWAY U 8"},
            {8, "SUB ROADWAY D 1"},
            {4, "SUB ROADWAY D 2"},
            {0, "SUB ROADWAY D 3"},
            {9, "SUB ROADWAY L"},
            {28, "SUB ROADWAY R"},
            {13, "ROADWAY L"},
            {24, "ROADWAY R"},
            {29, "WINDOW D 1"},
            {30, "WINDOW L 1"},
            {31, "WINDOW R 1"},
            {32, "WINDOW U 1"},
            {46, "WINDOW D 2"},
            {47, "WINDOW L 2"},
            {48, "WINDOW R 2"},
            {49, "WINDOW U 2"},
            {11, "PILE L 1"},
            {35, "PILE R 1"},
            {37, "PILE L 2"},
            {44, "PILE R 2"},
            {43, "PILE L 3"},
            {41, "PILE R 3"},
            {51, "PILE L 4"},
            {55, "PILE R 4"},
            {54, "PILE L 5"},
            {26, "PILE R 5"},
            {12, "ARC 1"},
            {34, "ARC 2"},
            {36, "ARC 3"},
            {42, "ARC 4"},
            {39, "ARC 5"},
            {50, "ARC 6"},
            {53, "ARC 7"},
            {25, "ARC 8"},
            {15, "STAYCABLES L U 1"},
            {102, "STAYCABLES L U 2"},
            {98, "STAYCABLES L U 3"},
            {94, "STAYCABLES L U 4"},
            {90, "STAYCABLES L U 5"},
            {85, "STAYCABLES L U 6"},
            {82, "STAYCABLES L U 7"},
            {66, "STAYCABLES L U 8"},
            {63, "STAYCABLES L U 9"},
            {101, "STAYCABLES L D 1"},
            {97, "STAYCABLES L D 2"},
            {93, "STAYCABLES L D 3"},
            {87, "STAYCABLES L D 4"},
            {84, "STAYCABLES L D 5"},
            {79, "STAYCABLES L D 6"},
            {65, "STAYCABLES L D 7"},
            {60, "STAYCABLES L D 8"},
            {59, "STAYCABLES L D 9"},
            {22, "STAYCABLES R U 1"},
            {127, "STAYCABLES R U 2"},
            {123, "STAYCABLES R U 3"},
            {119, "STAYCABLES R U 4"},
            {113, "STAYCABLES R U 5"},
            {110, "STAYCABLES R U 6"},
            {106, "STAYCABLES R U 7"},
            {77, "STAYCABLES R U 8"},
            {73, "STAYCABLES R U 9"},
            {126, "STAYCABLES R D 1"},
            {122, "STAYCABLES R D 2"},
            {118, "STAYCABLES R D 3"},
            {116, "STAYCABLES R D 4"},
            {109, "STAYCABLES R D 5"},
            {105, "STAYCABLES R D 6"},
            {76, "STAYCABLES R D 7"},
            {72, "STAYCABLES R D 8"},
            {70, "STAYCABLES R D 9"},
            {58, "PYLON 1 L 1"},
            {62, "PYLON 1 L 2"},
            {64, "PYLON 1 L 3"},
            {81, "PYLON 1 L 4"},
            {83, "PYLON 1 L 5"},
            {68, "PYLON 1 R 1"},
            {74, "PYLON 1 R 2"},
            {78, "PYLON 1 R 3"},
            {107, "PYLON 1 R 4"},
            {111, "PYLON 1 R 5"},
            {89, "PYLON 2 L 1"},
            {92, "PYLON 2 L 2"},
            {96, "PYLON 2 L 3"},
            {100, "PYLON 2 L 4"},
            {114, "PYLON 2 R 1"},
            {120, "PYLON 2 R 2"},
            {124, "PYLON 2 R 3"},
            {128, "PYLON 2 R 4"},
            {17, "PYLON 3 L 1"},
            {20, "PYLON 3 R 1"},
            {91, "PYLON 1 U L"},
            {112, "PYLON 1 U R"},
            {16, "PYLON 2 U L"},
            {21, "PYLON 2 U R"},
            {18, "TOPBALL 1"},
            {19, "TOPBALL 2"},
            {14, "ROADWAY U 1"},
            {103, "ROADWAY U 2"},
            {99, "ROADWAY U 3"},
            {95, "ROADWAY U 4"},
            {88, "ROADWAY U 5"},
            {86, "ROADWAY U 6"},
            {80, "ROADWAY U 7"},
            {67, "ROADWAY U 8"},
            {61, "ROADWAY U 9"},
            {57, "ROADWAY U 10"},
            {69, "ROADWAY U 11"},
            {71, "ROADWAY U 12"},
            {75, "ROADWAY U 13"},
            {104, "ROADWAY U 14"},
            {108, "ROADWAY U 15"},
            {115, "ROADWAY U 16"},
            {117, "ROADWAY U 17"},
            {121, "ROADWAY U 18"},
            {125, "ROADWAY U 19"},
            {23, "ROADWAY U 20"}
        };

        const int NB_DOMAINS = sizeof(domain_mapping) / sizeof(domain_mapping[0]);
        for (int i = 0; i < NB_DOMAINS; ++i)
        {
            domain_Mapping_t domain = domain_mapping[i];
            geoSetDomainName(domain.id, domain.name);
        }
    }
}

void createBoundaryConditions(femProblem *theProblem, int bridgeSimplified)
{
    typedef struct domainBoundaryMapping {
        char *name;
        femBoundaryType type;
        double value1;
        double value2;
    } domainBoundaryMapping_t;

    const double length_camion = 16.5;
    const double width_camion = 2.55;
    const double mass_camion = 32000.0;
    const double mass_pedestrian = 70.0;

    if (bridgeSimplified == TRUE)
    {
        // Define constants
        const int nb_camion = 2;
        const int nb_pedestrian = 1;
        double weightPedestrianDensityBridge = 9.81 * mass_pedestrian * nb_pedestrian / 0.2; 
        double weightCamionDensityBridge = nb_camion * 9.81 * (mass_camion * width_camion) / (length_camion * width_camion);

        domainBoundaryMapping_t mapping[] = {
            {"PILAR R 1", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR D 1", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR L 1", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR R 2", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR D 2", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR L 2", DIRICHLET_XY, 0.0, 0.0},
            {"SUB ROADWAY U 1", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 2", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 3", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 4", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 5", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 6", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 7", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 8", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            // {"SUB ROADWAY D 1", NEUMANN_Y, 0.0, NAN},
            // {"SUB ROADWAY D 2", NEUMANN_Y, 0.0, NAN},
            // {"SUB ROADWAY D 3", NEUMANN_Y, 0.0, NAN},
            {"SUB ROADWAY L", DIRICHLET_XY, 0.0, 0.0},
            {"SUB ROADWAY R", DIRICHLET_XY, 0.0, 0.0},
            {"ROADWAY L", DIRICHLET_XY, 0.0, 0.0},
            {"ROADWAY R", DIRICHLET_X, 0.0, NAN},
            {"ROADWAY U", NEUMANN_Y, -weightCamionDensityBridge, NAN}
        };

        const int NB_DOMAINS = sizeof(mapping) / sizeof(mapping[0]);

        for (int i = 0; i < NB_DOMAINS; ++i)
        {
            domainBoundaryMapping_t domain = mapping[i];
            femElasticityAddBoundaryCondition(theProblem, domain.name, domain.type, domain.value1, domain.value2);
        }
    }
    else
    {
        int nb_camion = 2 * 30;
        int nb_pedestrian = 1;

        double weightCamionDensityBridge = nb_camion * 9.81 * (mass_camion * width_camion) / (length_camion * width_camion);
        double weightPedestrianDensityBridge = 9.81 * mass_pedestrian * nb_pedestrian / 0.2;

        domainBoundaryMapping_t mapping[] = {
            {"PILAR R 1", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR D 1", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR L 1", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR R 2", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR D 2", DIRICHLET_XY, 0.0, 0.0},
            {"PILAR L 2", DIRICHLET_XY, 0.0, 0.0},
            {"SUB ROADWAY U 1", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 2", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 3", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 4", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 5", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 6", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 7", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY U 8", NEUMANN_Y, -weightPedestrianDensityBridge, NAN},
            {"SUB ROADWAY D 1", NEUMANN_Y, 0.0, NAN},
            {"SUB ROADWAY D 2", NEUMANN_Y, 0.0, NAN},
            {"SUB ROADWAY D 3", NEUMANN_Y, 0.0, NAN},
            {"SUB ROADWAY L", DIRICHLET_XY, 0.0, 0.0},
            {"SUB ROADWAY R", DIRICHLET_XY, 0.0, 0.0},
            {"ROADWAY L", DIRICHLET_XY, 0.0, 0.0},
            {"ROADWAY R", DIRICHLET_X, 0.0, NAN},
            {"WINDOW D 1", NEUMANN_X, 0.0, NAN},
            {"WINDOW L 1", NEUMANN_Y, 0.0, NAN},
            {"WINDOW R 1", NEUMANN_Y, 0.0, NAN},
            {"WINDOW U 1", NEUMANN_X, 0.0, NAN},
            {"WINDOW D 2", NEUMANN_X, 0.0, NAN},
            {"WINDOW L 2", NEUMANN_Y, 0.0, NAN},
            {"WINDOW R 2", NEUMANN_Y, 0.0, NAN},
            {"WINDOW U 2", NEUMANN_X, 0.0, NAN},
            {"PILE L 1", NEUMANN_X, 0.0, NAN},
            {"PILE R 1", NEUMANN_X, 0.0, NAN},
            {"PILE L 2", NEUMANN_X, 0.0, NAN},
            {"PILE R 2", NEUMANN_X, 0.0, NAN},
            {"PILE L 3", NEUMANN_X, 0.0, NAN},
            {"PILE R 3", NEUMANN_X, 0.0, NAN},
            {"PILE L 4", NEUMANN_X, 0.0, NAN},
            {"PILE R 4", NEUMANN_X, 0.0, NAN},
            {"PILE L 5", NEUMANN_X, 0.0, NAN},
            {"PILE R 5", NEUMANN_X, 0.0, NAN},
            {"ARC 1", NEUMANN_X, 0.0, NAN},
            {"ARC 2", NEUMANN_X, 0.0, NAN},
            {"ARC 3", NEUMANN_X, 0.0, NAN},
            {"ARC 4", NEUMANN_X, 0.0, NAN},
            {"ARC 5", NEUMANN_X, 0.0, NAN},
            {"ARC 6", NEUMANN_X, 0.0, NAN},
            {"ARC 7", NEUMANN_X, 0.0, NAN},
            {"ARC 8", NEUMANN_X, 0.0, NAN},
            {"STAYCABLES L U 1", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 2", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 3", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 4", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 5", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 6", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 7", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 8", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L U 9", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 1", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 2", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 3", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 4", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 5", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 6", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 7", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 8", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES L D 9", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 1", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 2", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 3", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 4", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 5", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 6", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 7", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 8", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R U 9", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 1", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 2", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 3", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 4", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 5", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 6", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 7", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 8", DIRICHLET_T, 0.0, NAN},
            {"STAYCABLES R D 9", DIRICHLET_T, 0.0, NAN},
            {"PYLON 1 R 1", NEUMANN_X, 0.0, NAN},
            {"PYLON 1 R 2", NEUMANN_X, 0.0, NAN},
            {"PYLON 1 R 3", NEUMANN_X, 0.0, NAN},
            {"PYLON 1 R 4", NEUMANN_X, 0.0, NAN},
            {"PYLON 1 R 5", NEUMANN_X, 0.0, NAN},
            {"PYLON 2 R 1", NEUMANN_X, 0.0, NAN},
            {"PYLON 2 R 2", NEUMANN_X, 0.0, NAN},
            {"PYLON 2 R 3", NEUMANN_X, 0.0, NAN},
            {"PYLON 2 R 4", NEUMANN_X, 0.0, NAN},
            {"PYLON 3 R 1", NEUMANN_X, 0.0, NAN},
            {"PYLON 1 U L", NEUMANN_X, 0.0, NAN},
            {"PYLON 1 U R", NEUMANN_X, 0.0, NAN},
            {"PYLON 2 U L", NEUMANN_X, 0.0, NAN},
            {"PYLON 2 U R", NEUMANN_X, 0.0, NAN},
            {"TOPBALL 1", DIRICHLET_XY, 0.0, 0.0},
            {"TOPBALL 2", DIRICHLET_XY, 0.0, 0.0},
            {"ROADWAY U 1", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 2", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 3", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 4", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 5", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 6", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 7", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 8", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 9", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 10", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 11", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 13", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 14", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 15", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 16", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 17", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 18", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 19", NEUMANN_Y, -weightCamionDensityBridge, NAN},
            {"ROADWAY U 20", NEUMANN_Y, -weightCamionDensityBridge, NAN}
        };

        const int NB_DOMAINS = sizeof(mapping) / sizeof(mapping[0]);

        for (int i = 0; i < NB_DOMAINS; ++i)
        {
            domainBoundaryMapping_t domain = mapping[i];
            femElasticityAddBoundaryCondition(theProblem, domain.name, domain.type, domain.value1, domain.value2);
        }
    }
}


/****************************/
/********* GEO SIZE *********/
/****************************/

double hermiteInterpolation(double x, double h_max, double h_min, double dist_interp)
{
    if (x >= dist_interp) { return h_max; }
    if (x <= 0)           { return h_min; }

    return (1 / (dist_interp * dist_interp)) * (h_max - h_min) * x * x * (- 2 * x / dist_interp + 3) + h_min;
}

double getEuclidianDistance(double x1, double y1, double x2, double y2) { return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)); }

double getVerticalDistance(double y1, double y2) { return fabs(y1 - y2); }

double getHorizontalDistance(double x1, double x2) { return fabs(x1 - x2); }

double geoSizePillars(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    double h_max = theGeometry->defaultSize;
    double h_min = h_max / 7;

    double d_Vertical     = getVerticalDistance(y, 0.0);
    double d_Horizontal_L_1 = getHorizontalDistance(x, -theGeometry->rxLongArc - 2 * theGeometry->rxArc - 2 * theGeometry->widthPillars);
    double d_Horizontal_R_1 = getHorizontalDistance(x, -theGeometry->rxLongArc - 2 * theGeometry->rxArc - theGeometry->widthPillars);
    double d_Horizontal_L_2 = getHorizontalDistance(x, -theGeometry->rxLongArc - theGeometry->widthPillars);
    double d_Horizontal_R_2 = getHorizontalDistance(x, -theGeometry->rxLongArc);
    double d_interp      = theGeometry->heightPillars;
    double d_interp_side = theGeometry->widthPillars / 4;

    double h1 = hermiteInterpolation(d_Vertical, h_max, h_min, d_interp);
    double h2 = hermiteInterpolation(d_Horizontal_L_1, h_max, h_min, d_interp_side);
    double h3 = hermiteInterpolation(d_Horizontal_R_1, h_max, h_min, d_interp_side);
    double h4 = hermiteInterpolation(d_Horizontal_L_2, h_max, h_min, d_interp_side);
    double h5 = hermiteInterpolation(d_Horizontal_R_2, h_max, h_min, d_interp_side);

    return fmin(h1, fmin(h2, fmin(h3, fmin(h4, h5))));
}

double geoSizeSubRoadWay(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    double h_max = theGeometry->defaultSize;
    double h_min = h_max / 10;

    if ((y <= theGeometry->heightPillars + 0.1) && ((x >= 0.1 -theGeometry->rxLongArc - theGeometry->widthPillars && x <= -0.1 -theGeometry->rxLongArc) || (x >= 0.1 -theGeometry->rxLongArc - 2 * theGeometry->rxArc - 2 * theGeometry->widthPillars && x <= -0.1 -theGeometry->rxLongArc - 2 * theGeometry->rxArc - theGeometry->widthPillars))) { return 0.7 * h_max; }

    double d_Vertical_U  = getVerticalDistance(y, theGeometry->heightPillars + theGeometry->heightSubRoadWay);
    double d_Horizontal  = getHorizontalDistance(x, -theGeometry->widthSpanBridge / 2);
    double d_interp      = 3 * theGeometry->heightSubRoadWay / 4;
    double d_interp_side = 1.0;

    double d_Vertical_B = getVerticalDistance(y, theGeometry->heightPillars);

    double h1 = hermiteInterpolation(d_Vertical_U, h_max, h_min, 2 * d_interp);
    double h2 = hermiteInterpolation(d_Horizontal, h_max, h_min, d_interp_side);
    double h3 = hermiteInterpolation(d_Vertical_B, h_max, h_min, d_interp / 2);

    return fmin(h1, fmin(h2, h3));
}

double geoSizeTopBall(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    double h_max = theGeometry->defaultSize;
    double h_min = h_max / 7;

    double d_EuclTopBall = getEuclidianDistance(x, y, -theGeometry->rxLongArc - theGeometry->widthPillars / 2, theGeometry->heightPillars + theGeometry->heightBridge + 2 * theGeometry->heightPylons + 0.55);
    double rTopBall = 0.75;
    double d_interp = rTopBall - 0.1;

    return h_max + h_min - hermiteInterpolation(d_EuclTopBall, h_max, h_min, d_interp);
}

double geoSizePylons(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    double h_max = theGeometry->defaultSize;
    double h_min = h_max / 7;

    const double LIMIT_HEIGHT_PYLONS_STAGE1 = theGeometry->heightPillars + theGeometry->heightBridge + theGeometry->heightPylons;
    const double LIMIT_HEIGHT_PYLONS_STAGE2 = LIMIT_HEIGHT_PYLONS_STAGE1 + 2 * theGeometry->heightPylons / 3;
    const double LIMIT_HEIGHT_PYLONS_STAGE3 = LIMIT_HEIGHT_PYLONS_STAGE2 + theGeometry->heightPylons / 3;

    if (y <= LIMIT_HEIGHT_PYLONS_STAGE1)
    {
        double d_Vertical     = getVerticalDistance(y, theGeometry->heightPillars + theGeometry->heightBridge);
        double d_Hozizontal_L = getHorizontalDistance(x, -theGeometry->rxLongArc - theGeometry->widthPillars / 2 - theGeometry->widthPylons / 2);
        double d_Horizontal_R = getHorizontalDistance(x, -theGeometry->rxLongArc - theGeometry->widthPillars / 2 + theGeometry->widthPylons / 2);
        double d_interp_side = theGeometry->widthPylons / 4;
        double d_interp      = theGeometry->heightPylons / 2;

        double h1  = hermiteInterpolation(d_Vertical, h_max, h_min, d_interp);
        double h2  = hermiteInterpolation(d_Hozizontal_L, h_max, h_min, d_interp_side);
        double h3  = hermiteInterpolation(d_Horizontal_R, h_max, h_min, d_interp_side);
        return fmin(h1, fmin(h2, h3));
    }
    else if (y <= LIMIT_HEIGHT_PYLONS_STAGE2)
    {
        double d_Hozizontal_L = getHorizontalDistance(x, -theGeometry->rxLongArc - theGeometry->widthPillars / 2 - theGeometry->widthPylons / 4);
        double d_Horizontal_R = getHorizontalDistance(x, -theGeometry->rxLongArc - theGeometry->widthPillars / 2 + theGeometry->widthPylons / 4);
        double d_interp_side = theGeometry->widthPylons / 4;

        double h1  = hermiteInterpolation(d_Hozizontal_L, h_max, h_min, d_interp_side);
        double h2  = hermiteInterpolation(d_Horizontal_R, h_max, h_min, d_interp_side);
        return fmin(h1, h2);
    }
    else if (y <= LIMIT_HEIGHT_PYLONS_STAGE3)
    {
        double d_Hozizontal_L = getHorizontalDistance(x, -theGeometry->rxLongArc - theGeometry->widthPillars / 2 - theGeometry->widthPylons / 8);
        double d_Horizontal_R = getHorizontalDistance(x, -theGeometry->rxLongArc - theGeometry->widthPillars / 2 + theGeometry->widthPylons / 8);
        double d_interp_side = theGeometry->widthPylons / 6;

        double h1  = hermiteInterpolation(d_Hozizontal_L, h_max, h_min, d_interp_side);
        double h2  = hermiteInterpolation(d_Horizontal_R, h_max, h_min, d_interp_side);
        return fmin(h1, h2);
    }
    else { return geoSizeTopBall(x, y); }
}

double getX_Ellipse(double a_ellipse, double b_ellipse, double xc, double yc, double x, double y)
{
    if (fabs(x - xc) < 1e-6 || fabs(y - yc) < 1e-6) { return 0.0; }

    double m = (y - yc) / (x - xc);

    double a = 1 / (a_ellipse * a_ellipse) + (m * m) / (b_ellipse * b_ellipse);
    double b = - (2 * xc * m * m) / (b_ellipse * b_ellipse) - (2 * xc) / (a_ellipse * a_ellipse);
    double c = - 1 + (m * m * xc * xc) / (b_ellipse * b_ellipse) + (xc * xc) / (a_ellipse * a_ellipse);

    double discriminant = b * b - 4 * a * c;
    if (discriminant <= 0) { Error("The discriminant must be positive !\nThere is always two solutions to this problem.\n"); }
    return (-b + sqrt(discriminant)) / (2 * a);
}

double geoSizeBridge(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    double h_max = theGeometry->defaultSize;
    double h_min = h_max / 10;

    double d_Horizontal_RoadWay  = getHorizontalDistance(x, -theGeometry->widthSpanBridge / 2);
    double d_Vertical_RoadWay    = getVerticalDistance(y, theGeometry->heightPillars + theGeometry->heightBridge);
    double d_Vertical_SubRoadWay = getVerticalDistance(y, theGeometry->heightPillars + theGeometry->heightSubRoadWay);

    double d_interp_RoadWay      = theGeometry->heightBridge - theGeometry->ryArc - 0.2;
    double d_interp_SubRoadWay   = 2 * theGeometry->heightSubRoadWay;
    double d_interp_RoadWay_Side = 0.8;

    double h1 = hermiteInterpolation(d_Vertical_RoadWay, h_max, h_min, d_interp_RoadWay);
    double h2 = hermiteInterpolation(d_Horizontal_RoadWay, h_max, h_min, d_interp_RoadWay_Side);
    double h3 = hermiteInterpolation(d_Vertical_SubRoadWay, h_max, h_min, d_interp_SubRoadWay);

    h_min = h_max / 7;

    double r_min = fmin(fmin(fmin(theGeometry->rxArc, theGeometry->ryArc), theGeometry->ryLongArc), theGeometry->rxLongArc);
    double d_interp_Arcs = 1.5;
    double yc = theGeometry->heightPillars + theGeometry->heightSubRoadWay;

    double xc = 0.0;
    double a = theGeometry->rxLongArc;
    double b = theGeometry->ryLongArc;
    double xEllipse = getX_Ellipse(a, b, xc, yc, x, y);
    double yEllipse = yc + sqrt(b * b * (1 - ((xEllipse - xc) * (xEllipse - xc)) / (a * a)));
    double rELlipse = getEuclidianDistance(xEllipse, yEllipse, xc, yc);
    double d_EuclArc = getEuclidianDistance(x, y, xc, yc);

    double h4 = (d_EuclArc < rELlipse) ? h_min : hermiteInterpolation(d_EuclArc - rELlipse, h_max, h_min, d_interp_Arcs);

    a = theGeometry->rxArc;
    b = theGeometry->ryArc;

    xc = - theGeometry->rxLongArc - theGeometry->widthPillars - theGeometry->rxArc;
    xEllipse  = getX_Ellipse(a, b, xc, yc, x, y);
    yEllipse  = yc + sqrt(b * b * (1 - ((xEllipse - xc) * (xEllipse - xc)) / (a * a)));
    rELlipse  = getEuclidianDistance(xEllipse, yEllipse, xc, yc);
    d_EuclArc = getEuclidianDistance(x, y, xc, yc);

    double h5 = (d_EuclArc < rELlipse) ? h_min : hermiteInterpolation(d_EuclArc - rELlipse, h_max, h_min, d_interp_Arcs);

    xc = - theGeometry->widthSpanBridge / 2;
    xEllipse  = getX_Ellipse(a, b, xc, yc, x, y);
    yEllipse  = yc + sqrt(b * b * (1 - ((xEllipse - xc) * (xEllipse - xc)) / (a * a)));
    rELlipse  = getEuclidianDistance(xEllipse, yEllipse, xc, yc);
    d_EuclArc = getEuclidianDistance(x, y, xc, yc);

    double h6 = (d_EuclArc < rELlipse) ? h_min : hermiteInterpolation(d_EuclArc - rELlipse, h_max, h_min, d_interp_Arcs);

    return fmin(h1, fmin(h2, fmin(h3, fmin(h4, fmin(h5, h6)))));
}

double geoSize(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    double h_max = theGeometry->defaultSize;
    double h_stayCables_min = h_max / 15;

    double factor_size = 3.0;
    if (x > 0) { return factor_size * h_max; }

    if (y < theGeometry->heightPillars)                                       { return factor_size * geoSizePillars(x, y); }
    else if (y <= theGeometry->heightPillars + theGeometry->heightSubRoadWay) { return factor_size * geoSizeSubRoadWay(x, y); }
    else if (y <= theGeometry->heightPillars + theGeometry->heightBridge)     { return factor_size * geoSizeBridge(x, y); }
    else                                                                      { return (isStayCables(x, y) == TRUE) ? factor_size * h_stayCables_min : factor_size * geoSizePylons(x, y); }
}

/***************************************************/
/********* MESH + GEO SIZE (ON AN EXAMPLE) *********/
/***************************************************/

double geoSizeExample(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    return theGeometry->defaultSize;
}

void geoMeshGenerateExample(femDiscreteType discreteType, int beam_example)
{
    if (beam_example == TRUE) { geoMeshGenerate_Beam(discreteType); }
    else                      { geoMeshGenerate_UForm(discreteType); }
}

void geoMeshGenerate_UForm(femDiscreteType discreteType)
{
    /*
    4 ------------------ 3
    |                    |
    |                    |
    5 ------- 6          |
                \        |
                )        |
                /        |
    8 ------- 7          |
    |                    |
    |                    |
    1 ------------------ 2
    */

    femGeometry *theGeometry = geoGetGeometry();
    double Lx = 1.0;
    double Ly = 1.0;
    theGeometry->LxPlate = Lx;
    theGeometry->LyPlate = Ly;
    theGeometry->defaultSize = Lx * 0.02;
    theGeometry->geoSize = geoSizeExample;

    geoSetSizeCallback(theGeometry->geoSize);

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;

    int ierr;
    double r = w / 4;
    int idRect = createRectangle(0.0, 0.0, w, h);
    int idDisk = createDisk(w / 2.0, h / 2.0, r, r);
    int idSlit = createRectangle(w / 2.0, h / 2.0 - r, w, 2.0 * r);
    int rect[] = {2, idRect};
    int disk[] = {2, idDisk};
    int slit[] = {2, idSlit};

    cutElement(rect, disk);
    cutElement(rect, slit);

    gmshModelOccSynchronize(&ierr);

    if (theGeometry->elementType == FEM_QUAD)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    if (theGeometry->elementType == FEM_TRIANGLE)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    
    return;
}

void geoMeshGenerate_Beam(femDiscreteType discreteType)
{
    femGeometry *theGeometry = geoGetGeometry();
    double Lx = 8.0;
    double Ly = 1.0;
    theGeometry->LxPlate = Lx;
    theGeometry->LyPlate = Ly;
    theGeometry->defaultSize = Lx * 0.01;
    theGeometry->geoSize = geoSizeExample;

    geoSetSizeCallback(theGeometry->geoSize);

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;

    int ierr;
    double r = w / 4;
    int idRect = createRectangle(0.0, 0.0, w, h);

    gmshModelOccSynchronize(&ierr);

    if (theGeometry->elementType == FEM_QUAD)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    if (theGeometry->elementType == FEM_TRIANGLE)
    {
        // Generate quadratic elements
        if (discreteType == FEM_DISCRETE_TYPE_QUADRATIC) { gmshOptionSetNumber("Mesh.ElementOrder", 2, &ierr); }
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    
    return;
}