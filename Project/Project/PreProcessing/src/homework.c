#include "../../../fem_library/include/fem_geometry.h"
#include "../../../fem_library/include/fem_gmsh.h"


// TODO : Raffiner intelligemment


/*********************************************************************************************************************/
/************** POSITION GEOMETRY ***** POSITION GEOMETRY ***** POSITION GEOMETRY ***** POSITION GEOMETRY ************/
/*********************************************************************************************************************/

int isTablier(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();
    const double HEIGHT_TABLIER = 1.0;

    if ((y >= theGeometry->heightPillier + theGeometry->heightPlate - HEIGHT_TABLIER) && (y <= theGeometry->heightPillier + theGeometry->heightPlate)) { return TRUE; }
    return FALSE;
}

int isSubTablier(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();

    if ((y < theGeometry->heightPillier) || (y > theGeometry->heightSubPlate + theGeometry->heightPillier)) { return FALSE; }
    
    if ((x > theGeometry->rxLongArc) && (x < theGeometry->rxLongArc + theGeometry->widthPillier)) { return FALSE; }
    if ((x > theGeometry->rxLongArc + theGeometry->widthPillier + 2 * theGeometry->rxArc) &&
        (x < theGeometry->rxLongArc + 2 * theGeometry->widthPillier + 2 * theGeometry->rxArc))
        { return FALSE; }

    if ((x < -theGeometry->rxLongArc) && (x > -theGeometry->rxLongArc - theGeometry->widthPillier)) { return FALSE; }
    if ((x < -theGeometry->rxLongArc - theGeometry->widthPillier - 2 * theGeometry->rxArc) &&
        (x > -theGeometry->rxLongArc - 2 * theGeometry->widthPillier - 2 * theGeometry->rxArc))
        { return FALSE; }

    return TRUE;
}

int isCables(double x, double y)
{
    femGeometry *theGeometry = geoGetGeometry();

    if ((y <= theGeometry->heightPlate + theGeometry->heightPillier) || (y >= theGeometry->heightPlate + theGeometry->heightPillier + 5 * theGeometry->heightBigColumn / 3)) { return FALSE; }

    if ((x < theGeometry->rxLongArc + theGeometry->widthPillier / 2 - theGeometry->widthBigColumn / 2) &&
        (x > - theGeometry->rxLongArc - theGeometry->widthPillier / 2 + theGeometry->widthBigColumn / 2))
        { return TRUE; }

    if (x > theGeometry->rxLongArc + theGeometry->widthPillier / 2 + theGeometry->widthBigColumn / 2)  { return TRUE; }
    if (x < -theGeometry->rxLongArc - theGeometry->widthPillier / 2 - theGeometry->widthBigColumn / 2) { return TRUE; }

    if (y > theGeometry->heightPlate + theGeometry->heightPillier + theGeometry->heightBigColumn)
    {
        if ((x < theGeometry->rxLongArc + theGeometry->widthPillier / 2 - theGeometry->widthBigColumn / 4) &&
            (x > -theGeometry->rxLongArc - theGeometry->widthPillier / 2 + theGeometry->widthBigColumn / 4))
            { return TRUE; }

        if (x > theGeometry->rxLongArc + theGeometry->widthPillier / 2 + theGeometry->widthBigColumn / 4)  { return TRUE; }
        if (x < -theGeometry->rxLongArc - theGeometry->widthPillier / 2 - theGeometry->widthBigColumn / 4) { return TRUE; }
    }
    return FALSE;
}


/*************************************************************************************************************/
/********* MATERIALS ********* MATERIALS ********* MATERIALS ********* MATERIALS ********* MATERIALS *********/
/**************************************************************************************************************/

double *getMaterialProperties(char *material)
{
    double *properties = malloc(3 * sizeof(double));
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

char *getMaterials(double x, double y)
{
    if (isSubTablier(x, y) == TRUE || isTablier(x, y) == TRUE || isCables(x, y) == TRUE) { return "steel"; }
    return "reinforced_concrete";
}


/************************************************************************************/
/********* MESH ********* MESH ********* MESH ********* MESH ********* MESH *********/
/************************************************************************************/

double geoSize(double x, double y)
{  
    if (isSubTablier(x, y) == TRUE) { return 0.6; }
    if (isTablier(x, y) == TRUE)    { return 0.6; }
    if (isCables(x, y) == TRUE)     { return 0.1; }
    else                            { return 0.6; }
}

void createWindows(femGeometry *theGeometry, int *idWindows, int ierr)
{
    double offset = 0.5;
    double y = theGeometry->heightPillier + theGeometry->heightSubPlate + offset;
    double x1 = theGeometry->rxLongArc + theGeometry->widthPillier / 4;
    double x2 = theGeometry->rxLongArc + 2 * theGeometry->rxArc + 5 * theGeometry->widthPillier / 4;

    int id1 = gmshModelOccAddRectangle(- x1 - theGeometry->widthWindow, y, 0.0, theGeometry->widthWindow, theGeometry->heightWindow, -1, 0, &ierr);
    int id2 = gmshModelOccAddRectangle(- x2 - theGeometry->widthWindow, y, 0.0, theGeometry->widthWindow, theGeometry->heightWindow, -1, 0, &ierr);
    int id3 = gmshModelOccAddRectangle(x1, y, 0.0, theGeometry->widthWindow, theGeometry->heightWindow, -1, 0, &ierr);
    int id4 = gmshModelOccAddRectangle(x2, y, 0.0, theGeometry->widthWindow, theGeometry->heightWindow, -1, 0, &ierr);

    int tempId[] = {id1, id2, id3, id4};
    memcpy(idWindows, tempId, 4 * sizeof(int));
}

void createColumns(femGeometry *theGeometry, int *idColumns, int ierr)
{
    double offset = theGeometry->widthColumn / 2;
    double y = theGeometry->heightPillier + theGeometry->heightSubPlate;

    double x1 = theGeometry->rxLongArc / 4;
    double x2 = 3 * theGeometry->rxLongArc / 4;
    double x3 = theGeometry->rxLongArc + theGeometry->widthPillier + theGeometry->rxArc / 2;
    double x4 = theGeometry->rxLongArc + theGeometry->widthPillier + 3 * theGeometry->rxArc / 2;
    double x5 = theGeometry->widthPlate / 2 - theGeometry->rxArc / 2;

    int id1  = gmshModelOccAddRectangle(x1 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id2  = gmshModelOccAddRectangle(- x1 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id3  = gmshModelOccAddRectangle(x2 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id4  = gmshModelOccAddRectangle(- x2 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id5  = gmshModelOccAddRectangle(x3 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id6  = gmshModelOccAddRectangle(- x3 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id7  = gmshModelOccAddRectangle(x4 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id8  = gmshModelOccAddRectangle(- x4 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id9  = gmshModelOccAddRectangle(x5 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);
    int id10 = gmshModelOccAddRectangle(- x5 - offset, y, 0.0, theGeometry->widthColumn, theGeometry->ryArc, -1, 0, &ierr);

    int tempId[] = {id1, id2, id3, id4, id5, id6, id7, id8, id9, id10};
    memcpy(idColumns, tempId, 10 * sizeof(int));
}

void createPilliers(femGeometry *theGeometry, int *idPilliers, int ierr)
{
    double width  = theGeometry->widthPillier;
    double height = theGeometry->heightPillier;

    double y = 0.0;
    double x1 = theGeometry->rxLongArc;
    double x2 = theGeometry->widthPlate / 2 - theGeometry->rxArc;

    int id1 = gmshModelOccAddRectangle(- x2, y, 0.0, width, height, -1, 0, &ierr);
    int id2 = gmshModelOccAddRectangle(- x1 - width, y, 0.0, width, height, -1, 0, &ierr);
    int id3 = gmshModelOccAddRectangle(x1, y, 0.0, width, height, -1, 0, &ierr);
    int id4 = gmshModelOccAddRectangle(x2 - width, y, 0.0, width, height, -1, 0, &ierr);

    int tempId[] = {id1, id2, id3, id4};
    memcpy(idPilliers, tempId, 4 * sizeof(int));
}

void createArcs(femGeometry *theGeometry, int *idArcs, int ierr)
{
    double rxArc = theGeometry->rxArc;
    double ryArc = theGeometry->ryArc;
    double rxLongArc = theGeometry->rxLongArc;
    double ryLongArc = theGeometry->ryLongArc;

    double y = theGeometry->heightPillier + theGeometry->heightSubPlate;
    double x1 = theGeometry->rxLongArc + 2 * theGeometry->widthPillier + 3 * theGeometry->rxArc;
    double x2 = theGeometry->rxLongArc + theGeometry->widthPillier;

    int id1 = gmshModelOccAddDisk(- x1, y, 0.0, rxArc, ryArc, -1, NULL, 0, NULL, 0, &ierr);
    int id2 = gmshModelOccAddDisk(- x2 - rxArc, y, 0.0, rxArc, ryArc, -1, NULL, 0, NULL, 0, &ierr);
    int id3 = gmshModelOccAddDisk(0.0, y, 0.0, rxLongArc, ryLongArc, -1, NULL, 0, NULL, 0, &ierr);
    int id4 = gmshModelOccAddDisk(x2 + rxArc, y, 0.0, rxArc, ryArc, -1, NULL, 0, NULL, 0, &ierr);
    int id5 = gmshModelOccAddDisk(x1, y, 0.0, rxArc, ryArc, -1, NULL, 0, NULL, 0, &ierr);

    int tempId[] = {id1, id2, id3, id4, id5};
    memcpy(idArcs, tempId, 5 * sizeof(int));
}

void createBigColumns(femGeometry *theGeometry, int *idBigColumns, int ierr)
{
    double width = theGeometry->widthBigColumn;
    double height = theGeometry->heightBigColumn;

    double y = theGeometry->heightPillier + theGeometry->heightPlate;

    // To change if the column position change
    double x = theGeometry->rxLongArc + theGeometry->widthPillier / 2;

    int id1 = gmshModelOccAddRectangle(-x - width / 2, y, 0.0, width, height, -1, 0, &ierr);
    int id2 = gmshModelOccAddRectangle(-x - width / 4, y + height, 0.0, width / 2, 2 * height / 3, -1, 0, &ierr);
    int id3 = gmshModelOccAddRectangle(-x - width / 8, y + 5 * height / 3, 0.0, width / 4, height / 3, -1, 0, &ierr);
    int id4 = gmshModelOccAddRectangle(x - width / 2, y, 0.0, width, height, -1, 0, &ierr);
    int id5 = gmshModelOccAddRectangle(x - width / 4, y + height, 0.0, width / 2, 2 * height / 3, -1, 0, &ierr);
    int id6 = gmshModelOccAddRectangle(x - width / 8, y + 5 * height / 3, 0.0, width / 4, height / 3, -1, 0, &ierr);

    int tempId[] = {id1, id2, id3, id4, id5, id6};
    memcpy(idBigColumns, tempId, 6 * sizeof(int));
}

void createDisk(femGeometry *theGeometry, int *idDisks, int ierr)
{
    double r = theGeometry->widthBigColumn / 3;
    // To change if the column position change
    double x = theGeometry->rxLongArc + theGeometry->widthPillier / 2;
    double y = theGeometry->heightPillier + theGeometry->heightPlate + 2 * theGeometry->heightBigColumn + r / 3;

    int id1 = gmshModelOccAddDisk(- x, y, 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr);
    int id2 = gmshModelOccAddDisk(x, y, 0.0, r, r, -1, NULL, 0, NULL, 0, &ierr);

    int tempId[] = {id1, id2};
    memcpy(idDisks, tempId, 2 * sizeof(int));
}

void createCables(femGeometry *theGeometry, int *idCables, double *positionX, double *positionY, int ierr)
{
    double width = theGeometry->widthCable;
    double height = theGeometry->heightCable;
    double distance = theGeometry->distanceBetweenCable;

    // To change if the column position change
    double centerColumn = theGeometry->rxLongArc + theGeometry->widthPillier / 2;

    const double y  = theGeometry->heightPillier + theGeometry->heightPlate + 5 * theGeometry->heightBigColumn / 3;
    double x1 = centerColumn + theGeometry->widthBigColumn / 4;
    double x2 = centerColumn + theGeometry->widthBigColumn / 2;

    // Fils de droite sur la colonne de droite
    int idR1 = gmshModelOccAddRectangle(x1, y, 0.0, width, height, -1, 0, &ierr);
    int idR2 = gmshModelOccAddRectangle(x1, y - 1 * distance, 0.0, width, height - 1, -1, 0, &ierr);
    int idR3 = gmshModelOccAddRectangle(x1, y - 2 * distance, 0.0, width, height - 2, -1, 0, &ierr);
    int idR4 = gmshModelOccAddRectangle(x1, y - 3 * distance, 0.0, width, height - 3, -1, 0, &ierr);
    int idR5 = gmshModelOccAddRectangle(x2, y - 5 * distance, 0.0, width, height - 5, -1, 0, &ierr);
    int idR6 = gmshModelOccAddRectangle(x2, y - 6 * distance, 0.0, width, height - 6, -1, 0, &ierr);
    int idR7 = gmshModelOccAddRectangle(x2, y - 7 * distance, 0.0, width, height - 7, -1, 0, &ierr);
    int idR8 = gmshModelOccAddRectangle(x2, y - 8 * distance, 0.0, width, height - 8, -1, 0, &ierr);
    int idR9 = gmshModelOccAddRectangle(x2, y - 9 * distance, 0.0, width, height - 9, -1, 0, &ierr);

    for (int i = 0; i < 4; i++) { positionX[i] = x1; positionY[i] = y - i * distance; }
    for (int i = 4; i < 9; i++) { positionX[i] = x2; positionY[i] = y - (i + 1) * distance; }

    x1 += width;
    x2 += width;

    // Fils de gauche sur la colonne de gauche
    int idL1 = gmshModelOccAddRectangle(-x1, y, 0.0, width, height, -1, 0, &ierr);
    int idL2 = gmshModelOccAddRectangle(-x1, y - 1 * distance, 0.0, width, height - 1, -1, 0, &ierr);
    int idL3 = gmshModelOccAddRectangle(-x1, y - 2 * distance, 0.0, width, height - 2, -1, 0, &ierr);
    int idL4 = gmshModelOccAddRectangle(-x1, y - 3 * distance, 0.0, width, height - 3, -1, 0, &ierr);
    int idL5 = gmshModelOccAddRectangle(-x2, y - 5 * distance, 0.0, width, height - 5, -1, 0, &ierr);
    int idL6 = gmshModelOccAddRectangle(-x2, y - 6 * distance, 0.0, width, height - 6, -1, 0, &ierr);
    int idL7 = gmshModelOccAddRectangle(-x2, y - 7 * distance, 0.0, width, height - 7, -1, 0, &ierr);
    int idL8 = gmshModelOccAddRectangle(-x2, y - 8 * distance, 0.0, width, height - 8, -1, 0, &ierr);
    int idL9 = gmshModelOccAddRectangle(-x2, y - 9 * distance, 0.0, width, height - 9, -1, 0, &ierr);

    for (int i = 18; i < 22; i++) { positionX[i] = - x1 + width; positionY[i] = y - (i - 18) * distance; }
    for (int i = 22; i < 27; i++) { positionX[i] = - x2 + width; positionY[i] = y - (i + 1 - 18) * distance; }

    x1 = centerColumn - theGeometry->widthBigColumn / 4 - width;
    x2 = centerColumn - theGeometry->widthBigColumn / 2 - width;

    // Fils de gauche sur la colonne de droite
    int idL10 = gmshModelOccAddRectangle(x1, y, 0.0, width, height, -1, 0, &ierr);
    int idL11 = gmshModelOccAddRectangle(x1, y - 1 * distance, 0.0, width, height - 1, -1, 0, &ierr);
    int idL12 = gmshModelOccAddRectangle(x1, y - 2 * distance, 0.0, width, height - 2, -1, 0, &ierr);
    int idL13 = gmshModelOccAddRectangle(x1, y - 3 * distance, 0.0, width, height - 3, -1, 0, &ierr);
    int idL14 = gmshModelOccAddRectangle(x2, y - 5 * distance, 0.0, width, height - 5, -1, 0, &ierr);
    int idL15 = gmshModelOccAddRectangle(x2, y - 6 * distance, 0.0, width, height - 6, -1, 0, &ierr);
    int idL16 = gmshModelOccAddRectangle(x2, y - 7 * distance, 0.0, width, height - 7, -1, 0, &ierr);
    int idL17 = gmshModelOccAddRectangle(x2, y - 8 * distance, 0.0, width, height - 8, -1, 0, &ierr);
    int idL18 = gmshModelOccAddRectangle(x2, y - 9 * distance, 0.0, width, height - 9, -1, 0, &ierr);

    for (int i = 27; i < 31; i++) { positionX[i] = x1 + width; positionY[i] = y - (i - 27) * distance; }
    for (int i = 31; i < 36; i++) { positionX[i] = x2 + width; positionY[i] = y - (i + 1 - 27) * distance; }
    
    x1 += width;
    x2 += width;

    // Fils de droite sur la colonne de gauche
    int idR10 = gmshModelOccAddRectangle(-x1, y, 0.0, width, height, -1, 0, &ierr);
    int idR11 = gmshModelOccAddRectangle(-x1, y - 1 * distance, 0.0, width, height - 1, -1, 0, &ierr);
    int idR12 = gmshModelOccAddRectangle(-x1, y - 2 * distance, 0.0, width, height - 2, -1, 0, &ierr);
    int idR13 = gmshModelOccAddRectangle(-x1, y - 3 * distance, 0.0, width, height - 3, -1, 0, &ierr);
    int idR14 = gmshModelOccAddRectangle(-x2, y - 5 * distance, 0.0, width, height - 5, -1, 0, &ierr);
    int idR15 = gmshModelOccAddRectangle(-x2, y - 6 * distance, 0.0, width, height - 6, -1, 0, &ierr);
    int idR16 = gmshModelOccAddRectangle(-x2, y - 7 * distance, 0.0, width, height - 7, -1, 0, &ierr);
    int idR17 = gmshModelOccAddRectangle(-x2, y - 8 * distance, 0.0, width, height - 8, -1, 0, &ierr);
    int idR18 = gmshModelOccAddRectangle(-x2, y - 9 * distance, 0.0, width, height - 9, -1, 0, &ierr);

    for (int i = 9; i < 13; i++) { positionX[i] = - x1; positionY[i] = y - (i - 9) * distance; }
    for (int i = 13; i < 18; i++) { positionX[i] = - x2; positionY[i] = y - (i + 1 - 9) * distance; }

    int tempId[] = {idR1, idR2, idR3, idR4, idR5, idR6, idR7, idR8, idR9,
                    idR10, idR11, idR12, idR13, idR14, idR15, idR16, idR17, idR18,
                    idL1, idL2, idL3, idL4, idL5, idL6, idL7, idL8, idL9,
                    idL10, idL11, idL12, idL13, idL14, idL15, idL16, idL17, idL18};

    memcpy(idCables, tempId, 36 * sizeof(int));
}

void createPlate(femGeometry *theGeometry, int *idPlate, int ierr)
{
    *idPlate = gmshModelOccAddRectangle(- theGeometry->widthPlate / 2, theGeometry->heightPillier, 0.0, theGeometry->widthPlate, theGeometry->heightPlate, -1, 0, &ierr);
}

void createSubPlate(femGeometry *theGeometry, int *idSubPlate, int ierr)
{
    *idSubPlate = gmshModelOccAddRectangle(- theGeometry->widthPlate / 2, theGeometry->heightPillier, 0.0, theGeometry->widthSubPlate, theGeometry->heightSubPlate, -1, 0, &ierr);
}

void geoMeshGenerate()
{
    femGeometry *theGeometry = geoGetGeometry();

    // Define the geometry parameters
    theGeometry->widthPlate = 62.0;
    theGeometry->heightPlate = 6.5; // Minimum 6.0
    theGeometry->widthWindow = 1.5;
    theGeometry->heightWindow = 0.6;
    theGeometry->widthSubPlate = 62.0;
    theGeometry->heightSubPlate = 1.0;
    theGeometry->rxArc = 5.0;
    theGeometry->ryArc = 4.0;
    theGeometry->rxLongArc = 10.0;
    theGeometry->ryLongArc = 3.0;
    theGeometry->widthColumn = 0.5;
    theGeometry->widthPillier = 3.0;
    theGeometry->heightPillier = 5.0;
    theGeometry->widthBigColumn = 2;
    theGeometry->heightBigColumn = 6.0;
    theGeometry->angleCable = 135 * M_PI / 180;
    theGeometry->widthCable = 0.2;
    theGeometry->heightCable = 14.5;
    theGeometry->distanceBetweenCable = 0.8;
    theGeometry->h = 0.6;

    theGeometry->geoSize = geoSize;
    theGeometry->getMaterialProperties = getMaterialProperties;
    theGeometry->getMaterials = getMaterials;

    theGeometry->elementType = FEM_TRIANGLE;
    
    int ierr;
    int idPlate, idSubPlate;

    int *idWindows          = malloc(4 * sizeof(int));
    int *idColumns          = malloc(10 * sizeof(int));
    int *idPilliers         = malloc(4 * sizeof(int));
    int *idArcs             = malloc(5 * sizeof(int));
    int *idBigColumns       = malloc(6 * sizeof(int));
    int *idDisks            = malloc(2 * sizeof(int));
    int *idCables           = malloc(36 * sizeof(int));
    double *positionCablesX = malloc(36 * sizeof(double));
    double *positionCablesY = malloc(36 * sizeof(double));

    if (idWindows == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idColumns == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idPilliers == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idArcs == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idBigColumns == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idDisks == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (idCables == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (positionCablesX == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    if (positionCablesY == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }

    int plate[2], subPlate[2], pillier[4][2], bigColumn[6][2], cable[36][2], arc[5][2], window[4][2], disk[2][2], column[10][2];

    createPlate(theGeometry, &idPlate, ierr);
    createSubPlate(theGeometry, &idSubPlate, ierr);
    createPilliers(theGeometry, idPilliers, ierr);
    createBigColumns(theGeometry, idBigColumns, ierr);

    subPlate[0] = 2;
    plate[0]    = 2;
    plate[1]    = idPlate;
    subPlate[1] = idSubPlate;

    for (int i = 0; i < 4; i++)  { pillier[i][0] = 2; pillier[i][1] = idPilliers[i]; }
    for (int i = 0; i < 6; i++)  { bigColumn[i][0] = 2; bigColumn[i][1] = idBigColumns[i]; }
    for (int i = 0; i < 36; i++) { cable[i][0] = 2; cable[i][1] = idCables[i]; }

    createArcs(theGeometry, idArcs, ierr);
    createWindows(theGeometry, idWindows, ierr);

    for (int i = 0; i < 5; i++)  { arc[i][0] = 2; arc[i][1] = idArcs[i]; }
    for (int i = 0; i < 4; i++)  { window[i][0] = 2; window[i][1] = idWindows[i]; }

    for (int i = 0; i < 5; i++) { gmshModelOccCut(plate, 2, arc[i], 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); }
    for (int i = 0; i < 4; i++) { gmshModelOccCut(plate, 2, window[i], 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); }
    
    createColumns(theGeometry, idColumns, ierr);
    createDisk(theGeometry, idDisks, ierr);
    createCables(theGeometry, idCables, positionCablesX, positionCablesY, ierr);

    for (int i = 0; i < 2; i++)  { disk[i][0] = 2; disk[i][1] = idDisks[i]; }
    for (int i = 0; i < 10; i++) { column[i][0] = 2; column[i][1] = idColumns[i]; }
    for (int i = 0; i < 36; i++) { cable[i][0] = 2; cable[i][1] = idCables[i]; }

    // Rotate the cables
    for (int i = 0; i < 18; i++)  { gmshModelOccRotate(cable[i], 2, positionCablesX[i], positionCablesY[i], 0.0, 0.0, 0.0, 1.0, -theGeometry->angleCable, &ierr); }
    for (int i = 18; i < 36; i++) { gmshModelOccRotate(cable[i], 2, positionCablesX[i], positionCablesY[i], 0.0, 0.0, 0.0, 1.0, theGeometry->angleCable, &ierr); }
    
    // Fuse the elements
    gmshModelOccFuse(plate, 2, subPlate, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    for (int i = 0; i < 4; i++)  { gmshModelOccFuse(plate, 2, pillier[i], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); }
    for (int i = 0; i < 10; i++) { gmshModelOccFuse(plate, 2, column[i], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); }
    for (int i = 1; i < 3; i++)  { gmshModelOccFuse(bigColumn[0], 2, bigColumn[i], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);}
    for (int i = 4; i < 6; i++)  { gmshModelOccFuse(bigColumn[3], 2, bigColumn[i], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);}
    gmshModelOccFuse(bigColumn[0], 2, disk[0], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(bigColumn[3], 2, disk[1], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    for (int i = 0; i < 9; i++)   { gmshModelOccFuse(bigColumn[3], 2, cable[i], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); }
    for (int i = 27; i < 36; i++) { gmshModelOccFuse(bigColumn[3], 2, cable[i], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); }
    for (int i = 9; i < 27; i++)  { gmshModelOccFuse(bigColumn[0], 2, cable[i], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); }
    gmshModelOccFuse(plate, 2, bigColumn[0], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    gmshModelOccFuse(plate, 2, bigColumn[3], 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    
    // Begin : Cut the half of the bridge for the symmetry
    int idRectFilterLeft = gmshModelOccAddRectangle(0.0, 0.0, 0.0, theGeometry->widthPlate / 2, theGeometry->heightPillier + theGeometry->heightPlate + 3 * theGeometry->heightBigColumn, -1, 0, &ierr);
    int *filterLeft = malloc(2 * sizeof(int));
    if (filterLeft == NULL) { Error("Memory Allocation Failed."); exit(EXIT_FAILURE); return; }
    filterLeft[0] = 2;
    filterLeft[1] = idRectFilterLeft;
    gmshModelOccCut(plate, 2, filterLeft, 2, NULL ,NULL, NULL, NULL, NULL, -1, 1, 1, &ierr);
    // End : Cut the half of the bridge for the symmetry

    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD)
    {
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
        gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);
        gmshModelGeoMeshSetRecombine(2, 1, 45, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }
    else if (theGeometry->elementType == FEM_TRIANGLE)
    {
        gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
        gmshModelMeshGenerate(2, &ierr);
    }

    //  Plot of Fltk
    gmshFltkInitialize(&ierr);

    free(idWindows);
    idWindows = NULL;
    free(idColumns);
    idColumns = NULL;
    free(idPilliers);
    idPilliers = NULL;
    free(idArcs);
    idArcs = NULL;
    free(idBigColumns);
    idBigColumns = NULL;
    free(idDisks);
    idDisks = NULL;
    free(idCables);
    idCables = NULL;
    free(positionCablesX);
    positionCablesX = NULL;
    free(positionCablesY);
    positionCablesY = NULL;
    free(filterLeft);
    filterLeft = NULL;
}