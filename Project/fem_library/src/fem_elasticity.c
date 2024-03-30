#include "../include/fem_elasticity.h"

femProblem *femElasticityCreate(femGeometry *theGeometry, double E, double nu, double rho, double gx, double gy, femElasticCase iCase)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    if (theProblem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
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

    theProblem->planarStrainStress = iCase;
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;

    int nNodes = theGeometry->theNodes->nNodes;
    int size = 2 * nNodes;
    theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
    if (theProblem->constrainedNodes == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    for (int i = 0; i < nNodes; i++)
    {
        theProblem->constrainedNodes[i].type = UNDEFINED;
        theProblem->constrainedNodes[i].nx = NAN;
        theProblem->constrainedNodes[i].ny = NAN;
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
    theProblem->system = femFullSystemCreate(size);
    return theProblem;
}

void femElasticityFree(femProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    for (int i = 0; i < theProblem->nBoundaryConditions; i++) { free(theProblem->conditions[i]); }
    free(theProblem->conditions);
    free(theProblem->soluce);
    free(theProblem->residuals);
    free(theProblem->constrainedNodes);
    free(theProblem);
}

/*
`value2` is only used for `DIRICHLET_XY` and `DIRICHLET_NT` boundary conditions. Otherwise it is ignored and set to NAN.
*/
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value1, double value2)
{
    int iDomain = geoGetDomain(nameDomain);
    if (iDomain == -1) { Error("Undefined domain :-("); }
    value2 = ((type != DIRICHLET_XY) && (type != DIRICHLET_NT)) ? NAN : value2;

    femBoundaryCondition *theBoundary = malloc(sizeof(femBoundaryCondition));
    if (theBoundary == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    theBoundary->domain = theProblem->geometry->theDomains[iDomain];
    theBoundary->value1 = value1;
    theBoundary->value2 = value2;
    theBoundary->type = type;
    theProblem->nBoundaryConditions++;
    int nBoundaryConditions = theProblem->nBoundaryConditions;

    if (theProblem->conditions == NULL)
    {
        theProblem->conditions = malloc(nBoundaryConditions * sizeof(femBoundaryCondition *));
        if (theProblem->conditions == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
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
        else // Need to compute normals
        {
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

femProblem *femElasticityRead(femGeometry *theGeometry, const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) { printf("Error at %s:%d\nUnable to open file %s\n", __FILE__, __LINE__, filename); exit(-1); }
    femProblem *theProblem = malloc(sizeof(femProblem));
    if (theProblem == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;

    int nNodes = theGeometry->theNodes->nNodes;
    int size = 2 * nNodes;
    theProblem->soluce = malloc(size * sizeof(double));
    if (theProblem->soluce == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    theProblem->residuals = malloc(size * sizeof(double));
    if (theProblem->residuals == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
    for (int i = 0; i < size; i++) { theProblem->soluce[i] = 0.0; theProblem->residuals[i] = 0.0; }

    theProblem->constrainedNodes = malloc(nNodes * sizeof(femConstrainedNode));
    if (theProblem->constrainedNodes == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return NULL; }
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
    theProblem->system = femFullSystemCreate(size);

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