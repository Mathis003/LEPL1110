/*
*  main.c
*  Projet 2023-2024
*  Elasticite lineaire plane
*
*  Pre-processeur
*
*  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
*  All rights reserved.
*
*/

#include "../../Project/src/fem.h"
#include "gmsh.h"
#include "glfem.h"

#include <getopt.h>

femElementType elementType  = FEM_TRIANGLE; // FEM_QUAD or FEM_TRIANGLE
femDiscreteType discretType = FEM_DISCRETE_TYPE_LINEAR; // FEM_DISCRETE_TYPE_LINEAR or FEM_DISCRETE_TYPE_QUADRATIC     // Not used yet
femElasticCase theCase      = PLANAR_STRESS; // PLANAR_STRESS or PLANAR_STRAIN or AXISYM (PLANAR_STRESS for our bridge problem)

int main(int argc, char *argv[])
{
    int opt;
    int meshVisualizer       = TRUE;
    int exampleBeam_Example  = FALSE;
    int exampleUForm_Example = FALSE;
    int beam_mesh            = FALSE;
    int bridgeSimplified     = FALSE;
    while ((opt = getopt(argc, argv, "msubh")) != -1)
    {
        switch (opt)
        {
            case 'm':
                meshVisualizer = FALSE;
                break;
            case 's':
                bridgeSimplified = TRUE;
                break;
            case 'u':
                exampleUForm_Example = TRUE;
                break;
            case 'b':
                exampleBeam_Example = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-m] [-s] [-u] [-b] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -m : Disable the mesh visualizer\n");
                printf("  -s : Start the program with the bridge without stay cables and pylon\n");
                printf("  -u : Start the program with the example mesh\n");
                printf("  -b : Start the program with the beam mesh\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-m] [-s] [-u] [-b] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    /************************************/
    /* 1 : Construction de la geometrie */ 
    /************************************/

    geoInitialize();
    femGeometry *theGeometry = geoGetGeometry();
    theGeometry->elementType = elementType;

    if (exampleUForm_Example == TRUE)
    {
        geoMeshGenerateExample(discretType, FALSE);
        geoMeshImport(discretType);

        geoSetDomainName(0, "Symmetry");
        geoSetDomainName(7, "Bottom");
        geoSetDomainName(1, "Top");
        geoMeshWrite("../../Rapport/data/mesh_example.txt", discretType);
    }
    else if (exampleBeam_Example == TRUE)
    {
        geoMeshGenerateExample(discretType, TRUE);
        geoMeshImport(discretType);

        geoSetDomainName(0, "Bottom");
        geoSetDomainName(1, "Right");
        geoSetDomainName(2, "Top");
        geoSetDomainName(3, "Left");
        geoMeshWrite("../../Rapport/data/mesh_beam.txt", discretType);
    }
    else if (bridgeSimplified == TRUE)
    {
        geoMeshGenerate(discretType, bridgeSimplified);
        geoMeshImport(discretType);
        setDomainsName(bridgeSimplified);
        geoMeshWrite("../../Rapport/data/mesh_simplified.txt", discretType);
    }
    else
    {
        geoMeshGenerate(discretType, bridgeSimplified);
        geoMeshImport(discretType);
        setDomainsName(bridgeSimplified);
        geoMeshWrite("../../Rapport/data/mesh.txt", discretType);
    }

    /******************************/
    /* 2 : Definition du probleme */ 
    /******************************/

    femProblem *theProblem;
    double gx = 0.0;
    double gy = -9.81;

    double E = 211.e9;
    double nu = 0.3;
    double rho = 7.85e1;

    theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, theCase, discretType);

    if (exampleUForm_Example == TRUE)
    {
        femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_XY, 0.0, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_XY, 0.0, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_Y, -1e2, NAN);
        femElasticityWrite(theProblem, "../../Rapport/data/problem_example.txt");
    }
    else if (exampleBeam_Example)
    {
        double w = 1e7 / theGeometry->LxPlate; // 1250 kN/m
        femElasticityAddBoundaryCondition(theProblem, "Left", DIRICHLET_XY, 0.0, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_Y, -w, NAN);
        femElasticityWrite(theProblem, "../../Rapport/data/problem_beam.txt");
    }
    else if (bridgeSimplified == TRUE)
    {
        E = 200.e9;
        nu = 0.3;
        rho = 7.85e3;
        theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, theCase, discretType);
        createBoundaryConditions(theProblem, bridgeSimplified);
        femElasticityPrint(theProblem);

        femElasticityWrite(theProblem, "../../Rapport/data/problem_simplified.txt");
    }
    else
    {
        E = 200.e9;
        nu = 0.3;
        rho = 7.85e3;
        theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, theCase, discretType);
        createBoundaryConditions(theProblem, bridgeSimplified);
        femElasticityPrint(theProblem);
        femElasticityWrite(theProblem, "../../Rapport/data/problem.txt");
    }

    /***************************************************/
    /* 3 : Champ de la taille de référence du maillage */
    /*     (uniquement pour la visualisation)          */
    /***************************************************/

    double *meshSizeField = (double *) malloc(theGeometry->theNodes->nNodes * sizeof(double));
    if (meshSizeField == NULL) { printf("Memory allocation error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    
    femNodes *theNodes = theGeometry->theNodes;
    for (int i = 0; i < theNodes->nNodes; ++i) { meshSizeField[i] = theGeometry->geoSize(theNodes->X[i], theNodes->Y[i]); }

    double hMin = femMin(meshSizeField, theNodes->nNodes);
    double hMax = femMax(meshSizeField, theNodes->nNodes);

    printf(" ==== Global requested h : %14.7e \n", theGeometry->defaultSize);
    printf(" ==== Minimum h          : %14.7e \n", hMin);
    printf(" ==== Maximum h          : %14.7e \n", hMax);

    /*********************/
    /* 4 : Visualisation */ 
    /*********************/

    if (!meshVisualizer) { return EXIT_SUCCESS; }

    int mode = 1;
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME+5];

    GLFWwindow *window = glfemInit("EPL1110 : Project 2023-24 ");
    glfwMakeContextCurrent(window);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    do {
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glfemReshapeWindows(window, theGeometry->theNodes, w, h);

        t = glfwGetTime();
        if (glfwGetKey(window, 'D') == GLFW_PRESS) { mode = 0; }
        if (glfwGetKey(window, 'V') == GLFW_PRESS) { mode = 1; }
        if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t; }

        if (t - told > 0.5) { freezingButton = FALSE; }
        if (mode == 1)
        {
            glfemPlotField(theGeometry->theElements, meshSizeField);
            glfemPlotMesh(theGeometry->theElements);
            sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        }
        else if (mode == 0)
        {
            domain %= theGeometry->nDomains;
            glfemPlotDomain(theGeometry->theDomains[domain]);
            sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
    // Check if the ESC key was pressed or the window was closed

    free(meshSizeField); meshSizeField = NULL;
    femElasticityFree(theProblem); theProblem = NULL;
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}