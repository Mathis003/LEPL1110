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

#include "../../../fem_library/include/fem.h"
#include "../../../fem_library/include/gmsh.h"
#include "../../../glfem_library/glfem.h"
#include "../../var.h"

#include <getopt.h>

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int meshVisualizer = TRUE;
    int exampleUsage   = FALSE;
    while ((opt = getopt(argc, argv, "meh")) != -1)
    {
        switch (opt)
        {
            case 'm':
                meshVisualizer = FALSE;
                break;
            case 'e':
                exampleUsage = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-m] [-e] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -e : Start the program with the example mesh\n");
                printf("  -m : Disable the mesh visualizer\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-m] [-e] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    /************************************/
    /* 1 : Construction de la geometrie */ 
    /************************************/

    geoInitialize();
    femGeometry *theGeometry = geoGetGeometry();
    theGeometry->elementType = elementType;

    if (exampleUsage == TRUE)
    {
        geoMeshGenerateExample(discretType);
        geoMeshImport(discretType);
        geoSetDomainName(0, "Symmetry");
        geoSetDomainName(7, "Bottom");
        geoSetDomainName(1, "Top");
        geoMeshWrite("../../Processing/data/mesh_example.txt", discretType);
    }
    else
    {
        geoMeshGenerate(discretType);
        geoMeshImport(discretType);
        setDomainsName();
        geoMeshWrite("../../Processing/data/mesh.txt", discretType);
    }

    /******************************/
    /* 2 : Definition du probleme */ 
    /******************************/

    femProblem *theProblem;
    double gx = 0.0;
    double gy = -9.81;

    if (exampleUsage == TRUE)
    {
        double E_example = 211.e9;
        double nu_example = 0.3;
        double rho_example = 7.85e3;
        theProblem = femElasticityCreate(theGeometry, E_example, nu_example, rho_example, gx, gy, theCase, discretType);
        // femElasticityAddBoundaryCondition(theProblem, "Symmetry", NEUMANN_N, 1e3, NAN);
        // femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_NT, 0.0, 0.0);
        // femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_N, 1e4, NAN);
        // femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_NT, 0.0, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Symmetry", DIRICHLET_XY, 0.0, 0.0);
        femElasticityAddBoundaryCondition(theProblem, "Bottom", DIRICHLET_XY, 0.0, 0.0);
        //femElasticityAddBoundaryCondition(theProblem, "Top", NEUMANN_Y, -5e1, NAN);

        femElasticityPrint(theProblem);
        femElasticityWrite(theProblem, "../../Processing/data/problem_example.txt");
    }
    else
    {
        double E_steel = 200.e9;
        double E_reinforced_concrete = 35e9;
        double nu_steel = 0.3;
        double nu_reinforced_concrete = 0.2;
        double rho_steel = 7.85e3;
        double rho_reinforced_concrete = 2.5e3;
        theProblem = femElasticityCreate(theGeometry, E_steel, nu_steel, rho_steel, gx, gy, theCase, discretType);
        createBoundaryConditions(theProblem);
        femElasticityPrint(theProblem);
        femElasticityWrite(theProblem, "../../Processing/data/problem.txt");
    }

    /***************************************************/
    /* 3 : Champ de la taille de référence du maillage */
    /*     (uniquement pour la visualisation)          */
    /***************************************************/

    double *meshSizeField = malloc(theGeometry->theNodes->nNodes * sizeof(double));
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

    free(meshSizeField);
    meshSizeField = NULL;
    femElasticityFree(theProblem);
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}