/*
*  main.c
*  Projet 2023-2024
*  Elasticite lineaire plane
*
*  Preprocesseur
*
*  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
*  All rights reserved.
*
*/

#include <getopt.h>
#include <stdbool.h>
#include "glfem.h"

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    bool meshVisualizer = true;
    while ((opt = getopt(argc, argv, "mh")) != -1)
    {
        switch (opt)
        {
            case 'm':
                meshVisualizer = false;
                break;
            case 'h':
                printf("Usage: %s [-m] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -m : Disable the mesh visualizer\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-m] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    /************************************/
    /* 1 : Construction de la geometrie */ 
    /************************************/

    geoInitialize();
    femGeo *theGeometry = geoGetGeometry();

    /* Option 1 : Utilisation de GMSH avec OpenCascade */

    // theGeometry->h = 0.05;
    // geoMeshGenerate();

    /* Option 2 : Utilisation de GMSH directement */

    // theGeometry->h = 0.05;
    // geoMeshGenerateGeo();

    /* Option 3 : Lecture d'un fichier .geo */

    theGeometry->h = 0.05;
    geoMeshGenerateGeoFile("../data/mesh.geo");

    /* Option 4 : Lecture d'un fichier .msh */

    // geoMeshGenerateMshFile("../data/mesh.msh");

    geoMeshImport();

    geoSetDomainName(0, "Something");
    geoSetDomainName(1, "SomethingElse");

    // Write the mesh to a file in this directory and in the Project directory
    geoMeshWrite("../data/mesh.txt");
    geoMeshWrite("../../Project/data/mesh.txt");

    /******************************/
    /* 2 : Definition du probleme */ 
    /******************************/

    double E = 211.e9;
    double nu = 0.3;
    double rho = 7.85e3;
    double gx = 0;
    double gy = -9.81;

    femProblem *theProblem = femElasticityCreate(theGeometry, E, nu, rho, gx, gy, PLANAR_STRAIN);

    femElasticityAddBoundaryCondition(theProblem, "Something", DIRICHLET_XY, 0.0, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "SomethingElse", DIRICHLET_Y, 0.0, NAN);

    femElasticityPrint(theProblem);

    // Write the problem to a file in this directory and in the Project directory
    femElasticityWrite(theProblem, "../data/problem.txt");
    femElasticityWrite(theProblem, "../../Project/data/problem.txt");

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

    printf(" ==== Global requested h : %14.7e \n", theGeometry->h);
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
    char theMessage[MAXNAME];

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
    femElasticityFree(theProblem);
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}