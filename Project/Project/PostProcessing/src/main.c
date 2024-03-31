/*
*  main.c
*  Projet 2023-2024
*  Elasticite lineaire plane
*
*  Postprocesseur
*
*  Copyright (C) 2023 UCL-IMMC : Vincent Legat, Miguel De Le Court
*  All rights reserved.
*
*/

#include <getopt.h>

#include "../../glfem_library/glfem.h"
#include "../../../fem_library/include/fem_geometry.h"

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int resultVisualizer = TRUE;
    while ((opt = getopt(argc, argv, "rh")) != -1)
    {
        switch (opt)
        {
            case 'r':
                resultVisualizer = FALSE;
                break;
            case 'h':
                printf("Usage: %s [-r] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -r : Disable the result visualizer\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-r] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n\n\n");

    /***************************/
    /* 1 : Lecture des donnees */ 
    /***************************/

    femGeometry *theGeometry = geoGetGeometry();
    geoMeshRead("../../../data/mesh.txt");

    femSolverType typeSolver = FEM_FULL;
    femProblem *theProblem = femElasticityRead(theGeometry, typeSolver, "../../../data/problem.txt");
    double *theSoluce = theProblem->soluce;
    int n = theGeometry->theNodes->nNodes;
    femSolutiondRead(2 * n, theSoluce, "../../../data/UV.txt");
    femElasticityPrint(theProblem);

    /****************************************************/
    /* 2 : Deformation du maillage pour le plot final   */ 
    /*     Creation du champ de la norme du deplacement */
    /****************************************************/

    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    if (normDisplacement == NULL) { printf("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }

    for (int i = 0; i < n; i++)
    {
        theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
        theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] + theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
    }

    double hMin = femMin(normDisplacement, n);
    double hMax = femMax(normDisplacement, n);

    printf(" ==== Minimum displacement          : %14.7e \n", hMin);
    printf(" ==== Maximum displacement          : %14.7e \n", hMax);

    /*********************/
    /* 3 : Visualisation */ 
    /*********************/

    if (!resultVisualizer) { return EXIT_SUCCESS; }

    int mode = 1;
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];

    GLFWwindow *window = glfemInit("EPL1110 : Project 2022-23 ");
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
        if (glfwGetKey(window, 'N') == GLFW_PRESS && freezingButton == FALSE)
        {
            domain++;
            freezingButton = TRUE;
            told = t;
        }

        if (t - told > 0.5) { freezingButton = FALSE; }
        if (mode == 1)
        {
            glfemPlotField(theGeometry->theElements, normDisplacement);
            glfemPlotMesh(theGeometry->theElements);
            sprintf(theMessage, "Number of elements : %d ", theGeometry->theElements->nElem);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        }
        else if (mode == 0)
        {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain(theGeometry->theDomains[domain]);
            sprintf(theMessage, "%s : %d ", theGeometry->theDomains[domain]->name, domain);
            glColor3f(1.0, 0.0, 0.0);
            glfemMessage(theMessage);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    femElasticityFree(theProblem);
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}