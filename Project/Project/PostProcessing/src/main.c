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

#include "../../../fem_library/include/fem.h"
#include "../../../glfem_library/glfem.h"

double constFunct(double x, double y) { return 1.0; }

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int resultVisualizer = TRUE;
    int exampleUsage     = FALSE;
    while ((opt = getopt(argc, argv, "reh")) != -1)
    {
        switch (opt)
        {
            case 'r':
                resultVisualizer = FALSE;
                break;
            case 'e':
                exampleUsage = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-r] [-e] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -r : Disable the result visualizer\n");
                printf("  -e : Start the program with the example mesh\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-r] [-e] [-h]\n", argv[0]);
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

    int n;
    double *theSoluce;
    femProblem *theProblem;
    if (exampleUsage == TRUE)
    {
        geoMeshRead("../../Processing/data/mesh_example.txt");
        femSolverType typeSolver = FEM_FULL;
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Processing/data/problem_example.txt");
        theSoluce = theProblem->soluce;
        n = theGeometry->theNodes->nNodes;
        femSolutiondRead(2 * n, theSoluce, "../../Processing/data/UV_example.txt");
        femElasticityPrint(theProblem);
    }
    else
    {
        geoMeshRead("../../Processing/data/mesh.txt");
        femSolverType typeSolver = FEM_FULL;
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Processing/data/problem.txt");
        theSoluce = theProblem->soluce;
        n = theGeometry->theNodes->nNodes;
        femSolutiondRead(2 * n, theSoluce, "../../Processing/data/UV.txt");
        femElasticityPrint(theProblem);
    }

    double rho = theProblem->rho;
    double gy = theProblem->gy;

    /*************************/
    /* 2 : Calcul des forces */
    /*************************/

    double *theForces = femElasticityForces(theProblem);
    double area       = femElasticityIntegrate(theProblem, constFunct);

    /****************************************************/
    /* 3 : Deformation du maillage pour le plot final   */ 
    /*     Creation du champ de la norme du deplacement */
    /****************************************************/

    femNodes *theNodes = theGeometry->theNodes;

    double deformationFactor;
    if (exampleUsage == TRUE) { deformationFactor = 1e5; }
    else                      { deformationFactor = 1e4; }

    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    if (normDisplacement == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    if (forcesX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    if (forcesY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }

    for (int i = 0; i < n; i++)
    {
        theNodes->X[i] += theSoluce[2 * i + 0] * deformationFactor;
        theNodes->Y[i] += theSoluce[2 * i + 1] * deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2 * i + 0] * theSoluce[2 * i + 0] + theSoluce[2 * i + 1] * theSoluce[2 * i + 1]);
        forcesX[i] = theForces[2 * i + 0];
        forcesY[i] = theForces[2 * i + 1];
    }

    double hMin = femMin(normDisplacement, n);
    double hMax = femMax(normDisplacement, n);

    printf(" ==== Minimum displacement          : %14.7e \n", hMin);
    printf(" ==== Maximum displacement          : %14.7e \n", hMax);

    /*********************************************/
    /* 4 : Calcul de la force globale resultante */ 
    /*********************************************/

    double theGlobalForce[2] = {0, 0};
    for (int i = 0; i < theProblem->geometry->theNodes->nNodes; i++)
    {
        theGlobalForce[0] += theForces[2*i+0];
        theGlobalForce[1] += theForces[2*i+1];
    }

    printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * gy);

    /*********************/
    /* 5 : Visualisation */ 
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