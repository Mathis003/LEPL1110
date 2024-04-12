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

double *femElasticityForces(femProblem *theProblem)
{
    double *residuals = theProblem->residuals;
    double *soluce    = theProblem->soluce;

    double **A_copy;
    double *B_copy;
    int size;
    
    // Read the system from the file
    femSystemRead(&A_copy, &B_copy, &size, "../../Processing/data/dirichletUnconstrainedSystem.txt");

    for (int i = 0; i < size; i++) { residuals[i] = 0.0; }

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { residuals[i] += A_copy[i][j] * soluce[j]; }
        residuals[i] -= B_copy[i];
    }
    return residuals;
}

int main(int argc, char *argv[])
{
    // Deal with the options arguments
    int opt;
    int resultVisualizer = TRUE;
    int exampleUsage     = FALSE;
    int animation        = FALSE;
    int plot             = TRUE;
    while ((opt = getopt(argc, argv, "areph")) != -1)
    {
        switch (opt)
        {
            case 'a':
                animation = TRUE;
            case 'r':
                resultVisualizer = FALSE;
                break;
            case 'e':
                exampleUsage = TRUE;
                break;
            case 'p':
                plot = FALSE;
                break;
            case 'h':
                printf("Usage: %s [-r] [-e] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -r : Disable the result visualizer\n");
                printf("  -e : Start the program with the example mesh\n");
                printf("  -p : Disable the plot\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-r] [-e] [-h]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }
    
    printf("\n\n");
    printf("    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n");
    printf("    S : Display the matrix \n");
    printf("    X : Display the horizontal forces \n");
    printf("    Y : Display the vertical forces \n");
    printf("\n\n");

    /***************************/
    /* 1 : Lecture des donnees */ 
    /***************************/

    femGeometry *theGeometry = geoGetGeometry();

    femSolverType typeSolver = FEM_BAND; // FEM_FULL or FEM_BAND
    femRenumType renumType   = FEM_XNUM; // FEM_NO or FEM_XNUM or FEM_YNUM (or FEM_RCMK)
    int n;
    double *theSoluce;
    femProblem *theProblem;

    if (exampleUsage == TRUE)
    {
        geoMeshRead("../../Processing/data/mesh_example.txt");
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Processing/data/problem_example.txt", renumType);
        theSoluce = theProblem->soluce;
        n = theGeometry->theNodes->nNodes;
        femSolutiondRead(2 * n, theSoluce, "../../Processing/data/UV_example.txt");
        femElasticityPrint(theProblem);
    }
    else
    {
        geoMeshRead("../../Processing/data/mesh.txt");
        theProblem = femElasticityRead(theGeometry, typeSolver, "../../Processing/data/problem.txt", renumType);
        theSoluce = theProblem->soluce;
        n = theGeometry->theNodes->nNodes;
        femSolutiondRead(2 * n, theSoluce, "../../Processing/data/UV.txt");
        femElasticityPrint(theProblem);
    }


    // Create the solver with the final system to visualize the matrix by pressing 'S'
    femSolver *theSolver = theProblem->solver;
    double **A = getMatrixA(theSolver);
    double *B  = getVectorB(theSolver);
    int size   = theSolver->size;
    femSystemRead(&A, &B, &size, "../../Processing/data/finalSystem.txt");

    femSolverSetSystem(theSolver, A, B);

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
    if (exampleUsage == TRUE) { deformationFactor = 5e4; } // To change the deformation factor
    else                      { deformationFactor = 5e4; } // To change the deformation factor

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

    /******************************/
    /* 5 : Lance le script python */ 
    /******************************/

    if (plot == TRUE)
    {
        char command[100] = "python3 ../src/plot.py";
        if (exampleUsage == TRUE) { strcat(command, " -e"); }
        if (animation == TRUE)    { strcat(command, " -a"); }
        system(command);
    }

    /*********************/
    /* 6 : Visualisation */ 
    /*********************/

    if (!resultVisualizer) { return EXIT_SUCCESS; }

    int mode = 1;
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME + 5];

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
        if (glfwGetKey(window, 'S') == GLFW_PRESS) { mode = 2; }
        if (glfwGetKey(window,'X')  == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window,'Y')  == GLFW_PRESS) { mode = 4;}
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
        else if (mode == 2)
        {
            glColor3f(1.0, 0.0, 0.0);
            glfemPlotSolver(theSolver, theSolver->size, w, h);
        } 
        else if (mode==3)
        {
            glfemPlotField(theGeometry->theElements,forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode==4)
        {
            glfemPlotField(theGeometry->theElements,forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        } else {
            printf("Mode %d not implemented\n", mode);
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