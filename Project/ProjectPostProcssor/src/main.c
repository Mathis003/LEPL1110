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

#include "../../Project/src/fem.h"
#include "../../Project/src/homework.c"

#include "../../ProjectPreProcessor/src/glfem.h"

femDiscreteType discretType = FEM_DISCRETE_TYPE_LINEAR; // FEM_DISCRETE_TYPE_LINEAR or FEM_DISCRETE_TYPE_QUADRATIC // Not used yet
femSolverType typeSolver    = FEM_BAND;  // FEM_FULL or FEM_BAND
femRenumType renumType      = FEM_RCMK;  // FEM_NO or FEM_XNUM or FEM_YNUM or FEM_RCMK

double constFunct(double x, double y) { return 1.0; }

double *femElasticityVonMises(femProblem *theProblem, double *sigmaXX, double *sigmaYY, double *sigmaXY, int nNodes)
{
    double *vonMises = (double *) malloc(theProblem->geometry->theNodes->nNodes * sizeof(double));
    if (vonMises == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }

    for (int i = 0; i < nNodes; i++)
    {
        double trace = sigmaXX[i] + sigmaYY[i];
        sigmaXX[i] -= trace / 3.0;
        sigmaYY[i] -= trace / 3.0;
        double final_value = sigmaXX[i] * sigmaXX[i] + sigmaYY[i] * sigmaYY[i] + 2 * sigmaXY[i] * sigmaXY[i];
        final_value = 3.0 * final_value / 2.0;
        vonMises[i] = sqrt(final_value);
    }
    return vonMises;
}

void femElasticitySigma(femProblem *theProblem, double *sigmaXX, double *sigmaYY, double *sigmaXY)
{
    femIntegration *theRule  = theProblem->rule;
    femDiscrete *theSpace    = theProblem->space;
    femGeometry *theGeometry = theProblem->geometry;
    femNodes *theNodes       = theGeometry->theNodes;
    femMesh *theMesh         = theGeometry->theElements;

    int nLocal = theSpace->n;

    femSolver *theSolver;
    double **M1, **M2, **M3, *B1, *B2, *B3, *epsilonXX, *epsilonYY, *epsilonXY, *theSoluce;
    double x[nLocal], y[nLocal], phi[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    double u[nLocal], v[nLocal], xsi, eta, weight, jac, weightedJac, dxdxsi, dxdeta, dydxsi, dydeta, a, b, c;
    int iElem, iInteg, iEdge, i, j, map[nLocal], nNodes;

    nNodes    = theNodes->nNodes;
    theSoluce = theProblem->soluce;    
    theSolver = femSolverFullCreate(nNodes);

    M1 = femSolverGetA(theSolver);
    B1 = femSolverGetB(theSolver);

    B2 = (double *) malloc(nNodes * sizeof(double));
    B3 = (double *) malloc(nNodes * sizeof(double));

    if (B2 == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (B3 == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }

    for (i = 0; i < nNodes; i++) { B2[i] = 0.0; B3[i] = 0.0; }
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (i = 0; i < theSpace->n; i++)
        {
            map[i] = theMesh->elem[iElem * nLocal + i];
            x[i] = theNodes->X[map[i]];
            y[i] = theNodes->Y[map[i]];
            u[i] = theSoluce[2 * map[i]];
            v[i] = theSoluce[2 * map[i] + 1];
        }

        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {
            xsi    = theRule->xsi[iInteg];
            eta    = theRule->eta[iInteg];
            weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            dxdxsi = 0.0; dydxsi = 0.0;
            dxdeta = 0.0; dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }

            jac = dxdxsi * dydeta - dxdeta * dydxsi;
            if (jac < 0.0) { printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n"); }
            jac = fabs(jac);

            for (i = 0; i < theSpace->n; i++)
            {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            weightedJac = jac * weight;

            for (i = 0; i < theSpace->n; i++)
            {
                for (j = 0; j < theSpace->n; j++) { M1[map[i]][map[j]] += phi[i] * phi[j] * weightedJac; }
                B1[map[i]] += phi[i] * dphidx[i] * u[i] * weightedJac;
                B2[map[i]] += phi[i] * dphidy[i] * v[i] * weightedJac;
                B3[map[i]] += phi[i] * (dphidy[i] * u[i] + dphidx[i] * v[i]) * weightedJac / 2.0;
            }
        }
    }

    M2 = (double **) malloc(sizeof(double *) * nNodes);
    M3 = (double **) malloc(sizeof(double *) * nNodes);

    if (M2 == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
    if (M3 == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

    for (i = 0; i < nNodes; i++)
    {
        M2[i] = (double *) malloc(nNodes * sizeof(double));
        M3[i] = (double *) malloc(nNodes * sizeof(double));

        if (M2[i] == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }
        if (M3[i] == NULL) { Error("Memory allocation error\n"); exit(EXIT_FAILURE); return; }

        memcpy(M2[i], M1[i], nNodes * sizeof(double));
        memcpy(M3[i], M1[i], nNodes * sizeof(double));
    }

    epsilonXX = (double *) malloc(nNodes * sizeof(double));
    epsilonYY = (double *) malloc(nNodes * sizeof(double));
    epsilonXY = (double *) malloc(nNodes * sizeof(double));

    if (epsilonXX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (epsilonYY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }
    if (epsilonXY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return; }

    femSolverEliminate(theSolver);
    memcpy(epsilonXX, femSolverGetB(theSolver), nNodes * sizeof(double));

    femSolverSet(theSolver, M2, B2);
    femSolverEliminate(theSolver);
    memcpy(epsilonYY, femSolverGetB(theSolver), nNodes * sizeof(double));

    femSolverSet(theSolver, M3, B3);
    femSolverEliminate(theSolver);
    memcpy(epsilonXY, femSolverGetB(theSolver), nNodes * sizeof(double));

    a = theProblem->A;
    b = theProblem->B;
    c = theProblem->C;

    for (i = 0; i < nNodes; i++)
    {
        sigmaXX[i] = a * epsilonXX[i] + b * epsilonYY[i];
        sigmaXY[i] = 2 * c * epsilonXY[i];
        sigmaYY[i] = b * epsilonXX[i] + a * epsilonYY[i];
    }

    free(epsilonXX); epsilonXX = NULL;
    free(epsilonYY); epsilonYY = NULL;
    free(epsilonXY); epsilonXY = NULL;
    for (i = 0; i < nNodes; i++) { free(M2[i]); M2[i] = NULL; free(M3[i]); M3[i] = NULL; }
    free(M2); M2 = NULL;
    free(M3); M3 = NULL;
    free(B2); B2 = NULL;
    free(B3); B3 = NULL;
}

double *femElasticityForces(femProblem *theProblem)
{
    femNodes *theNodes   = theProblem->geometry->theNodes;
    double *theSoluce    = theProblem->soluce;
    femSolver *theSolver = theProblem->solver;
    int size = theSolver->size;

    double *soluce       = (double *) malloc(size * sizeof(double));
    double *residuals    = (double *) malloc(size * sizeof(double));
    int *inverted_number = (int *) malloc(theNodes->nNodes * sizeof(int));
    if (soluce == NULL)          { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }
    if (residuals == NULL)       { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }
    if (inverted_number == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return NULL; }

    for (int i = 0; i < size; i++) { residuals[i] = 0.0; }

    // Numerote the nodes of the solution
    for (int i = 0; i < theNodes->nNodes; i++)
    {
        soluce[2 * theNodes->number[i]]     = theSoluce[2 * i];
        soluce[2 * theNodes->number[i] + 1] = theSoluce[2 * i + 1];
    }

    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem, 1.0, -1);

    double **A = femSolverGetA(theSolver);
    double *B  = femSolverGetB(theSolver);

    femSolverGetResidual(theSolver, residuals, soluce);

    // Remove the numerotation of the nodes of the residuals
    for (int i = 0; i < theNodes->nNodes; i++) { inverted_number[theNodes->number[i]] = i; }

    for (int i = 0; i < theNodes->nNodes; i++)
    {
        theProblem->residuals[2 * inverted_number[i]]     = residuals[2 * i];
        theProblem->residuals[2 * inverted_number[i] + 1] = residuals[2 * i + 1];
    }

    // Apply the Dirichlet boundary conditions to the system to plot the final system
    femElasticityApplyDirichlet(theProblem, 1.0);

    free(residuals); residuals = NULL;
    free(soluce); soluce = NULL;
    free(inverted_number); inverted_number = NULL;
    return theProblem->residuals;
}

int main(int argc, char *argv[])
{
    int opt;
    int resultVisualizer   = TRUE;
    int example_UForm      = FALSE;
    int example_beam       = FALSE;
    int example_simplified = FALSE;
    int animation          = FALSE;
    int plot               = TRUE;
    while ((opt = getopt(argc, argv, "rpaubsh")) != -1)
    {
        switch (opt)
        {
            case 'r':
                resultVisualizer = FALSE;
                break;
            case 'p':
                plot = FALSE;
                break;
            case 'a':
                animation = TRUE;
                break;
            case 'u':
                example_UForm = TRUE;
                break;
            case 'b':
                example_beam = TRUE;
                break;
            case 's':
                example_simplified = TRUE;
                break;
            case 'h':
                printf("Usage: %s [-r] [-p] [-a] [-u] [-b] [-s] [-h]\n", argv[0]);
                printf("Options:\n");
                printf("  -r : Disable the result visualizer\n");
                printf("  -p : Disable the plot\n");
                printf("  -a : Enable the animation\n");
                printf("  -u : Start the program with the example UForm\n");
                printf("  -b : Start the program with the example beam\n");
                printf("  -s : Start the program with the example simplified\n");
                printf("  -h : Display this help message\n");
                return EXIT_SUCCESS;
            default:
                fprintf(stderr, "Usage: %s [-r] [-p] [-a] [-u] [-b] [-s] [-h]\n", argv[0]);
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
    printf("    U : Display the sigmaXX \n");
    printf("    I : Display the sigmaYY \n");
    printf("    O : Display the sigmaXY \n");
    printf("    P : Display the von Mises \n");

    printf("\n\n");

    /***************************/
    /* 1 : Lecture des donnees */ 
    /***************************/

    femGeometry *theGeometry = geoGetGeometry();

    char *meshPath;
    char *problemPath;
    char *resultPath;
    if      (example_UForm)      { meshPath = "../../Rapport/data/mesh_example.txt";     problemPath = "../../Rapport/data/problem_example.txt";     resultPath = "../../Rapport/data/UV_example.txt"; }
    else if (example_beam)       { meshPath = "../../Rapport/data/mesh_beam.txt";        problemPath = "../../Rapport/data/problem_beam.txt";        resultPath = "../../Rapport/data/UV_beam.txt"; }
    else if (example_simplified) { meshPath = "../../Rapport/data/mesh_simplified.txt";  problemPath = "../../Rapport/data/problem_simplified.txt";  resultPath = "../../Rapport/data/UV_simplified.txt"; }
    else                         { meshPath = "../../Rapport/data/mesh.txt";             problemPath = "../../Rapport/data/problem.txt";             resultPath = "../../Rapport/data/UV.txt"; }

    geoMeshRead(meshPath, discretType);
    femMeshRenumber(theGeometry->theElements, renumType);
    femProblem *theProblem = femElasticityRead(theGeometry, typeSolver, problemPath, renumType, discretType);
    double *theSoluce = theProblem->soluce;
    int n = theGeometry->theNodes->nNodes;
    femSolutiondRead(2 * n, theSoluce, resultPath);
    femElasticityPrint(theProblem);

    double rho = theProblem->rho;
    double gy = theProblem->gy;

    /*********************************/
    /* 2 : Calcul des forces         */
    /*     Création du système final */
    /*********************************/
    
    double *theForces = femElasticityForces(theProblem);
    double area       = femElasticityIntegrate(theProblem, constFunct);

    femNodes *theNodes = theGeometry->theNodes;
    int nNodes = theNodes->nNodes;

    double *sigmaXX = (double *) malloc(nNodes * sizeof(double));
    double *sigmaYY = (double *) malloc(nNodes * sizeof(double));
    double *sigmaXY = (double *) malloc(nNodes * sizeof(double));

    if (sigmaXX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    if (sigmaYY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    if (sigmaXY == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }

    double *vonMises = NULL;

    // Only for the complete bridge
    if (example_UForm == FALSE && example_beam == FALSE && example_simplified == FALSE)
    {
        femElasticitySigma(theProblem, sigmaXX, sigmaYY, sigmaXY);

        vonMises = femElasticityVonMises(theProblem, sigmaXX, sigmaYY, sigmaXY, nNodes);

        FILE *file = fopen("../../Rapport/data/maxConstrainBridge.txt", "w");
        if (file == NULL) { Error("File opening error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
        for (int i = 0; i < nNodes; i++) { fprintf(file, "%f\n", vonMises[i]); }
        fclose(file);
    }

    /****************************************************/
    /* 3 : Deformation du maillage pour le plot final   */ 
    /*     Creation du champ de la norme du deplacement */
    /****************************************************/

    double deformationFactor;
    if      (example_UForm == TRUE)      { deformationFactor = 1e6; }
    else if (example_beam == TRUE)       { deformationFactor = 1e2; }
    else if (example_simplified == TRUE) { deformationFactor = 1e4; }
    else                                 { deformationFactor = 1e4; }

    double *normDisplacement = (double *) malloc(nNodes * sizeof(double));
    double *forcesX = (double *) malloc(nNodes * sizeof(double));
    double *forcesY = (double *) malloc(nNodes * sizeof(double));

    if (normDisplacement == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
    if (forcesX == NULL) { Error("Allocation Error\n"); exit(EXIT_FAILURE); return EXIT_FAILURE; }
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
    for (int i = 0; i < nNodes; i++)
    {
        theGlobalForce[0] += theForces[2 * i];
        theGlobalForce[1] += theForces[2 * i + 1];
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

        if (animation) { strcat(command, " -a"); }

        if      (example_UForm)      { strcat(command, " -u"); }
        else if (example_beam)       { strcat(command, " -b"); }
        else if (example_simplified) { strcat(command, " -s"); }
        
        system(command);
    }

    /*********************/
    /* 6 : Visualisation */ 
    /*********************/

    if (!resultVisualizer)
    {
        free(normDisplacement); normDisplacement = NULL;
        femElasticityFree(theProblem); theProblem = NULL;
        geoFree();
        exit(EXIT_SUCCESS); return EXIT_SUCCESS;
    }

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
        if (glfwGetKey(window, 'X') == GLFW_PRESS) { mode = 3; }
        if (glfwGetKey(window, 'Y') == GLFW_PRESS) { mode = 4; }
        if (glfwGetKey(window, 'U') == GLFW_PRESS) { mode = 5; }
        if (glfwGetKey(window, 'I') == GLFW_PRESS) { mode = 6; }
        if (glfwGetKey(window, 'O') == GLFW_PRESS) { mode = 7; }
        if (glfwGetKey(window, 'P') == GLFW_PRESS) { mode = 8; }
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
            glfemPlotSolver(theProblem->solver, theProblem->solver->size, w, h);
        }
        else if (mode == 3)
        {
            glfemPlotField(theGeometry->theElements, forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 4)
        {
            glfemPlotField(theGeometry->theElements, forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 5)
        {
            glfemPlotField(theGeometry->theElements, sigmaXX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 6)
        {
            glfemPlotField(theGeometry->theElements, sigmaYY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 7)
        {
            glfemPlotField(theGeometry->theElements, sigmaXY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else if (mode == 8)
        {
            glfemPlotField(theGeometry->theElements, vonMises);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);
        }
        else { printf("Mode %d not implemented\n", mode); }  

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1);
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement); normDisplacement = NULL;
    femElasticityFree(theProblem); theProblem = NULL;
    geoFree();
    glfwTerminate();

    exit(EXIT_SUCCESS);
    return EXIT_SUCCESS;
}