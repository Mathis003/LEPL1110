/*
 *  main.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2024 UCL-EPL : Vincent Legat    
 *  All rights reserved.
 *
 */


#include "glfem.h"
#include <time.h>



int main(void)
{  
    printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    F-B-I : Full solver - Band solver - Iterative solver \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_YNUM;
    char meshFileName[] = "../data/triangles101.txt";  
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles101.txt") 
    // par :
    // ("..\\data\\triangles101.txt") 
    //
    // Sorry for the inconvenience :-)
    // On reflechit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
    
    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
    clock_t tic = clock();
    femDiffusionCompute(theProblem);  
    femSolverPrintInfos(theProblem->solver);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    

    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;
    femMesh *theMesh = theProblem->geo->theElements;

    GLFWwindow* window = glfemInit("LEPL1110 : Band Solver ");
    glfwMakeContextCurrent(window);
     
    do 
    {
        int testConvergence,w,h;
        char theMessage[256];
        sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        glfwGetFramebufferSize(window,&w,&h);
      
        if (option == 1) {
            glfemReshapeWindows(theMesh,w,h);
            glfemPlotField(theMesh,theProblem->soluce);   }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); }
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
      
        if (solverType != newSolverType || renumType != newRenumType) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            theMesh = theProblem->geo->theElements;
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }
        if (glfwGetKey(window,'V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey(window,'B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        if (glfwGetKey(window,'I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
       
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
}