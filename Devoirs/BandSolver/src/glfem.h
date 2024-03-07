/*
 *  glfem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2022 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 3.1)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "fem.h"


void        glfemDrawColorElement(float *x, float *y, double *u, int n);
void 		glfemDrawElement(float *x, float *y, int n);
void 		glfemDrawNodes(double* x, double* y,int n);

void 		glfemReshapeWindows(femMesh *theMesh, int width, int heigh);
void 		glfemPlotField(femMesh *theMesh, double *u);
void 		glfemPlotMesh(femMesh *theMesh);

// void 		glfemPlotEdges(femEdges *theEdges);
// void 		glfemPlotBnd(femEdges *theEdges);

void        glfemMatrix(double **A, int size, int width, int heigh);
void        glfemPlotSolver(femSolver *theSolver, int size, int width, int heigh);

void 		glfemMessage(char *message);
void 		glfemDrawMessage(int h, int v, char *message);
void 		glfemSetRasterSize(int width, int height);
GLFWwindow* glfemInit(char *windowName);

#endif