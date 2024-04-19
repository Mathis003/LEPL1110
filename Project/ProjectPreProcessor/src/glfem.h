/*
*  glfem.h
*  Library for LEPL1110 : Finite Elements for dummies
*
*  Copyright (C) 2017 UCL-EPL : Vincent Legat
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
#include "../../Project/src/fem.h"

void glfemDrawColorElement(float *x, float *y, double *u, int n);
void glfemDrawElement(float *x, float *y, int n);
void glfemDrawNodes(double *x, double *y, int n);
int glfemGetAction(void);

void glfemMatrix(double **A, int n, int w, int h);
void glfemPlotSolver(femSolver *mySolver, int n, int w, int h);
void glfemReshapeWindows(GLFWwindow *window, femNodes *theNodes, int width, int height);
void glfemPlotField(femMesh *theMesh, double *u);
void glfemPlotMesh(femMesh *theMesh);
void glfemPlotDomain(femDomain *theDomain);

void glfemMessage(char *message);
void glfemDrawMessage(int h, int v, char *message);
void glfemSetRasterSize(int width, int height);

GLFWwindow *glfemInit(char *windowName);
static void glfemKeyCallback(GLFWwindow *self, int key, int scancode, int action, int mods);

void scroll_callback(GLFWwindow *window, double xoffset, double yoffset);
void mouse_button_callback(GLFWwindow *window, int button, int action, int mods);

#endif