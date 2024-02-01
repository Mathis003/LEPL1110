/*
 *  glfem.h
 *  Library for EPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  GLFW  http://www.glfw.org/ (version utilisée 3.3.2)
 *  BOV   https://git.immc.ucl.ac.be/hextreme/NGP/-/tree/master/deps/BOV
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#include <stdlib.h>
#if __has_include("BOV.h")
  #include "BOV.h"
#endif
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include "fem.h"


static float GLFEM_BLACK[4] = {0.0,0.0,0.0,1.0};
static float GLFEM_BLUE[4]  = {0.0,0.0,1.0,1.0};
static float GLFEM_RED[4]   = {1.0,0.0,0.0,1.0};


void          glfemWindowCreate(const char *windowName,int w,int h,int n,double *x,double *y); 
void          glfemWindowFree(); 
void          glfemWindowUpdate();
int           glfemWindowShouldClose();
void          glfemReshape(double x[3],double y[3], int n);

void          glfemSetColor(float color[4]);
void          glfemSetTextColor(float color[4]);
void          glfemDrawMessage(char *message, double pos[2]);
void          glfemDrawNodes(double *x, double *y, int n);
void          glfemDrawElement(double *x, double *y, int n);
void          glfemDrawSolution(double *x, double *y, double* u, int n);

static void   glfemKeyCallback(GLFWwindow* self,int key,int scancode,int action,int mods);


#endif


