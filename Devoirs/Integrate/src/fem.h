/*
 *  fem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2016 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



double       femMin(double *x, int n);
double       femMax(double *x, int n);
void         femError(char *text, int line, char *file);
void         femErrorScan(int test, int line, char *file);
void         femWarning(char *text, int line, char *file);



#endif
