// Made by Vincent Theeuwes
// Last modified: 05.08.2014

#ifndef _DEFINITIONS_H_INCLUDED_
#define _DEFINITIONS_H_INCLUDED_

#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <iomanip>
#include <cstring>
using namespace std;

extern const double nl, Tf, CA, CF, pi, GammaEuler, Zeta3, Zeta2;

extern int LLorder, fixedorder;
extern bool Laplace;

extern bool end_point, F_HypGeo;

extern int f_gluon, exp_coef;

extern double alpha_Obs;

extern int test;

extern double Grid[10000];

extern double ppar;

extern int DeltaPS_flag;

// extern char* wd;

#endif
