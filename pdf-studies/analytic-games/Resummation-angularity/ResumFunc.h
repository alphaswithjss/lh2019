#ifndef _RESUMFUNC_H_INCLUDED_
#define _RESUMFUNC_H_INCLUDED_

#include "Definitions.h"

double KLL(double mu2, double mu02, double Gamma, double alphas);
double KNLL(double mu2, double mu02, double Gamma, double gamma0, double alphas);
double KNNLL(double mu2, double mu02, double Gamma, double gamma0, double gamma1, double alphas);

double wLL(double mu2, double mu02, double Gamma, double alphas);
double wNLL(double mu2, double mu02, double Gamma, double alphas);
double wNNLL(double mu2, double mu02, double Gamma, double alphas);

double RLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas);
double RNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double alphas);
double RNNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double gamma1, double alphas);

double RpLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas);
double RppLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas);
double RpppLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas);
double RpNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double alphas);
double RppNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double alphas);
double RpNNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double gamma1, double alphas);

double exp1p(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0);
double exp2p(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0, double gamma1, double C1=0.);
double exp1(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0);
double exp2(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0, double gamma1, double C1=0.);
double exp3(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0, double gamma1, double C1=0.);

#endif
