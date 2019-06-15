#ifndef _INTEGRAND_H_INCLUDED_
#define _INTEGRAND_H_INCLUDED_

#include "Definitions.h"

typedef struct{
	double alphas, rscaleR, rscaleQ, tau;
} para;

typedef struct{
	double alphas, rscaleR, rscaleQ, tau, zcut, beta;
} paraSD;

double Gamma(double z);


double dQCDSD_cumu(void *params);
double dQCD_cumu(void *params);

double ExpansionSD_cumu(void *params);
double Expansion_cumu(void *params);

#endif
