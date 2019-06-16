#ifndef __QCD_CONSTANTS_HH__
#define __QCD_CONSTANTS_HH__

#include <cmath>

//----------------------------------------------------------------------
// physics constants
//----------------------------------------------------------------------
// const double CF = 1.33333333333333333333333;
// const double CA = 3.0;
// const double TR = 0.5;
// const unsigned int nf = 5;
// 
// const double b0 = (11*CA-4*nf*TR)/(12*M_PI);
// const double b1 = (17*CA*CA-5*CA*nf-3*CF*nf)/(24*M_PI*M_PI);
// const double K  = (67.0/18.0-M_PI*M_PI/6.0)*CA - 5.0/9.0*nf;
// 
// const double Bq = -3.0/4.0;
// const double Bg = (-11*CA+4*nf*TR)/(12*CA);

extern double CF;       // default: 4.0/3
extern double CA;       // default: 3.0
extern double TR;       // default: 0.5
extern unsigned int nf; // default 5

extern double b0; // = (11*CA-4*nf*TR)/(12*M_PI);
extern double b1; // = (17*CA*CA-5*CA*nf-3*CF*nf)/(24*M_PI*M_PI);
extern double b2; // = (2857./54.*CA*CA*CA+TR*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+TR*TR*nl*nl*(44./9.*CF+158./27.*CA))/(64*M_PI*M_PI*M_PI);
extern double K;  // = (67.0/18.0-M_PI*M_PI/6.0)*CA - 5.0/9.0*nf;

extern double Bq; // = -3.0/4.0;
extern double Bg; // = (-11*CA+4*nf*TR)/(12*CA);

const double MZ = 91.1876;

// update the main parameters
// negative values are not updated
// the other ones will automatically be updated
void set_fundamental_qcd_parameters(double _CF, double _CA,
                                    double _TR, int _nf);

// more detailed tweaks (-ve values below -0.5 are unchanged)
void set_qcd_running_parameters(double _b0=-1.0,
                                double _b1=-1.0,
                                double _K =-1.0);
void set_qcd_Bi(double _Bq=-1.0, double _Bg=-1.0);
                                
//----------------------------------------------------------------------
// running coupling
//----------------------------------------------------------------------
// this is the expression for alphas implementin the exact two-loop
// beta function if needed
//
// Note that this is a MSbar running, i.e. it does not include K terms!
double alphas_rg(double kt, double muR, double alphas_muR,
                 bool two_loops = true);


#endif  // __QCD_CONSTANTS_HH__
