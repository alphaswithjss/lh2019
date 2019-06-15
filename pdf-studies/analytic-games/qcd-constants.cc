#include "qcd-constants.hh"

//----------------------------------------------------------------------
// physics constants
//----------------------------------------------------------------------
double CF = 1.33333333333333333333333;
double CA = 3.0;
double TR = 0.5;
unsigned int nf = 5;

double b0 = (11*CA-4*nf*TR)/(12*M_PI);
double b1 = (17*CA*CA-5*CA*nf-3*CF*nf)/(24*M_PI*M_PI);
double K  = (67.0/18.0-M_PI*M_PI/6.0)*CA - 5.0/9.0*nf;

double Bq = -3.0/4.0;
double Bg = (-11*CA+4*nf*TR)/(12*CA);

void set_fundamental_qcd_parameters(double _CF, double _CA, double _TR, int _nf){
  if (_CF>-0.5) CF=_CF;
  if (_CA>-0.5) CA=_CA;
  if (_TR>-0.5) TR=_TR;
  if (_nf>=0) nf=(unsigned int) _nf;

  b0 = (11*CA-4*nf*TR)/(12*M_PI);
  b1 = (17*CA*CA-5*CA*nf-3*CF*nf)/(24*M_PI*M_PI);
  K  = (67.0/18.0-M_PI*M_PI/6.0)*CA - 5.0/9.0*nf;
  
  Bq = -3.0/4.0;
  Bg = (-11*CA+4*nf*TR)/(12*CA);
}

void set_qcd_running_parameters(double _b0,
                                double _b1,
                                double _K){
  if (_b0>-0.5){ b0 = _b0;}
  if (_b1>-0.5){ b1 = _b1;}
  if (_K >-0.5){ K  = _K;}
}

void set_qcd_Bi(double _Bq, double _Bg){
  if (_Bq>-0.5){ Bq = _Bq;}
  if (_Bg>-0.5){ Bg = _Bg;}
}

//----------------------------------------------------------------------
// running copling
//----------------------------------------------------------------------


// this is the expression for alphas implementin the exact two-loop
// beta function if needed
//
// Note that this is a MSbar running, i.e. it does not include K terms!
double alphas_rg(double kt, double muR, double alphas_muR, bool two_loops){
  if (!two_loops) return alphas_muR/(1+2*b0*alphas_muR*log(kt/muR));

  double as = alphas_muR;
  double asprev;
  double L = log(kt/muR);
  do{
    asprev = as;
    as = 1.0/(2*b0*L + 1/alphas_muR
              -b1/b0*log((as/alphas_muR)
                         *((1+b1/b0*alphas_muR)/(1+b1/b0*as))));
  } while (std::abs(as-asprev)>1e-6);

  return as;
}
