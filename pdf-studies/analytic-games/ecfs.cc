#include "qcd-constants.hh"
#include "ecfs.hh"
#include "resum-blocks-vincent.hh"
#include "hypergeometric.hh"

#include <cstring>
#include <algorithm>
#include <cfloat>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

double Gamma(double z){
  double g = 7.0;
  double p[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                 771.32342877765313, -176.61502916214059, 12.507343278686905,
                 -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
  if(z < 0.5) return(pi/(sin(pi*z)*Gamma(1.0-z)));
  else{
    z = z-1.0;
    double x = p[0];
    for(double i=1.0; i<g+2.0;i++){
      x = x+p[(int)i]/(z+i);
    }
    double t = z+g+0.5;
    return(sqrt(2.0*pi)*pow(t,z+0.5)*exp(-t)*x);
  }
}


double ECF::fraction_below(double lambda, double pt, bool gluon){
  if (_zcut<=0) return fraction_below_plain(lambda, pt, gluon);
                    
  double alphas = alphas_rg(_muR*pt*_R, MZ, _alphas_MZ, _two_loops);
  int old_LLorder = LLorder;
    
  double rscaleR = 1/(_muR*_muR);
  double rscaleQ = 1/(_muQ*_muQ);
  double O_alpha = lambda;
  if(O_alpha <= 0.) return 0.0;
    
  double zc = _zcut;
  double beta = _beta;
  double mu2 = 1./rscaleR;
  double O_max = (_endpoint_order==0)
    ? sqrt(rscaleQ)
    : ((_endpoint_order==0) ? 0.25 : 1.0/3.0);
    
  double alpha = _alpha;
    
  if(O_alpha>O_max)
    return fraction_below(O_max, pt, gluon);
    
  LLorder = _two_loops ? 2 : 1;

  O_alpha = O_alpha/sqrt(rscaleQ);
  O_max = O_max/sqrt(rscaleQ);
  zc = zc/sqrt(rscaleQ);
    
  double O_L, zcL;
  O_L = O_alpha*O_max/pow(pow(O_max,_p)-pow(O_alpha,_p)+pow(O_alpha*O_max,_p),1./_p);
  zcL = zc;
    
  double pJnu = -1./alpha;
  double pSCnu = -(beta+1.)/(beta+alpha);
  double pSGnu = 0.;
  double pJzc = 0.;
  double pSCzc = (alpha-1.)/(beta+alpha);
  double pSGzc = 1.;
  double pJ2 = 0.;
  double pSC2 = 0.;
  double pSG2 = 0.;
    
  double pSnu = -1.;
  double pSzc = 0.;
  double pS2 = 0.;
    
  pJ2 += log(rscaleQ)/log(2.)/2.*(-pJnu);
  pSC2 += log(rscaleQ)/log(2.)/2.*(pSCzc-pSCnu);
  pSG2 += log(rscaleQ)/log(2.)/2.*(pSGzc);
  pS2 += log(rscaleQ)/log(2.)/2.*(-pSnu);
    
  double C1;
  if(gluon){
    C1 = CA;
  } else {
    C1 = CF;
  }
    
  double beta0 = 4*pi*b0;
    
  double GammaJ, gammaJ01;
  double GammaSC, gammaSC01;
  double GammaSG, gammaSG01;
  //double GammaH, gammaH01;
  double GammaS, gammaS01;
  double JLL = 1., JNLL = 1., JNNLL = 1., J0 = 0.;
  double SCLL = 1., SCNLL = 1., SCNNLL = 1., SC0 = 0.;
  double SGLL = 1., SGNLL = 1., SGNNLL = 1., SG0 = 0.;
  //double HLL = 1., HNLL = 1., HNNLL = 1., H0 = 0.;
  //double HLL = 1., H0 = 0.;
  double H0 = 0.;
  //double S0 = 0., S0p = 0., S0pp = 0.;
    
  double InvNLL = 1., InvNNLL = 1.;
  //double InvpNLL = 0., InvpNNLL = 0.;
    
  //double RpJLL1 = 0., RppJLL1 = 0; 
  //double RpSCLL1 = 0., RppSCLL1 = 0.; 
  //double RpSGLL1 = 0., RppSGLL1 = 0.; 
  double RpTLL1 = 0., RppTLL1 = 0.; 
    
  double RpJALL1  = 0., RppJALL1 = 0.;
  double RpSCALL1 = 0., RppSCALL1 = 0.;   
  double RpSGALL1 = 0., RppSGALL1 = 0.;   
  // double RpJALL1 = 0., RpJANLL1 = 0., RppJALL1 = 0.;
  // double RpSCALL1 = 0., RpSCANLL1 = 0., RppSCALL1 = 0.;   
  // double RpSGALL1 = 0., RpSGANLL1 = 0., RppSGALL1 = 0.;   
    
  //double TranspNLL = 0., TransNLL = 1.;
  double TransNLL = 1.;
  double RppTLLzc1 = 0., RppSCLLzc1 = 0., RppSGLLzc1 = 0., RppJLLzc1 = 0.;
    
  GammaJ = alpha/(alpha-1.);
  GammaSC = -(beta+alpha)/(beta+1.)/(alpha-1.);
  GammaSG = 2./(1.+beta);
  //GammaH = -2.;    
    
  GammaS = -2./(alpha-1.);    
    
  double transp = zc;
        
  if(O_L>transp){
    SCLL = 1.;
    SGLL = exp(RLL(mu2, 1., pSnu, pSzc, 1./O_L, zcL, GammaS*C1, alphas)/2.);
            
    RpSCALL1 = 0.;
    RpSGALL1 = RpLL(mu2, 1., pSnu, pSzc, 1./O_L, zcL, GammaS*C1, alphas)/2.;
  }
  else{
    SCLL = exp(RLL(mu2, 1., pSCnu, pSCzc, 1./O_L, zcL, GammaSC*C1, alphas)/2.);
    SGLL = exp(RLL(mu2, 1., pSGnu, pSGzc, 1./O_L, zcL, GammaSG*C1, alphas)/2.);
      
    RpSCALL1 = RpLL(mu2, 1., pSCnu, pSCzc, 1./O_L, zcL, GammaSC*C1, alphas)/2.;
    RpSGALL1 = 0.;
  }
    
  JLL = exp(RLL(mu2, 1., pJnu, pJzc, 1./O_L, O_L, GammaJ*C1, alphas)/2.);
  //HLL = exp(RLL(mu2, 1., 0., 0., 1./O_L, zcL, GammaH*C1, alphas)/2.);
    
  RpJALL1 = RpLL(mu2, 1., pJnu, pJzc, 1./O_L, zcL, GammaJ*C1, alphas)/2.;
    
  RpTLL1 = 2.*RpJALL1+2.*RpSCALL1+RpSGALL1;
     

  if(LLorder>1){
    if(gluon){
      gammaJ01 = 2.*beta0;
      //gammaH01 = -4.*beta0;
    } else {
      gammaJ01 = 6.*CF;
      //gammaH01 = -12.*CF;
    }
    gammaSC01 = 0.;
    gammaSG01 = 0.;
        
    gammaS01 = 0.;
        
    if(O_L>transp){
      SCNLL = 1.;
      SGNLL = exp(RNLL(mu2, 1., pSnu, pSzc, pS2, 1./O_L, zcL, GammaS*C1, gammaS01, alphas)/2.);
            
      SGNLL *= exp(pow(O_alpha/O_max,_p)*RpNLL(mu2, 1., pSnu, pSzc, pS2, 1., zc, GammaS*C1, gammaS01, alphas)*log(O_L)/2.);
            
      RppSCALL1 = 0.;
      RppSGALL1 = RppLL(mu2, 1., pSnu, pSzc, 1./O_L, zcL, GammaS*C1, alphas)/2.;
    } else {        
      SCNLL = exp(RNLL(mu2, 1., pSCnu, pSCzc, pSC2, 1./O_L, zcL, GammaSC*C1, gammaSC01, alphas)/2.);
      SGNLL = exp(RNLL(mu2, 1., pSGnu, pSGzc, pSG2, 1./O_L, zcL, GammaSG*C1, gammaSG01, alphas)/2.);
            
      SGNLL *= exp(pow(O_alpha/O_max,_p)*RpNLL(mu2, 1., pSnu, pSzc, pS2, 1., zc, GammaS*C1, gammaS01, alphas)*log(O_L)/2.);
            
      RppSCALL1 = RppLL(mu2, 1., pSCnu, pSCzc, 1./O_L, zcL, GammaSC*C1, alphas)/2.;
      RppSGALL1 = 0.;
    }
        
    JNLL = exp(RNLL(mu2, 1., pJnu, pJzc, pJ2, 1./O_L, zcL, GammaJ*C1, gammaJ01, alphas)/2.);
    JNLL *= exp(pow(O_alpha/O_max,_p)*RpNLL(mu2, 1., pJnu, pJzc, pJ2, 1., zc, GammaJ*C1, gammaJ01, alphas)*log(O_L)/2.);
        
    RppJALL1 = RppLL(mu2, 1., pJnu, pJzc, 1./O_L, zcL, GammaJ*C1, alphas)/2.;
    RppSCLLzc1 = RppLL(mu2, 1., pSCnu, pSCzc, 1./transp, zcL, GammaSC*C1, alphas)/2.;
    RppJLLzc1 = RppLL(mu2, 1., pJnu, pJzc, 1./transp, zcL, GammaJ*C1, alphas)/2.;
    RppTLL1 = 2.*RppJALL1+2.*RppSCALL1+RppSGALL1;
    RppTLLzc1 = 2.*RppJLLzc1+2.*RppSCLLzc1+RppSGLLzc1;
    InvNLL = exp(-log(Gamma(1.-RpTLL1))+GammaEuler*RpTLL1);
        
    gsl_sf_result realres;
    gsl_sf_result imagres;
    double rez, imz; //, *reres, *imres;
        
    double argLi2 = zc/1.;
        
    rez = argLi2;
    imz = 0.;
        
    //int i = gsl_sf_complex_dilog_xy_e(rez,imz,&realres,&imagres);
    gsl_sf_complex_dilog_xy_e(rez,imz,&realres,&imagres);
        
    double reLi2 = realres.val;
    //double imLi2 = imagres.val; 
    double Li20 = reLi2;
        
        
    if ((O_alpha != O_max) && (_transition_point_corrections)){
      if(O_L>transp) TransNLL *= exp(HypGeo_Py( -RpTLL1, transp/O_L)*log(O_L/transp)*(RppTLLzc1-RppTLL1));
      double Fzc = 1.;
      //double lambda = alphas/4./pi*beta0*log(zc);
      
      double RppTLL10 = 2.*RppLL(mu2, 1., pJnu, pJzc, 1., zcL, GammaJ*C1, alphas)/2.+RppLL(mu2, 1., pSnu, pSzc, 1., zcL, GammaS*C1, alphas)/2.;
      
      Fzc = (RppTLLzc1-RppTLL10);
      
      TransNLL *= exp(-pow(O_alpha/O_max,_p)*2.*alphas/pi*C1/alpha*Li20*log(1./zc)*Fzc*log(O_L));
    }
    
  }
      
  double J = JLL*JNLL*JNNLL;
  double SC = SCLL*SCNLL*SCNNLL;
  double SG = SGLL*SGNLL*SGNNLL;
  double Inv = InvNLL*InvNNLL;
  double T0 = 1.+2.*J0+2.*SC0+SG0+H0;
  T0 = 1.;
  
  double total = T0*J*J*SC*SC*SG*Inv*TransNLL;
  
  double finalresult = total;
  
  if(!(finalresult == finalresult && (finalresult <= DBL_MAX && finalresult>= -DBL_MAX))) finalresult = 0.;
  
  LLorder = old_LLorder;

  return(finalresult);
}

double ECF::fraction_below_plain(double lambda, double pt, bool gluon){
  double alphas = alphas_rg(_muR*pt*_R, MZ, _alphas_MZ, _two_loops);
  int old_LLorder = LLorder;
    
  double rscaleR = _muR;
  double rscaleQ = _muQ;
  double O_alpha = lambda;
  if(O_alpha <= 0.) return 0.0;
    
  double zc = 1.0;
  double mu2 = 1./rscaleR;
  double O_max = (_endpoint_order==0)
    ? sqrt(rscaleQ)
    : ((_endpoint_order==0) ? 0.25 : 1.0/3.0);
    
  double alpha = _alpha;
    
  if (O_alpha>O_max)
    return fraction_below_plain(O_max, pt, gluon);
    
  LLorder = _two_loops ? 2 : 1;

  double C1;
  if(gluon){
    C1 = CA;
  }
  else{
    C1 = CF;
  }
    
  O_alpha = O_alpha/sqrt(rscaleQ);
  O_max = O_max/sqrt(rscaleQ);
  zc = zc/sqrt(rscaleQ);
    
  double O_L;
  O_L = O_alpha*O_max/pow(pow(O_max,_p)-pow(O_alpha,_p)+pow(O_alpha*O_max,_p),1./_p);
    
  double pJnu = -1./alpha;
  double pJzc = 0.;
  double pJ2 = 0.;
    
  double pSnu = -1.;
  double pSzc = 0.;
  double pS2 = 0.;
    
  pJ2 += log(rscaleQ)/log(2.)/2.*(-pJnu);
  pS2 += log(rscaleQ)/log(2.)/2.*(-pSnu);
    
  double beta0 = 4*pi*b0;
    
  double GammaJ, gammaJ01;
  //double GammaH; //, gammaH01;
  double GammaS, gammaS01;
  double JLL = 1., JNLL = 1., JNNLL = 1., J0 = 0.;
  double SLL = 1., SNLL = 1., SNNLL = 1., S0 = 0.;
  double H0 = 0.;
    
  double InvNLL = 1., InvNNLL = 1.;

  double RpTLL1 = 0.; 
    
  double RpJALL1 = 0.;
  double RpSALL1 = 0.;
    
  GammaJ = (alpha/(alpha-1.));
  //GammaH = -2.;    
    
  GammaS = -2./(alpha-1.);    
    
  SLL = exp(RLL(mu2, 1., pSnu, pSzc, 1./O_L, zc, GammaS*C1, alphas)/2.);
   
  RpSALL1 = RpLL(mu2, 1., pSnu, pSzc, 1./O_L, zc, GammaS*C1, alphas)/2.;
    
  JLL = exp(RLL(mu2, 1., pJnu, pJzc, 1./O_L, zc, GammaJ*C1, alphas)/2.);
        
  RpJALL1 = RpLL(mu2, 1., pJnu, pJzc, 1./O_L, zc, GammaJ*C1, alphas)/2.;
    
  RpTLL1 = 2.*RpJALL1+RpSALL1;

  if(LLorder>1){
    if (gluon){
      gammaJ01 = 2.*beta0;
      //gammaH01 = -4.*beta0;
    }
    else{
      gammaJ01 = 6.*CF;
      //gammaH01 = -12.*CF;
    }
        
    gammaS01 = 0.;
            
            
    SNLL = exp(RNLL(mu2, 1., pSnu, pSzc, pS2, 1./O_L, zc, GammaS*C1, gammaS01, alphas)/2.);
    SNLL *= exp(pow(O_alpha/O_max,_p)*RpNLL(mu2, 1., pSnu, pSzc, pS2, 1., zc, GammaS*C1, gammaS01, alphas)*log(O_L)/2.);
                   
            
        
    JNLL = exp(RNLL(mu2, 1., pJnu, pJzc, pJ2, 1./O_L, zc, GammaJ*C1, gammaJ01, alphas)/2.);
    JNLL *= exp(pow(O_alpha/O_max,_p)*RpNLL(mu2, 1., pJnu, pJzc, pJ2, 1., zc, GammaJ*C1, gammaJ01, alphas)*log(O_L)/2.);
        
        
    InvNLL = exp(-log(Gamma(1.-RpTLL1))+GammaEuler*RpTLL1);
        
  }
  double J = JLL*JNLL*JNNLL;
  double S = SLL*SNLL*SNNLL;
  double H = 1.;
  double Inv = InvNLL*InvNNLL;
  double T0 = 1.+2.*J0+S0+H0;
  T0 = 1.;
    
  double total = T0*J*J*S*H*Inv;
    
  double finalresult = total;
    
  if(!(finalresult == finalresult && (finalresult <= DBL_MAX && finalresult>= -DBL_MAX))) finalresult = 0.;
    
  LLorder = old_LLorder;
    
  return(finalresult);
}

