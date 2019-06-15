#include "Definitions.h"
#include "Integrand.h"
#include "ResumFunc.h"
#include "HypGeo.h"

#include <cstring>
#include <algorithm>
#include <cfloat>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_psi.h>


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


double dQCDSD_cumu(void *params){
    paraSD *pa = (paraSD*)params;
    
    double alphas = pa->alphas;
    double rscaleR = pa->rscaleR;
    double rscaleQ = pa->rscaleQ;
    double O_alpha = pa->tau;
    if(O_alpha <= 0.) return(0.);
    
    double zc = pa->zcut;
    double beta = pa->beta;
    double mu2 = 1./rscaleR;
    double O_max = sqrt(rscaleQ);
    
    double alpha = alpha_Obs;
    
    if(fixedorder == 1) O_max = 1./4.;
    if(fixedorder == 2) O_max = 1./3.;
    
//     cout << "test res: " << tau << endl;
    
    if(O_alpha>O_max && end_point){
        paraSD temp_end = {alphas,rscaleR,rscaleQ,O_max,pa->zcut,beta};
//         paraSD temp_end_alt = {alphas,rscaleR,rscaleQ,taum,zc,beta,Q2,pa->R};
//         cout << "test res2: " << taum << " " << pa->zcut << " " << zc << " " << beta << " " << dQCDSD_cumu(&temp_end) << " " << dQCDSD_cumu(&temp_end_alt) << endl;
        return(dQCDSD_cumu(&temp_end));
    }
    
    O_alpha = O_alpha/sqrt(rscaleQ);
    O_max = O_max/sqrt(rscaleQ);
    zc = zc/sqrt(rscaleQ);
    
    double O_L, zcL;
    if(end_point) O_L = O_alpha*O_max/pow(pow(O_max,ppar)-pow(O_alpha,ppar)+pow(O_alpha*O_max,ppar),1./ppar);
    else O_L = O_alpha;
//     cout << "test res3: " << tauL << endl;
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
    if(f_gluon == 1){
        C1 = CA;
    }
    else{
        C1 = CF;
    }
    
//     cout << pJ2 << "\t" << pSC2 << "\t" << pSG2 << endl;
    
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double GammaJ, gammaJ01;
    double GammaSC, gammaSC01;
    double GammaSG, gammaSG01;
    double GammaH, gammaH01;
    double GammaS, gammaS01;
    double JLL = 1., JNLL = 1., JNNLL = 1., J0 = 0.;
    double SCLL = 1., SCNLL = 1., SCNNLL = 1., SC0 = 0.;
    double SGLL = 1., SGNLL = 1., SGNNLL = 1., SG0 = 0.;
    double HLL = 1., HNLL = 1., HNNLL = 1., H0 = 0.;
    double S0 = 0., S0p = 0., S0pp = 0.;
    
    double InvNLL = 1., InvNNLL = 1.;
    double InvpNLL = 0., InvpNNLL = 0.;
    
    double RpJLL1 = 0., RppJLL1 = 0; 
    double RpSCLL1 = 0., RppSCLL1 = 0.; 
    double RpSGLL1 = 0., RppSGLL1 = 0.; 
    double RpTLL1 = 0., RppTLL1 = 0.; 
    
    double RpJALL1 = 0., RpJANLL1 = 0., RppJALL1 = 0.;
    double RpSCALL1 = 0., RpSCANLL1 = 0., RppSCALL1 = 0.;   
    double RpSGALL1 = 0., RpSGANLL1 = 0., RppSGALL1 = 0.;   
    
    double TranspNLL = 0., TransNLL = 1.;
    double RppTLLzc1 = 0., RppSCLLzc1 = 0., RppSGLLzc1 = 0., RppJLLzc1 = 0.;
    
    GammaJ = alpha/(alpha-1.);
    GammaSC = -(beta+alpha)/(beta+1.)/(alpha-1.);
    GammaSG = 2./(1.+beta);
    GammaH = -2.;    
    
    GammaS = -2./(alpha-1.);    
    
    double transp = zc;
        
    if(LLorder>0){
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
        HLL = exp(RLL(mu2, 1., 0., 0., 1./O_L, zcL, GammaH*C1, alphas)/2.);
        
        RpJALL1 = RpLL(mu2, 1., pJnu, pJzc, 1./O_L, zcL, GammaJ*C1, alphas)/2.;
        
        RpTLL1 = 2.*RpJALL1+2.*RpSCALL1+RpSGALL1;
    }    
    if(LLorder>1){
        if(f_gluon == 1){
            gammaJ01 = 2.*beta0;
            gammaH01 = -4.*beta0;
        }
        else{
            gammaJ01 = 6.*CF;
            gammaH01 = -12.*CF;
        }
        gammaSC01 = 0.;
        gammaSG01 = 0.;
        
        gammaS01 = 0.;
        
        if(O_L>transp){
            SCNLL = 1.;
            SGNLL = exp(RNLL(mu2, 1., pSnu, pSzc, pS2, 1./O_L, zcL, GammaS*C1, gammaS01, alphas)/2.);
            
            if(end_point) SGNLL = SGNLL*exp(pow(O_alpha/O_max,ppar)*RpNLL(mu2, 1., pSnu, pSzc, pS2, 1., zc, GammaS*C1, gammaS01, alphas)*log(O_L)/2.);
            
            RppSCALL1 = 0.;
            RppSGALL1 = RppLL(mu2, 1., pSnu, pSzc, 1./O_L, zcL, GammaS*C1, alphas)/2.;
        }
        else{        
            SCNLL = exp(RNLL(mu2, 1., pSCnu, pSCzc, pSC2, 1./O_L, zcL, GammaSC*C1, gammaSC01, alphas)/2.);
            SGNLL = exp(RNLL(mu2, 1., pSGnu, pSGzc, pSG2, 1./O_L, zcL, GammaSG*C1, gammaSG01, alphas)/2.);
            
            if(end_point) SGNLL = SGNLL*exp(pow(O_alpha/O_max,ppar)*RpNLL(mu2, 1., pSnu, pSzc, pS2, 1., zc, GammaS*C1, gammaS01, alphas)*log(O_L)/2.);
            
            
            RppSCALL1 = RppLL(mu2, 1., pSCnu, pSCzc, 1./O_L, zcL, GammaSC*C1, alphas)/2.;
            RppSGALL1 = 0.;
        }
        
        JNLL = exp(RNLL(mu2, 1., pJnu, pJzc, pJ2, 1./O_L, zcL, GammaJ*C1, gammaJ01, alphas)/2.);
        if(end_point) JNLL = JNLL*exp(pow(O_alpha/O_max,ppar)*RpNLL(mu2, 1., pJnu, pJzc, pJ2, 1., zc, GammaJ*C1, gammaJ01, alphas)*log(O_L)/2.);
        
        
        RppJALL1 = RppLL(mu2, 1., pJnu, pJzc, 1./O_L, zcL, GammaJ*C1, alphas)/2.;
        
        RppSCLLzc1 = RppLL(mu2, 1., pSCnu, pSCzc, 1./transp, zcL, GammaSC*C1, alphas)/2.;
        
        RppJLLzc1 = RppLL(mu2, 1., pJnu, pJzc, 1./transp, zcL, GammaJ*C1, alphas)/2.;
        
        RppTLL1 = 2.*RppJALL1+2.*RppSCALL1+RppSGALL1;
        
        RppTLLzc1 = 2.*RppJLLzc1+2.*RppSCLLzc1+RppSGLLzc1;
        
        InvNLL = exp(-log(Gamma(1.-RpTLL1))+GammaEuler*RpTLL1);
        
        gsl_sf_result realres;
        gsl_sf_result imagres;
        double rez, imz, *reres, *imres;
        
//         double argLi2 = zc/tauL;
        double argLi2 = zc/1.;
        
        rez = argLi2;
        imz = 0.;
        
        int i = gsl_sf_complex_dilog_xy_e(rez,imz,&realres,&imagres);
        
        double reLi2 = realres.val;
        double imLi2 = imagres.val; 
        
//             complex<double> Li2 = (complex<double> (reLi2,imLi2));
        double Li20 = reLi2;
        
        
        if(O_alpha != O_max && F_HypGeo){
            if(O_L>transp) TransNLL = TransNLL*exp(HypGeo_Py( -RpTLL1, transp/O_L)*log(O_L/transp)*(RppTLLzc1-RppTLL1));
            double Fzc = 1.;
            if(end_point){
                double lambda = alphas/4./pi*beta0*log(zc);
                
                double RppTLL10 = 2.*RppLL(mu2, 1., pJnu, pJzc, 1., zcL, GammaJ*C1, alphas)/2.+RppLL(mu2, 1., pSnu, pSzc, 1., zcL, GammaS*C1, alphas)/2.;
                
                Fzc = (RppTLLzc1-RppTLL10);
                
                TransNLL = TransNLL*exp(-pow(O_alpha/O_max,ppar)*2.*alphas/pi*C1/alpha*Li20*log(1./zc)*Fzc*log(O_L));
            }
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
    
    return(finalresult);
}


double dQCD_cumu(void *params){
    paraSD *pa = (paraSD*)params;
    
    double alphas = pa->alphas;
    double rscaleR = pa->rscaleR;
    double rscaleQ = pa->rscaleQ;
    double O_alpha = pa->tau;
    if(O_alpha <= 0.) return(0.);
    double zc = 1.;
    double beta = pa->beta;
    
    double alpha = alpha_Obs;
    
    double mu2 = 1./rscaleR;
    
    double O_max = sqrt(rscaleQ);
    
    double C1;
    if(f_gluon == 1){
        C1 = CA;
    }
    else{
        C1 = CF;
    }
    
    
    
    if(fixedorder == 1) O_max = 1./4.;
    if(fixedorder == 2) O_max = 1./3.;
    
    if(O_alpha>O_max && end_point){
        paraSD temp_end = {alphas,rscaleR,rscaleQ,O_max,pa->zcut,beta};
        return(dQCD_cumu(&temp_end));
    }
    
    O_alpha = O_alpha/sqrt(rscaleQ);
    O_max = O_max/sqrt(rscaleQ);
    zc = zc/sqrt(rscaleQ);
    
    double O_L;
    if(end_point) O_L = O_alpha*O_max/pow(pow(O_max,ppar)-pow(O_alpha,ppar)+pow(O_alpha*O_max,ppar),1./ppar);
    else O_L = O_alpha;
    
    double pJnu = -1./alpha;
    double pJzc = 0.;
    double pJ2 = 0.;
    
    double pSnu = -1.;
    double pSzc = 0.;
    double pS2 = 0.;
    
    pJ2 += log(rscaleQ)/log(2.)/2.*(-pJnu);
    pS2 += log(rscaleQ)/log(2.)/2.*(-pSnu);
    
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double GammaJ, gammaJ01;
    double GammaH, gammaH01;
    double GammaS, gammaS01;
    double JLL = 1., JNLL = 1., JNNLL = 1., J0 = 0.;
    double SLL = 1., SNLL = 1., SNNLL = 1., S0 = 0.;
    double H0 = 0.;
    
    double InvNLL = 1., InvNNLL = 1.;
    double InvpNLL = 0., InvpNNLL = 0.;
    
    double RpJLL1 = 0., RppJLL1 = 0; 
    double RpSLL1 = 0., RppSLL1 = 0.; 
    double RpTLL1 = 0., RppTLL1 = 0.;
    
    double RpJLL2 = 0., RppJLL2 = 0; 
    double RpSLL2 = 0., RppSLL2 = 0.; 
    double RpTLL2 = 0., RppTLL2 = 0.;     
    
    double RpJALL1 = 0., RpJANLL1 = 0., RppJALL1 = 0., RpJANNLL1 = 0.;
    double RpSALL1 = 0., RpSANLL1 = 0., RppSALL1 = 0., RpSANNLL1 = 0.;
    double RpJALL2 = 0., RpJANLL2 = 0., RppJALL2 = 0., RpJANNLL2 = 0.;
    double RpSALL2 = 0., RpSANLL2 = 0., RppSALL2 = 0., RpSANNLL2 = 0.;   
    
    GammaJ = (alpha/(alpha-1.));
    GammaH = -2.;    
    
    GammaS = -2./(alpha-1.);    
    
    if(LLorder>0){
        SLL = exp(RLL(mu2, 1., pSnu, pSzc, 1./O_L, zc, GammaS*C1, alphas)/2.);
            
        RpSALL1 = RpLL(mu2, 1., pSnu, pSzc, 1./O_L, zc, GammaS*C1, alphas)/2.;
        
        
        JLL = exp(RLL(mu2, 1., pJnu, pJzc, 1./O_L, zc, GammaJ*C1, alphas)/2.);
        
        RpJALL1 = RpLL(mu2, 1., pJnu, pJzc, 1./O_L, zc, GammaJ*C1, alphas)/2.;
        
        RpTLL1 = 2.*RpJALL1+RpSALL1;
    }    
    if(LLorder>1){
        if(f_gluon == 1){
            gammaJ01 = 2.*beta0;
            gammaH01 = -4.*beta0;
        }
        else{
            gammaJ01 = 6.*CF;
            gammaH01 = -12.*CF;
        }
        
        gammaS01 = 0.;
            
            
        SNLL = exp(RNLL(mu2, 1., pSnu, pSzc, pS2, 1./O_L, zc, GammaS*C1, gammaS01, alphas)/2.);
        if(end_point) SNLL = SNLL*exp(pow(O_alpha/O_max,ppar)*RpNLL(mu2, 1., pSnu, pSzc, pS2, 1., zc, GammaS*C1, gammaS01, alphas)*log(O_L)/2.);
                   
            
        
        JNLL = exp(RNLL(mu2, 1., pJnu, pJzc, pJ2, 1./O_L, zc, GammaJ*C1, gammaJ01, alphas)/2.);
        if(end_point) JNLL = JNLL*exp(pow(O_alpha/O_max,ppar)*RpNLL(mu2, 1., pJnu, pJzc, pJ2, 1., zc, GammaJ*C1, gammaJ01, alphas)*log(O_L)/2.);
        
        
        InvNLL = exp(-log(Gamma(1.-RpTLL1))+GammaEuler*RpTLL1);
        
    }
    double J = JLL*JNLL*JNNLL;
    double S = SLL*SNLL*SNNLL;
    double H = 1.;
    double Inv = InvNLL*InvNNLL;
    double T0 = 1.+2.*J0+S0+H0;
    T0 = 1.;
    
    double total = T0*J*J*S*H*Inv;
    
    cout << T0 << " " << J << " " << S << " " << H << " " << Inv << endl;
    
    double finalresult = total;
    
    if(!(finalresult == finalresult && (finalresult <= DBL_MAX && finalresult>= -DBL_MAX))) finalresult = 0.;
    
    return(finalresult);
}



double ExpansionSD_cumu(void *params){
    paraSD *pa = (paraSD*)params;
    
    double alphas = pa->alphas;
    double rscale = pa->rscaleR;
    double rscaleQ = pa->rscaleQ;
    double O_alpha = pa->tau;
    if(O_alpha<=0.) return(0.);
    double zc = pa->zcut;
    double beta = pa->beta;
    
    double alpha = alpha_Obs;
    
    double mu2 = 1./rscale;
    double O_max = sqrt(rscaleQ);
    if(fixedorder == 1) O_max = 1./4.;
    if(fixedorder == 2) O_max = 1./3.;
    
    if(O_alpha>O_max && end_point){
        paraSD tempSD_end = {alphas,rscale,rscaleQ,O_max,pa->zcut,beta};
//         cout << "Above endpoint" << endl;
        return(ExpansionSD_cumu(&tempSD_end));
    }   
    
    
    double C1;
    if(f_gluon == 1){
        C1 = CA;
    }
    else{
        C1 = CF;
    }
    
    
    para pa2 = {alphas,rscale, rscaleQ,O_alpha};
    
    
    O_alpha = O_alpha/sqrt(rscaleQ);
    O_max = O_max/sqrt(rscaleQ);
    zc = zc/sqrt(rscaleQ);    
    
    double O_L, zcL;
    if(end_point) O_L = O_alpha*O_max/pow(pow(O_max,ppar)-pow(O_alpha,ppar)+pow(O_alpha*O_max,ppar),1./ppar);
    else O_L = O_alpha;
    
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    if(O_L>zc){
        double total = 0;
        total = Expansion_cumu(&pa2);
        
        if(fixedorder>1 && F_HypGeo && exp_coef != 1){    
            gsl_sf_result realres;
            gsl_sf_result imagres;
            double rez, imz, *reres, *imres;
    
            double transp = zc;
        
            double argLi2 = transp/O_L;
        
            rez = argLi2;
            imz = 0.;
        
            int i = gsl_sf_complex_dilog_xy_e(rez,imz,&realres,&imagres);
        
            double reLi2 = realres.val;
            double imLi2 = imagres.val; 
        
            double Li2 = reLi2;
            
            double argLi20 = transp/1.;
        
            rez = argLi20;
            imz = 0.;
        
            i = gsl_sf_complex_dilog_xy_e(rez,imz,&realres,&imagres);
        
            double reLi20 = realres.val;
            double imLi20 = imagres.val; 
        
            double Li20 = reLi20;
            
            double trans_above_as2 = 8.*alphas*alphas/pi/pi/alpha/(beta+alpha)*C1*C1*Li2*log(O_L/transp)*log(O_L);
            
            total += trans_above_as2;
            
            double endpoint_cor2 = 0.;
            if(end_point){
                endpoint_cor2 = -pow(O_alpha/O_max,ppar)*8.*alphas*alphas/pi/pi/alpha/(beta+alpha)*C1*C1*Li20*log(1./transp)*log(O_L);
                total += endpoint_cor2;
            }
        }
        return(total);
    }
    

    
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
    
    pS2 += log(rscaleQ)/log(2.)/2.*(-pSnu);
    
    double GammaS, gammaS0;
    
    
    pJ2 += log(rscaleQ)/log(2.)/2.*(-pJnu);
    pSG2 += log(rscaleQ)/log(2.)/2.*(pSGzc);
    pSC2 += log(rscaleQ)/log(2.)/2.*(pSCzc-pSCnu);
    
    double GammaJ, gammaJ01, gammaJ11;
    double GammaSC, gammaSC01, gammaSC11;
    double GammaSG, gammaSG01, gammaSG11;
    double GammaH, gammaH01, gammaH11;
    
    GammaJ = alpha/(alpha-1.);
    GammaSC = -(beta+alpha)/(beta+1.)/(alpha-1.);
    GammaSG = 2./(1.+beta);
    GammaH = -2.;    
    
    GammaS = -2./(alpha-1.);
    
    if(f_gluon == 0){
        gammaJ01 = 2.*beta0;
        gammaH01 = -4.*beta0;
    }
    else{
        gammaJ01 = 6.*CF;
        gammaH01 = -12.*CF;
    }
    
    gammaSC01 = 0.;
    gammaSG01 = 0.;
    
    
    gammaJ11 = 0.;
    gammaJ11 = 0.;
    gammaSC11 = 0.;
    gammaSG11 = 0.;
    
   
    double LTJ  = pJnu*log(1./O_L)+pJzc*log(zc);
    double LTSC = pSCnu*log(1./O_L)+pSCzc*log(zc);
    
    double expJ1 = exp1(mu2, 1., pJnu, pJzc, pJ2,  1./O_L, zcL, GammaJ*C1, gammaJ01)/2.;
    double expSC1 = exp1(mu2, 1., pSCnu, pSCzc, pSC2,  1./O_L, zcL, GammaSC*C1, gammaSC01)/2.;
    double expSG1 = exp1(mu2, 1., pSGnu, pSGzc, pSG2,  1./O_L, zcL, GammaSG*C1, gammaSG01)/2.;
    
    double expJ2 = exp2(mu2, 1., pJnu, pJzc, pJ2,  1./O_L, zcL, GammaJ*C1, gammaJ01, gammaJ11)/2.;
    double expSC2 = exp2(mu2, 1., pSCnu, pSCzc, pSC2,  1./O_L, zcL, GammaSC*C1, gammaSC01, gammaSC11)/2.;
    double expSG2 = exp2(mu2, 1., pSGnu, pSGzc, pSG2,  1./O_L, zcL, GammaSG*C1, gammaSG01, gammaSG11)/2.;
    
    double expinv21 = -1./3.*pi*pi*(2.*GammaJ*C1*pJnu*LTJ+2.*GammaSC*C1*pSCnu*LTSC)*(2.*GammaJ*C1*pJnu*LTJ+2.*GammaSC*C1*pSCnu*LTSC)*16./2.;
    
    double expinv2 = expinv21;
    
    double total = 0.;
    
    double endpoint_corp = 0.;
    double endpoint_cor = 0.;
    double endpoint_cor2 = 0.;
    
    double J0 = 0.;
    double SC0 = 0.;
    double SG0 = 0.;
    double H0 = 0.;
    
    if(fixedorder>0){
      
      if(end_point){
        endpoint_cor += -pow(O_alpha/O_max,ppar)*(gammaS0*C1-8.*pS2*GammaS*C1*log(2.))*pSnu*log(O_L)/2.;
        endpoint_cor += -2.*pow(O_alpha/O_max,ppar)*(gammaJ01-8.*pJ2*GammaJ*C1*log(2.))*pJnu*log(O_L)/2.;
      }
      if(exp_coef == 0){
        total = 1.+alphas/4./pi*(2.*expJ1+2.*expSC1+expSG1+2.*J0+2.*SC0+SG0+H0);
        total += alphas/4./pi*endpoint_cor;
      }
      else if(exp_coef == 1){
          total = alphas/4./pi*(2.*expJ1+2.*expSC1+expSG1+2.*J0+2.*SC0+SG0+H0);
          total += alphas/4./pi*endpoint_cor;
      }
    }   
    double Li20 = 0.;
    double transp = zc;
    if(fixedorder>1){        
        gsl_sf_result realres;
        gsl_sf_result imagres;
        double rez, imz, *reres, *imres;
        
        double argLi20 = transp/1.;
        
        rez = argLi20;
        imz = 0.;
        
        int i = gsl_sf_complex_dilog_xy_e(rez,imz,&realres,&imagres);
        
        double reLi20 = realres.val;
        double imLi20 = imagres.val; 
        
        Li20 = reLi20;
            
        if(end_point && F_HypGeo) endpoint_cor2 += -pow(O_alpha/O_max,ppar)*16.*8./alpha/(beta+alpha)*C1*C1*Li20*log(1./transp)*log(O_L);
        
        if(exp_coef == 0){
            total += pow(alphas/4./pi,2.)*(2.*expJ2+2.*expSC2+expSG2+expinv2+((2.*expJ1+2.*expSC1+expSG1+endpoint_cor)/2.+2.*J0+2.*SC0+SG0+H0) *(2.*expJ1+2.*expSC1+expSG1+endpoint_cor));
            total += pow(alphas/4./pi,2.)*endpoint_cor2;
        }
        else if(exp_coef == 2){
            total = pow(alphas/4./pi,2.)*(2.*expJ2+2.*expSC2+expSG2+expinv2+((2.*expJ1+2.*expSC1+expSG1+endpoint_cor)/2.+2.*J0+2.*SC0+SG0+H0) *(2.*expJ1+2.*expSC1+expSG1+endpoint_cor));
            total += pow(alphas/4./pi,2.)*endpoint_cor2;
        }
    }
    
    double finalresult = total;    
    
    return(finalresult);
}

double Expansion_cumu(void *params){
    para *pa = (para*)params;
    
    double alphas = pa->alphas;
    double rscale = pa->rscaleR;
    double rscaleQ = pa->rscaleQ;
    double O_alpha = pa->tau;
    if(O_alpha<=0.) return(0.);
    double mu2 = 1./rscale;
    double O_max = sqrt(rscaleQ);
    
    double alpha = alpha_Obs;
    
    if(fixedorder == 1) O_max = 1./4.;
    if(fixedorder == 2) O_max = 1./3.;
    
    if(O_alpha>O_max && end_point){
        para temp_end = {alphas,rscale,rscaleQ,O_max};
        return(Expansion_cumu(&temp_end));
    }   
    
    
    double C1;
    if(f_gluon == 1){
        C1 = CA;
    }
    else{
        C1 = CF;
    }
    
    
    O_alpha = O_alpha/sqrt(rscaleQ);
    O_max = O_max/sqrt(rscaleQ);
    
    double O_L;
    if(end_point) O_L = O_alpha*O_max/pow(pow(O_max,ppar)-pow(O_alpha,ppar)+pow(O_alpha*O_max,ppar),1./ppar);
    else O_L = O_alpha;
    
    double pJnu = -1./alpha;
    double pJzc = 0.;
    double pJ2 = 0.;
    
    double pSnu = -1.;
    double pSzc = 0.;
    double pS2 = 0.;
    
    pJ2 += log(rscaleQ)/log(2.)/2.*(-pJnu);
    pS2 += log(rscaleQ)/log(2.)/2.*(-pSnu);
    
    double GammaJ, gammaJ01, gammaJ11;
    double GammaS, gammaS01, gammaS11;
    double GammaH, gammaH01, gammaH11;
    
    GammaJ = (alpha/(alpha-1.));
    GammaH = -2.;
    GammaS = -2./(alpha-1.);
    
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    if(f_gluon == 0){
        gammaJ01 = 2.*beta0;
        gammaH01 = -4.*beta0;
    }
    else{
        gammaJ01 = 6.*CF;
        gammaH01 = -12.*CF;
    }
    gammaS01 = 0.;
    
    gammaJ11 = 0.;
    gammaH11 = 0.;
    gammaS11 = 0.;
    
    
    double LTJ = pJnu*log(1./O_L)+pJzc*log(1.);//+pJ2*log(2.);
    double LTS = pSnu*log(1./O_L)+pSzc*log(1.);//+pS2*log(2.);
    
    double expJ1 = exp1(mu2, 1., pJnu, pJzc, pJ2,  1./O_L, 1., GammaJ*C1, gammaJ01)/2.;
    double expS1 = exp1(mu2, 1., pSnu, pSzc, pS2,  1./O_L, 1., GammaS*C1, gammaS01)/2.;
    double expH1 = 0.;
                             
    double expJ2 = exp2(mu2, 1., pJnu, pJzc, pJ2,  1./O_L, 1., GammaJ*C1, gammaJ01, gammaJ11)/2.;
    double expS2 = exp2(mu2, 1., pSnu, pSzc, pS2,  1./O_L, 1., GammaS*C1, gammaS01, gammaS11)/2.;
    
    double expinv21 = -1./3.*pi*pi*(2.*GammaJ*C1*pJnu*LTJ+GammaS*C1*pSnu*LTS)*(2.*GammaJ*C1*pJnu*LTJ+GammaS*C1*pSnu*LTS)*16./2.;
    
    double expinv2 = expinv21;
    
    double total;
    
    double endpoint_cor = 0., endpoint_corp = 0.;
    
    
    double J0 = 0.;
    double S0 = 0.;
    double H0 = 0.;
    
    if(fixedorder>0){
      if(end_point){
        endpoint_cor += -pow(O_alpha/O_max,ppar)*(gammaS01-8.*pS2*GammaS*C1*log(2.))*pSnu*log(O_L)/2.;
        endpoint_cor += -2.*pow(O_alpha/O_max,ppar)*(gammaJ01-8.*pJ2*GammaJ*C1*log(2.))*pJnu*log(O_L)/2.;
        
      }
      if(exp_coef == 0){
          total = 1.+alphas/4./pi*(2.*expJ1+expS1+2.*J0+S0+H0);
          total += alphas/4./pi*endpoint_cor;
      }
      else if(exp_coef == 1){
          total = alphas/4./pi*(2.*expJ1+expS1+2.*J0+S0+H0);
          total += alphas/4./pi*endpoint_cor;
      }
    }  
    if(fixedorder>1){ 
        if(exp_coef == 0){
            total += pow(alphas/4./pi,2.)*(2.*expJ2+expS2+expinv2+((2.*expJ1+expS1+expH1+endpoint_cor)/2.+2.*J0+S0+H0)*(2.*expJ1+expS1+expH1+endpoint_cor));
        }
        else if(exp_coef == 2){
            total = pow(alphas/4./pi,2.)*(2.*expJ2+expS2+expinv2+((2.*expJ1+expS1+expH1+endpoint_cor)/2.+2.*J0+S0+H0)*(2.*expJ1+expS1+expH1+endpoint_cor));
        }       
    }
    
    double finalresult = total;
    
    return(finalresult);
}
