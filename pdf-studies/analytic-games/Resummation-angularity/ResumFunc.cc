#include "LHAPDF/LHAPDF.h"
#include "Definitions.h"
#include "ResumFunc.h"

#include <cstring>
#include <algorithm>
#include <cfloat>

double KLL(double mu2, double mu02, double Gamma, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    double beta2 = 2857./54.*CA*CA*CA+Tf*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+Tf*Tf*nl*nl*(44./9.*CF+158./27.*CA);
    
    double Gamma0 = Gamma*A0;
    
    double l = 1.+beta0/4./pi*alphas*log(mu02/mu2);
    double alphasmu0;
    
    if(LLorder == 1.) alphasmu0 = alphas/l;
    else if(LLorder == 2.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l));
    else if(LLorder == 3.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l)+pow(alphas/4./pi,2.)*(beta1*beta1/beta0/beta0*(log(l)*log(l)-log(l)-1.)-beta2/beta0*(l-1.)));
    
    double r = alphas/alphasmu0;
    
    double result = Ci*Gamma0/2./beta0/beta0*4.*pi/alphasmu0*(log(r)+1./r-1.);
    
    return(result);
}

double KNLL(double mu2, double mu02, double Gamma, double gamma0, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    double beta2 = 2857./54.*CA*CA*CA+Tf*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+Tf*Tf*nl*nl*(44./9.*CF+158./27.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    
    double l = 1.+beta0/4./pi*alphas*log(mu02/mu2);
    double alphasmu0;
    
    if(LLorder == 1.) alphasmu0 = alphas/l;
    else if(LLorder == 2.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l));
    else if(LLorder == 3.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l)+pow(alphas/4./pi,2.)*(beta1*beta1/beta0/beta0*(log(l)*log(l)-log(l)-1.)-beta2/beta0*(l-1.)));
    
    double r = alphas/alphasmu0;
    
    double result = Ci*Gamma0/2./beta0/beta0*((Gamma1/Gamma0-beta1/beta0)*(r-1.-log(r))-beta1/2./beta0*log(r)*log(r));
    
    result += -gamma0/2./beta0*log(r);
    
    return(result);
}

double KNNLL(double mu2, double mu02, double Gamma, double gamma0, double gamma1, double alphas){
    
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double A2 = 4.*CA*CA*(245./6.-134.*pi*pi/27.+11.*pow(pi,4.)/45.+22.*Zeta3/3.)+32.*CA*Tf*nl*(-209./108.+5.*pi*pi/27.-7.*Zeta3/3.)+4.*CF*Tf*nl*(16*Zeta3-55./3.)-64./27.*Tf*Tf*nl*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    double beta2 = 2857./54.*CA*CA*CA+Tf*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+Tf*Tf*nl*nl*(44./9.*CF+158./27.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    double Gamma2 = Gamma*A2;
    
    double l = 1.+beta0/4./pi*alphas*log(mu02/mu2);
    double alphasmu0;
    
    if(LLorder == 1.) alphasmu0 = alphas/l;
    else if(LLorder == 2.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l));
    else if(LLorder == 3.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l)+pow(alphas/4./pi,2.)*(beta1*beta1/beta0/beta0*(log(l)*log(l)-log(l)-1.)-beta2/beta0*(l-1.)));
    
    double r = alphas/alphasmu0;
    
    double result = Ci*Gamma0/2./beta0/beta0*alphasmu0/4./pi*((Gamma1/Gamma0-beta1/beta0)*beta1/beta0*(r-1.-r*log(r))-(beta1*beta1/beta0/beta0-beta2/beta0)*log(r)+(Gamma2/Gamma0-Gamma1*beta1/Gamma0/beta0+beta1*beta1/beta0/beta0-beta2/beta0)*(r*r-1.)/2.+(Gamma2/Gamma0-Gamma1*beta1/Gamma0/beta0)*(1.-r));
    
    result += -1./2./beta0*alphasmu0/4./pi*(gamma1-beta1*gamma0/beta0)*(r-1.);
    
    return(result);
}

double wLL(double mu2, double mu02, double Gamma, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    double beta2 = 2857./54.*CA*CA*CA+Tf*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+Tf*Tf*nl*nl*(44./9.*CF+158./27.*CA);
    
    double Gamma0 = Gamma*A0;
    
    double l = 1.+beta0/4./pi*alphas*log(mu02/mu2);
    double alphasmu0;
    
    if(LLorder == 1.) alphasmu0 = alphas/l;
    else if(LLorder == 2.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l));
    else if(LLorder == 3.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l)+pow(alphas/4./pi,2.)*(beta1*beta1/beta0/beta0*(log(l)*log(l)-log(l)-1.)-beta2/beta0*(l-1.)));
    
    double r = alphas/alphasmu0;
    
    double result = -Ci*Gamma0/2./beta0*(log(r));
    
    return(result);
}

double wNLL(double mu2, double mu02, double Gamma, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    double beta2 = 2857./54.*CA*CA*CA+Tf*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+Tf*Tf*nl*nl*(44./9.*CF+158./27.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    
    double l = 1.+beta0/4./pi*alphas*log(mu02/mu2);
    double alphasmu0;
    
    if(LLorder == 1.) alphasmu0 = alphas/l;
    else if(LLorder == 2.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l));
    else if(LLorder == 3.) alphasmu0 = alphas/l*(1.-alphas/4./pi/l*beta1/beta0*log(l)+pow(alphas/4./pi,2.)*(beta1*beta1/beta0/beta0*(log(l)*log(l)-log(l)-1.)-beta2/beta0*(l-1.)));
    
    double r = alphas/alphasmu0;
    
    double result = -Ci*Gamma0/2./beta0*((Gamma1/Gamma0-beta1/beta0)*(r-1.));
    
    return(result);
}

double wNNLL(double mu2, double mu02, double Gamma, double alphas){
    return(0.);
}

double RLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double Gamma0 = Gamma*A0;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = Ci*Gamma0/2./beta0/beta0*4.*pi/alphas*(lambdaT2p1*log(lambdaT2p1)-2.*lambdaT);
    
    if(lambdaT2p1<0.) return(0./0.);
    return(result);
}

double RNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = Ci/beta0/beta0*(Gamma1*(lambdaT-log(lambdaT2p1)/2.)
    +Gamma0*beta1/4./beta0*(log(lambdaT2p1)*(2.+log(lambdaT2p1))-4.*lambdaT)
    +Gamma0*beta0/2.*log(lambdaT2p1)*(log(Q2/mu2)+2.*p2*log(2.)));
    
    result += -gamma0/2./beta0*log(lambdaT2p1);
    
    return(result);
}

double RNNLL(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0, double gamma1, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double A2 = 4.*CA*CA*(245./6.-134.*pi*pi/27.+11.*pow(pi,4.)/45.+22.*Zeta3/3.)+32.*CA*Tf*nl*(-209./108.+5.*pi*pi/27.-7.*Zeta3/3.)+4.*CF*Tf*nl*(16*Zeta3-55./3.)-64./27.*Tf*Tf*nl*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    double beta2 = 2857./54.*CA*CA*CA+Tf*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+Tf*Tf*nl*nl*(44./9.*CF+158./27.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    double Gamma2 = Gamma*A2;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = Gamma2*lambdaT*lambdaT/(lambdaT2p1);
    
    result += -Gamma1*beta1/beta0*(log(lambdaT2p1)-2.*(1.-lambdaT)*lambdaT)/2./lambdaT2p1;
    
    result += Gamma0*beta2/beta0*(lambdaT2p1*log(lambdaT2p1)-2.*lambdaT*(1.+lambdaT))/2./lambdaT2p1;
    
    result += Gamma0*beta1*beta1/beta0/beta0*pow(lambdaT-log(lambdaT2p1)/2.,2.)/lambdaT2p1;
    
    result += 2.*p2*Gamma1*beta0*log(2.)*lambdaT/lambdaT2p1;
    
    result += 2.*p2*Gamma0*beta1*log(2.)*(log(lambdaT2p1)-2.*lambdaT)/2./lambdaT2p1;
    
    result += Gamma1*beta0*log(Q2/mu2)*lambdaT/(lambdaT2p1);
    
    result += Gamma0*beta1*log(Q2/mu2)*(log(lambdaT2p1)-2.*lambdaT)/2./(lambdaT2p1);
    
    result += -2.*p2*Gamma0*beta0*beta0*log(Q2/mu2)*log(2.)*lambdaT/lambdaT2p1;
    
    result += -2.*p2*p2*Gamma0*beta0*beta0*log(2.)*log(2.)*lambdaT/(lambdaT2p1);
    
    result += -Gamma0*beta0*beta0*log(Q2/mu2)*log(Q2/mu2)/2.*lambdaT/(lambdaT2p1);
    
    result += -gamma1*beta0*lambdaT/lambdaT2p1;
    
    result += gamma0*beta1*(lambdaT-log(lambdaT2p1)/2.)/lambdaT2p1;
    
    result += 2.*p2*gamma0*beta0*beta0*log(2.)*lambdaT/lambdaT2p1;
    
    result += gamma0*beta0*beta0*log(Q2/mu2)*lambdaT/lambdaT2p1;
    
    result = result*Ci*alphas/4./pi/beta0/beta0;
    
    return(result);
    
}

double RpLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double Gamma0 = Gamma*A0;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = Ci*Gamma0/2./beta0*log(lambdaT2p1)*2.*pnu;
    
    return(result);
}

double RppLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double Gamma0 = Gamma*A0;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = alphas/4./pi*Ci*Gamma0*2./lambdaT2p1*pnu*pnu;
    
    return(result);
}

double RpppLL(double mu2, double Q2, double pnu, double pzc, double nub, double zc, double Gamma, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double Gamma0 = Gamma*A0;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = -pow(alphas/4./pi,2)*Ci*Gamma0*4.*beta0/lambdaT2p1/lambdaT2p1*pnu*pnu*pnu;
    
    return(result);
}

double RpNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = alphas/4./pi*Ci/beta0*(Gamma1*2.*lambdaT/lambdaT2p1
    +Gamma0*beta1/beta0*(log(lambdaT2p1)-2.*lambdaT)/lambdaT2p1
    +Gamma0*beta0/lambdaT2p1*(log(Q2/mu2)+2.*p2*log(2.)))*pnu;
    
    result += -alphas/4./pi*gamma0/2.*(2./lambdaT2p1)*pnu;
    
    return(result);
}

double RppNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = pow(alphas/4./pi,2)*Ci*(Gamma1*2./lambdaT2p1/lambdaT2p1
    -2.*Gamma0*beta1/beta0*log(lambdaT2p1)/lambdaT2p1/lambdaT2p1
    -2.*Gamma0*beta0/lambdaT2p1/lambdaT2p1*(log(Q2/mu2)+2.*p2*log(2.)))*pnu*pnu;
    
    result += pow(alphas/4./pi,2)*beta0*gamma0/2.*(4./lambdaT2p1/lambdaT2p1)*pnu*pnu;
    
    return(result);
}

double RpNNLL(double mu2, double Q2, double pnu, double pzc, double p2, double nub, double zc, double Gamma, double gamma0, double gamma1, double alphas){
    double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double A2 = 4.*CA*CA*(245./6.-134.*pi*pi/27.+11.*pow(pi,4.)/45.+22.*Zeta3/3.)+32.*CA*Tf*nl*(-209./108.+5.*pi*pi/27.-7.*Zeta3/3.)+4.*CF*Tf*nl*(16*Zeta3-55./3.)-64./27.*Tf*Tf*nl*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    double beta2 = 2857./54.*CA*CA*CA+Tf*nl*(2.*CF*CF-205./9.*CF*CA-1415./27.*CA*CA)+Tf*Tf*nl*nl*(44./9.*CF+158./27.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    double Gamma2 = Gamma*A2;
    
    double lambdanub = beta0/4./pi*alphas*log(nub);
    double lambdazc = beta0/4./pi*alphas*log(zc);
    
    double lambdaT = pnu*lambdanub+pzc*lambdazc;
    
    double lambdaT2p1 = 1.+2.*lambdaT;
    
    double result = Gamma2*2.*lambdaT*(1.+lambdaT)/lambdaT2p1/lambdaT2p1;
    
    result += Gamma1*beta1/beta0*(log(lambdaT2p1)-2.*(1.+lambdaT)*lambdaT)/lambdaT2p1/lambdaT2p1;
    
    result += Gamma0*beta2/beta0*(-2.*lambdaT*lambdaT)/lambdaT2p1/lambdaT2p1;
    
    result += Gamma0*beta1*beta1/beta0/beta0*(log(lambdaT2p1)*log(lambdaT2p1)-4.*lambdaT*lambdaT)/2./lambdaT2p1/lambdaT2p1;
    
    result += 2.*p2*Gamma1*beta0*log(2.)/lambdaT2p1/lambdaT2p1;
    
    result += -2.*p2*Gamma0*beta1*log(2.)*log(lambdaT2p1)/lambdaT2p1/lambdaT2p1;
    
    result += Gamma1*beta0*log(Q2/mu2)/lambdaT2p1/lambdaT2p1;
    
    result += -Gamma0*beta1*log(Q2/mu2)*log(lambdaT2p1)/lambdaT2p1/lambdaT2p1;
    
    result += -2.*p2*Gamma0*beta0*beta0*log(Q2/mu2)*log(2.)/lambdaT2p1/lambdaT2p1;
    
    result += -2.*p2*p2*Gamma0*beta0*beta0*log(2.)*log(2.)/lambdaT2p1/lambdaT2p1;
    
    result += -Gamma0*beta0*beta0*log(Q2/mu2)*log(Q2/mu2)/2./lambdaT2p1/lambdaT2p1;
    
    result += -gamma1*beta0/lambdaT2p1/lambdaT2p1;
    
    result += gamma0*beta1*log(lambdaT2p1)/lambdaT2p1/lambdaT2p1;
    
    result += 2.*p2*gamma0*beta0*beta0*log(2.)/lambdaT2p1/lambdaT2p1;
    
    result += gamma0*beta0*beta0*log(Q2/mu2)/lambdaT2p1/lambdaT2p1;
    
    result = result*Ci*pow(alphas/4./pi,2)/beta0*pnu;
    
    return(result);
}

double exp1p(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0){
//     double Ci = 1.;
    double A0 = 4.;
//     double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double Gamma0 = Gamma*A0;
    
    double LT = pnu*log(nub)+pzc*log(zc);
    
    double expansion = 0.;
    if(LLorder > 0.) expansion += 2.*Gamma0*LT;
    if(LLorder > 1.) expansion += -gamma0+Gamma0*(log(Q2/mu2)+2.*p2*log(2.));
    
    expansion = expansion*(-pnu);
    
    return(expansion);
}

double exp2p(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0, double gamma1, double C1){
//     double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
//     double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    
    double LT = pnu*log(nub)+pzc*log(zc);
    
    double expansion = 0.;
    if(LLorder > 0.) expansion += -2.*beta0*Gamma0*LT*LT;
    if(LLorder > 1.) expansion += 2.*(Gamma1+beta0*gamma0-beta0*Gamma0*(log(Q2/mu2)+2.*p2*log(2.)))*LT;
    if(LLorder > 2.) expansion += -(gamma1+2.*beta0*C1)+(Gamma1+beta0*gamma0-beta0*Gamma0*(log(Q2/mu2)+2.*p2*log(2.))/2.)*(log(Q2/mu2)+2.*p2*log(2.));
    
    expansion = expansion*(-pnu);
    
    return(expansion);
}

double exp1(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0){
//     double Ci = 1.;
    double A0 = 4.;
//     double beta0 = 11./3.*CA-4./3.*Tf*nl;
    
    double Gamma0 = Gamma*A0;
    
    double LT = pnu*log(nub)+pzc*log(zc);
    
    double expansion = 0.;
    if(LLorder > 0.) expansion += Gamma0*LT*LT;
    if(LLorder > 1.) expansion += (-gamma0+Gamma0*(log(Q2/mu2)+2.*p2*log(2.)))*LT;
    
    return(expansion);
}

double exp2(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0, double gamma1, double C1){
//     double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
//     double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    
    double LT = pnu*log(nub)+pzc*log(zc);
    
    double expansion = 0.;
    if(LLorder > 0.) expansion += -2./3.*beta0*Gamma0*LT*LT*LT;
    if(LLorder > 1.) expansion += (Gamma1+beta0*gamma0-beta0*Gamma0*(log(Q2/mu2)+2.*p2*log(2.)))*LT*LT;
    if(LLorder > 2.) expansion += (-(gamma1+2.*beta0*C1)+(Gamma1+beta0*gamma0-beta0*Gamma0*(log(Q2/mu2)+2.*p2*log(2.))/2.)*(log(Q2/mu2)+2.*p2*log(2.)))*LT;
    
    return(expansion);
}

double exp3(double mu2, double Q2, double pnu, double pzc, double p2,  double nub, double zc, double Gamma, double gamma0, double gamma1, double C1){
//     double Ci = 1.;
    double A0 = 4.;
    double A1 = 4.*CA*(67./9.-pi*pi/3.)-80./9.*Tf*nl;
    double A2 = 4.*CA*CA*(245./6.-134.*pi*pi/27.+11.*pow(pi,4.)/45.+22.*Zeta3/3.)+32.*CA*Tf*nl*(-209./108.+5.*pi*pi/27.-7.*Zeta3/3.)+4.*CF*Tf*nl*(16*Zeta3-55./3.)-64./27.*Tf*Tf*nl*nl;
    double beta0 = 11./3.*CA-4./3.*Tf*nl;
    double beta1 = 34./3.*CA*CA-4.*Tf*nl*(CF+5./3.*CA);
    
    double Gamma0 = Gamma*A0;
    double Gamma1 = Gamma*A1;
    double Gamma2 = Gamma*A2;
    
    double LT = pnu*log(nub)+pzc*log(zc);
    double y = log(Q2/mu2)+2.*p2*log(2.);
    
    double expansion = 0.;
    if(LLorder > 0.) expansion += 2./3.*beta0*beta0*Gamma0*LT*LT*LT*LT;
    if(LLorder > 1.) expansion += -2./3.*(2.*beta0*(Gamma1+beta0*gamma0-beta0*Gamma0*y)+beta1*Gamma0)*LT*LT*LT;
    if(LLorder > 2.) expansion += (Gamma2+2.*beta0*(gamma1+2.*beta0*C1)-2.*beta0*Gamma1*y+(beta1-2.*beta0*beta0*y)*gamma0-(beta1-beta0*beta0*y)*y*Gamma0)*LT*LT;
    
    return(expansion);
}
