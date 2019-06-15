#include "resum-blocks.hh"

#include <cmath>
#include <cassert>
#include <limits>
#include <iostream>

using namespace std;


// sudakov frequent bits
//
// Note: to obtain fixed coupling (at scale kt), set Lfr=0
// i.e. mufr=1.0 and the appropriate alphas value in the ConfigBase.


//----------------------------------------------------------------------
// internal flags and helpers
//----------------------------------------------------------------------
const double go_below_freeze = 1.0; // 1 for continuous freezing, 0 for a sharp cut

#define lam(X) (2.0*alphas_ref*b0*(X))

//========================================================================
// running coupling
//========================================================================

// this is the expression for alphas suited for NLL resummation
double ConfigBase::alphas(double log_kt, bool allow_freezing, bool force_no_2loops) const{
  double lkt_eff = (allow_freezing && (log_kt>Lfr)) ? Lfr : log_kt;
  double den = 1.0-lam(lkt_eff);

  // if we're sitting below the lamdau pole, return "infinity"
  if (den<=0) return numeric_limits<double>::max();

  return alphas_ref/den
    * ((two_loops && (!force_no_2loops))
       ? 1.0+alphas_ref*(-b1/b0*log(den) + 0.5*Keff/M_PI)/den
       : 1.0);
}


// this is the expression for alphas at 1 loop suited for NLL resummation
double ConfigBase::alphas_1loop(double log_kt) const{
  double den = 1.0-lam(log_kt);

  // if we're sitting below the lamdau pole, return "infinity"
  return (den<=0) ? numeric_limits<double>::max() : alphas_ref/den;
}


// this is the exact derivative.
// it is trivial to show that it is equal to
//    2 b0 alphas^2(log_kt) + 2 b1 alphas^2(log_kt) + O(alphas^4)
// This bit is used in double-derivatives of the radiators 
double ConfigBase::alphas_derivative(double log_kt, bool allow_freezing, bool force_no_2loops) const{
  if ((allow_freezing) && (log_kt>Lfr)) return 0.0;

  double den = 1.0-lam(log_kt);
  // if we're sitting below the lamdau pole, return "infinity"
  if (den<=0) return numeric_limits<double>::max();

  return 2*b0*alphas_ref*alphas_ref/(den*den)
    * ((two_loops && (!force_no_2loops))
       ? 1.0+alphas_ref*(-b1/b0*(2*log(den)-1) + Keff/M_PI)/den
       : 1.0);
}


//========================================================================
// radiators
//========================================================================

// a triangle bound by an upper theta limit, a kt limit and a z
// theta^alpha limit
//
// For alpha>1, the kt limit is an upper one, the z theta^alpha limit
// is a lower one
// For alpha<1, the kt limit is a lower one, the z theta^alpha limit
// is an upper one
//
// For alpha=0 and B!=0, a B term is included. It is the user's
// responsibility to included it only when needed
//
// constraints: - alpha != 1
//              - Ltop<Lbot
//
// assumptions: the full triangle is in the kinematic acceptance
// [note however that we can trick it if we want to take
// differences of two such triangles where th epart outside of the
// soundary cancels]
double ConfigBase::kt_triangle(double alpha, double Ltop, double Lbot, double B) const{
  // consistency checks
  assert(alpha!=1.0);
  if (Lbot-Ltop<=-1e-8) cerr << "---> " << Ltop << " " << Lbot << endl;
  assert(Lbot-Ltop>=-1e-8);
  assert((std::abs(alpha)<1e-5) || (std::abs(B)<1e-5));

  // handle the Landau pole if we do not freeze the coupling
  if ((!freeze) && (2*alphas_ref*b0*Lbot>=1.0)) return numeric_limits<double>::max();

  if (alpha>1){
    // both limits above the freezing scale
    if (Lbot<Lfr){
      double lbot = lam(Lbot);
      double ltop = lam(Ltop);
      double lgbot = log(1-lbot);
      double lgtop = log(1-ltop);
      double lg = lgbot-lgtop;
      return CR/(2*M_PI*alphas_ref*b0*b0*(alpha-1))
        *((1-lbot)*lg + (lbot-ltop)
          + (two_loops
             ? 0.5*alphas_ref*Keff/M_PI*((ltop-lbot)/(1-ltop)-lg)
             - alphas_ref*b1/b0*(0.5*lgtop*lgtop-0.5*lgbot*lgbot-lgbot+(1-lbot)/(1-ltop)*lgtop+(ltop-lbot)/(1-ltop))
             : 0.0));
    }
    
    // sitting over the freezing scale
    if (Ltop<Lfr){
      double lbot = lam(Lbot);
      double lfr  = lam(Lfr);
      double ltop = lam(Ltop);
      double lgfr  = log(1-lfr);
      double lgtop = log(1-ltop);
      return CR/(2*M_PI*alphas_ref*b0*b0*(alpha-1))
        *((1-lbot)*(lgfr-lgtop) + (lfr-ltop)
        + (two_loops
           ? 0.5*alphas_ref*Keff/M_PI*((ltop-lfr)*(1-lbot)/((1-ltop)*(1-lfr))+(lgtop-lgfr))
           - alphas_ref*b1/b0*(0.5*lgtop*lgtop-0.5*lgfr*lgfr-(1-lbot)/(1-lfr)*lgfr+(1-lbot)/(1-ltop)*lgtop+(ltop-lfr)*(1-lbot)/((1-ltop)*(1-lfr)))
           : 0.0))
        +go_below_freeze*alphas(Lfr)*CR/(M_PI*(alpha-1))*(Lbot-Lfr)*(Lbot-Lfr);
    }
    
    // both limits below the freezing scale
    return go_below_freeze*alphas(Lfr)*CR/(M_PI*(alpha-1))*(Lbot-Ltop)*(Lbot-Ltop);
  }

  // now triangles pointing upwards
  if (Lbot<Lfr){
    double lbot = lam(Lbot);
    double ltop = lam(Ltop);
    double lgbot = log(1-lbot);
    double lgtop = log(1-ltop);
    double lg = lgtop-lgbot;
    double Bfact = 1.0+lam(B)/(1-ltop);
    return CR/(2*M_PI*alphas_ref*b0*b0*(1-alpha))
      *((1-ltop)*Bfact*lg + (ltop-lbot)
        + (two_loops
           ? 0.5*alphas_ref*Keff/M_PI*((lbot-ltop)/(1-lbot)-lg)
           - alphas_ref*b1/b0*(0.5*lgbot*lgbot-0.5*lgtop*lgtop
                         -lgtop+(1-ltop)/(1-lbot)*lgbot
                         +(lbot-ltop)/(1-lbot))
           : 0.0));
  }
    
  // sitting over the freezing scale
  if (Ltop<Lfr){
    double lbot = lam(Lfr);  // misnamed so we can reuse the same final expression as above
    double ltop = lam(Ltop);
    double lgbot = log(1-lbot);
    double lgtop = log(1-ltop);
    double lg = lgtop-lgbot;
    double Bfact = 1.0+lam(B)/(1-ltop);
    return CR/(2*M_PI*alphas_ref*b0*b0*(1-alpha))
      *((1-ltop)*Bfact*lg + (ltop-lbot)
        + (two_loops
           ? 0.5*alphas_ref*Keff/M_PI*((lbot-ltop)/(1-lbot)-lg)
           - alphas_ref*b1/b0*(0.5*lgbot*lgbot-0.5*lgtop*lgtop
                         -lgtop+(1-ltop)/(1-lbot)*lgbot
                         +(lbot-ltop)/(1-lbot))
           : 0.0))
      +go_below_freeze*CR/(M_PI*(1-alpha))*(alphas(Lfr)*(Lbot+Lfr-2*Ltop)+2*alphas_1loop(Lfr)*B)*(Lbot-Lfr);
  }


  // both limits below the freezing scale
  return go_below_freeze*CR/(M_PI*(1-alpha))*(Lbot-Ltop)*(alphas(Lfr)*(Lbot-Ltop)+2*alphas_1loop(Lfr)*B);
}

//------------------------------------------------------------------------
double ConfigBase::triangle(double alpha, double beta, double Ltop, double Lbot, double B) const{
  assert(beta>alpha);

  // start dealing with cases with simplifications
  if (std::abs(alpha-1)<1e-4) return kt_triangle(beta ,Ltop,Lbot,0.0);
  if (std::abs(beta -1)<1e-4) return kt_triangle(alpha,Ltop,Lbot,  B);

  // the intermediate point
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);

  if (alpha<1){
    if (beta>1)
      return kt_triangle(alpha,Ltop,Lmed,B)+kt_triangle(beta,Lmed,Lbot,0.0);

    // note that we force the 2nd B to 0. This could probably be done a bit more elegantly
    return kt_triangle(alpha,Ltop,Lmed,B)-kt_triangle(beta,Lbot,Lmed,0.0);
  }
  
  return kt_triangle(beta,Lmed,Lbot,0.0)-kt_triangle(alpha,Lmed,Ltop,0.0);
}

//========================================================================
// derivatives of triangles
//========================================================================

double ConfigBase::delta_kt_line(double alpha, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(alpha!=1);

  // if the line is upside-down, just revert it (this is merely a convenience)
  if (Lbot<Ltop) return delta_kt_line(alpha, Lbot, Ltop, B, force_no_2loops);
  assert(Lbot-Ltop>=-1e-8);

  // for the Bi term, we need alpha>0
  if (std::abs(B)>1e-4) assert(alpha>0);
  
  // handle the Landau pole if we do not freeze the coupling
  if ((!freeze) && (2*alphas_ref*b0*Lbot>=1.0)) return 0.0;

  bool local_two_loops = two_loops && (!force_no_2loops);

  // determine the scale that optionally enters in the B term
  double LB = ((alpha>0) && (alpha<1)) ? Lbot : Ltop;
  if (LB>Lfr) LB=Lfr;

  // now the results
  if (Lbot<Lfr){
    double lbot = lam(Lbot);
    double ltop = lam(Ltop);
    return CR/(M_PI*b0*std::abs(alpha-1))*
      (log((1-ltop)/(1-lbot))
       +(local_two_loops
         ? alphas_ref*Keff/(2*M_PI)*(lbot-ltop)/((1-ltop)*(1-lbot))
         - alphas_ref*b1/b0*( log(1-lbot)/(1-lbot) - log(1-ltop)/(1-ltop)
                          + (lbot-ltop)/((1-ltop)*(1-lbot)) )
         : 0.0))
      +2*alphas_1loop(LB)*CR/(alpha*M_PI)*B;
  }

  if (Ltop<Lfr){
    // split into a contribution with RC from Lfr to Ltop and a frozen
    // contribution between Lbot and Lfr
    //
    // Note that we need to be careful to discard te B term if it
    // falls below the freezing region and we're actually not going
    // there!
    double lbot = lam(Lfr);  // trick to get the same RC expression as above
    double ltop = lam(Ltop);
    return CR/(M_PI*b0*std::abs(alpha-1))*
      (log((1-ltop)/(1-lbot))
       +(local_two_loops
         ? alphas_ref*Keff/(2*M_PI)*(lbot-ltop)/((1-ltop)*(1-lbot))
         - alphas_ref*b1/b0*( log(1-lbot)/(1-lbot) - log(1-ltop)/(1-ltop)
                        + (lbot-ltop)/((1-ltop)*(1-lbot)) )
         : 0.0))
      +go_below_freeze*2*alphas(Lfr,false,force_no_2loops)*CR/(std::abs(alpha-1)*M_PI)*(Lbot-Lfr)
      +(LB<Lfr || go_below_freeze ? 1.0 : 0.0)*2*alphas_1loop(LB)*CR/(alpha*M_PI)*B;
  }
  
  return go_below_freeze*2*CR/M_PI
    *( alphas(Lfr,false,force_no_2loops)*(Lbot-Ltop)/std::abs(alpha-1) + alphas_1loop(Lfr)*B/alpha );
}


// same as above, specified on a basis of a line given by the
// coordinates of its endpoints. This also supports a constant kt
// line.
//
// The line starts from (theta0, kt0) to (theta1, kt1)
// Note that this is theta, not theta^2
//
// Note also that there is no check performed if z=1 is included,
// i.e. the B term is included blindly and it it the end-user's
// responsibility to set B to 0 when needed.
double ConfigBase::line(double Ltheta0, double Lkt0, double Ltheta1, double Lkt1, double B, bool force_no_2loops) const{
  
  // figure out the slope
  double alpha = 1-(Lkt1-Lkt0)/(Ltheta1-Ltheta0);
  if (std::abs(alpha-1)>1e-4)
    return delta_kt_line(alpha, Lkt0, Lkt1, B, force_no_2loops);

  double Las = Lkt0 < Lfr ? Lkt0 : Lfr;
  return 2.0*CR/M_PI*(alphas(Las,true,force_no_2loops)*std::abs(Ltheta0-Ltheta1)+alphas_1loop(Las)*B);
}

//------------------------------------------------------------------------
// the R' corresponding to the bottom-line of the triangle
double ConfigBase::triangle_prime(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);

  // the kt of the intermediate point
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);
  if (std::abs(beta-1)>1e-4)
    return (beta>1) ? delta_kt_line(beta, Lmed, Lbot, B, force_no_2loops)
                    : delta_kt_line(beta, Lbot, Lmed, B, force_no_2loops);
  
  // angle of the intermediate point
  //
  // The top line has z theta^alpha = cst
  // The top can be taken to have Ltheta=0, Lz=Ltop
  //  => the theta of the intermediate point is
  //     (alpha-1) Ltheta + Lmed = Ltop
  //  -> Ltheta = (Ltop-Lmed)/(alpha-1)
  //            = (Lbot-Ltop))/(beta-alpha)
  double Ltheta = (Lbot-Ltop)/(beta-alpha);
  return line(0.0, Lbot, Ltheta, Lmed, B, force_no_2loops);
}

//========================================================================
// double derivatives of triangles
//========================================================================

double ConfigBase::triangle_double_prime(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);
  //assert(alpha>=0.0);

  // note in what follows that the B term uses the 1-loop alphas
  
  if (std::abs(beta-1)>1e-4){
    double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);
    return 2*CR/(M_PI*(beta-1))*(alphas(Lbot,true,force_no_2loops)
                                -(1-alpha)/(beta-alpha)*alphas(Lmed,true,force_no_2loops))
      +2*B*CR/(beta*M_PI)*(1-alpha)/(beta-alpha)*alphas_derivative(Lmed,true,true);
  }
  
  // note; bottom is likely wrong  
  return 2.0*CR/(M_PI*(1-alpha))*alphas(Lbot,true,force_no_2loops)
    +2*CR/M_PI*((Lbot-Ltop)/(1-alpha)*alphas_derivative(Lbot,true,force_no_2loops)+B*alphas_derivative(Lbot,true,true));
}

//========================================================================
// LO expansions
//========================================================================
double ConfigBase::kt_triangle_LO(double alpha, double Ltop, double Lbot, double B) const{
  // consistency checks
  assert(alpha!=1.0);
  if (Lbot-Ltop<=-1e-8) cerr << "---> " << Ltop << " " << Lbot << endl;
  assert(Lbot-Ltop>=-1e-8);
  assert((std::abs(alpha)<1e-5) || (std::abs(B)<1e-5));
  return CR*alphas_ref/(M_PI*std::abs(alpha-1))*(Lbot-Ltop)*(Lbot-Ltop-2*(alpha-1)*B);
}

//------------------------------------------------------------------------
double ConfigBase::triangle_LO(double alpha, double beta, double Ltop, double Lbot, double B) const{
  assert(beta>alpha);
  if (std::abs(alpha-1)<1e-4) return kt_triangle_LO(beta ,Ltop,Lbot,0.0);
  if (std::abs(beta -1)<1e-4) return kt_triangle_LO(alpha,Ltop,Lbot,  B);
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);

  if (alpha<1)
    return (beta>1) ? kt_triangle_LO(alpha,Ltop,Lmed,B)+kt_triangle_LO(beta,Lmed,Lbot,0.0)
                    : kt_triangle_LO(alpha,Ltop,Lmed,B)-kt_triangle_LO(beta,Lbot,Lmed,0.0);
  return kt_triangle_LO(beta,Lmed,Lbot,0.0)-kt_triangle_LO(alpha,Lmed,Ltop,0.0);
}

double ConfigBase::delta_kt_line_LO(double alpha, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(alpha!=1);
  if (Lbot<Ltop) return delta_kt_line_LO(alpha, Lbot, Ltop, B, force_no_2loops);
  if (std::abs(B)>1e-4) assert(alpha>0);
  return 2*CR*alphas_ref/M_PI*((Lbot-Ltop)/std::abs(alpha-1)+B/alpha);
}

double ConfigBase::line_LO(double Ltheta0, double Lkt0, double Ltheta1, double Lkt1, double B, bool force_no_2loops) const{
  
  // figure out the slope
  double alpha = 1-(Lkt1-Lkt0)/(Ltheta1-Ltheta0);
  if (std::abs(alpha-1)>1e-4)
    return delta_kt_line_LO(alpha, Lkt0, Lkt1, B, force_no_2loops);

  return 2*CR/M_PI*(std::abs(Ltheta0-Ltheta1)+B)*alphas_ref;
}

double ConfigBase::triangle_prime_LO(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);
  if (std::abs(beta-1)>1e-4)
    return (beta>1) ? delta_kt_line_LO(beta, Lmed, Lbot, B, force_no_2loops)
                    : delta_kt_line_LO(beta, Lbot, Lmed, B, force_no_2loops);
  double Ltheta = (Ltop-Lmed)/(alpha-1);
  return line_LO(0.0, Lbot, Ltheta, Lmed, B, force_no_2loops);
}

double ConfigBase::triangle_double_prime_LO(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);
  //assert(alpha>=0.0);
  return 2*CR*alphas_ref/(M_PI*(beta-alpha));
}


//========================================================================
// NLO expansions
//========================================================================
double ConfigBase::kt_triangle_NLO(double alpha, double Ltop, double Lbot, double B) const{
  // consistency checks
  assert(alpha!=1.0);
  if (Lbot-Ltop<=-1e-8) cerr << "---> " << Ltop << " " << Lbot << endl;
  assert(Lbot-Ltop>=-1e-8);
  assert((std::abs(alpha)<1e-5) || (std::abs(B)<1e-5));

  if (alpha<1){ // possibly a B term
    return CR*alphas_ref*alphas_ref/(M_PI*(1-alpha))
      *(Lbot-Ltop)
      *(2.0/3.0*b0*((Lbot-Ltop)*(2*Lbot+Ltop)+3*(1-alpha)*(Ltop+Lbot)*B)
        + (two_loops ? Keff/(2*M_PI)*(Lbot-Ltop) : 0.0));
  }

  // definitely no B term in this setup
  return CR*alphas_ref*alphas_ref/(M_PI*(alpha-1))
    *(Lbot-Ltop)*(Lbot-Ltop) *(2.0/3.0*b0*(2*Ltop+Lbot) + (two_loops ? Keff/(2*M_PI) : 0.0));
}

//------------------------------------------------------------------------
double ConfigBase::triangle_NLO(double alpha, double beta, double Ltop, double Lbot, double B) const{
  assert(beta>alpha);
  if (std::abs(alpha-1)<1e-4) return kt_triangle_NLO(beta ,Ltop,Lbot,0.0);
  if (std::abs(beta -1)<1e-4) return kt_triangle_NLO(alpha,Ltop,Lbot,  B);
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);

  if (alpha<1)
    return (beta>1) ? kt_triangle_NLO(alpha,Ltop,Lmed,B)+kt_triangle_NLO(beta,Lmed,Lbot,0.0)
                    : kt_triangle_NLO(alpha,Ltop,Lmed,B)-kt_triangle_NLO(beta,Lbot,Lmed,0.0);
  return kt_triangle_NLO(beta,Lmed,Lbot,0.0)-kt_triangle_NLO(alpha,Lmed,Ltop,0.0);
}

double ConfigBase::delta_kt_line_NLO(double alpha, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(alpha!=1);
  if (Lbot<Ltop) return delta_kt_line_NLO(alpha, Lbot, Ltop, B, force_no_2loops);
  if (std::abs(B)>1e-4) assert(alpha>0);
  if (alpha<1){
    return 2*(Lbot-Ltop)*CR/(M_PI*(1-alpha))*alphas_ref*alphas_ref
      *((Lbot+Ltop)*b0 + ((two_loops && (!force_no_2loops)) ? Keff/(2*M_PI) : 0.0))
      +4*b0*CR*B*alphas_ref*alphas_ref/(alpha*M_PI)*Lbot;
  }
  return 2*(Lbot-Ltop)*CR/(M_PI*(alpha-1))*alphas_ref*alphas_ref
    *((Lbot+Ltop)*b0 + ((two_loops && (!force_no_2loops)) ? Keff/(2*M_PI) : 0.0))
    +4*b0*CR*B*alphas_ref*alphas_ref/(alpha*M_PI)*Ltop;
}

double ConfigBase::line_NLO(double Ltheta0, double Lkt0, double Ltheta1, double Lkt1, double B, bool force_no_2loops) const{
  
  // figure out the slope
  double alpha = 1-(Lkt1-Lkt0)/(Ltheta1-Ltheta0);
  if (std::abs(alpha-1)>1e-4)
    return delta_kt_line_NLO(alpha, Lkt0, Lkt1, B, force_no_2loops);

  return 2*CR*alphas_ref*alphas_ref/M_PI
    *(2*(std::abs(Ltheta0-Ltheta1)+B)*Lkt0*b0
      + ((two_loops && (!force_no_2loops)) ? Keff/(2*M_PI)*std::abs(Ltheta0-Ltheta1) : 0.0));
}

double ConfigBase::triangle_prime_NLO(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);
  if (std::abs(beta-1)>1e-4)
    return (beta>1) ? delta_kt_line_NLO(beta, Lmed, Lbot, B, force_no_2loops)
                    : delta_kt_line_NLO(beta, Lbot, Lmed, B, force_no_2loops);
  double Ltheta = (Ltop-Lmed)/(alpha-1);
  return line_NLO(0.0, Lbot, Ltheta, Lmed, B, force_no_2loops);
}


double ConfigBase::triangle_double_prime_NLO(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);
  //assert(alpha>=0.0);
  return CR*alphas_ref*alphas_ref/((beta-alpha)*M_PI)
    *(4*b0*((1-2*alpha+beta)*Lbot-(1-alpha)*Ltop)/(beta-alpha)
      +4*(1-alpha)/beta*b0*B
      +((two_loops && (!force_no_2loops)) ? Keff/M_PI : 0.0));
}


//========================================================================
// NNLO expansions
//========================================================================
double ConfigBase::kt_triangle_NNLO(double alpha, double Ltop, double Lbot, double B) const{
  // consistency checks
  assert(alpha!=1.0);
  if (Lbot-Ltop<=-1e-8) cerr << "---> " << Ltop << " " << Lbot << endl;
  assert(Lbot-Ltop>=-1e-8);
  assert((std::abs(alpha)<1e-5) || (std::abs(B)<1e-5));

  if (alpha<1)
    return 2*alphas_ref*alphas_ref*alphas_ref*CR*(Lbot-Ltop)/(3*(1-alpha)*M_PI)
      *(((3*Lbot*Lbot+2*Lbot*Ltop+Ltop*Ltop)*(Lbot-Ltop)
         +(Lbot*Lbot+Lbot*Ltop+Ltop*Ltop)*4*(1-alpha)*B)*b0*b0
        + (two_loops ? (2*Lbot+Ltop)*(Lbot-Ltop)*(K*b0/M_PI+b1) : 0.0));
  
  return 2*alphas_ref*alphas_ref*alphas_ref*CR*(Lbot-Ltop)*(Lbot-Ltop)/(3*(alpha-1)*M_PI)
    *((Lbot*Lbot+2*Lbot*Ltop+3*Ltop*Ltop)*b0*b0
      + (two_loops ? (Lbot+2*Ltop)*(K*b0/M_PI+b1) : 0.0));
}

//------------------------------------------------------------------------
double ConfigBase::triangle_NNLO(double alpha, double beta, double Ltop, double Lbot, double B) const{
  assert(beta>alpha);
  if (std::abs(alpha-1)<1e-4) return kt_triangle_NNLO(beta ,Ltop,Lbot,0.0);
  if (std::abs(beta -1)<1e-4) return kt_triangle_NNLO(alpha,Ltop,Lbot,  B);
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);

  if (alpha<1)
    return (beta>1) ? kt_triangle_NNLO(alpha,Ltop,Lmed,B)+kt_triangle_NNLO(beta,Lmed,Lbot,0.0)
                    : kt_triangle_NNLO(alpha,Ltop,Lmed,B)-kt_triangle_NNLO(beta,Lbot,Lmed,0.0);
  return kt_triangle_NNLO(beta,Lmed,Lbot,0.0)-kt_triangle_NNLO(alpha,Lmed,Ltop,0.0);
}

double ConfigBase::delta_kt_line_NNLO(double alpha, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(alpha!=1);
  if (Lbot<Ltop) return delta_kt_line_NNLO(alpha, Lbot, Ltop, B, force_no_2loops);
  if (std::abs(B)>1e-4) assert(alpha>0);
  if (alpha<1){
    return 4*(Lbot-Ltop)*CR/(M_PI*(1-alpha))*alphas_ref*alphas_ref*alphas_ref
      *(2.0/3.0*(Lbot*Lbot+Lbot*Ltop+Ltop*Ltop)*b0*b0
        + ((two_loops && (!force_no_2loops)) ? (Lbot+Ltop)*(Keff/(2*M_PI)*b0+0.5*b1) : 0.0))
      +8*b0*b0*CR*B*alphas_ref*alphas_ref*alphas_ref/(alpha*M_PI)*Lbot*Lbot;
  }
  return 4*(Lbot-Ltop)*CR/(M_PI*(alpha-1))*alphas_ref*alphas_ref*alphas_ref
    *(2.0/3.0*(Lbot*Lbot+Lbot*Ltop+Ltop*Ltop)*b0*b0
      + ((two_loops && (!force_no_2loops)) ? (Lbot+Ltop)*(Keff/(2*M_PI)*b0+0.5*b1) : 0.0))
    +8*b0*b0*CR*B*alphas_ref*alphas_ref*alphas_ref/(alpha*M_PI)*Ltop*Ltop;
}

double ConfigBase::line_NNLO(double Ltheta0, double Lkt0, double Ltheta1, double Lkt1, double B, bool force_no_2loops) const{
  
  // figure out the slope
  double alpha = 1-(Lkt1-Lkt0)/(Ltheta1-Ltheta0);
  if (std::abs(alpha-1)>1e-4)
    return delta_kt_line_NNLO(alpha, Lkt0, Lkt1, B, force_no_2loops);

  return 8*CR*alphas_ref*alphas_ref*alphas_ref/M_PI*Lkt0
    *(Lkt0*b0*b0*(std::abs(Ltheta0-Ltheta1)+B) +
      ((two_loops && (!force_no_2loops)) ? (Keff/(2*M_PI)*b0+0.5*b1)*std::abs(Ltheta0-Ltheta1) : 0.0));
}

double ConfigBase::triangle_prime_NNLO(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);
  double Lmed = ((beta-1)*Ltop + (1-alpha)*Lbot)/(beta-alpha);
  if (std::abs(beta-1)>1e-4)
    return (beta>1) ? delta_kt_line_NNLO(beta, Lmed, Lbot, B, force_no_2loops)
                    : delta_kt_line_NNLO(beta, Lbot, Lmed, B, force_no_2loops);
  double Ltheta = (Ltop-Lmed)/(alpha-1);
  return line_NNLO(0.0, Lbot, Ltheta, Lmed, B, force_no_2loops);
}

double ConfigBase::triangle_double_prime_NNLO(double alpha, double beta, double Ltop, double Lbot, double B, bool force_no_2loops) const{
  assert(beta>alpha);
  //assert(alpha>=0.0);
  if (std::abs(beta-1)>1e-4){
    double Lmed = ((1-alpha)*Lbot+(beta-1)*Ltop)/(beta-alpha);
    return 4*CR*alphas_ref*alphas_ref*alphas_ref/M_PI
      *(2*b0*b0/(beta-1)*(Lbot*Lbot-(1-alpha)/(beta-alpha)*Lmed*Lmed)
        +4*(1-alpha)/(beta*(beta-alpha))*b0*b0*B*Lmed
        +((two_loops && (!force_no_2loops)) ? (Keff*b0/M_PI+b1)/(beta-1)*(Lbot-(1-alpha)/(beta-alpha)*Lmed) : 0.0));
  }

  return 4*CR*alphas_ref*alphas_ref*alphas_ref/M_PI
    *(2*(3*Lbot-2*Ltop)*Lbot*b0*b0/(1-alpha)+4*Lbot*b0*b0*B
      +((two_loops && (!force_no_2loops)) ? (Keff*b0/M_PI+b1)*(2*Lbot-Ltop)/(1-alpha) : 0.0));
}
