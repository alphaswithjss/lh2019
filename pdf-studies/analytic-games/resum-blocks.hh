#ifndef __BUILDING_BLOCKS_HH__
#define __BUILDING_BLOCKS_HH__

#include <cmath>
#include "qcd-constants.hh"

//----------------------------------------------------------------------
// put all the parameters in a structure (easier to pass as an argument)
//----------------------------------------------------------------------
/// \class ConfigBase
/// Basic class for Sudakov building blocks
///
/// In the code below, we will refer to a "line with slope alpha" a
/// line such that z \theta^\alpha = constant
class ConfigBase{
public:
  /// default ctor
  /// WATCH OUT: leaves things in an unphysical state
  ConfigBase()
    : CR(0.0), alphas_ref(0.1), two_loops(false){
    Bi = 0.0;
    Keff = 0.0;
    set_freeze(-1.0);
  }

  /// ctor
  ///  CR      colour factor (CF for quarks, CA for gluons)
  ///  alphas  alphas(ptR) [optionally includes scale variations]
  ///  mufr    alphas freezing scale (in units of the alphas scale)
  ///
  /// By default, a 2-loop running coupling is used
  ConfigBase(double CR_in, double alphas_in, double mufr=-1.0, double muR_factor=1.0)
    : CR(CR_in), alphas_ref(alphas_in),
      two_loops(true){
    Bi = (std::abs(CR-CA)<1e-4) ? Bg : Bq;
    Keff = K+4*M_PI*b0*log(muR_factor);
    set_freeze(mufr);
  }

  /// set the freezing scale (in units of the scale used for alphas)
  void set_freeze(double mufr=-1.0){
    if (mufr<0){
      freeze=false;
      Lfr=-1.0;
    }
    
    freeze=true;
    Lfr = log(1/mufr);
  }

  /// specify if we should use a 2-loop running coupling
  void set_two_loops(bool value){ two_loops = value; }

  /// some internal parameters (kept public for simplicity)
  double CR, Bi, Keff;
  double alphas_ref;
  double Lfr;
  bool freeze, two_loops;

  //----------------------------------------------------------------------
  // running copling
  //----------------------------------------------------------------------
  /// this is the expression for alphas suited for NLL resummation
  double alphas(double log_kt, bool allow_freezing=false, bool force_no_2loops=false) const;

  /// this is the expression for alphas at one loop  suited for NLL resummation (Bi terms)
  double alphas_1loop(double log_kt) const;

    /// the (exact) derivative of alphas wrt log_kt
  double alphas_derivative(double log_kt, bool allow_freezing=false, bool force_no_2loops=false) const;

  //========================================================================
  // sudakov frequent bits
  //
  // Note: to obtain fixed coupling (at scale kt), set Lfr=0
  // i.e. mufr=1.0 and the appropriate alphas value in the ConfigBase.
  //========================================================================

  
  /// a triangle bound by an upper theta limit, a kt limit and a z
  /// theta^alpha limit
  ///
  /// For alpha>1, the kt limit is an upper one, the z theta^alpha limit
  /// is a lower one
  /// For alpha<1, the kt limit is a lower one, the z theta^alpha limit
  /// is an upper one
  ///
  /// For alpha=0 and B!=0, a B term is included. It is the user's
  /// responsibility to included it only when needed
  ///
  /// constraints: - alpha != 1
  ///              - Ltop<Lbot
  ///
  /// assumptions: the full triangle is in the kinematic acceptance
  /// [note however that we can trick it if we want to take
  /// differences of two such triangles where the part outside of the
  /// soundary cancels. This would not work without freezing the coupling]
  double kt_triangle(double alpha, double Ltop, double Lbot, double B=0.0) const;

  /// a triangle with 
  ///  - an upper bound set as z\theta^alpha = cst
  ///  - a lower  bound set as z\theta^beta  = cst (beta>alpha)
  ///  - a maximal angle (unspecified since not needed)
  ///  - the log of (1/) the max and min kt (Ltop and Lbot,
  ///    respectively are computed at the maximal angle
  ///  - it is the user's responsibility to see if B needs to be set
  ///    to 0 when not needed
  double triangle(double alpha, double beta, double Ltop, double Lbot, double B=0.0) const;

  
  /// a line with a given slope between 2 different values of kt
  /// (should probably not been used alone, use line below insstead)
  double delta_kt_line(double alpha, double Ltop, double Lbot, double B=0.0, bool force_no_2loops=false) const;

  /// a line, given by the coordinates of its endpoints. Compared to the
  /// above, this also supports a constant kt line and should be
  /// preferred.
  ///
  /// The line starts from (theta0, kt0) to (theta1, kt1)
  /// Note that this is theta, not theta^2
  ///
  /// Note also that there is no check performed if z=1 is included,
  /// i.e. the B term is included blindly and it it the end-user's
  /// responsibility to set B to 0 when needed.
  double line(double Ltheta0, double Lkt0, 
              double Ltheta1, double Lkt1,
              double B=0.0, bool force_no_2loops=false) const;

  /// the R' corresponding to the bottom-line of the triangle
  double triangle_prime(double alpha, double beta,
                        double Ltop, double Lbot,
                        double B=0.0, bool force_no_2loops=false) const;

  /// the R" corresponding to the bottom-line of the triangle
  ///currently only for B=0 and beta!=1
  double triangle_double_prime(double alpha, double beta,
                               double Ltop, double Lbot,
                               double B=0.0, bool force_no_2loops=false) const;


  //========================================================================
  // expansions at LO
  //========================================================================
  double kt_triangle_LO(double alpha, double Ltop, double Lbot,
                        double B=0.0) const;
  double triangle_LO(double alpha, double beta,
                     double Ltop, double Lbot, double B=0.0) const;
  double delta_kt_line_LO(double alpha, double Ltop, double Lbot,
                          double B=0.0, bool force_no_2loops=false) const;
  double line_LO(double Ltheta0, double Lkt0, 
                 double Ltheta1, double Lkt1,
                 double B=0.0, bool force_no_2loops=false) const;
  double triangle_prime_LO(double alpha, double beta,
                           double Ltop, double Lbot,
                           double B=0.0, bool force_no_2loops=false) const;
  double triangle_double_prime_LO(double alpha, double beta,
                                  double Ltop, double Lbot,
                                  double B=0.0, bool force_no_2loops=false) const;

  //========================================================================
  // expansions at NLO
  //========================================================================
  double kt_triangle_NLO(double alpha, double Ltop, double Lbot,
                         double B=0.0) const;
  double triangle_NLO(double alpha, double beta,
                      double Ltop, double Lbot, double B=0.0) const;
  double delta_kt_line_NLO(double alpha, double Ltop, double Lbot,
                           double B=0.0, bool force_no_2loops=false) const;
  double line_NLO(double Ltheta0, double Lkt0, 
                  double Ltheta1, double Lkt1,
                  double B=0.0, bool force_no_2loops=false) const;
  double triangle_prime_NLO(double alpha, double beta,
                            double Ltop, double Lbot,
                            double B=0.0, bool force_no_2loops=false) const;
  double triangle_double_prime_NLO(double alpha, double beta,
                                   double Ltop, double Lbot,
                                   double B=0.0, bool force_no_2loops=false) const;

  //========================================================================
  // expansions at NNLO
  //========================================================================
  double kt_triangle_NNLO(double alpha, double Ltop, double Lbot,
                          double B=0.0) const;
  double triangle_NNLO(double alpha, double beta,
                       double Ltop, double Lbot, double B=0.0) const;
  double delta_kt_line_NNLO(double alpha, double Ltop, double Lbot,
                            double B=0.0, bool force_no_2loops=false) const;
  double line_NNLO(double Ltheta0, double Lkt0, 
                   double Ltheta1, double Lkt1,
                   double B=0.0, bool force_no_2loops=false) const;
  double triangle_prime_NNLO(double alpha, double beta,
                             double Ltop, double Lbot,
                             double B=0.0, bool force_no_2loops=false) const;
  double triangle_double_prime_NNLO(double alpha, double beta,
                                    double Ltop, double Lbot,
                                    double B=0.0, bool force_no_2loops=false) const;
  
};

#endif
