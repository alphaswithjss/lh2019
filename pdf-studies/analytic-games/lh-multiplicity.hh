#ifndef __LH_MULTIPLICITY_HH__
#define  __LH_MULTIPLICITY_HH__

#include "resum-blocks.hh"
#include <cmath>

/// \class LHMultiplicity
/// LL implementation of the Les-Houches multiplicity distribution
///
/// This is a rough implementation which is correct at LL. It (can)
/// include(s) some of the NLL corrections to the density but no NLL
/// corrections to the Poisson distribution itself
class LHMultiplicity{
public:
  /// ctor
  ///  \param ktcut      kt scale down to which the emissions are included
  ///  \param R          jet radius
  ///  \param zcut       optional SD zcut parameter (set to -ve to disable)
  ///  \param beta       optional SD beta parameter (set to -ve to disable)
  ///  \param alphas_MZ  strong coupling at the Z mass
  ///  \param two_loops  include 2-loop corrections in the emission probability
  ///  \param muR        choice of renormalisation scale (factor of ptR)
  ///  \param muQ        choice of resummation scale (factor of ptR)
  LHMultiplicity(double ktcut, double R,
                 double zcut, double beta,
                 double alphas_MZ,
                 bool two_loops,
                 double muR=1.0, double muQ=1.0)
    : _ktcut(ktcut), _R(R),
      _zcut(zcut), _beta(beta),
      _alphas_MZ(alphas_MZ),
      _two_loops(two_loops),
      _muR(muR), _muQ(muQ){}

  /// probability to get a multipliccity of n
  double Pn(unsigned int n, double pt, bool gluon){
    double nu = _get_density(pt, gluon);
    
    // get the probability
    double res=1.0;
    for (unsigned int i=0; i<n; ++i){
      res *= nu/(i+1);
    }
    return res*exp(-nu);
  }

  /// probability to get a multipliccity of n
  double fraction_below(double n, double pt, bool gluon){
    double nu = _get_density(pt, gluon);
    
    // get the probability
    double loc=1.0, sum=0.0;
    unsigned int i=0;
    while (i<=n){
      sum += loc;
      loc *= nu/(i+1);
      ++i; 
    }
    return (sum+(n-(i-1))*loc)*exp(-nu);
  }

protected:
  /// get the emission density
  double _get_density(double pt, bool gluon){
    double alphas = alphas_rg(_muR*pt*_R, MZ, _alphas_MZ, _two_loops);
    ConfigBase cfg((gluon ? CA : CF), alphas, 1.0/(pt*_R), _muR);
    cfg.set_two_loops(_two_loops);
  
    double nu = cfg.kt_triangle(0.0, -cfg.Bi, log(_muQ*pt*_R/_ktcut), 0.0);
    if (_zcut>0)
      nu -= cfg.kt_triangle(-_beta, -log(_zcut), log(_muQ*pt*_R/_ktcut), 0.0);
    
    return nu;
  }

  double _ktcut;  ///< cut on kt for each primary emission
  double _R;      ///< jet radius
  double _zcut;   ///< optional SD cut: zcut parameter
  double _beta;   ///< optional SD cut: beta parameter

  double _alphas_MZ;  ///< strong coupling at the Z mass
  bool _two_loops;    ///< whether to include the 2-loop corrections
  
  double _muR;  ///< renormalisation scale
  double _muQ;  ///< resummation scale
};

#endif
