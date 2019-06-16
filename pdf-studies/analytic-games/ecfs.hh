#ifndef __ECFS_HH__
#define __ECFS_HH__

/// \class ECF
/// NLL implementation of angularity distributions
///
/// This is a rough implementation which is correct at LL. It (can)
/// include(s) some of the NLL corrections to the density but no NLL
/// corrections to the Poisson distribution itself
class ECF{
public:
  /// ctor
  ///  \param alpha      exponent of the angles
  ///  \param R          jet radius
  ///  \param zcut       optional SD zcut parameter (set to -ve to disable)
  ///  \param beta       optional SD beta parameter (set to -ve to disable)
  ///  \param alphas_MZ  strong coupling at the Z mass
  ///  \param two_loops  include 2-loop corrections in the emission probability
  ///  \param muR        choice of renormalisation scale (factor of ptR)
  ///  \param muQ        choice of resummation scale (factor of ptR)
  ECF(double alpha, double R,
      double zcut, double beta,
      double alphas_MZ,
      bool two_loops,
      double muR=1.0, double muQ=1.0,
      unsigned int endpoint_order=2, double p=1.0,
      bool transition_point_corrections=false)
    : _alpha(alpha), _R(R),
      _zcut(zcut), _beta(beta),
      _alphas_MZ(alphas_MZ),
      _two_loops(two_loops),
      _muR(muR), _muQ(muQ),
      _endpoint_order(endpoint_order), _p(p),
      _transition_point_corrections(transition_point_corrections){}
  
  double fraction_below(double lambda, double pt, bool gluon);
  double fraction_below_plain(double lambda, double pt, bool gluon);

protected:
  double _alpha;  ///< angular exponent
  double _R;      ///< jet radius
  double _zcut;   ///< optional SD cut: zcut parameter
  double _beta;   ///< optional SD cut: beta parameter

  double _alphas_MZ;  ///< strong coupling at the Z mass
  bool _two_loops;    ///< whether to include the 2-loop corrections
  
  double _muR;  ///< renormalisation scale
  double _muQ;  ///< resummation scale
  unsigned int _endpoint_order;
  double _p;    ///< endpoint approach
  bool _transition_point_corrections; ///< enable transition point corrections
};

#endif


