#ifndef __LO_HH__
#define __LO_HH__

#include <gsl/gsl_integration.h>

const unsigned int PARTON_BOTH  = 0;
const unsigned int PARTON_QUARK = 1;
const unsigned int PARTON_GLUON = 2;

/*
 * compute the LO jet cross-section
 */
class LO{
public:
  /// ctor
  LO(double _sqrts, double _mur=1.0, double _muf=1.0, bool gg2qq=true, bool gg2gg=true,
     bool qq2gg=true, bool qq2qq=true, bool qg2qg=true);

  /// dtor
  ~LO();

  // differential x-sect in pt and y [nb/GeV]
  double dsigma_dpt_dy(double pt, double y, unsigned int parton_type = PARTON_BOTH);

  // differential x-sect in pt [nb/GeV]
  // Note: use ymax<ymin to integrate over the complete domain
  double dsigma_dpt(double pt, unsigned int parton_type = PARTON_BOTH, double ymin=0.1, double ymax=-0.1);

  // internal vars
  double sqrts;                    ///< the com energy
  
  double mur, muf;                 ///< renorm and factorisation scales
  bool _gg2qq, _gg2gg, _qq2gg, _qq2qq, _qg2qg;///< determine which channels are used

  gsl_integration_workspace *_w1;  ///< GSL integration workspace for the Y integration
  gsl_integration_workspace *_w2;  ///< GSL integration workspace for the y integration
  double _pt, _y;                  ///< for transition to the "static" GSL integrands

  double _q, _g;
  unsigned int _parton_type;
};

#endif  // __LO_HH__
