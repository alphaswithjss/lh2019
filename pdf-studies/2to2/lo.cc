#include "lo.hh"
#include "constants.hh"
#include "me.hh"

#include <iostream>
#include <cassert>
#include <cmath>
#include <LHAPDF/LHAPDF.h>

using namespace std;

bool protonA=true;
bool protonB=true;

inline unsigned int A(unsigned int i){return protonA ? i : 12-i;}
inline unsigned int B(unsigned int i){return protonB ? i : 12-i;}

// LO *fortranwrapper;
// 
// extern "C" {
//   void initspectrum_(double* sqrts, double* mur, double* muf, bool* gg2qq, bool* gg2gg,
// 		     bool* qq2gg, bool* qq2qq, bool* qg2qg) {
//     fortranwrapper=new LO(*sqrts,*mur,*muf,*gg2qq,*gg2gg,*qq2gg,*qq2qq,*qg2qg);
//   }
//   double quarkspectrum_(double* pt, double* ymin, double* ymax) {
//     return fortranwrapper->dsigma_dpt(*pt,1,*ymin,*ymax)*(*pt);
//   }
//   double gluonspectrum_(double* pt, double* ymin, double* ymax) {
//     return fortranwrapper->dsigma_dpt(*pt,2,*ymin,*ymax)*(*pt);
//   }
// }



//-------------------------------------
// ctor
//-------------------------------------
LO::LO(double _sqrts, double _mur, double _muf, const std::string &pdfset, unsigned int subset,
       bool gg2qq, bool gg2gg, bool qq2gg, bool qq2qq, bool qg2qg){
  sqrts = _sqrts;

  muf = _muf;
  mur = _mur;
  _gg2qq = gg2qq;
  _gg2gg = gg2gg;
  _qq2gg = qq2gg;
  _qq2qq = qq2qq;
  _qg2qg = qg2qg;

  _q = 1.0;
  _g = 1.0;

  // init the PDF set from LHAPDF
  cout << "LO: load PDF set" << endl;
  LHAPDF::initPDFSet(pdfset);
  //LHAPDF::initPDFSet("cteq6l1.LHgrid");
  //LHAPDF::initPDFSet("/home/frederic/Work/Software/local/share/lhapdf/PDFsets/MSTW2008nlo90cl.LHgrid");
  //LHAPDF::initPDFSet("/home/frederic/Work/Software/local/share/lhapdf/PDFsets/MSTW2008nnlo90cl.LHgrid");
  LHAPDF::initPDF(subset);
  //LHAPDF::initPDFSet("/afs/cern.ch/user/s/soyez/jets/jets-hadronisation/cross-sections/pp_ratio/analytic/cteq66.LHgrid");
  //cout << "  alphas(MZ) = " << LHAPDF::alphasPDF(91.2) << endl;
  //cout << "    done" << endl;

  // init GSL integration workspaces  
  _w1 = gsl_integration_workspace_alloc(1000000);
  _w2 = gsl_integration_workspace_alloc(1000000);
}

//-------------------------------------
// dtor
//-------------------------------------
LO::~LO(){
  gsl_integration_workspace_free(_w1);
  gsl_integration_workspace_free(_w2);
}


//-------------------------------------
// dsigma/(dpt dy)  in nb/GeV
//-------------------------------------
// Note: the factor 16 \pi^2 \alphas^2 is taken out of the amplitudes
// and put in here
double integrand_lo_dsigma_dpt_dy(double Y, void* params){
  const LO & lo = * ((LO*) params);
  const double & pt = lo._pt;
  const double & y  = lo._y;

  double yhat = y-Y;
  double ch = 0.5*(exp(yhat)+exp(-yhat));

  double x1 = 2.0*pt/lo.sqrts*ch*exp(+Y);
  double x2 = 2.0*pt/lo.sqrts*ch*exp(-Y);

  if (x1>1) cerr << "x1 = " << x1 << endl;
  if (x2>1) cerr << "x2 = " << x2 << endl;

  double t = -pt*pt*(1+exp(-2*yhat));
  double u = -pt*pt*(1+exp(+2*yhat));
  double s = -t-u;

  double xpdf1[13]; LHAPDF::xfx(x1, lo.muf * pt, xpdf1);
  double xpdf2[13]; LHAPDF::xfx(x2, lo.muf * pt, xpdf2);

  double sum=0.0;

  // qq' (qqbar', qbar q' and qbarqbar') contrib
  double pdf=0.0;
  for (unsigned int i=0;i<6;i++){
    double q1 = xpdf1[i]+xpdf1[12-i];
    for (unsigned int j=0;j<6;j++){
      if (i==j) continue;
      pdf += q1*(xpdf2[j]+xpdf2[12-j]);
    }
  }
  if (lo._qq2qq)
    sum += lo._q * pdf * (Mqqp2qqp(s,t,u) + Mqqp2qqp(s,u,t));

  // qq (and qbarqbar) contrib
  pdf=0.0;
  for (unsigned int i=0;i<6;i++)
    pdf += xpdf1[A(i)]*xpdf2[B(i)] + xpdf1[A(12-i)]*xpdf2[B(12-i)];
  if (lo._qq2qq)
    sum += lo._q * pdf * 2.0 * Mqq2qq(s,t,u);
 

  // qqbar (and qbarq) contrib
  pdf=0.0;
  for (unsigned int i=0;i<6;i++)
    pdf += xpdf1[A(i)]*xpdf2[B(12-i)] + xpdf1[A(12-i)]*xpdf2[B(i)];
  if (lo._qq2qq)
    sum += pdf * lo._q * (Mqqbar2qpqbarp(s,t,u) + Mqqbar2qpqbarp(s,u,t) + Mqqbar2qqbar(s,t,u) + Mqqbar2qqbar(s,u,t));
  if (lo._qq2gg)
    sum += pdf * lo._g * 2 * Mqqbar2gg(s,t,u);

  // // qg (gq, qbqrg and gqbqr) contrib
  // pdf=0.0;
  // for (unsigned int i=0;i<6;i++)
  //   pdf += (xpdf1[i]+xpdf1[12-i])*xpdf2[6] + xpdf1[6]*(xpdf2[12-i]+xpdf2[i]);
  // if (lo._qg2qg)
  //   sum += 0.5*(lo._q+lo._g) * pdf * (Mqg2qg(s,t,u) + Mqg2qg(s,u,t));

  // qg (gq, qbqrg and gqbqr) contrib
  pdf=0.0;
  for (unsigned int i=0;i<6;i++)
    pdf += (xpdf1[i]+xpdf1[12-i])*xpdf2[6];
  if (lo._qg2qg)
    sum += pdf * (lo._q*Mqg2qg(s,t,u) + lo._g*Mqg2qg(s,u,t));
    // sum += 0.5*(lo._q+lo._g) * pdf * (Mqg2qg(s,t,u) + Mqg2qg(s,u,t));
  pdf=0.0;
  for (unsigned int i=0;i<6;i++)
    pdf += xpdf1[6]*(xpdf2[12-i]+xpdf2[i]);
  if (lo._qg2qg)
    sum += pdf * (lo._q*Mqg2qg(s,u,t) + lo._g*Mqg2qg(s,t,u));
    // sum += 0.5*(lo._q+lo._g) * pdf * (Mqg2qg(s,t,u) + Mqg2qg(s,u,t));
    
  // gg contrib
  pdf = xpdf1[6] * xpdf2[6];
  if (lo._gg2gg)
    sum += pdf*lo._g*2.0*Mgg2gg(s,t,u);
  if (lo._gg2qq)
    sum += pdf*lo._q*(Mgg2qqbar(s,t,u) + Mgg2qqbar(s,u,t));
 
  return sum/ch/ch/ch/ch;
}


double LO::dsigma_dpt_dy(double pt, double y, unsigned int parton_type){
  // impose boundaries on y
  double tmp = sqrts/2.0/pt;
  double maxrap = log(tmp + sqrt(tmp*tmp-1));
  if (abs(y) > maxrap){
    cerr << "Going above the rapidity limit" << endl;
    return 0.0;
  }


  //cout << "  LO with pt = " << pt << endl;

  _pt = pt;
  _y  = y;

  _q = ((parton_type == PARTON_BOTH) || (parton_type == PARTON_QUARK));
  _g = ((parton_type == PARTON_BOTH) || (parton_type == PARTON_GLUON));

  gsl_function f;
  f.function = &integrand_lo_dsigma_dpt_dy;
  f.params = (void*) this;

  double Ymin = y-0.5*log(sqrts*exp(+y)/pt-1);
  double Ymax = y+0.5*log(sqrts*exp(-y)/pt-1);

  //cerr << y << " " << pt << " " << Ymin << " " << Ymax << endl;
  assert(Ymin<=Ymax);

  double res, err;
  //double pts[2] = {Ymin, Ymax};
  gsl_integration_qag(&f, Ymin, Ymax, 0.0, 1e-3, 1000000, GSL_INTEG_GAUSS15, _w1, &res, &err);

  double scale = mur * pt;
  double alphas = LHAPDF::alphasPDF(scale);
  return M_PI/4.0*alphas*alphas*res/pt/pt/pt * GeV2nb;
}


//-------------------------------------
// dsigma/dpt
//-------------------------------------
double integrand_lo_dsigma_dpt(double y, void* params){
  //cout << "   y= " << y << endl;
  LO * lo = ((LO*) params);
  double res = lo->dsigma_dpt_dy(lo->_pt, y, lo->_parton_type);
  //cout << "   done" << endl;
  return res;
}

// use ymax<ymin to integrate over the complete domain
double LO::dsigma_dpt(double pt, unsigned int parton_type, double ymin, double ymax){
  if (pt>=0.99*0.5*sqrts) return 0.0;

  //cout << "pt is " << pt << endl;

  _pt = pt;
  _parton_type = parton_type;

  gsl_function f;
  f.function = &integrand_lo_dsigma_dpt;
  f.params = this;

  double tmp = sqrts/2.0/pt;
  double maxrap = log(tmp + sqrt(tmp*tmp-1));
  double lower = ((ymax<ymin) || (ymin<-maxrap)) ? -maxrap : ymin;
  double upper = ((ymax<ymin) || (ymax>+maxrap)) ? +maxrap : ymax;
  double res, err;
  gsl_integration_qag(&f, lower, upper, 0.0, 1e-3, 
		      1000000, GSL_INTEG_GAUSS15, _w2, &res, &err);  

  return res;
}
