#include "lo.hh"
#include "constants.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <LHAPDF/LHAPDF.h> 
#include <CmdLine.hh>

//////////
// Simple program to write out the inclusive quark and gluon jet spectrum
// obtained from the LO code in the file ../parton-spectrum.dat
// 
// The output file can then be read by inclusive_spectrum.f90,
// using it as input spectrum. For example:
// ./inclusive_spectrum -out output/lo-code-R01.dat -R 0.1 -in spectrum-lo-code.dat
/////////

using namespace std;

int main(int argc, char *argv[]){
  CmdLine cmd(argc, argv);

  double sqrts = cmd.value("-rts", 13000.0);
  double mur   = cmd.value("-mur", 1.0);
  double muf   = cmd.value("-muf", 1.0);

  double ptmin = cmd.value("-ptmin",100.0);
  double ptmax = cmd.value("-ptmax",0.5*sqrts);
  unsigned int npt = cmd.value<unsigned int>("-npt", 100);

  double ymin=cmd.value("-ymin", -0.5);
  double ymax=cmd.value("-ymax",  0.5);

  bool gg2qq = cmd.present("-gg2qq");
  bool gg2gg = cmd.present("-gg2gg");
  bool qq2gg = cmd.present("-qq2gg");
  bool qq2qq = cmd.present("-qq2qq");
  bool qg2qg = cmd.present("-qg2qg");

  string out = cmd.value<string>("-out", "");
  shared_ptr<ofstream> ostr_ptr;
  ostream *ostr;
  if (out.empty()){
    ostr = &cout;
  } else {
    ostr_ptr.reset(new ofstream(out.c_str()));
    ostr = ostr_ptr.get();
  }
  
  LO lo(sqrts, mur, muf, gg2qq, gg2gg, qq2gg, qq2qq, qg2qg);

  (*ostr) << "# ran: " << cmd.command_line() << endl;
  (*ostr) << "# pt  gluon_jet  quark_jet" << endl;
  for (unsigned int ipt=0; ipt<npt; ++ipt){
    double f = ipt/(1.0*npt);
    double pt = exp((1-f)*log(ptmin) + f*log(ptmax));
    double sigma_q, sigma_g;
    sigma_g=lo.dsigma_dpt(pt, PARTON_GLUON, ymin, ymax);
    sigma_q=lo.dsigma_dpt(pt, PARTON_QUARK, ymin, ymax);
    
    (*ostr) << setw(15) << pt << setw(15) << sigma_g << setw(15) << sigma_q << endl;
  }

  return 0;
}
