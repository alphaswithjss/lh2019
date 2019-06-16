// multiplicity distribution for ISD
//
// currently just for -ve beta

#include <iostream>
#include <cassert>
#include <CmdLine.hh>
#include "ecfs.hh"

using namespace std;

// main code
int main(int argc, char *argv[]){

  //--------------------------------------------------
  // parse the command line
  CmdLine cmd(argc, argv);

  double alpha = cmd.value("-alpha", 0.5 );
  double zcut  = cmd.value("-zcut", -1.0);
  double beta  = cmd.value("-beta",  0.0);
  
  double alphas_MZ = cmd.value("-alphasMZ", 0.1265);
  
  double pt = cmd.value("-pt", 500.0);
  double R  = cmd.value("-R", 0.5);

  bool one_loops  = cmd.present("-one-loops");
  
  double muR = cmd.value("-muR", 1.0);
  double muQ = cmd.value("-muQ", 1.0);

  unsigned int eporder = cmd.value<unsigned int>("-eporder", 2);
  
  // use distribution instead of cumulative
  bool cumul = cmd.present("-cumul");
  assert(cumul);

  cmd.assert_all_options_used();
  
  cout << "# Ran: " << cmd.command_line() << endl;

  //--------------------------------------------------
  // build objects needed
  ECF ecf(alpha, R, zcut, beta, alphas_MZ, !one_loops, muR, muQ, eporder);

  double vmax = (eporder==0)
    ? 1.0 : ((eporder==1) ? 0.25 : 0.35);
  
  if (cumul){
    // get the cumulative probability
    cout << "# lambda P_quark(<=lambda) P_gluon(<=lambda)" << endl;
    for (double x=0.01*vmax;x<=vmax*1.0001;x+=0.01*vmax)
      cout << x << " "
           << ecf.fraction_below(x, pt, false) << " "
           << ecf.fraction_below(x, pt, true)  << endl;
  }
  // else {
  //   // get the probability
  //   cout << "# n P_quark(=n) P_gluon(=n)" << endl;
  //   for (unsigned int n=0; n<20;++n)
  //     cout << n << " "
  //          << lhm.Pn(n, pt, false) << " "
  //          << lhm.Pn(n, pt, true)  << endl;
  // }
  
  return 0;
}
