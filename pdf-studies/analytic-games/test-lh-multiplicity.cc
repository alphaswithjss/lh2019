// multiplicity distribution for ISD
//
// currently just for -ve beta

#include <iostream>
#include <CmdLine.hh>
#include "lh-multiplicity.hh"

using namespace std;

// main code
int main(int argc, char *argv[]){

  //--------------------------------------------------
  // parse the command line
  CmdLine cmd(argc, argv);

  double ktcut = cmd.value("-ktcut", 1.0);
  double zcut  = cmd.value("-zcut", -1.0);
  double beta  = cmd.value("-beta",  0.0);
  
  double alphas_MZ = cmd.value("-alphasMZ", 0.1265);
  
  double pt = cmd.value("-pt", 500.0);
  double R  = cmd.value("-R", 0.5);

  bool two_loops  = cmd.present("-two-loops");
  
  double muR = cmd.value("-muR", 1.0);
  double muQ = cmd.value("-muQ", 1.0);

  // use distribution instead of cumulative
  bool cumul = cmd.present("-cumul");

  cmd.assert_all_options_used();
  
  cout << "# Ran: " << cmd.command_line() << endl;

  //--------------------------------------------------
  // build objects needed
  LHMultiplicity lhm(ktcut, R, zcut, beta, alphas_MZ, two_loops, muR, muQ);

  if (cumul){
    // get the cumulative probability
    cout << "# n P_quark(<=n) P_gluon(<=n)" << endl;
    for (double n=0; n<20;n+=0.2)
      cout << n << " "
           << lhm.fraction_below(n, pt, false) << " "
           << lhm.fraction_below(n, pt, true)  << endl;
  } else {
    // get the probability
    cout << "# n P_quark(=n) P_gluon(=n)" << endl;
    for (unsigned int n=0; n<20;++n)
      cout << n << " "
           << lhm.Pn(n, pt, false) << " "
           << lhm.Pn(n, pt, true)  << endl;
  }
  
  return 0;
}
