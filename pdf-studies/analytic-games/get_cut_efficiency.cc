// multiplicity distribution for ISD
//
// currently just for -ve beta

#include <iostream>
#include <cassert>
#include <CmdLine.hh>
#include "shape.hh"
#include "lh-multiplicity.hh"
#include "ecfs.hh"
#include <memory>

using namespace std;

// main code
int main(int argc, char *argv[]){

  //--------------------------------------------------
  // parse the command line
  CmdLine cmd(argc, argv);

  double vcut = cmd.value("-vcut", 0.1);
  
  double ktcut = cmd.value("-ktcut", 1.0);
  double alpha = cmd.value("-alpha", 0.5 );
  double zcut  = cmd.value("-zcut", -1.0);
  double beta  = cmd.value("-beta",  0.0);
  
  double alphas_MZ = cmd.value("-alphasMZ", 0.1265);
  
  double pt = cmd.value("-pt", 500.0);
  double R  = cmd.value("-R", 0.5);

  bool one_loop = cmd.present("-one-loop");
  
  double muR = cmd.value("-muR", 1.0);
  double muQ = cmd.value("-muQ", 1.0);

  bool gluon = cmd.present("-gluon");

  unsigned int eporder = cmd.value<unsigned int>("-eporder", 2);

  shared_ptr<Shape> v;
  if (cmd.present("-ecf")){
    v.reset(new ECF(alpha, R, zcut, beta, alphas_MZ, !one_loop, muR, muQ, eporder));
  } else if (cmd.present("-nlh")){
    v.reset(new LHMultiplicity(ktcut, R, zcut, beta, alphas_MZ, !one_loop, muR, muQ));
  } else {
    cerr << "A shape (-ecf or -nlh) should be specified" << endl;
    exit(1);
  }

  // use distribution instead of cumulative
  cmd.assert_all_options_used();
  
  cout << 1-v->fraction_below(vcut, pt, gluon) << endl;
  
  return 0;
}
