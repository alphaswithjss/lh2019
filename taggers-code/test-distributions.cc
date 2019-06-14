#include "AnalysisFramework.hh"
#include <fastjet/contrib/ModifiedMassDropTagger.hh>
#include <fastjet/contrib/SoftDrop.hh>
#include "qg-taggers.hh"

using namespace std;
using namespace fastjet;

class TaggerWithName{
public:
  TaggerWithName(FunctionOfPseudoJet<double> *tagger_in, const string &name_in, bool use_log_in = false)
    : tagger(tagger_in), name(name_in), use_log(use_log_in){}

  shared_ptr<FunctionOfPseudoJet<double> > tagger;
  string name;
  bool use_log;
};

/// Example class derived from AnalysisFramework that will help in evaluation of
/// some basic cross sections. 
class XSctAnalysis : public AnalysisFramework {
public:
  XSctAnalysis(CmdLine & cmdline) : AnalysisFramework(cmdline) {}
  void user_startup() {
    // extra user parameters
    param["jet.ptmin"]  = 500.0;
    param["jet.rapmax"] =   4.0;
    param["R"]          =   0.5;

    param["mmdt.zcut"] = 0.1;
    param["sd.zcut"] = 0.05;
    param["sd.beta"] = 2.0;
    
    jet_def = JetDefinition(antikt_algorithm, cmdline.value("-R",param["R"]));

    DefaultHist::set_defaults(0.0, 15.0, 0.1);
  }

  void user_post_startup() {
    _sel_no_neutrino = !(SelectorAbsPDGId(12) || SelectorAbsPDGId(14) || SelectorAbsPDGId(16));
    _sel_jets = SelectorNHardest(2) * SelectorAbsRapMax(param["jet.rapmax"]) * SelectorPtMin(param["jet.ptmin"]);

    _mmdt.reset(new contrib::ModifiedMassDropTagger(param["mmdt.zcut"]));
    _sd.reset(new contrib::SoftDrop(param["sd.beta"], param["sd.zcut"], param["R"]));

    // angularities
    _taggers.push_back(TaggerWithName(new Angularity(0.5, param["R"]), "plain_ang_0.5", true));
    _taggers.push_back(TaggerWithName(new Angularity(1.0, param["R"]), "plain_ang_1.0", true));
    _taggers.push_back(TaggerWithName(new Angularity(2.0, param["R"]), "plain_ang_2.0", true));
    _taggers.push_back(TaggerWithName(new Angularity(0.5, param["R"], _sd.get()), "sd_ang_0.5", true));
    _taggers.push_back(TaggerWithName(new Angularity(1.0, param["R"], _sd.get()), "sd_ang_1.0", true));
    _taggers.push_back(TaggerWithName(new Angularity(2.0, param["R"], _sd.get()), "sd_ang_2.0", true));
    _taggers.push_back(TaggerWithName(new Angularity(0.5, param["R"], _mmdt.get()), "mmdt_ang_0.5", true));
    _taggers.push_back(TaggerWithName(new Angularity(1.0, param["R"], _mmdt.get()), "mmdt_ang_1.0", true));
    _taggers.push_back(TaggerWithName(new Angularity(2.0, param["R"], _mmdt.get()), "mmdt_ang_2.0", true));

    // ECFs
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(0.5, param["R"]), "plain_ecf_0.5", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(1.0, param["R"]), "plain_ecf_1.0", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(2.0, param["R"]), "plain_ecf_2.0", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(0.5, param["R"], _sd.get()), "sd_ecf_0.5", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(1.0, param["R"], _sd.get()), "sd_ecf_1.0", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(2.0, param["R"], _sd.get()), "sd_ecf_2.0", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(0.5, param["R"], _mmdt.get()), "mmdt_ecf_0.5", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(1.0, param["R"], _mmdt.get()), "mmdt_ecf_1.0", true));
    _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(2.0, param["R"], _mmdt.get()), "mmdt_ecf_2.0", true));

    // ISD muly
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(1.0, param["R"]), "plain_isd_1.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(2.0, param["R"]), "plain_isd_2.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(5.0, param["R"]), "plain_isd_5.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(1.0, param["R"], param["sd.zcut"], param["sd.beta"]), "sd_isd_1.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(2.0, param["R"], param["sd.zcut"], param["sd.beta"]), "sd_isd_2.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(5.0, param["R"], param["sd.zcut"], param["sd.beta"]), "sd_isd_5.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(1.0, param["R"], param["sd.zcut"]), "mmdt_isd_1.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(2.0, param["R"], param["sd.zcut"]), "mmdt_isd_2.0"));
    _taggers.push_back(TaggerWithName(new IteratedSoftDropKtMultiplicity(5.0, param["R"], param["sd.zcut"]), "mmdt_isd_5.0"));

    declare_hists("plain");
    declare_hists("mmdt");
    declare_hists("sd");
  }

  void declare_hists(const string &tag){
    // angularities
    norm_hists[tag+"_ang_0.5"].declare(-15.0, 0.0, 0.025);
    norm_hists[tag+"_ang_1.0"].declare(-15.0, 0.0, 0.025);
    norm_hists[tag+"_ang_2.0"].declare(-15.0, 0.0, 0.025);

    // ECFs
    norm_hists[tag+"_ecf_0.5"].declare(-15.0, 0.0, 0.025);
    norm_hists[tag+"_ecf_1.0"].declare(-15.0, 0.0, 0.025);
    norm_hists[tag+"_ecf_2.0"].declare(-15.0, 0.0, 0.025);

    // ISD multiplicities
    norm_hists[tag+"_isd_1.0"].declare(-0.5, 40.5, 1.0);
    norm_hists[tag+"_isd_2.0"].declare(-0.5, 40.5, 1.0);
    norm_hists[tag+"_isd_5.0"].declare(-0.5, 40.5, 1.0);
    
  }

  void analyse_event() {
    double evwgt = driver->generator->hadron_level().weight();
    xsections["total_cross_section"] += evwgt;

    vector<PseudoJet> jets = _sel_jets(jet_def(_sel_no_neutrino(driver->generator->hadron_level().particles())));
    
    for (const auto &jet : jets){
      xsections["jet_cross_section"] += evwgt;

      for (const auto &tagger : _taggers){
        double res = (*(tagger.tagger))(jet);
        if (tagger.use_log) res = (res>0) ? log(res) : -14.999;
        norm_hists[tagger.name].add_entry(res);
      }
    }
  }


protected:
  shared_ptr<contrib::ModifiedMassDropTagger> _mmdt;
  shared_ptr<contrib::SoftDrop> _sd;
  vector<TaggerWithName> _taggers;
  Selector _sel_no_neutrino, _sel_jets;
};

//----------------------------------------------------------------------
int main (int argc, char ** argv) {
  
  CmdLine cmdline(argc,argv);
  XSctAnalysis analysis(cmdline);
  analysis.run();
}
