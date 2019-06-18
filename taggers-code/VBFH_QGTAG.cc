// -*- C++ -*-
// author: y. Haddad (yhaddad@cern.ch)

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "qg-taggers.hh"

#include <iostream>
#include <fstream>
#include <ctime>
#include <locale>


using namespace fastjet;
using namespace std;

class TaggerWithName{
public:
  TaggerWithName(FunctionOfPseudoJet<double> *tagger_in, const string &name_in, bool use_log_in = false)
    : tagger(tagger_in), name(name_in), use_log(use_log_in){}

  shared_ptr<FunctionOfPseudoJet<double> > tagger;
  string name;
  bool use_log;
};

namespace Rivet {
  /// @brief Add a short analysis description here
  class VBFH_QGTAG : public Analysis {
  public:

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(VBFH_QGTAG);

  public:

    /// @brief Whether particle p originate from any of the ptcls
    bool originateFrom(const Particle& p, const Particles& ptcls ) {
      const GenVertex* prodVtx = p.genParticle()->production_vertex();
      if (prodVtx == nullptr) return false;
      // for each ancestor, check if it matches any of the input particles
      for (auto ancestor:particles(prodVtx, HepMC::ancestors)){
	for ( auto part:ptcls )
	  if ( ancestor==part.genParticle() ) return true;
      }
      // if we get here, no ancetor matched any input particle
      return false;
    }
    /// @brief Return flavour of jet given collection of partons, -1 means no maching found
    int jetFlavour(const Jet &jet, const Particles& partons){
      double delta_r = 99999;
      int tag = -1;
      for(const Particle &p: partons){
	double dr = deltaR(jet, p);
	if ((dr < delta_r) && dr < 0.4){
	  tag = p.pid();
	  delta_r = dr;
	}
      }
      return tag;
    }
    /// @brief Return true is particle decays to quarks
    bool quarkDecay(const Particle &p) {
      for (auto child:p.children())
        if (PID::isQuark(child.pdgId())) return true;
      return false;
    }

    /// @brief Sets the Higgs production mode
    // Book histograms and initialise projections before the run
    void init() {
      std::cout << "==============================================================" << std::endl;
      std::cout << "========          QG TAG VBF ROUTINE                 =========" << std::endl;
      std::cout << "==============================================================" << std::endl;

      // Projections for final state particles
      const FinalState FS;//(-5.0, 5.0, 0.0*GeV);
      addProjection(FS,"FS");
      m_sumw = 0.0;

      // Projection for dressed electrons and muons
      IdentifiedFinalState higgs(FS);
      higgs.acceptIdPair(PID::HIGGS);
      declare(higgs,"Higgs");

      // Jets, anti-kt 0.4
      VetoedFinalState fsJets; //final state for jet finding: veto Higgs
      fsJets.addVetoPairId(PID::HIGGS);
      fsJets.addVetoPairId(6);
      fsJets.addVetoPairId(11);
      fsJets.addVetoPairId(12);
      fsJets.addVetoPairId(13);
      fsJets.addVetoPairId(14);
      fsJets.addVetoPairId(15);
      fsJets.addVetoPairId(16);

      declare(FastJets(fsJets, FastJets::ANTIKT, 0.4), "AK4JETS");
      // declare qg taggers
      _taggers.push_back(TaggerWithName(new Angularity(0.5, 0.4), "plain_ang_0.5", true));
      _taggers.push_back(TaggerWithName(new Angularity(1.0, 0.4), "plain_ang_1.0", true));
      _taggers.push_back(TaggerWithName(new Angularity(2.0, 0.4), "plain_ang_2.0", true));

      _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(0.5, 0.4), "plain_ecf_0.5", true));
      _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(1.0, 0.4), "plain_ecf_1.0", true));
      _taggers.push_back(TaggerWithName(new EnergyCorrelationFunction(2.0, 0.4), "plain_ecf_2.0", true));

      // creating a CSV file
      csv_file_name = "./QGTAG_variables.csv" ;
      csv_file.open ( csv_file_name );
      csv_file << "njets"   << ","
	       << "mjj"     << ","
	       << "detajj"  << ","
               << "dphijj"  << ","
	       << "ptj1"    << ","
	       << "ptj2"    << ","
	       << "ptj3"    << ","
	       << "etaj1"   << ","
	       << "etaj2"   << ","
	       << "etaj3"   << ","
	       << "lead_plain_ang_0.5"   << ","
	       << "subl_plain_ang_0.5"   << ","
	       << "lead_plain_ang_1.0"   << ","
	       << "subl_plain_ang_1.0"   << ","
	       << "lead_plain_ang_2.0"   << ","
	       << "subl_plain_ang_2.0"   << ","
	       << "lead_plain_ecf_0.5"   << ","
	       << "subl_plain_ecf_0.5"   << ","
	       << "lead_plain_ecf_1.0"   << ","
	       << "subl_plain_ecf_1.0"   << ","
	       << "lead_plain_ecf_2.0"   << ","
	       << "subl_plain_ecf_2.0"   << ","
	       << "Npart1"    << ","
	       << "Npart2"    << ","
	       << "Npart3"    << ","
	       << "tag1"    << ","
	       << "tag2"    << ","
	       << "tag3"    << ","
	       << "weight";
      csv_file << std::endl;
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      //const GenEvent* weightevent = event.genEvent();
      const double central_weight = event.weight();

      // Looping on jets above 30 GeV
      const Jets jets = applyProjection<FastJets>(event, "AK4JETS").jetsByPt(30.0);
      if(jets.size() < 2) vetoEvent;

      double mjj =  (jets[0].momentum() + jets[1].momentum()).mass();
      double etajj = deltaEta(jets[0],jets[1]);
      double phijj = deltaPhi(jets[0],jets[1]);

      // checking the parton level particles
      Particles partons;
      for (auto ptcl: Rivet::particles(event.genEvent())) {
      	Particle part =  Particle(ptcl);
	if ( Rivet::PID::isQuark  ( part.pid() )  && ptcl->status() == 23 ){
	  //std::cout << "Quark ? : " << ptcl->pdg_id() << " : " << ptcl->status() << " is parton : " << PID::isParton(part.pid()) << std::endl;
	  partons.push_back(part);
	}
	if ( (part.pid() == PID::GLUON)  && ptcl->status() == 23 ){
	  //std::cout << "Gluon ? : " << ptcl->pdg_id() << " : " << ptcl->status() << " is parton : " << PID::isParton(part.pid()) << std::endl;
	  partons.push_back(part);
	}
      }

      // do matching to partons
      int tag1 =  jetFlavour(jets[0], partons);
      int tag2 =  jetFlavour(jets[1], partons);
      int tag3 =  (jets.size() > 2 ? jetFlavour(jets[2], partons) : -1.0);

      // ----------
      m_sumw += central_weight;

      // fill the row in the csv file
      csv_file << jets.size()    << ","
	       << mjj            << ","
	       << etajj          << ","
	       << phijj          << ","
	       << jets[0].pt()   << ","
	       << jets[1].pt()   << ","
	       << (jets.size() > 2 ? jets[2].pt() : -1.0)  << ","
	       << jets[0].eta()  << ","
	       << jets[1].eta()  << ","
	       << (jets.size() > 2 ? jets[2].eta() : -99.0) << ",";

      for (const auto &tagger : _taggers){
	double lead_res = (*(tagger.tagger))(jets[0].pseudojet());
        double subl_res = (*(tagger.tagger))(jets[1].pseudojet());
        if (tagger.use_log) lead_res = (lead_res>0) ? log(lead_res) : -14.999;
        if (tagger.use_log) subl_res = (subl_res>0) ? log(subl_res) : -14.999;
	std::cout << "values : "<< tagger.name  << " : " << lead_res  << " : " << subl_res << std::endl;
        csv_file << lead_res   << ","
		 << subl_res   << ",";
      }
      csv_file << jets[0].size()<< ","
	       << jets[1].size()<< ","
	       << (jets.size() > 2 ? jets[2].size() : -1.0) << ","
	       << tag1    << ","
	       << tag2    << ","
	       << tag3    << ","
	       << central_weight * 1e24; // the factor 1e24 because pythia8 output weird numbers, to be normalised later
      csv_file << std::endl;
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      csv_file.close();
      std::cout << " crossSection() : " << crossSection() << std::endl;
      std::cout << " sumOfWeights() : " << sumOfWeights() << std::endl;
      std::cout << " m_sumw         : " << m_sumw         << std::endl;
      std::cout << " numEvents ()   : " << numEvents () << std::endl;
    }

    const std::string currentDateTime() {
      std::locale::global(std::locale("ja_JP.utf8"));
      std::time_t t = std::time(nullptr);

      char mbstr[100];
      std::strftime(mbstr, sizeof(mbstr), "%Y-%m-%d %H:%M", std::localtime(&t));
      return mbstr;
    }

    void insertCSVMetaData(){
      // ofstream

      std::fstream processed_csv(csv_file_name.c_str());
      std::stringstream csv_metadata;
      csv_metadata << "# user:"              << getenv("USER")   << std::endl;
      csv_metadata << "# date:"              << currentDateTime()<< std::endl;
      csv_metadata << "# crossSection [pb]:" << crossSection()   << std::endl;
      csv_metadata << "# sumOfWeights:"      << sumOfWeights()   << std::endl;
      csv_metadata << "# scale_factor:"      << crossSection()/sumOfWeights() << std::endl;

      csv_metadata << processed_csv.rdbuf();
      processed_csv.close();

      processed_csv.open(csv_file_name.c_str(), std::fstream::out | std::fstream::trunc);
      processed_csv << csv_metadata.rdbuf();

    }
    /// Define histograms
    void initializeHistos(){}


  private:
    double m_sumw;
    vector<TaggerWithName> _taggers;
    std::string csv_file_name;
    std::ofstream csv_file;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(VBFH_QGTAG);


}
