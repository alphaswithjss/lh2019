//Based on Andy Buckley's Boost 2014 Analysis

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "qg-taggers.hh"

using namespace fastjet;


namespace Rivet {

  using namespace Cuts;

  /// \class QGClassifiedXS
  /// a helper class for computing the inclusive x-section with Q/G tagging.
  /// It holds
  ///  - a discriminating variable (v)
  ///  - the cut on v (v<vcut => q; v>vcut => g)
  ///  - associated pt histograms for the qq, qg, gg channels
  class QGClassifiedXS{
  public:
    /// ctor w proper initialisation
    QGClassifiedXS(FunctionOfPseudoJet<double> *v_ptr, double vcut,
                   Histo1DPtr h_qq, Histo1DPtr h_qg, Histo1DPtr h_gg)
      : _v(v_ptr), _vcut(vcut), _h_qq(h_qq), _h_qg(h_qg), _h_gg(h_gg){}

    /// process an event (given by a set of jets and a weight)
    void process(const vector<PseudoJet> & jets, double weight){
      // classify the event
      Histo1DPtr h;
      if (jets.size()<2){
        // let's say it is qq
        h = _h_qq;
      } else {
        unsigned int nquarks = 0;
        if ((*_v)(jets[0])<_vcut) ++nquarks;
        if ((*_v)(jets[1])<_vcut) ++nquarks;
        h = (nquarks==0) ? _h_gg : ((nquarks==1) ? _h_qg : _h_qq);
      }

      // bin the jets
      for (const auto &jet : jets){
        h->fill(jet.pt(), weight);
      }
    }
    
  protected:
    shared_ptr<FunctionOfPseudoJet<double> > _v;
    double _vcut;
    Histo1DPtr _h_qq, _h_qg, _h_gg;
  };

  /// Standard jet radius used in this analysis (for both kT and anti-kT)

  class MC_LH19_QG4PDF : public Analysis {
  public:

    /// parameters
    const double PARTICLE_RAPMAX; ///< maximal rapidity for particles

    const double JET_RADIUS;      ///< R for jet clustering
    const double JET_PTMIN;       ///< min pt for binning the jets
    const double JET_PTMAX;       ///< max pt for binning the jets
    const unsigned int JET_NPT;   ///< numbber of bins
    const double JET_RAPMAX;      ///< maximal rapidity for the jets

    const double MMDT_ZCUT;       ///< zcut for mMDT (tight grooming)
    const double SD_BETA;         ///< beta for SoftDrop (loose grooming)
    const double SD_ZCUT;         ///< zcut for SoftDrop (loose grooming)

    const double CUT_PLAIN_ECF;
    const double CUT_LOOSE_ECF;
    const double CUT_TIGHT_ECF;
    const double CUT_PLAIN_ISD;
    const double CUT_LOOSE_ISD;
    const double CUT_TIGHT_ISD;
    
    /// Constructor
    MC_LH19_QG4PDF()
      : Analysis("MC_LH19_QG4PDF"),
        PARTICLE_RAPMAX(5.0),
        JET_RADIUS(0.4),
        JET_PTMIN (100.0),
        JET_PTMAX (5000.0),
        JET_NPT   (50),
        JET_RAPMAX(4.5),
        MMDT_ZCUT(0.1),
        SD_BETA  (2.0),
        SD_ZCUT  (0.05),
        CUT_PLAIN_ECF(0.1),
        CUT_LOOSE_ECF(0.1),
        CUT_TIGHT_ECF(0.1),
        CUT_PLAIN_ISD(2.5),
        CUT_LOOSE_ISD(2.5),
        CUT_TIGHT_ISD(2.5)
    {}
    
    /// Book histograms and initialise projections before the run
    void init() {
      
      FinalState fs(-PARTICLE_RAPMAX, PARTICLE_RAPMAX, 0.0*GeV);
      
      // for the jets
      VetoedFinalState jet_input(fs);
      jet_input.vetoNeutrinos();
      jet_input.addVetoPairId(PID::MUON);
      addProjection(jet_input, "JET_INPUT");
      _jet_def = JetDefinition(antikt_algorithm, JET_RADIUS);
      _sel_jets = SelectorAbsRapMax(JET_RAPMAX); //Yes, this is rapidity. Do NOT use pseudo-rapidity here!

      // groomers
      _mmdt.reset(new contrib::ModifiedMassDropTagger(MMDT_ZCUT));
      _sd.reset(new contrib::SoftDrop(SD_BETA, SD_ZCUT, JET_RADIUS));

      // q/g taggers and histogram bookings
      _qg_xs.push_back(QGClassifiedXS(new EnergyCorrelationFunction(0.5), CUT_PLAIN_ECF,
                                      bookHisto1D("plain_ecf_qq",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("plain_ecf_qg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("plain_ecf_gg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX))));
      _qg_xs.push_back(QGClassifiedXS(new EnergyCorrelationFunction(0.5, _sd.get()), CUT_LOOSE_ECF,
                                      bookHisto1D("loose_ecf_qq",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("loose_ecf_qg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("loose_ecf_gg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX))));
      _qg_xs.push_back(QGClassifiedXS(new EnergyCorrelationFunction(0.5, _mmdt.get()), CUT_TIGHT_ECF,
                                      bookHisto1D("tight_ecf_qq",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("tight_ecf_qg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("tight_ecf_gg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX))));

      _qg_xs.push_back(QGClassifiedXS(new IteratedSoftDropKtMultiplicity(2.0), CUT_PLAIN_ISD,
                                      bookHisto1D("plain_isd_qq",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("plain_isd_qg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("plain_isd_gg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX))));
      _qg_xs.push_back(QGClassifiedXS(new IteratedSoftDropKtMultiplicity(2.0, SD_ZCUT, SD_BETA), CUT_LOOSE_ISD,
                                      bookHisto1D("loose_isd_qq",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("loose_isd_qg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("loose_isd_gg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX))));
      _qg_xs.push_back(QGClassifiedXS(new IteratedSoftDropKtMultiplicity(2.0, MMDT_ZCUT), CUT_TIGHT_ISD,
                                      bookHisto1D("tight_isd_qq",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("tight_isd_qg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)),
                                      bookHisto1D("tight_isd_gg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX))));
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {
      // a few shortcuts
      const double weight = e.weight();

      // get the jet inputs
      const VetoedFinalState &fs = applyProjection<VetoedFinalState>(e, "JET_INPUT");
      vector<PseudoJet> particles;
      particles.reserve(fs.particles().size());
      foreach (const Particle &p, fs.particles()){
        particles.push_back(p.pseudojet());
      }

      // cluster the event and get the jets once and for all
      vector<PseudoJet> jets = _sel_jets(_jet_def(particles));

      // loop over all classifieers
      for (auto &classifier : _qg_xs)
        classifier.process(jets, weight);      
    }
     
    /// Normalise histograms etc., after the run
    void finalize() { }

  private:
    vector<QGClassifiedXS> _qg_xs;
    
    shared_ptr<contrib::ModifiedMassDropTagger> _mmdt;
    shared_ptr<contrib::SoftDrop> _sd;

    JetDefinition _jet_def;
    Selector _sel_jets;
  };
  

  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_LH19_QG4PDF);

}
