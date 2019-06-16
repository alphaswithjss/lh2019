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
    QGClassifiedXS(FunctionOfPseudoJet<double> *v_ptr,
                   const vector<double> &vcuts,
                   vector<Histo1DPtr> &hs_qq, vector<Histo1DPtr> &hs_qg, vector<Histo1DPtr> &hs_gg,
                   bool integer_valued=false)
      : _v(v_ptr), _vcuts(vcuts), _hs_qq(hs_qq), _hs_qg(hs_qg), _hs_gg(hs_gg), _integer_valued(integer_valued){}

    /// process an event (given by a set of jets and a weight)
    void process(const vector<PseudoJet> & jets, double weight){
      // compute the discriminant on the 2 leading jets
      vector<double> vs = {(jets.size()>0) ? (*_v)(jets[0]) : -1.0,
                           (jets.size()>1) ? (*_v)(jets[1]) : -1.0};

      // compute the jet pts
      vector<double> pts;
      pts.reserve(jets.size());
      for (const auto &jet : jets) pts.push_back(jet.pt());

      if (_integer_valued){
        // for multiplicities, we interpolate between integer cuts
        // for each cut, classify and bin
        for (unsigned int icut=0; icut<_vcuts.size(); ++icut){
          vector<double> gluon_weights;
          for (auto const v : vs){
            if (v>=_vcuts[icut]) gluon_weights.push_back(1.0);
            else{
              if (v<_vcuts[icut]-1) gluon_weights.push_back(0.0);
              else gluon_weights.push_back(1-(_vcuts[icut]-v));
            }
          }

          double w_qq = (1-gluon_weights[0])*(1-gluon_weights[1]);
          double w_gg = gluon_weights[0]*gluon_weights[1];
          double w_qg = 1-w_qq-w_gg;
        
          // bin the jets
          for (const auto &pt : pts){
            _hs_qq[icut]->fill(pt, weight*w_qq);
            _hs_qg[icut]->fill(pt, weight*w_qg);
            _hs_gg[icut]->fill(pt, weight*w_gg);
          }
        } // loop over cuts
      } else {
        // for each cut, classify and bin
        for (unsigned int icut=0; icut<_vcuts.size(); ++icut){
          // classify the event
          unsigned int nquarks = 0;
          for (auto const v : vs) if (v<_vcuts[icut]) ++nquarks;

          Histo1DPtr h = (nquarks==0) ? _hs_gg[icut] : ((nquarks==1) ? _hs_qg[icut] : _hs_qq[icut]);

          // bin the jets
          for (const auto &pt : pts) h->fill(pt, weight);
        }
      } // loop over cuts
    }

    vector<Histo1DPtr> & hs_qq(){ return _hs_qq;}
    vector<Histo1DPtr> & hs_qg(){ return _hs_qg;}
    vector<Histo1DPtr> & hs_gg(){ return _hs_gg;}
    
  protected:
    shared_ptr<FunctionOfPseudoJet<double> > _v;
    vector<double> _vcuts;
    vector<Histo1DPtr> _hs_qq, _hs_qg, _hs_gg;
    bool _integer_valued;
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

    const vector<double> CUTS_PLAIN_ECF; ///< ./get-cuts-and-stats.py -pt 2000 -log -fstep 0.05 -shape plain_ecf_0.5
    const vector<double> CUTS_LOOSE_ECF; ///< ./get-cuts-and-stats.py -pt 2000 -log -fstep 0.05 -shape sd_ecf_0.5
    const vector<double> CUTS_TIGHT_ECF; ///< ./get-cuts-and-stats.py -pt 2000 -log -fstep 0.05 -shape mmdt_ecf_0.5
    const vector<double> CUTS_PLAIN_ISD; ///< ./get-cuts-and-stats.py -pt 2000 -fstep 0.1 -shape plain_isd_1.0
    const vector<double> CUTS_LOOSE_ISD; ///< ./get-cuts-and-stats.py -pt 2000 -fstep 0.1 -shape sd_isd_1.0
    const vector<double> CUTS_TIGHT_ISD; ///< ./get-cuts-and-stats.py -pt 2000 -fstep 0.1 -shape mmdt_isd_1.0
    
    /// Constructor
    MC_LH19_QG4PDF()
      : Analysis("MC_LH19_QG4PDF"),
        PARTICLE_RAPMAX(5.0),
        JET_RADIUS(0.4),
        JET_PTMIN (100.0),
        JET_PTMAX (3925.0),
        JET_NPT   (50),
        JET_RAPMAX(4.5),
        MMDT_ZCUT(0.1),
        SD_BETA  (2.0),
        SD_ZCUT  (0.05),
        CUTS_PLAIN_ECF({0.07168, 0.11203, 0.16660, 0.23204}),
        CUTS_LOOSE_ECF({0.06832, 0.11205, 0.16871, 0.23468}),
        CUTS_TIGHT_ECF({0.05072, 0.11762, 0.18024, 0.24325}),
        CUTS_PLAIN_ISD({4.29675, 6.10855, 7.31715, 8.41689}),
        CUTS_LOOSE_ISD({3.36962, 5.06116, 6.26175, 7.41097}),
        CUTS_TIGHT_ISD({1.95088, 3.61970, 4.95275, 6.59159})
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
      _declare_classifiers(new EnergyCorrelationFunction(0.5, JET_RADIUS),              CUTS_PLAIN_ECF, "plain_ecf");
      _declare_classifiers(new EnergyCorrelationFunction(0.5, JET_RADIUS, _sd.get()),   CUTS_LOOSE_ECF, "loose_ecf");
      _declare_classifiers(new EnergyCorrelationFunction(0.5, JET_RADIUS, _mmdt.get()), CUTS_TIGHT_ECF, "tight_ecf");
      
      _declare_classifiers(new IteratedSoftDropKtMultiplicity(1.0, JET_RADIUS),                   CUTS_PLAIN_ISD, "plain_isd", true);
      _declare_classifiers(new IteratedSoftDropKtMultiplicity(1.0, JET_RADIUS, SD_ZCUT, SD_BETA), CUTS_LOOSE_ISD, "loose_isd", true);
      _declare_classifiers(new IteratedSoftDropKtMultiplicity(1.0, JET_RADIUS, MMDT_ZCUT),        CUTS_TIGHT_ISD, "tight_isd", true);
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
    void finalize() {
      double norm = 1.0/sumOfWeights();
      for (auto &classifier : _qg_xs){
        for (auto &h : classifier.hs_qq()) scale(h, norm); 
        for (auto &h : classifier.hs_qg()) scale(h, norm); 
        for (auto &h : classifier.hs_gg()) scale(h, norm); 
      }
    }

  protected:
    void _declare_classifiers(FunctionOfPseudoJet<double> *v, const vector<double> &vcuts, const std::string & vname, bool integer_valued=false){
      assert(vcuts.size() == 4);
      vector<Histo1DPtr> hs_qq, hs_qg, hs_gg;
      for (const auto cutname : {"very_loose", "loose", "tight", "very_tight"}){
        hs_qq.push_back(bookHisto1D(vname+"_"+cutname+"_qq",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)));
        hs_qg.push_back(bookHisto1D(vname+"_"+cutname+"_qg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)));
        hs_gg.push_back(bookHisto1D(vname+"_"+cutname+"_gg",logspace(JET_NPT, JET_PTMIN, JET_PTMAX)));
      }
      _qg_xs.push_back(QGClassifiedXS(v, vcuts, hs_qq, hs_qg, hs_gg, integer_valued));        
    }


    
    vector<QGClassifiedXS> _qg_xs;
    
    shared_ptr<contrib::ModifiedMassDropTagger> _mmdt;
    shared_ptr<contrib::SoftDrop> _sd;

    JetDefinition _jet_def;
    Selector _sel_jets;
  };
  

  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_LH19_QG4PDF);

}
