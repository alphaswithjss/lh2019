#ifndef __QG_TAGGERS_HH__
#define __QG_TAGGERS_HH__

#include <fastjet/tools/Recluster.hh>
#include <sstream>
#include <fastjet/contrib/LundGenerator.hh>

/// \class Angularity
/// compute the angularity for a given jet
///
/// This is defined as
/// \f[
///   \frac{1}{\sum_i p_{t,i}R^\alpha} \sum_{i} p_{t,i} \theta_{i,jet}^\alpha
/// \f]
/// where the jet axis is defined by reclustering the jet with the
/// antikt algorithm using the WTA recombination scheme
///
class Angularity : public fastjet::FunctionOfPseudoJet<double>{
public:
  /// ctor
  ///  - alpha    is the exponent
  ///  - groomer  is an optional (pointer to) a groomer to be applied before the computation
  Angularity(double alpha, double R, const fastjet::Transformer *groomer = nullptr)
    : _alpha(alpha), _R(R), _groomer(groomer),
      _recluster(fastjet::JetDefinition(fastjet::antikt_algorithm, 999.0, fastjet::WTA_pt_scheme)){}

  /// dummy dtor
  virtual ~Angularity(){}

  /// description
  std::string description() const override{
    std::ostringstream oss;
    oss << "Angularity with alpha=" << _alpha << " (using WTA axes)" << ", R=" << _R;
    if (_groomer)
      oss << " and grooming using [" << _groomer->description() << "]";
    return oss.str();
  }
  
  /// actual angularity computation
  double result(const fastjet::PseudoJet &jet) const override{
    fastjet::PseudoJet processed_jet = jet;
      
    // first groom the jet (if needed)
    if (_groomer){
      processed_jet = (*_groomer)(jet);
      if (processed_jet == 0) return 0.0; //< this is a convention
    }

    // recluster to get the axis right
    processed_jet = _recluster(processed_jet);

    // compute the angularity
    double num=0.0, den=0.0;
    for (const auto &c : processed_jet.constituents()){
      double pt = c.pt();
      den += pt;
      num += pt * pow(c.squared_distance(processed_jet), 0.5*_alpha);
    }
    if (den<=0.0) return 0.0;
    return num/(den*pow(_R,_alpha));
  }

protected:
  double _alpha;         ///< the angular exponent
  double _R;             ///< the jet radius
  const fastjet::Transformer *_groomer; ///< an optional pre-grooming step
  fastjet::Recluster _recluster;  ///< recluster the jet with C/A
};


/// \class EnergyCorrelationFunction
/// compute the energy correlation function for a given jet
///
/// This is defined as
/// \f[
///   \frac{1}{(\sum_i p_{t,i})^2R^\alpha} \sum_{i<j} p_{t,i} p_{t,j} \theta_{ij}^\alpha
/// \f]
class EnergyCorrelationFunction : public fastjet::FunctionOfPseudoJet<double>{
public:
  /// ctor
  ///  - alpha    is the exponent
  ///  - groomer  is an optional (pointer to) a groomer to be applied before the computation
  EnergyCorrelationFunction(double alpha, double R, const fastjet::Transformer *groomer = nullptr)
    : _alpha(alpha), _R(R), _groomer(groomer){}

  /// dummy dtor
  virtual ~EnergyCorrelationFunction(){}

  /// description
  std::string description() const override{
    std::ostringstream oss;
    oss << "EnergyCorrelationFunction with alpha=" << _alpha << ", R=" << _R;
    if (_groomer)
      oss << " and grooming using [" << _groomer->description() << "]" ;
    return oss.str();
  }
  
  /// actual angularity computation
  double result(const fastjet::PseudoJet &jet) const override{
    fastjet::PseudoJet processed_jet = jet;
      
    // first groom the jet (if needed)
    if (_groomer){
      processed_jet = (*_groomer)(jet);
      if (processed_jet == 0) return 0.0; //< this is a convention
    }

    // compute the shape
    double num=0.0, den=0.0;

    std::vector<fastjet::PseudoJet> constits = processed_jet.constituents();
    std::vector<double> pts;
    pts.reserve(constits.size());
    for (const auto &c : constits){
      double pt = c.pt();
      pts.push_back(pt);
      den += pt;
    }
    if (den<=0.0) return 0.0;
    
    for (unsigned int i=0;i<constits.size();++i){
      for (unsigned int j=i+1;j<constits.size();++j){
        num += pts[i] * pts[j] * pow(constits[i].squared_distance(constits[j]), 0.5*_alpha);
      }
    }
    return num/(den*den*pow(_R,_alpha));
  }

protected:
  double _alpha;         ///< the angular exponent
  double _R;             ///< the jet radius
  const fastjet::Transformer *_groomer; ///< an optional pre-grooming step
};


/// \class IteratedSoftDropKtMultiplicity
/// compute the multuplicity of emissions from the "leaing" branch aove a given kt
/// (with an optional SD condition)
///
/// Here, we look at all 
class IteratedSoftDropKtMultiplicity : public fastjet::FunctionOfPseudoJet<double>{
public:
  /// ctor
  ///  - ktcut    the minimal kt for emissions (same energy dimension as the input)
  ///  - zcut     (optional) SoftDrop condition: zcut parameter (-ve means no SD)
  ///  - beta     (optional) SoftDrop condition: beta parameter
  IteratedSoftDropKtMultiplicity(double ktcut, double R, double zcut = -1.0, double beta =0.0)
    : _ktcut(ktcut), _R(R), _zcut(zcut), _beta(beta){}

  /// dummy dtor
  virtual ~IteratedSoftDropKtMultiplicity(){}

  /// description
  std::string description() const override{
    std::ostringstream oss;
    oss << "IteratedSoftDropKtMultiplicity with kt>=" << _ktcut << ", R=" << _R;
    if (_zcut>0)
      oss << " and condition z>" << _zcut << " theta^" << _beta;
    return oss.str();
  }
  
  /// actual angularity computation
  double result(const fastjet::PseudoJet &jet) const override{
    // get the list of primary Lund de-clusterings
    std::vector<fastjet::contrib::LundDeclustering> declusterings = _lund(jet);

    // looop over those and keep the ones above the ktcut and (optionally) parssing the SD condition
    unsigned int multiplicity = 0;
    for (const auto &d : declusterings){
      if (d.kt()<_ktcut) continue;
      if ((_zcut>0) && (d.z() < _zcut*pow(d.Delta()/_R, _beta))) continue;
      ++multiplicity;
    }

    return 1.0*multiplicity;
  }

protected:
  double _ktcut;  ///< ktcut imposed on emissions
  double _R;      ///< the jet radius
  double _zcut;   ///< (optional SoftDrop condition: zcut parameter
  double _beta;   ///< (optional SoftDrop condition: beta parameter
  fastjet::contrib::LundGenerator _lund;
};

#endif
