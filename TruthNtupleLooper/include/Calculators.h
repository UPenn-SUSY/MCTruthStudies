#ifndef CALCULATORS_H
#define CALCULATORS_H

#include <vector>
#include "TruthNtupleLooper/include/ObjectDefs.h"

// =============================================================================
// namespace TruthNtuple
// {
//   class Particle;
//   class Lepton;
//   class Electron;
//   class Muon;
//   class Jet;
// }

// =============================================================================
namespace TruthNtuple
{
  double comEnergy( const TruthNtuple::Particle*
                  , const TruthNtuple::Particle*
                  );
  double invariantMass( const TruthNtuple::Particle*
                      , const TruthNtuple::Particle*
                      );
  double ptDiObject( const TruthNtuple::Particle*
                   , const TruthNtuple::Particle*
                   );
  double emmaMt( const TruthNtuple::Particle*
               , const TruthNtuple::Particle*
               );
  double deltaEta( const TruthNtuple::Particle*
                 , const TruthNtuple::Particle*
                 );
  double deltaEta( double eta1
                 , double eta2
                 );
  double deltaPhi( const TruthNtuple::Particle*
                 , const TruthNtuple::Particle*
                 );
  double deltaPhi( double phi1
                 , double phi2
                 );
  double deltaR( const TruthNtuple::Particle*
               , const TruthNtuple::Particle*
               );
  double deltaR( double eta1
               , double phi1
               , double eta2
               , double phi2
               );
  double ht(  const std::vector<TruthNtuple::Particle*>&
	    , const std::vector<TruthNtuple::Particle*>&
	    , const std::vector<TruthNtuple::Particle*>&
	    );

  double thetaFromEta(double eta);
  double pFromEM(double e, double m);
  double ptFromPTheta(double p, double theta);
  double ptFromPEta(double p, double eta);
  double pxFromPtPhi(double pt, double phi);
  double pyFromPtPhi(double pt, double phi);
  double pzFromPtTheta(double pt, double theta);
  double pzFromPtEta(double pt, double eta);
  double eFromPxPyPzM(double px, double py, double pz, double m);

  double rapidity(const TruthNtuple::Particle*);

  // get list of invariant masses by combining two lists of particles.
  // The method is defined as follows
  //    0: truth matching to parent barcode
  //    1: max dphi matching
  //    2: min dr matching
  template <class T1, class T2>
    std::vector<double> getInvariantMassList( const std::vector<T1*>&
                                            , const std::vector<T2*>&
                                            , unsigned int method = 0
                                            );
}

#include "TruthNtupleLooper/include/Calculators.icc"

#endif
