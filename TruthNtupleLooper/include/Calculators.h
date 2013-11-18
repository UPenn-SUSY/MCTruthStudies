#ifndef CALCULATORS_H
#define CALCULATORS_H

// =============================================================================
namespace TruthNtuple
{
  class Particle;
  class Lepton;
  class Electron;
  class Muon;
  class Jet;
}

// =============================================================================
namespace TruthNtupleLooper
{
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
}

#endif
