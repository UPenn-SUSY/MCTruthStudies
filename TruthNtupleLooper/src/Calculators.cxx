#include "include/Calculators.h"
#include <math.h>

#include "include/ObjectDefs.h"

// -----------------------------------------------------------------------------
double TruthNtuple::invariantMass( const TruthNtuple::Particle* p1
                                 , const TruthNtuple::Particle* p2
                                 )
{
  double px = p1->getPx() + p2->getPx();
  double py = p1->getPy() + p2->getPy();
  double pz = p1->getPz() + p2->getPz();
  double e  = p1->getE()  + p2->getE();

  double m2 = e*e - px*px - py*py - pz*pz;
  return m2*m2/fabs(m2);
}

// -----------------------------------------------------------------------------
double TruthNtuple::ptDiObject( const TruthNtuple::Particle* p1
                              , const TruthNtuple::Particle* p2
                              )
{
  double px = p1->getPx() + p2->getPx();
  double py = p1->getPy() + p2->getPy();

  double pt2 = px*px + py*py;
  return pt2*pt2/fabs(pt2);
}

// -----------------------------------------------------------------------------
double TruthNtuple::emmaMt( const TruthNtuple::Particle* p1
                          , const TruthNtuple::Particle* p2
                          )
{
  double mll = invariantMass(p1, p2);
  double ptll = ptDiObject(p1, p2);

  double emma_mt_2 = mll*mll + ptll*ptll;
  return emma_mt_2*emma_mt_2/fabs(emma_mt_2);
}

// -----------------------------------------------------------------------------
double TruthNtuple::deltaEta( const TruthNtuple::Particle* p1
                            , const TruthNtuple::Particle* p2
                            )
{
  return deltaEta(p1->getEta(), p2->getEta());
}

// -----------------------------------------------------------------------------
double TruthNtuple::deltaEta(double eta1, double eta2)
{
  return fabs( fabs(eta1) - fabs(eta2) );
}

// -----------------------------------------------------------------------------
double TruthNtuple::deltaPhi( const TruthNtuple::Particle* p1
                            , const TruthNtuple::Particle* p2
                            )
{
  return deltaPhi(p1->getPhi(), p2->getPhi());
}

// -----------------------------------------------------------------------------
double TruthNtuple::deltaPhi(double phi1, double phi2)
{
  double delta_phi = fabs( phi1 - phi2 );
  while (delta_phi > 3.14159265359) {
    delta_phi -= 3.14159265359;
  }
  return delta_phi;
}

// -----------------------------------------------------------------------------
double TruthNtuple::deltaR( const TruthNtuple::Particle* p1
                          , const TruthNtuple::Particle* p2
                          )
{
  // return 1.;
  return deltaR(p1->getEta(), p1->getPhi(), p2->getEta(), p2->getPhi());
}

// -----------------------------------------------------------------------------
double TruthNtuple::deltaR( double eta1
                          , double phi1
                          , double eta2
                          , double phi2
                          )
{
  double delta_phi = deltaPhi(phi1, phi2);
  double delta_eta = deltaEta(eta1, eta2);
  return sqrt( delta_phi*delta_phi + delta_eta*delta_eta );
}
