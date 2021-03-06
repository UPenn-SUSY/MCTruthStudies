#include "TruthNtupleLooper/include/Calculators.h"
#include <math.h>
#include <stdlib.h>

#include "TruthNtupleLooper/include/ObjectDefs.h"

// =============================================================================
static const double PI = 3.14159265359;

// -----------------------------------------------------------------------------
double TruthNtuple::comEnergy( const TruthNtuple::Particle* p1
                             , const TruthNtuple::Particle* p2
                             )
{
  double px = p1->getPx() + p2->getPx();
  double py = p1->getPy() + p2->getPy();
  double pz = p1->getPz() + p2->getPz();
  double m  = p1->getM()  + p2->getM();

  double e2 = m*m + px*px + py*py + pz*pz;
  return e2/sqrt(fabs(e2));
}

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
  return m2/sqrt(fabs(m2));
}

// -----------------------------------------------------------------------------
double TruthNtuple::ptDiObject( const TruthNtuple::Particle* p1
                              , const TruthNtuple::Particle* p2
                              )
{
  double px = p1->getPx() + p2->getPx();
  double py = p1->getPy() + p2->getPy();

  double pt2 = px*px + py*py;
  return sqrt(pt2)*pt2/fabs(pt2);
}

// -----------------------------------------------------------------------------
double TruthNtuple::emmaMt( const TruthNtuple::Particle* p1
                          , const TruthNtuple::Particle* p2
                          )
{
  double mll = invariantMass(p1, p2);
  double ptll = ptDiObject(p1, p2);

  double emma_mt_2 = mll*mll + ptll*ptll;
  return emma_mt_2/sqrt(fabs(emma_mt_2));
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
  // if (eta2 > eta1) {
  //   double tmp = eta1;
  //   eta1 = eta2;
  //   eta2 = eta1;
  // }
  return fabs( eta1 - eta2 );
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
  while (delta_phi > +PI) delta_phi -= 2*PI;
  while (delta_phi < -PI) delta_phi += 2*PI;

  return fabs(delta_phi);
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

// -----------------------------------------------------------------------------
double TruthNtuple::ht(  const std::vector<TruthNtuple::Particle*>& el_list
	   ,const std::vector<TruthNtuple::Particle*>& mu_list
	   ,const std::vector<TruthNtuple::Particle*>& jet_list
	   )
{
  double ht = 0.;

  size_t num_el = el_list.size();
  for (size_t el_it = 0; el_it != num_el; ++el_it) {
    ht += el_list.at(el_it)->getPt();
  }
  size_t num_mu = mu_list.size();
  for (size_t mu_it = 0; mu_it != num_mu; ++mu_it) {
    ht += mu_list.at(mu_it)->getPt();
  }
  size_t num_jet = jet_list.size();
  for (size_t jet_it = 0; jet_it != num_jet; ++jet_it) {
    ht += jet_list.at(jet_it)->getPt();
  }

  return ht;
}

// -----------------------------------------------------------------------------
double TruthNtuple::thetaFromEta(double eta)
{
  return 2*atan( exp( -eta));
}

// -----------------------------------------------------------------------------
double TruthNtuple::pFromEM(double e, double m)
{
  return sqrt(e*e - m*m);
}

// -----------------------------------------------------------------------------
double TruthNtuple::ptFromPTheta(double p, double theta)
{
  return p*sin(theta);
}

// -----------------------------------------------------------------------------
double TruthNtuple::ptFromPEta(double p, double eta)
{
  return p/cosh(eta);
}

// -----------------------------------------------------------------------------
double TruthNtuple::pxFromPtPhi(double pt, double phi)
{
  return pt*cos(phi);
}

// -----------------------------------------------------------------------------
double TruthNtuple::pyFromPtPhi(double pt, double phi)
{
  return pt*sin(phi);
}

// -----------------------------------------------------------------------------
double TruthNtuple::pzFromPtTheta(double pt, double theta)
{
  return pt/tan(theta);
}

// -----------------------------------------------------------------------------
double TruthNtuple::pzFromPtEta(double pt, double eta)
{
  return pzFromPtTheta(pt, thetaFromEta(eta));
}

// -----------------------------------------------------------------------------
double TruthNtuple::eFromPxPyPzM(double px, double py, double pz, double m)
{
  return sqrt(px*px + py*py + pz*pz + m*m);
}

// -----------------------------------------------------------------------------
double TruthNtuple::rapidity(const TruthNtuple::Particle* p)
{
  double e  = p->getE();
  double pz = p->getPz();
  return (0.5*log( (e + pz)/(e-pz) ));
}
