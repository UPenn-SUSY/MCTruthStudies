#include "TruthNtupleLooper/include/Calculators.h"
#include <math.h>

#include "TruthNtupleLooper/include/ObjectDefs.h"

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

// -----------------------------------------------------------------------------
std::vector<double> TruthNtuple::getMbl( const std::vector<Lepton*>& leptons
                                       , const std::vector<Jet*>& jets
                                       , unsigned int mbl_method
                                       )
{
  const size_t num_lep = leptons.size();
  const size_t num_jet = jets.size();

  std::vector<double> mbl_list;

  if (mbl_method == 0) {
    // TODO implement truth matching method
  }
  else if (mbl_method == 1) {
    // TODO implement dphi matching method
  }
  else if (mbl_method == 2) {
    for (size_t lep_it = 0; lep_it != num_lep; ++lep_it) {
      float min_dr = 999;
      float dr = 999;
      int chosen_jet_it = -1;
      for (size_t jet_it = 0; jet_it != num_jet; ++jet_it) {
        dr = TruthNtuple::deltaR(leptons.at(lep_it), jets.at(jet_it));
        if (dr < min_dr) {
          min_dr = dr;
          chosen_jet_it = jet_it;
        }
      }
      if (chosen_jet_it >= 0) {
        mbl_list.push_back( TruthNtuple::invariantMass( leptons.at(lep_it)
                                                      , jets.at(chosen_jet_it)
                                                      )
                          );
      }
    }
  }

  return mbl_list;
}

