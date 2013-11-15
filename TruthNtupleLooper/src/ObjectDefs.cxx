#include "include/ObjectDefs.h"
#include "include/TruthNtupleLooper.h"
#include "TruthRecordHelpers/include/ParentFinder.h"
// #include "TruthRecordHelpers/include/JetFlavorFinder.h"
#include "TruthRecordHelpers/src/JetFlavorFinder.cxx"

#include <iostream>
#include <math.h>

// -----------------------------------------------------------------------------
TruthNtuple::Particle::Particle()
{
  m_index = 0;
  m_pt = 0;
  m_eta = 0;
  m_phi = 0;
  m_e = 0;
  m_px = 0;
  m_py = 0;
  m_pz = 0;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setIndex(unsigned int val)
{
  m_index = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setPt(double val)
{
  m_pt = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setEta(double val)
{
  m_eta = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setPhi(double val)
{
  m_phi = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setE(double val)
{
  m_e = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setPx(double val)
{
  m_px = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setPy(double val)
{
  m_py = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setPz(double val)
{
  m_pz = val;
}


// -----------------------------------------------------------------------------
unsigned int TruthNtuple::Particle::getIndex() const
{
  return m_index;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getPt() const
{
  return m_pt;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getEta() const
{
  return m_eta;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getPhi() const
{
  return m_phi;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getE() const
{
  return m_e;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getPx() const
{
  return m_px;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getPy() const
{
  return m_py;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getPz() const
{
  return m_pz;
}

// -----------------------------------------------------------------------------
TruthNtuple::Electron::Electron()
{
}

// -----------------------------------------------------------------------------
TruthNtuple::Electron::Electron( const TruthNtuple::TruthNtupleLooper* tnl
                               , unsigned int el_index
                               )
{
  setIndex(el_index);

  setPt(    tnl->el_pt->at(el_index));
  setEta(   tnl->el_eta->at(el_index));
  setPhi(   tnl->el_phi->at(el_index));
  setE(     tnl->el_E->at(el_index));
  setPx(    tnl->el_px->at(el_index));
  setPy(    tnl->el_py->at(el_index));
  setPz(    tnl->el_pz->at(el_index));
  setCharge(tnl->el_charge->at(el_index));

  setParentPdgid(TruthRecordHelpers::getParentPdgIdFromBarcode( tnl->el_barcode->at(el_index)
                                                              , tnl->mc_barcode
                                                              , tnl->mc_pdgId
                                                              , tnl->mc_parent_index
                                                              )
                );
}

// -----------------------------------------------------------------------------
void TruthNtuple::Electron::setCharge(double val)
{
  m_charge = val;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Electron::getCharge() const
{
  return m_charge;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Electron::setParentPdgid(int val)
{
  m_parent_pdgid = val;
}
//
// -----------------------------------------------------------------------------
int TruthNtuple::Electron::getParentPdgid() const
{
  return m_parent_pdgid;
}

// -----------------------------------------------------------------------------
TruthNtuple::Muon::Muon()
{
}

// -----------------------------------------------------------------------------
TruthNtuple::Muon::Muon( const TruthNtuple::TruthNtupleLooper* tnl
                       , unsigned int mu_index
                       )
{
  setIndex(mu_index);

  setPt(    tnl->mu_staco_pt->at(mu_index));
  setEta(   tnl->mu_staco_eta->at(mu_index));
  setPhi(   tnl->mu_staco_phi->at(mu_index));
  setE(     tnl->mu_staco_E->at(mu_index));
  setPx(    tnl->mu_staco_px->at(mu_index));
  setPy(    tnl->mu_staco_py->at(mu_index));
  setPz(    tnl->mu_staco_pz->at(mu_index));
  setCharge(tnl->mu_staco_charge->at(mu_index));

  setParentPdgid(TruthRecordHelpers::getParentPdgIdFromBarcode( tnl->mu_staco_barcode->at(mu_index)
                                                              , tnl->mc_barcode
                                                              , tnl->mc_pdgId
                                                              , tnl->mc_parent_index
                                                              )
                );
}

// -----------------------------------------------------------------------------
void TruthNtuple::Muon::setCharge(double val)
{
  m_charge = val;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Muon::getCharge() const
{
  return m_charge;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Muon::setParentPdgid(int val)
{
  m_parent_pdgid = val;
}

// -----------------------------------------------------------------------------
int TruthNtuple::Muon::getParentPdgid() const
{
  return m_parent_pdgid;
}

// -----------------------------------------------------------------------------
TruthNtuple::Jet::Jet()
{
}

// -----------------------------------------------------------------------------
TruthNtuple::Jet::Jet( const TruthNtuple::TruthNtupleLooper* tnl
                     , unsigned int jet_index
                     )
{
  setIndex(jet_index);

  setPt( tnl->jet_AntiKt4TruthJets_pt->at(jet_index));
  setEta(tnl->jet_AntiKt4TruthJets_eta->at(jet_index));
  setPhi(tnl->jet_AntiKt4TruthJets_phi->at(jet_index));
  setE(  tnl->jet_AntiKt4TruthJets_E->at(jet_index));

  setTheta( 2*atan( exp( -fabs(m_eta))) * m_eta/fabs(m_eta));

  setPx( m_pt*cos(m_phi));
  setPy( m_pt*sin(m_phi));
  setPz( m_pt*sin(m_theta));

  setIsBJet( TruthRecordHelpers::isBJet( m_eta
                                       , m_phi
                                       , tnl->mc_pdgId
                                       , tnl->mc_status
                                       , tnl->mc_barcode
                                       , tnl->mc_pt
                                       , tnl->mc_eta
                                       , tnl->mc_phi
                                       // , false
                                       )
           );
}

// -----------------------------------------------------------------------------
void TruthNtuple::Jet::setTheta(double val)
{
  m_theta = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Jet::setIsBJet(bool val)
{
  m_is_b_jet = val;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Jet::getTheta() const
{
  return m_theta;
}

// -----------------------------------------------------------------------------
bool TruthNtuple::Jet::getIsBJet() const
{
  return m_is_b_jet;
}
