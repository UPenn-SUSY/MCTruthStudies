#include "TruthNtupleLooper/include/ObjectDefs.h"
#include "TruthNtupleLooper/include/TruthNtupleLooper.h"
#include "TruthNtupleLooper/include/Calculators.h"

#include "TruthRecordHelpers/include/ParentFinder.h"
#include "TruthRecordHelpers/include/JetFlavorFinder.h"

#include <vector>
#include <iostream>
#include <math.h>

// =============================================================================
static const double PI = 3.14159265359;

// =============================================================================
// = Particle
// =============================================================================
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
void TruthNtuple::Particle::setParentPdgid(int val)
{
  m_parent_pdgid = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setParentIndex(int val)
{
  m_parent_index = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setParentBarcode(int val)
{
  m_parent_barcode = val;
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
int TruthNtuple::Particle::getParentPdgid() const
{
  return m_parent_pdgid;
}

// -----------------------------------------------------------------------------
int TruthNtuple::Particle::getParentIndex() const
{
  return m_parent_index;
}

// -----------------------------------------------------------------------------
int TruthNtuple::Particle::getParentBarcode() const
{
  return m_parent_barcode;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::printGeneralInfo() const
{
  std::cout << "\tindex: " << m_index
            << "\n"
            << "\tpt: " << m_pt
            << "\teta: " << m_eta
            << "\tphi: " << m_phi
            << "\tE: " << m_e
            << "\n"
            << "\tpx: " << m_px
            << "\tpy: " << m_py
            << "\tpz: " << m_pz
            << "\n"
            << "\tparent index: "   << m_parent_index
            << "\tparent pdg id: "  << m_parent_pdgid
            << "\tparent barcode: " << m_parent_barcode
            << "\n";
}

// =============================================================================
// = Lepton
// =============================================================================
// -----------------------------------------------------------------------------
TruthNtuple::Lepton::Lepton()
{
}

// -----------------------------------------------------------------------------
void TruthNtuple::Lepton::setIsElectron(bool val)
{
  m_is_electron = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Lepton::setCharge(double val)
{
  m_charge = val;
}

// -----------------------------------------------------------------------------
bool TruthNtuple::Lepton::isElectron() const
{
  return m_is_electron;
}

// -----------------------------------------------------------------------------
bool TruthNtuple::Lepton::isMuon() const
{
  return !m_is_electron;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Lepton::getCharge() const
{
  return m_charge;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Lepton::print(TruthNtuple::TruthNtupleLooper* /*tnl*/) const
{
  std::cout << "lepton (" << (m_is_electron ? "electron" : "muon") << ")"
            << "\tcharge: " << m_charge
            << "\n";
  printGeneralInfo();
}

// =============================================================================
// = Electron
// =============================================================================
// -----------------------------------------------------------------------------
TruthNtuple::Electron::Electron()
{
  setIsElectron(true);
}

// -----------------------------------------------------------------------------
TruthNtuple::Electron::Electron( const TruthNtuple::TruthNtupleLooper* tnl
                               , unsigned int el_index
                               )
{
  setIsElectron(true);

  setIndex(el_index);

  setPt(    tnl->el_pt->at(el_index));
  setEta(   tnl->el_eta->at(el_index));
  setPhi(   tnl->el_phi->at(el_index));
  setE(     tnl->el_E->at(el_index));
  setPx(    tnl->el_px->at(el_index));
  setPy(    tnl->el_py->at(el_index));
  setPz(    tnl->el_pz->at(el_index));
  setCharge(tnl->el_charge->at(el_index));

  setParentIndex(TruthRecordHelpers::getParentIndexFromBarcode( tnl->el_barcode->at(el_index)
                                                              , tnl->mc_barcode
                                                              , tnl->mc_pdgId
                                                              , tnl->mc_parent_index
                                                              )
                );
  if (m_parent_index >= 0) {
    setParentPdgid(tnl->mc_pdgId->at(m_parent_index));
    setParentBarcode(tnl->mc_barcode->at(m_parent_index));
  }
  else {
    setParentPdgid(0);
    setParentBarcode(0);
  }
}

// =============================================================================
// = Muon
// =============================================================================
// -----------------------------------------------------------------------------
TruthNtuple::Muon::Muon()
{
  setIsElectron(false);
}

// -----------------------------------------------------------------------------
TruthNtuple::Muon::Muon( const TruthNtuple::TruthNtupleLooper* tnl
                       , unsigned int mu_index
                       )
{
  setIsElectron(false);

  setIndex(mu_index);

  setPt(    tnl->mu_staco_pt->at(mu_index));
  setEta(   tnl->mu_staco_eta->at(mu_index));
  setPhi(   tnl->mu_staco_phi->at(mu_index));
  setE(     tnl->mu_staco_E->at(mu_index));
  setPx(    tnl->mu_staco_px->at(mu_index));
  setPy(    tnl->mu_staco_py->at(mu_index));
  setPz(    tnl->mu_staco_pz->at(mu_index));
  setCharge(tnl->mu_staco_charge->at(mu_index));

  setParentIndex(TruthRecordHelpers::getParentIndexFromBarcode( tnl->mu_staco_barcode->at(mu_index)
                                                              , tnl->mc_barcode
                                                              , tnl->mc_pdgId
                                                              , tnl->mc_parent_index
                                                              )
                );
  if (m_parent_index >= 0) {
    setParentPdgid(tnl->mc_pdgId->at(m_parent_index));
    setParentBarcode(tnl->mc_barcode->at(m_parent_index));
  }
  else {
    setParentPdgid(0);
    setParentBarcode(0);
  }
}

// =============================================================================
// = Jet
// =============================================================================
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

  setBQuarkIndex(TruthRecordHelpers::isBJet( m_eta
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
  if (m_b_quark_index >= 0) {
    setIsBJet(true);

    setParentIndex(TruthRecordHelpers::getParentIndex( m_b_quark_index
                                                     , tnl->mc_pdgId
                                                     , tnl->mc_parent_index
                                                     )
                  );
    if (m_parent_index >= 0) {
      setParentPdgid(tnl->mc_pdgId->at(m_parent_index));
      setParentBarcode(tnl->mc_barcode->at(m_parent_index));
    }
    else {
      setParentPdgid(0);
      setParentBarcode(0);
    }
  }
  else {
    setIsBJet(false);
    setParentPdgid(0);
  }
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
void TruthNtuple::Jet::setBQuarkIndex(int val)
{
  m_b_quark_index = val;
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

// -----------------------------------------------------------------------------
int TruthNtuple::Jet::getBQuarkIndex() const
{
  return m_b_quark_index;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Jet::print(TruthNtuple::TruthNtupleLooper* tnl) const
{
  std::cout << "Jet (b jet: " << m_is_b_jet << ")" << "\n";
  printGeneralInfo();
  if (m_is_b_jet && tnl != NULL) {
    int immediate_parent_index   = tnl->mc_parent_index->at(m_b_quark_index).at(0);
    int immediate_parent_pdgid   = tnl->mc_pdgId->at(immediate_parent_index);
    int immediate_parent_barcode = tnl->mc_barcode->at(immediate_parent_index);

    std::cout << "\tb quark index: "    << m_b_quark_index
              << "\tb quark pdg id: " << tnl->mc_pdgId->at(m_b_quark_index)
              << "\tb quark barcode: " << tnl->mc_barcode->at(m_b_quark_index)
              << "\tquark status: " << tnl->mc_status->at(m_b_quark_index)
              << "\n"
              << "\tb quark pt: "     << tnl->mc_pt->at(   m_b_quark_index)
              << "\tb quark eta: "    << tnl->mc_eta->at(  m_b_quark_index)
              << "\tb quark phi: "    << tnl->mc_phi->at(  m_b_quark_index)
              << "\n"
              << "\tdetaR(jet, quark): " << TruthNtuple::deltaR( m_eta
                                                               , m_phi
                                                               , tnl->mc_eta->at(m_b_quark_index)
                                                               , tnl->mc_phi->at(m_b_quark_index)
                                                               )
              << "\tdetaEta(jet, quark): " << TruthNtuple::deltaEta(m_eta, tnl->mc_eta->at(m_b_quark_index))
              << "\tdetaPhi(jet, quark): " << TruthNtuple::deltaPhi(m_phi, tnl->mc_phi->at(m_b_quark_index))
              << "\n";
    std::cout << "\timmediate parent index: "   << immediate_parent_index
              << "\timmediate parent pdg id: "  << immediate_parent_pdgid
              << "\timmediate parent barcode: " << immediate_parent_barcode
              << "\n";
  }
}

// =============================================================================
// = Met
// =============================================================================
// -----------------------------------------------------------------------------
TruthNtuple::Met::Met() : m_met_etx_noint(0.)
                        , m_met_ety_noint(0.)
                        , m_met_et_noint(0.)
{
}

// -----------------------------------------------------------------------------
TruthNtuple::Met::Met(double met_etx_noint, double met_ety_noint) : m_met_etx_noint(met_etx_noint)
                                                                  , m_met_ety_noint(met_etx_noint)
{
  setMetNoint(met_etx_noint, met_ety_noint);
}

// -----------------------------------------------------------------------------
void TruthNtuple::Met::clear()
{
  m_met_etx_noint = 0;
  m_met_ety_noint = 0;
  m_met_et_noint = 0;
  m_met_phi_noint = 0;
  m_met_rel_noint = 0;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Met::setMetNoint(double met_etx, double met_ety)
{
  m_met_etx_noint = met_etx;
  m_met_ety_noint = met_ety;
  m_met_et_noint = sqrt( met_etx*met_etx + met_ety*met_ety );
  m_met_phi_noint = atan2(m_met_ety_noint, m_met_etx_noint);
  m_met_rel_noint = 0;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Met::calculateMetRelNoint(const std::vector<TruthNtuple::Particle*>& particles)
{
  double min_dphi = 999.;
  double dphi = 999.;

  size_t num_particles = particles.size();
  for (size_t p_itr = 0; p_itr != num_particles; ++p_itr) {
    dphi = TruthNtuple::deltaPhi( particles.at(p_itr)->getPhi()
                                , m_met_phi_noint
                                );
    if (dphi < min_dphi)
      min_dphi = dphi;
  }

  m_met_rel_noint = m_met_et_noint;
  if (min_dphi < PI) {
    m_met_rel_noint *= sin(min_dphi);
  }
}

// -----------------------------------------------------------------------------
void TruthNtuple::Met::calculateMetRelNoint( const std::vector<TruthNtuple::Electron*>& electrons
                                           , const std::vector<TruthNtuple::Muon*>& muons
                                           , const std::vector<TruthNtuple::Jet*>& jets
                                           )
{
  std::vector<TruthNtuple::Particle*> particles;
  particles.reserve(electrons.size() + muons.size() + jets.size());
  particles.insert(particles.end(), electrons.begin(), electrons.end());
  particles.insert(particles.end(), muons.begin(), muons.end());
  particles.insert(particles.end(), jets.begin(), jets.end());

  return calculateMetRelNoint(particles);
}

// -----------------------------------------------------------------------------
double TruthNtuple::Met::getMetNoint() const
{
  return m_met_et_noint;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Met::getMetPhiNoint() const
{
  return m_met_phi_noint;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Met::getMetRelNoint() const
{
  return m_met_rel_noint;
}
