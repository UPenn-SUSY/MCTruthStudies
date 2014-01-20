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
  m_mc_index = 0;
  m_pdgid = 0;
  m_pt = 0;
  m_eta = 0;
  m_phi = 0;
  m_e = 0;
  m_m = 0;
  m_px = 0;
  m_py = 0;
  m_pz = 0;
}

// -----------------------------------------------------------------------------
TruthNtuple::Particle::Particle( const TruthNtuple::TruthNtupleLooper* tnl
                               , int mc_index
                               , bool verbose
                               )
{
  if (verbose) {
    std::cout << "Initializing particle - mc_index: " << mc_index << "\n";
  }

  // set quantities directly from d3pd
  setMCIndex(mc_index);
  setPdgid(tnl->mc_pdgId->at(mc_index));
  setStatus(tnl->mc_status->at(mc_index));
  setPt(tnl->mc_pt->at(mc_index));
  setEta(tnl->mc_eta->at(mc_index));
  setPhi(tnl->mc_phi->at(mc_index));
  setM(tnl->mc_m->at(mc_index));

  // set derived quantities
  setTheta(TruthNtuple::thetaFromEta(m_eta));
  setPx(TruthNtuple::pxFromPtPhi(m_pt, m_phi));
  setPy(TruthNtuple::pyFromPtPhi(m_pt, m_phi));
  setPz(TruthNtuple::pzFromPtTheta(m_pt, m_theta));
  setE(eFromPxPyPzM(m_px, m_py, m_pz, m_m));

  // set info about parents
  if (verbose) {
    std::cout << "about to get parent index\n";
  }
  setParentMCIndex( TruthRecordHelpers::getParentIndex( m_mc_index
                                                      , tnl->mc_pdgId
                                                      , tnl->mc_parent_index
                                                      // , verbose
                                                      // , true, true
                                                      )
                  );
  if (verbose) {
    std::cout << "parent index: " << m_parent_index << "\n";
  }
  if (m_parent_index >= 0) {
    setParentPdgid(tnl->mc_pdgId->at(m_parent_index));
    setParentBarcode(tnl->mc_barcode->at(m_parent_index));
  }
  else {
    setParentPdgid(0);
    setParentBarcode(0);
  }
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setMCIndex(int val)
{
  m_mc_index = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setPdgid(int val)
{
  m_pdgid = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setStatus(int val)
{
  m_status = val;
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
void TruthNtuple::Particle::setTheta(double val)
{
  m_theta = val;
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
void TruthNtuple::Particle::setM(double val)
{
  m_m = val;
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
void TruthNtuple::Particle::setParentMCIndex(int val)
{
  m_parent_index = val;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Particle::setParentBarcode(int val)
{
  m_parent_barcode = val;
}


// -----------------------------------------------------------------------------
int TruthNtuple::Particle::getMCIndex() const
{
  return m_mc_index;
}

// -----------------------------------------------------------------------------
int TruthNtuple::Particle::getPdgid() const
{
  return m_pdgid;
}

// -----------------------------------------------------------------------------
int TruthNtuple::Particle::getStatus() const
{
  return m_status;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getPt() const
{
  return m_pt;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getP() const
{
  return sqrt(m_px*m_px + m_py*m_py + m_pz*m_pz);
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getEta() const
{
  return m_eta;
}

// -----------------------------------------------------------------------------
double TruthNtuple::Particle::getTheta() const
{
  return m_theta;
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
double TruthNtuple::Particle::getM() const
{
  return m_m;
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
int TruthNtuple::Particle::getParentMCIndex() const
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
  std::cout << "\tmc index: " << m_mc_index
            << "\tpdg id: "  << m_pdgid
            // << "\tbarcode: "  << m_barcode
            << "\n"
            << "\t\tpt: " << m_pt
            << "\t\teta: " << m_eta
            << "\t\tphi: " << m_phi
            << "\t\tE: " << m_e
            << "\t\tm: " << m_m
            << "\n"
            << "\t\tpx: " << m_px
            << "\t\tpy: " << m_py
            << "\t\tpz: " << m_pz
            << "\n"
            << "\t\tparent index: "   << m_parent_index
            << "\t\tparent pdg id: "  << m_parent_pdgid
            << "\t\tparent barcode: " << m_parent_barcode
            << "\n";
}

// =============================================================================
// = Lepton
// =============================================================================
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
                               , int el_index
                               , bool get_final_state
                               , bool verbose
                               )
{
  if (verbose) {
    std::cout << "Initializing electron\n";
  }

  setIsElectron(true);
  setElIndex(el_index);

  if (verbose) {
    std::cout << "about to look for index from barcode\n";
  }
  setMCIndex( TruthRecordHelpers::getIndexFromBarcode( tnl->el_barcode->at(el_index)
                                                     , tnl->mc_barcode
                                                     // , true
                                                     )
            );

  // set quantities directly from d3pd
  if (get_final_state) {
    setStatus(tnl->el_status->at(el_index));
    setPt(    tnl->el_pt->at(el_index));
    setEta(   tnl->el_eta->at(el_index));
    setPhi(   tnl->el_phi->at(el_index));
    setM(     tnl->el_m->at(el_index));
    setCharge(tnl->el_charge->at(el_index));
  }
  else {
    if (verbose) std::cout << "about to get initial index for electron - mc index: " << m_mc_index << "\n";
    setMCIndex( TruthRecordHelpers::getInitialIndex( m_mc_index
                                                   , tnl->mc_pdgId
                                                   , tnl->mc_parent_index
                                                   , verbose
                                                   )
              );
    if (m_mc_index >= 0) {
      if (verbose) std::cout << "getting electron quantities from MC truth block - mc index: " << m_mc_index << "\n";
      setStatus(tnl->mc_status->at(m_mc_index));
      setPt(    tnl->mc_pt->at(m_mc_index));
      setEta(   tnl->mc_eta->at(m_mc_index));
      setPhi(   tnl->mc_phi->at(m_mc_index));
      setM(     tnl->mc_m->at(m_mc_index));
      setCharge(tnl->el_charge->at(el_index));
    }
    else {
      if (verbose) std::cout << "getting electron quantities from electron block - mc index: " << m_mc_index << "\n";
      setStatus(tnl->el_status->at(el_index));
      setPt(    tnl->el_pt->at(el_index));
      setEta(   tnl->el_eta->at(el_index));
      setPhi(   tnl->el_phi->at(el_index));
      setM(     tnl->el_m->at(el_index));
      setCharge(tnl->el_charge->at(el_index));
    }
  }

  // set derived quantities
  setPdgid(11*m_charge);
  if (fabs(m_pdgid) != 11) std::cout << "ERROR pdg id for electron != 11\n";
  setTheta(TruthNtuple::thetaFromEta(m_eta));
  setPx(TruthNtuple::pxFromPtPhi(m_pt, m_phi));
  setPy(TruthNtuple::pyFromPtPhi(m_pt, m_phi));
  setPz(TruthNtuple::pzFromPtTheta(m_pt, m_theta));
  setE(eFromPxPyPzM(m_px, m_py, m_pz, m_m));

  // set info about parents
  if (verbose) std::cout << "about to get parent index for electron\n";

  setParentMCIndex( TruthRecordHelpers::getParentIndex( m_mc_index
                                                      , tnl->mc_pdgId
                                                      , tnl->mc_parent_index
                                                      // , true
                                                      // , true, true
                                                      )
                  );
  if (verbose) std::cout << "parent index: " << m_parent_index << "\n";
  if (m_parent_index >= 0) {
    setParentPdgid(tnl->mc_pdgId->at(m_parent_index));
    setParentBarcode(tnl->mc_barcode->at(m_parent_index));
  }
  else {
    setParentPdgid(0);
    setParentBarcode(0);
  }
}

// -----------------------------------------------------------------------------
void TruthNtuple::Electron::setElIndex(int val)
{
  m_el_index = val;
}

// -----------------------------------------------------------------------------
int TruthNtuple::Electron::getElIndex() const
{
  return m_el_index;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Electron::print(TruthNtuple::TruthNtupleLooper* tnl) const
{
  Lepton::print(tnl);

  size_t num_parents = tnl->el_parent_index->at(m_el_index).size();
  std::cout << "parents - " << num_parents << "\n";
  for (size_t parent_it = 0; parent_it != num_parents; ++parent_it) {
    std::cout << "\tparent it: " << parent_it << " - " << tnl->el_parent_index->at(m_el_index).at(parent_it) << "\n";
  }
}

// =============================================================================
// = Muon
// =============================================================================
TruthNtuple::Muon::Muon()
{
  setIsElectron(false);
}

// -----------------------------------------------------------------------------
TruthNtuple::Muon::Muon( const TruthNtuple::TruthNtupleLooper* tnl
                       , int mu_index
                       , bool get_final_state
                       , bool verbose
                       )
{
  if (verbose) {
    std::cout << "Initializing muon\n";
  }

  setIsElectron(false);
  setMuIndex(mu_index);

  setMCIndex( TruthRecordHelpers::getIndexFromBarcode( tnl->mu_barcode->at(mu_index)
                                                     , tnl->mc_barcode
                                                     // , true
                                                     )
            );

  // set quantities directly from d3pd
  if (get_final_state) {
    setStatus(tnl->mu_status->at(mu_index));
    setPt(    tnl->mu_pt->at(mu_index));
    setEta(   tnl->mu_eta->at(mu_index));
    setPhi(   tnl->mu_phi->at(mu_index));
    setM(     tnl->mu_m->at(mu_index));
    setCharge(tnl->mu_charge->at(mu_index));
  }
  else {
    if (verbose) std::cout << "about to get initial index for muon\n";
    setMCIndex( TruthRecordHelpers::getInitialIndex( m_mc_index
                                                  , tnl->mc_pdgId
                                                  , tnl->mc_parent_index
                                                  , verbose
                                                  )
              );
    if (m_mc_index >= 0) {
      setStatus(tnl->mc_status->at(m_mc_index));
      setPt(    tnl->mc_pt->at(m_mc_index));
      setEta(   tnl->mc_eta->at(m_mc_index));
      setPhi(   tnl->mc_phi->at(m_mc_index));
      setM(     tnl->mc_m->at(m_mc_index));
      setCharge(tnl->mu_charge->at(mu_index));
    }
    else {
      setStatus(tnl->mu_status->at(mu_index));
      setPt(    tnl->mu_pt->at(mu_index));
      setEta(   tnl->mu_eta->at(mu_index));
      setPhi(   tnl->mu_phi->at(mu_index));
      setM(     tnl->mu_m->at(mu_index));
      setCharge(tnl->mu_charge->at(mu_index));
    }
  }

  // set derived quantities
  setPdgid(13*m_charge);
  if (fabs(m_pdgid) != 13) std::cout << "ERROR pdg id for muon != 13\n";
  setTheta(TruthNtuple::thetaFromEta(m_eta));
  setPx(TruthNtuple::pxFromPtPhi(m_pt, m_phi));
  setPy(TruthNtuple::pyFromPtPhi(m_pt, m_phi));
  setPz(TruthNtuple::pzFromPtTheta(m_pt, m_theta));
  setE(eFromPxPyPzM(m_px, m_py, m_pz, m_m));

  // set info about parents
  if (verbose) {
    std::cout << "about to get parent index\n";
  }
  setParentMCIndex(TruthRecordHelpers::getParentIndex( m_mc_index
                                                     , tnl->mc_pdgId
                                                     , tnl->mc_parent_index
                                                     // , true
                                                     // , true, true
                                                     )
                  );
  if (verbose) {
    std::cout << "muon parent index: " << m_parent_index << "\n";
  }
  if (m_parent_index >= 0) {
    setParentPdgid(tnl->mc_pdgId->at(m_parent_index));
    setParentBarcode(tnl->mc_barcode->at(m_parent_index));
  }
  else {
    setParentPdgid(0);
    setParentBarcode(0);
  }
}

// -----------------------------------------------------------------------------
void TruthNtuple::Muon::setMuIndex(int val)
{
  m_mu_index = val;
}

// -----------------------------------------------------------------------------
int TruthNtuple::Muon::getMuIndex() const
{
  return m_mu_index;
}

// -----------------------------------------------------------------------------
void TruthNtuple::Muon::print(TruthNtuple::TruthNtupleLooper* tnl) const
{
  Lepton::print(tnl);

  size_t num_parents = tnl->mu_parent_index->at(m_mu_index).size();
  std::cout << "parents - " << num_parents << "\n";
  for (size_t parent_it = 0; parent_it != num_parents; ++parent_it) {
    std::cout << "\tparent it: " << parent_it << " - " << tnl->mu_parent_index->at(m_mu_index).at(parent_it) << "\n";
  }
}

// =============================================================================
// = Jet
// =============================================================================
TruthNtuple::Jet::Jet()
{
}

// -----------------------------------------------------------------------------
TruthNtuple::Jet::Jet( const TruthNtuple::TruthNtupleLooper* tnl
                     , int jet_index
                     , bool verbose
                     )
{
  if (verbose) {
    std::cout << "Initializing jet - jet index: " << jet_index << "\n";
  }

  setJetIndex(jet_index);

  // set quantities directly from d3pd
  setPt( tnl->jet_AntiKt4TruthJets_WZ_pt->at(jet_index));
  setEta(tnl->jet_AntiKt4TruthJets_WZ_eta->at(jet_index));
  setPhi(tnl->jet_AntiKt4TruthJets_WZ_phi->at(jet_index));
  setM(  tnl->jet_AntiKt4TruthJets_WZ_m->at(jet_index));

  // set derived quantities
  setTheta(TruthNtuple::thetaFromEta(m_eta));
  setPx(TruthNtuple::pxFromPtPhi(m_pt, m_phi));
  setPy(TruthNtuple::pyFromPtPhi(m_pt, m_phi));
  setPz(TruthNtuple::pzFromPtTheta(m_pt, m_theta));
  setE(eFromPxPyPzM(m_px, m_py, m_pz, m_m));

  // set info about parents
  setIsBJet(fabs(tnl->jet_AntiKt4TruthJets_WZ_flavor_partonFlavor->at(jet_index)) == 5);
  if (m_is_b_jet) {
    if (verbose) {
      std::cout << "about to find b quark index\n";
    }
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
    setMCIndex(m_b_quark_index);

    if (verbose) {
      std::cout << "about to find parent index\n";
    }
    setParentMCIndex(TruthRecordHelpers::getParentIndex( m_b_quark_index
                                                       , tnl->mc_pdgId
                                                       , tnl->mc_parent_index
                                                       // , true, true
                                                       )
                    );
    if (verbose) {
      std::cout << "jet parent index: " << m_parent_index << "\n";
    }
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
    setMCIndex(0);
    setParentPdgid(0);
  }
}

// -----------------------------------------------------------------------------
void TruthNtuple::Jet::setJetIndex(int val)
{
  m_jet_index = val;
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
int TruthNtuple::Jet::getJetIndex() const
{
  return m_jet_index;
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
TruthNtuple::Met::Met() : m_met_etx_noint(0.)
                        , m_met_ety_noint(0.)
                        , m_met_et_noint(0.)
{
}

// -----------------------------------------------------------------------------
TruthNtuple::Met::Met( double met_etx_noint
                     , double met_ety_noint
                     , bool verbose
                     ) : m_met_etx_noint(met_etx_noint)
                       , m_met_ety_noint(met_etx_noint)
{
  if (verbose) {
    std::cout << "initializing met\n";
  }
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
