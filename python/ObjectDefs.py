#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import ROOT
import rootlogon
import metaroot

import TruthRecordHelpers as truth_helpers

# ==============================================================================
class Electron(object):
    # ------------------------------------------------------------------------------
    def __init__(self, event, electron_index):
        self.index        = electron_index
        self.pt           = event.el_pt.at(electron_index)
        self.eta          = event.el_eta.at(electron_index)
        self.phi          = event.el_phi.at(electron_index)
        self.E            = event.el_E.at(electron_index)
        self.charge       = event.el_charge.at(electron_index)
        self.px           = event.el_px.at(electron_index)
        self.py           = event.el_py.at(electron_index)
        self.pz           = event.el_pz.at(electron_index)
        self.parent_pdgid = truth_helpers.getParentPdgIDFromBarcode( event
                                                                   , event.el_barcode.at(electron_index)
                                                                   )

# ==============================================================================
class Muon(object):
    # ------------------------------------------------------------------------------
    def __init__(self, event, muon_index):
        self.index        = muon_index
        self.pt           = event.mu_staco_pt.at(muon_index)
        self.eta          = event.mu_staco_eta.at(muon_index)
        self.phi          = event.mu_staco_phi.at(muon_index)
        self.E            = event.mu_staco_E.at(muon_index)
        self.charge       = event.mu_staco_charge.at(muon_index)
        self.px           = event.mu_staco_px.at(muon_index)
        self.py           = event.mu_staco_py.at(muon_index)
        self.pz           = event.mu_staco_pz.at(muon_index)
        self.parent_pdgid = truth_helpers.getParentPdgIDFromBarcode( event
                                                                   , event.mu_staco_barcode.at(muon_index)
                                                                   )

# ==============================================================================
class Jet(object):
    # ------------------------------------------------------------------------------
    def __init__(self, event, jet_index):
        self.index        = jet_index
        self.pt           = event.jet_AntiKt4TruthJets_pt.at(jet_index)
        self.eta          = event.jet_AntiKt4TruthJets_eta.at(jet_index)
        self.phi          = event.jet_AntiKt4TruthJets_phi.at(jet_index)
        self.E            = event.jet_AntiKt4TruthJets_E.at(jet_index)
        self.theta        = math.copysign( 2*math.atan(math.exp(-abs(self.eta)))
                                         , self.eta
                                         )
        self.px           = self.pt*math.cos(self.phi)
        self.py           = self.pt*math.sin(self.phi)
        self.pz           = self.pt*math.sin(self.theta)
        self.is_b_jet     = truth_helpers.isBJet( self
                                                , event.mc_pdgId
                                                , event.mc_pt
                                                , event.mc_eta
                                                , event.mc_phi
                                                )
