#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import ROOT
import rootlogon
import metaroot

import mt2 as mt2_calc
import TruthRecordHelpers as truth_helpers
import Calculators as calc
import ObjectSelection as object_sel
import ObjectDefs as object_defs

# ==============================================================================
ROOT.gInterpreter.GenerateDictionary("vector<vector<int> >","vector");

# ==============================================================================
# helper definitions
mother_type_list = [ 'st'
                   , 'none'
                   ]
decay_categories = [ 'dc_all'
                   , 'dc_st_st'
                   , 'dc_none'
                   ]
flavor_channels = [ 'fc_all'
                  , 'fc_ee'
                  , 'fc_mm'
                  , 'fc_em'
                  , 'fc_me'
                  , 'fc_none'
                  ]

min_num_leptons = 2
max_num_leptons = 2
max_num_jets = 2

el_prefix  = 'el_'
mu_prefix  = 'mu_staco_'
jet_prefix = 'jet_AntiKt4TruthJets_'

# # ------------------------------------------------------------------------------
# baseline_el_pt_cut  = 10.e3
# baseline_mu_pt_cut  = 10.e3
# baseline_jet_pt_cut = 20.e3
# 
# baseline_el_eta_cut  = 2.47
# baseline_mu_eta_cut  = 2.4
# baseline_jet_eta_cut = 5
# 
# 
# signal_el_pt_cut  = 10.e3
# signal_mu_pt_cut  = 10.e3
# signal_jet_pt_cut = 20.e3
# 
# signal_el_eta_cut  = 2.4
# signal_mu_eta_cut  = 2.4
# signal_jet_eta_cut = 2.4

# # ==============================================================================
# class Electron(object):
#     # ------------------------------------------------------------------------------
#     def __init__(self, event, electron_index):
#         self.index        = electron_index
#         self.pt           = event.el_pt.at(electron_index)
#         self.eta          = event.el_eta.at(electron_index)
#         self.phi          = event.el_phi.at(electron_index)
#         self.E            = event.el_E.at(electron_index)
#         self.charge       = event.el_charge.at(electron_index)
#         self.px           = event.el_px.at(electron_index)
#         self.py           = event.el_py.at(electron_index)
#         self.pz           = event.el_pz.at(electron_index)
#         self.parent_pdgid = truth_helpers.getParentPdgIDFromBarcode( event
#                                                                    , event.el_barcode.at(electron_index)
#                                                                    )
# 
# # ==============================================================================
# class Muon(object):
#     # ------------------------------------------------------------------------------
#     def __init__(self, event, muon_index):
#         self.index        = muon_index
#         self.pt           = event.mu_staco_pt.at(muon_index)
#         self.eta          = event.mu_staco_eta.at(muon_index)
#         self.phi          = event.mu_staco_phi.at(muon_index)
#         self.E            = event.mu_staco_E.at(muon_index)
#         self.charge       = event.mu_staco_charge.at(muon_index)
#         self.px           = event.mu_staco_px.at(muon_index)
#         self.py           = event.mu_staco_py.at(muon_index)
#         self.pz           = event.mu_staco_pz.at(muon_index)
#         self.parent_pdgid = truth_helpers.getParentPdgIDFromBarcode( event
#                                                                    , event.mu_staco_barcode.at(muon_index)
#                                                                    )
# 
# # ==============================================================================
# class Jet(object):
#     # ------------------------------------------------------------------------------
#     def __init__(self, event, jet_index):
#         print 'jet init'
#         self.index        = jet_index
#         self.pt           = event.jet_AntiKt4TruthJets_pt.at(jet_index)
#         self.eta          = event.jet_AntiKt4TruthJets_eta.at(jet_index)
#         self.phi          = event.jet_AntiKt4TruthJets_phi.at(jet_index)
#         self.E            = event.jet_AntiKt4TruthJets_E.at(jet_index)
#         self.theta        = math.copysign( 2*math.atan(math.exp(-abs(self.eta)))
#                                          , self.eta
#                                          )
#         self.px           = self.pt*math.cos(self.phi)
#         self.py           = self.pt*math.sin(self.phi)
#         self.pz           = self.pt*math.sin(self.theta)
#         self.is_b_jet     = truth_helpers.isBJet( self
#                                                 , event.mc_pdgId
#                                                 , event.mc_pt
#                                                 , event.mc_eta
#                                                 , event.mc_phi
#                                                 )

# ==============================================================================
class BMinusLCutFlow(object):
    # ------------------------------------------------------------------------------
    def __init__(self, event):
        self.event = event
        self.valid_cutflow = self.doObjectSelection(
                                                   # verbose = True
                                                   )

    # ------------------------------------------------------------------------------
    def doObjectSelection(self, verbose = False):
        if verbose:
            print '----------------------------------------'
            print 'doing object selection for event: %s' % self.event.EventNumber

        # get baseline objects
        self.baseline = object_sel.getBaselineObjects( self.event, verbose)
        self.num_baseline_el  = len(self.baseline['el'])
        self.num_baseline_mu  = len(self.baseline['mu'])
        self.num_baseline_jet = len(self.baseline['jet'])

        if verbose:
            print '  num baseline electrons: %s' % len(self.baseline['el'])
            print '  num baseline muons:     %s' % len(self.baseline['mu'])
            print '  num baseline jets:      %s' % len(self.baseline['jet'])

        if self.num_baseline_el + self.num_baseline_mu < min_num_leptons:
            return False

        # do overlap removal
        self.overlap_removed = object_sel.doOverlapRemoval(self.baseline, verbose)
        self.num_overlap_removed_el  = len(self.overlap_removed['el'])
        self.num_overlap_removed_mu  = len(self.overlap_removed['mu'])
        self.num_overlap_removed_jet = len(self.overlap_removed['jet'])

        if verbose:
            print '  num overlap removal electrons: %s' % len(self.overlap_removed['el'])
            print '  num overlap removal muons:     %s' % len(self.overlap_removed['mu'])
            print '  num overlap removal jets:      %s' % len(self.overlap_removed['jet'])

        if self.num_overlap_removed_el + self.num_overlap_removed_mu < min_num_leptons:
            return False

        # get signal objets
        self.signal = object_sel.getSignalObjects(self.overlap_removed, verbose)
        self.num_signal_el  = len(self.signal['el'])
        self.num_signal_mu  = len(self.signal['mu'])
        self.num_signal_jet = len(self.signal['jet'])

        if verbose:
            print '  num singal electrons: %s' % len(self.signal['el'])
            print '  num singal muons:     %s' % len(self.signal['mu'])
            print '  num singal jets:      %s' % len(self.signal['jet'])

        if self.num_signal_el + self.num_signal_mu < min_num_leptons:
            return False

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get flavor channel
        self.flavor_channel = getFlavorChannel(self.signal)
        if verbose:
            print 'flavor channel: %s' % self.flavor_channel

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get mll
        if verbose:
            print '    get mll'
        self.mll = calc.getMll(self.signal['el'], self.signal['mu'])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get mbl
        if verbose:
            print '    get mbl'
        self.mbl_list = calc.getMbl( self.signal['el']
                                   , self.signal['mu']
                                   , self.signal['jet']
                                   )
        self.mbl_truth_list = calc.getTruthMbl( self.signal['el']
                                              , self.signal['mu']
                                              , self.signal['jet']
                                              )

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get ptll
        if verbose:
            print '    get ptll'
        self.ptll = calc.getPtll(self.signal['el'], self.signal['mu'])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get emma_mt
        if verbose:
            print '    get emma_mt'
        self.emma_mt = calc.getEmmaMt(self.signal['el'], self.signal['mu'])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get met and met-related variables
        if verbose:
            print '    get met'
        self.getMet()

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose:
            print '    get mt2'
        self.mt2 = mt2_calc.getMT2( self.signal['el']
                                  , self.signal['mu']
                                  , self.event.MET_Truth_NonInt_etx
                                  , self.event.MET_Truth_NonInt_ety
                                  , minv = 0.
                                  , verbose = False
                                  )

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if verbose:
            print '    checking regions'
        self.regions = []
        # if self.isSRSS1():    self.regions.append('srss1')

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # checking decay categories
        self.decay_category = getDecayCategory(self.signal)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        return True

    # ------------------------------------------------------------------------------
    def getMet(self):
        met_etx_int = self.event.MET_Truth_Int_etx
        met_ety_int = self.event.MET_Truth_Int_ety
        for muon in self.signal['mu']:
            met_etx_int -= muon.px
            met_ety_int -= muon.py
        met_int = math.sqrt( met_etx_int*met_etx_int
                           + met_ety_int*met_ety_int
                           )/1000.
        metrel_int = calc.getMetRel( met_etx_int
                                   , met_ety_int
                                   , self.signal['el']
                                   , self.signal['mu']
                                   , self.signal['jet']
                                   )/1000.

        met_etx_noint = self.event.MET_Truth_NonInt_etx
        met_ety_noint = self.event.MET_Truth_NonInt_ety
        met_noint = math.sqrt( met_etx_noint*met_etx_noint
                             + met_ety_noint*met_ety_noint
                             )/1000.
        metrel_noint = calc.getMetRel( met_etx_noint
                                     , met_ety_noint
                                     , self.signal['el']
                                     , self.signal['mu']
                                     , self.signal['jet']
                                     )/1000.

        self.met = { 'int':met_int    , 'rel_int':metrel_int
                   , 'noint':met_noint, 'rel_noint':metrel_noint
                   }

# ------------------------------------------------------------------------------
def getFlavorChannel(signal_objects, verbose = False):
    if verbose:
        print 'getting flavor channel'

    num_el  = len(signal_objects['el'])
    num_mu  = len(signal_objects['mu'])
    num_lep = num_el+num_mu

    # 2-lepton events
    if num_lep == 2:
        if num_el == 2: return 'fc_ee'
        if num_el == 1: return 'fc_em' if signal_objects['el'][0].pt >= signal_objects['mu'][0].pt else 'fc_me'
        if num_el == 0: return 'fc_mm'

        print 'Oh no! Di-lepton event did not fall into any channel!!!'
        print '    num el: %s' % num_el
        print '    num mu: %s' % num_mu
        print '    charge product: %s' % charge_product
        assert False

    # don't classify events with num_lep != 2
    return 'fc_none'

# ------------------------------------------------------------------------------
def getDecayCategory(signal_objects, verbose = False):
    mother_pdgid_list = []
    lepton_pt = []
    for el in signal_objects['el']:
        mother_pdgid_list.append(el.parent_pdgid)
        lepton_pt.append(el.pt)
    for mu in signal_objects['mu']:
        mother_pdgid_list.append(mu.parent_pdgid)
        lepton_pt.append(mu.pt)

    # TODO do pt sorting

    if verbose:
        print mother_pdgid_list
    num_stop_mothers = 0

    for mpl in mother_pdgid_list:
        if mpl is None: return 'dc_none'
        if abs(mpl) == 1000006: num_stop_mothers += 1

    if verbose:
        print '  num stop mothers: %s' % num_stop_mothers

    # classify decays
    dc_classification = 'dc_none'
    if num_stop_mothers == 2: dc_classification = 'dc_st_st'

    if dc_classification not in decay_categories:
        dc_classification = 'dc_none'
    return dc_classification

