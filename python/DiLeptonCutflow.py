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

# ==============================================================================
ROOT.gInterpreter.GenerateDictionary("vector<vector<int> >","vector");

# ==============================================================================
# helper definitions
mother_type_list = [ 'c1'
                   , 'n2'
                   , 'sl'
                   , 'none'
                   ]
decay_categories = [ 'dc_all'
                   , 'dc_c1_c1'
                   , 'dc_sl_sl'
                   , 'dc_c1_sl'
                   , 'dc_sl_c1'
                   , 'dc_n2_n2'
                   , 'dc_n2_c1'
                   , 'dc_c1_n2'
                   , 'dc_n2_sl'
                   , 'dc_sl_n2'
                   , 'dc_none'
                   ]
flavor_channels = [ 'fc_all'
                  , 'fc_ee_os'
                  , 'fc_ee_ss'
                  , 'fc_mm_os'
                  , 'fc_mm_ss'
                  , 'fc_em_os'
                  , 'fc_em_ss'
                  , 'fc_me_os'
                  , 'fc_me_ss'
                  , 'fc_eee'
                  , 'fc_eem'
                  , 'fc_emm'
                  , 'fc_mmm'
                  , 'fc_multi'
                  , 'fc_none'
                  ]

min_num_leptons = 2
max_num_leptons = 4
max_num_jets = 2

el_prefix  = 'el_'
mu_prefix  = 'mu_staco_'
jet_prefix = 'jet_AntiKt4TruthJets_'

# ------------------------------------------------------------------------------
baseline_el_pt_cut  = 10.e3
baseline_mu_pt_cut  = 10.e3
baseline_jet_pt_cut = 20.e3

baseline_el_eta_cut  = 2.47
baseline_mu_eta_cut  = 2.4
baseline_jet_eta_cut = 5


signal_el_pt_cut  = 10.e3
signal_mu_pt_cut  = 10.e3
signal_jet_pt_cut = 20.e3

signal_el_eta_cut  = 2.4
signal_mu_eta_cut  = 2.4
signal_jet_eta_cut = 2.4

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
        self.parent_pdgid = getParentPdgIDFromBarcode( event
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
        self.parent_pdgid = getParentPdgIDFromBarcode( event
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

# ==============================================================================
class EwkCutFlow(object):
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
        self.baseline = getBaselineObjects( self.event, verbose)
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
        self.overlap_removed = doOverlapRemoval(self.baseline, verbose)
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
        self.signal = getSignalObjects(self.overlap_removed, verbose)
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
        self.mll = getMll(self.signal['el'], self.signal['mu'])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get ptll
        if verbose:
            print '    get ptll'
        self.ptll = getPtll(self.signal['el'], self.signal['mu'])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # get emma_mt
        if verbose:
            print '    get emma_mt'
        self.emma_mt = getEmmaMt(self.signal['el'], self.signal['mu'])

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
        if self.isSRSS1():    self.regions.append('srss1')
        if self.isSRSS2():    self.regions.append('srss2')
        if self.isSRSS3():    self.regions.append('srss3')
        if self.isSRSS4():    self.regions.append('srss4')
        if self.isSRSS5():    self.regions.append('srss5')
        if self.isSROSMT2a(): self.regions.append('srmt2a')
        if self.isSROSMT2b(): self.regions.append('srmt2b')
        if self.isSROSMT2c(): self.regions.append('srmt2c')

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
        metrel_int = getMetRel( met_etx_int
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
        metrel_noint = getMetRel( met_etx_noint
                                , met_ety_noint
                                , self.signal['el']
                                , self.signal['mu']
                                , self.signal['jet']
                                )/1000.

        self.met = { 'int':met_int    , 'rel_int':metrel_int
                   , 'noint':met_noint, 'rel_noint':metrel_noint
                   }

    # ------------------------------------------------------------------------------
    def isSRSS1(self):
        # num_el = len(self.signal['el'])
        # num_mu = len(self.signal['mu'])
        # num_lep = num_el+num_mu

        if "ss" not in self.flavor_channel: return False
        if len(self.signal['jet']) == 0:    return False
        if self.met['rel_noint'] < 50.:     return False
        if self.emma_mt/1000. > 40:         return False

        return True

    # ------------------------------------------------------------------------------
    def isSRSS2(self):
        # len(num_el = self.signal['el'])
        # len(num_mu = self.signal['mu'])
        # num_lep = num_el+num_mu

        if "ss" not in self.flavor_channel: return False
        if len(self.signal['jet']) == 0:    return False
        if self.met['rel_noint'] < 50.:     return False
        if self.mll/1000. > 100:            return False
        if self.ptll/1000. > 100:           return False

        return True

    # ------------------------------------------------------------------------------
    def isSRSS3(self):
        if "ss" not in self.flavor_channel: return False
        if len(self.signal['jet']) == 0:    return False
        if self.met['rel_noint'] < 50.:     return False
        if self.mll/1000. > 75:             return False
        if self.ptll/1000. > 75:            return False

        return True

    # ------------------------------------------------------------------------------
    def isSRSS4(self):
        if "ss" not in self.flavor_channel: return False
        if self.met['noint'] < 200.:        return False

        return True

    # ------------------------------------------------------------------------------
    def isSRSS5(self):
        if "ss" not in self.flavor_channel: return False
        if self.met['rel_noint'] < 200.:    return False

        return True

    # ------------------------------------------------------------------------------
    def isSROSMT2a(self):
        if 'os' not in self.flavor_channel: return False
        if self.met['rel_noint'] < 40.:     return False
        if self.mt2/1000. < 90.:            return False

        #Z veto for ee/mm
        if 'ee' in self.flavor_channel or 'mm' in self.flavor_channel:
            if abs(self.mll/1000. - 91.2) < 10: return False

        return True

    # ------------------------------------------------------------------------------
    def isSROSMT2b(self):
        if 'os' not in self.flavor_channel: return False
        if self.met['rel_noint'] < 40.:     return False
        if self.mt2/1000. < 120.:           return False

        #Z veto for ee/mm
        if 'ee' in self.flavor_channel or 'mm' in self.flavor_channel:
            if abs(self.mll/1000. - 91.2) < 10: return False

        return True

    # ------------------------------------------------------------------------------
    def isSROSMT2c(self):
        if 'os' not in self.flavor_channel: return False
        if self.met['rel_noint'] < 40.:     return False
        if self.mt2/1000. < 150.:           return False

        #Z veto for ee/mm
        if 'ee' in self.flavor_channel or 'mm' in self.flavor_channel:
            if abs(self.mll/1000. - 91.2) < 10: return False

        return True

# ------------------------------------------------------------------------------
def getBaselineObjects( event
                      , verbose = False
                      ):
    baseline_el  = []
    baseline_mu  = []
    baseline_jet = []

    # get baseline electrons
    el_index_order = getPtSortedIndices(event.el_n, event.el_pt)
    for el_index in el_index_order:
        this_el = Electron(event, el_index)
        if this_el.pt < baseline_el_pt_cut:
            if verbose:
                print '  electron %d failed pt cut (%f < %f)' % ( el_index
                                                                , this_el.pt
                                                                , baseline_el_pt_cut
                                                                )
            continue
        if abs(this_el.eta) > baseline_el_eta_cut:
            if verbose:
                print '  electron %d failed eta cut (%f < %f)' % ( el_index
                                                                 , this_el.eta
                                                                 , baseline_el_eta_cut
                                                                 )
            continue
        baseline_el.append(this_el)

    # get baseline muons
    mu_index_order = getPtSortedIndices(event.mu_staco_n, event.mu_staco_pt)
    for mu_index in mu_index_order:
        this_mu = Muon(event, mu_index)
        if this_mu.pt < baseline_mu_pt_cut:
            if verbose:
                print '  muon %d failed pt cut (%f < %f)' % ( mu_index
                                                            , this_mu.pt
                                                            , baseline_mu_pt_cut
                                                            )
            continue
        if abs(this_mu.eta) > baseline_mu_eta_cut:
            if verbose:
                print '  muon %d failed eta cut (%f < %f)' % ( mu_index
                                                             , this_mu.eta
                                                             , baseline_mu_eta_cut
                                                             )
            continue
        baseline_mu.append(this_mu)

    # get baseline jets
    jet_index_order = getPtSortedIndices( event.jet_AntiKt4TruthJets_n
                                        , event.jet_AntiKt4TruthJets_pt
                                        )
    for jet_index in jet_index_order:
        this_jet = Jet(event, jet_index)
        if this_jet.pt < baseline_jet_pt_cut:
            if verbose:
                print '  jet %d failed pt cut (%f < %f)' % ( jet_index
                                                           , this_jet.pt
                                                           , baseline_jet_pt_cut
                                                           )
            continue
        if abs(this_jet.eta) > baseline_jet_eta_cut:
            if verbose:
                print '  jet %d failed eta cut (%f < %f)' % ( jet_index
                                                            , this_jet.eta
                                                            , baseline_jet_eta_cut
                                                            )
            continue
        baseline_jet.append(this_jet)

    return {'el':baseline_el, 'mu':baseline_mu, 'jet':baseline_jet}

# ------------------------------------------------------------------------------
def getPtSortedIndices(n, pt_list):
    if n == 0: return []
    pt_order, index_order = (list(t) for t in zip(*sorted(zip(pt_list, range(n)), reverse=True)))
    return index_order

# ------------------------------------------------------------------------------
def doOverlapRemoval(baseline, verbose = False):
    overlap_removed_el  = baseline['el']
    overlap_removed_mu  = baseline['mu']
    overlap_removed_jet = baseline['jet']

    if verbose:
        print 'before overlap removal:'
        print '    el:  %d - %s' % (len(overlap_removed_el) , len(overlap_removed_el))
        print '    mu:  %d - %s' % (len(overlap_removed_mu) , len(overlap_removed_mu))
        print '    jet: %d - %s' % (len(overlap_removed_jet), len(overlap_removed_jet))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # do e-e overlap removal
    el_to_remove_ee = []
    dr_cut_ee = 0.1
    for i_ee in xrange(len(overlap_removed_el)):
        eta_i_ee = overlap_removed_el[i_ee].eta
        phi_i_ee = overlap_removed_el[i_ee].eta
        for j_ee in xrange(i_ee+1, len(overlap_removed_el)):
            eta_j_ee = overlap_removed_el[j_ee].eta
            phi_j_ee = overlap_removed_el[j_ee].eta
            if isOverlap(dr_cut_ee, eta_i_ee, phi_i_ee, eta_j_ee, phi_j_ee):
                pt_i_ee = overlap_removed_el[i_ee].pt
                pt_j_ee = overlap_removed_el[j_ee].pt
                this_removal_ee = j_ee if pt_i_ee > pt_j_ee else i_ee

                if verbose:
                    print 'e-e overlap: %d - %d' % (i_ee, j_ee)
                    print '  e %d -- pT: %f - eta: %f - phi: %f' % ( i_ee
                                                                   , pt_i_ee
                                                                   , eta_i_ee
                                                                   , phi_i_ee
                                                                   )
                    print '  e %d -- pT: %f - eta: %f - phi: %f' % ( j_ee
                                                                   , pt_j_ee
                                                                   , eta_j_ee
                                                                   , phi_j_ee
                                                                   )
                    print '  remove: %d' % (this_removal_ee)
                el_to_remove_ee.append(this_removal_ee)
    overlap_removed_el = removeElements(overlap_removed_el, el_to_remove_ee)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # do e-jet overlap removal
    jet_to_remove_ej = []
    dr_cut_ej = 0.2
    for jet_it_ej in xrange(len(overlap_removed_jet)):
        eta_jet_ej = overlap_removed_jet[jet_it_ej].eta
        phi_jet_ej = overlap_removed_jet[jet_it_ej].phi
        for el_it_ej in xrange(len(overlap_removed_el)):
            eta_el_ej = overlap_removed_el[el_it_ej].eta
            phi_el_ej = overlap_removed_el[el_it_ej].phi
            if isOverlap(dr_cut_ej, eta_jet_ej, phi_jet_ej, eta_el_ej, phi_el_ej):
                if verbose:
                    pt_jet_ej = overlap_removed_jet[jet_it_ej].pt
                    pt_el_ej  = overlap_removed_el[el_it_ej].pt

                    print 'e-jet overlap: el %d - jet %d' % (el_it_ej, jet_it_ej)
                    print '  el %d -- pT: %f - eta: %f - phi: %f' % ( el_it_ej
                                                                    , pt_el_ej
                                                                    , eta_el_ej
                                                                    , phi_el_ej
                                                                    )
                    print '  jet %d -- pT: %f - eta: %f - phi: %f' % ( jet_it_ej
                                                                     , pt_jet_ej
                                                                     , eta_jet_ej
                                                                     , phi_jet_ej
                                                                     )
                    print '  remove: jet %d' % jet_it_ej
                jet_to_remove_ej.append(jet_it_ej)
                break # no need to keep looping over electrons
    overlap_removed_jet = removeElements(overlap_removed_jet, jet_to_remove_ej)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # do jet-e overlap removal
    el_to_remove_je = []
    dr_cut_je = 0.4
    for el_it_je in xrange(len(overlap_removed_el)):
        eta_el_je = overlap_removed_el[el_it_je].eta
        phi_el_je = overlap_removed_el[el_it_je].phi
        for jet_it_je in xrange(len(overlap_removed_jet)):
            eta_jet_je = overlap_removed_jet[jet_it_je].eta
            phi_jet_je = overlap_removed_jet[jet_it_je].phi
            if isOverlap(dr_cut_je, eta_el_je, phi_el_je, eta_jet_je, phi_jet_je):
                if verbose:
                    pt_el_je  = overlap_removed_el[el_it_je].pt
                    pt_jet_je = overlap_removed_jet[jet_it_je].pt

                    print 'jet-e overlap: jet %d - el %d' % (jet_it_je, el_it_je)
                    print '  jet %d -- pT: %f - eta: %f - phi: %f' % ( jet_it_je
                                                                     , pt_jet_je
                                                                     , eta_jet_je
                                                                     , phi_jet_je
                                                                     )
                    print '  el %d -- pT: %f - eta: %f - phi: %f' % ( el_it_je
                                                                    , pt_el_je
                                                                    , eta_el_je
                                                                    , phi_el_je
                                                                    )
                    print '  remove: el %d' % el_it_je
                el_to_remove_je.append(el_it_je)
                break # no need to keep looping over jets
    overlap_removed_el = removeElements(overlap_removed_el, el_to_remove_je)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # do jet-mu overlap removal
    mu_to_remove_jm = []
    dr_cut_jm = 0.4
    for mu_it_jm in xrange(len(overlap_removed_mu)):
        eta_mu_jm = overlap_removed_mu[mu_it_jm].eta
        phi_mu_jm = overlap_removed_mu[mu_it_jm].phi
        for jet_it_jm in xrange(len(overlap_removed_jet)):
            eta_jet_jm = overlap_removed_jet[jet_it_jm].eta
            phi_jet_jm = overlap_removed_jet[jet_it_jm].phi
            if isOverlap(dr_cut_jm, eta_mu_jm, phi_mu_jm, eta_jet_jm, phi_jet_jm):
                if verbose:
                    pt_mu_jm  = overlap_removed_mu[mu_it_jm].pt
                    pt_jet_jm = overlap_removed_jet[jet_it_jm].pt

                    print 'jet-e overlap: jet %d - el %d' % (jet_it_jm, mu_it_jm)
                    print '  jet %d -- pT: %f - eta: %f - phi: %f' % ( jet_it_jm
                                                                     , pt_jet_jm
                                                                     , eta_jet_jm
                                                                     , phi_jet_jm
                                                                     )
                    print '  mu %d -- pT: %f - eta: %f - phi: %f' % ( mu_it_jm
                                                                    , pt_mu_jm
                                                                    , eta_mu_jm
                                                                    , phi_mu_jm
                                                                    )
                    print '  remove: el %d' % mu_it_jm
                mu_to_remove_jm.append(mu_it_jm)
                break # no need to keep looping over jets
    overlap_removed_mu = removeElements(overlap_removed_mu, mu_to_remove_jm)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # do e-mu overlap removal
    el_to_remove_em = []
    mu_to_remove_em = []
    dr_cut_em = 0.1
    for el_it_em in xrange(len(overlap_removed_el)):
        eta_el_em = overlap_removed_el[el_it_em].eta
        phi_el_em = overlap_removed_el[el_it_em].phi
        for mu_it_em in xrange(len(overlap_removed_mu)):
            eta_mu_em = overlap_removed_mu[mu_it_em].eta
            phi_mu_em = overlap_removed_mu[mu_it_em].phi
            if isOverlap(dr_cut_em, eta_el_em, phi_el_em, eta_mu_em, phi_mu_em):
                if verbose:
                    pt_el_em = overlap_removed_el[el_it_em].pt
                    pt_mu_em  = overlap_removed_mu[mu_it_em].pt

                    print 'e-mu overlap: el %d - mu %d' % (el_it_em, mu_it_em)
                    print '  el %d -- pT: %f - eta: %f - phi: %f' % ( el_it_em
                                                                    , pt_el_em
                                                                    , eta_el_em
                                                                    , phi_el_em
                                                                    )
                    print '  mu %d -- pT: %f - eta: %f - phi: %f' % ( mu_it_em
                                                                    , pt_mu_em
                                                                    , eta_mu_em
                                                                    , phi_mu_em
                                                                    )
                    print '  remove: el %d' % el_it_em
                    print '  remove: mu %d' % mu_it_em
                el_to_remove_em.append(el_it_em)
                mu_to_remove_em.append(mu_it_em)
    overlap_removed_el = removeElements(overlap_removed_el, el_to_remove_em)
    overlap_removed_mu = removeElements(overlap_removed_mu, mu_to_remove_em)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # do mu-mu overlap removal
    mu_to_remove_mm = []
    dr_cut_mm = 0.05
    for i_mm in xrange(len(overlap_removed_mu)):
        eta_i_mm = overlap_removed_mu[i_mm].eta
        phi_i_mm = overlap_removed_mu[i_mm].phi
        for j_mm in xrange(i_mm+1, len(overlap_removed_mu)):
            eta_j_mm = overlap_removed_mu[j_mm].eta
            phi_j_mm = overlap_removed_mu[j_mm].phi
            if isOverlap(dr_cut_mm, eta_i_mm, phi_i_mm, eta_j_mm, phi_j_mm):
                if verbose:
                    pt_i_mm = overlap_removed_mu[i_mm].pt
                    pt_j_mm = overlap_removed_mu[j_mm].pt

                    print 'mu-mu overlap: %d - %d' % (i_mm, j_mm)
                    print '  mu %d -- pT: %f - eta: %f - phi: %f' % ( i_mm
                                                                    , pt_i_mm
                                                                    , eta_i_mm
                                                                    , phi_i_mm
                                                                    )
                    print '  mu %d -- pT: %f - eta: %f - phi: %f' % ( j_mm
                                                                    , pt_j_mm
                                                                    , eta_j_mm
                                                                    , phi_j_mm
                                                                    )
                    print '  remove: %d' % (i_mm)
                    print '  remove: %d' % (j_mm)
                mu_to_remove_mm.append(i_mm)
                mu_to_remove_mm.append(j_mm)
    overlap_removed_mu = removeElements(overlap_removed_mu, mu_to_remove_mm)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if verbose:
        print 'after overlap removal:'
        print '    el:  %d - %s' % (len(overlap_removed_el) , len(overlap_removed_el))
        print '    mu:  %d - %s' % (len(overlap_removed_mu) , len(overlap_removed_mu))
        print '    jet: %d - %s' % (len(overlap_removed_jet), len(overlap_removed_jet))
    return { 'el':overlap_removed_el
           , 'mu':overlap_removed_mu
           , 'jet':overlap_removed_jet
           }

# ------------------------------------------------------------------------------
def removeElements(particle_lists, to_remove):
    to_remove = list(set(to_remove))
    to_remove.sort(reverse=True)
    for tr in to_remove:
        del particle_lists[tr]
        #~~~~ # particle_lists['num'] -= 1
        #~~~~ for pl in particle_lists:
        #~~~~     if isinstance(particle_lists[pl], list):
        #~~~~         del particle_lists[pl][tr]
    return particle_lists

# ------------------------------------------------------------------------------
def isOverlap(cut, eta_1, phi_1, eta_2, phi_2):
    deta = eta_1-eta_2
    dphi = abs(phi_1-phi_2)
    while dphi > 3.14159:
        dphi -= 3.14159
    dr = math.sqrt(deta*deta + dphi*dphi)
    return (dr < cut)

# ------------------------------------------------------------------------------
def getSignalObjects( baseline_objects
                    , verbose = False
                    ):
    signal_el  = baseline_objects['el']
    signal_mu  = baseline_objects['mu']
    signal_jet = baseline_objects['jet']

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get signal electrons
    if verbose:
        print '    get signal electrons'
    to_remove_el = []
    for el_it in xrange(len(signal_el)):
        el_pt  = signal_el[el_it].pt
        el_eta = signal_el[el_it].eta

        if el_pt < signal_el_pt_cut:
            if verbose:
                print '  electron %d failed pt cut (%f < %f)' % ( el_it
                                                                , el_pt
                                                                , signal_el_pt_cut
                                                                )
            to_remove_el.append(el_it)
            continue
        if abs(el_eta) > signal_el_eta_cut:
            if verbose:
                print '  electron %d failed eta cut (%f < %f)' % ( el_it
                                                                 , el_eta
                                                                 , signal_el_eta_cut
                                                                 )
            to_remove_el.append(el_it)
            continue
    if verbose:
        print '   removing electrons from signal selection: %s' % to_remove_el
    removeElements(signal_el, to_remove_el)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get signal muons
    if verbose:
        print '    get signal muons'
    to_remove_mu = []
    for mu_it in xrange(len(signal_mu)):
        mu_pt  = signal_mu[mu_it].pt
        mu_eta = signal_mu[mu_it].eta

        if mu_pt < signal_el_pt_cut:
            if verbose:
                print '  muon %d failed pt cut (%f < %f)' % (mu_it, mu_pt, lep_pt_cut)
            to_remove_mu.append(mu_it)
            continue
        if abs(mu_eta) > signal_mu_eta_cut:
            if verbose:
                print '  muon %d failed eta cut (|%f| < %f)' % (mu_it, mu_eta, lep_eta_cut)
            to_remove_mu.append(mu_it)
            continue
    if verbose:
        print '    removing muons from signal selection: %s' % to_remove_mu
    removeElements(signal_mu, to_remove_mu)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get signal jets
    if verbose:
        print '    get signal jets'
    to_remove_jet = []
    for jet_it in xrange(len(signal_jet)):
        jet_pt  = signal_jet[jet_it].pt
        jet_eta = signal_jet[jet_it].eta

        if jet_pt < signal_jet_pt_cut:
            if verbose:
                print '  jet %d failed pt cut (%f < %f)' % ( jet_it
                                                           , jet_pt
                                                           , signal_jet_pt_cut
                                                           )
            to_remove_jet.append(jet_it)
            continue
        if abs(jet_eta) > signal_jet_eta_cut:
            if verbose:
                print '  jet %d failed eta cut (|%f| < %f)' % ( jet_eta
                                                              , jet_eta
                                                              , signal_jet_eta_cut
                                                              )
            to_remove_jet.append(jet_it)
            continue
    if verbose:
        print '    removing jets from signal selection: %s' % to_remove_jet
    removeElements(signal_jet, to_remove_jet)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return { 'el':signal_el
           , 'mu':signal_mu
           , 'jet':signal_jet
           }

# ------------------------------------------------------------------------------
def getMetRel( met_etx, met_ety , signal_el, signal_mu, signal_jet):
    met_phi = math.atan2(met_ety, met_etx)

    min_dphi = 9999999999

    for el_phi in [el.phi for el in signal_el]:
        dphi = abs(met_phi - el_phi)
        while dphi > 3.14159: dphi -= 3.14159
        if dphi < min_dphi: min_dphi = dphi

    for mu_phi in [mu.phi for mu in signal_mu]:
        dphi = abs(met_phi - mu_phi)
        while dphi > 3.14159: dphi -= 3.14159
        if dphi < min_dphi: min_dphi = dphi

    for jet_phi in [jet.phi for jet in signal_jet]:
        dphi = abs(met_phi - jet_phi)
        while dphi > 3.14159: dphi -= 3.14159
        if dphi < min_dphi: min_dphi = dphi

    met_rel = math.sqrt(met_etx*met_etx + met_ety*met_ety)
    if min_dphi < 3.14159:
        met_rel *= math.cos(min_dphi)

    return met_rel

# ------------------------------------------------------------------------------
def getFlavorChannel(signal_objects, verbose = False):
    if verbose:
        print 'getting flavor channel'

    num_el  = len(signal_objects['el'])
    num_mu  = len(signal_objects['mu'])
    num_lep = num_el+num_mu

    # 2-lepton events
    if num_lep == 2:
        # get charge product
        charge_product = 1
        for el in signal_objects['el']:
            charge_product *= el.charge
        for mu in signal_objects['mu']:
            charge_product *= mu.charge

        if num_el == 2 and charge_product < 0: return 'fc_ee_os'
        if num_el == 1 and charge_product < 0:
            if signal_objects['el'][0].pt >= signal_objects['mu'][0].pt:
                return 'fc_em_os'
            return 'fc_me_os'
        if num_el == 0 and charge_product < 0: return 'fc_mm_os'
        if num_el == 2 and charge_product > 0: return 'fc_ee_ss'
        if num_el == 1 and charge_product > 0:
            if signal_objects['el'][0].pt >= signal_objects['mu'][0].pt:
                return 'fc_em_ss'
            return 'fc_me_ss'
        if num_el == 0 and charge_product > 0: return 'fc_mm_ss'

        print 'Oh no! Di-lepton event did not fall into any channel!!!'
        print '    num el: %s' % num_el
        print '    num mu: %s' % num_mu
        print '    charge product: %s' % charge_product
        assert False

    # 3-lepton events
    if num_lep == 3:
        if num_el == 3: return 'fc_eee'
        if num_el == 2: return 'fc_eem'
        if num_el == 1: return 'fc_emm'
        if num_el == 0: return 'fc_mmm'

        print 'Oh no! Tri-lepton event did not fall into any channel!!!'
        assert False

    # multi-lepton events
    if num_lep >= 4:
        return 'fc_multi'

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
    num_sl_mothers = 0
    num_n2_mothers = 0
    num_c1_mothers = 0

    for mpl in mother_pdgid_list:
        if mpl is None: return 'dc_none'
        if abs(mpl) >= 1000011 and abs(mpl) <= 1000016: num_sl_mothers += 1
        if abs(mpl) == 1000023: num_n2_mothers += 1
        if abs(mpl) == 1000024: num_c1_mothers += 1

    if verbose:
        print '  num sl mothers: %s' % num_sl_mothers
        print '  num n2 mothers: %s' % num_n2_mothers
        print '  num c1 mothers: %s' % num_c1_mothers

    if num_sl_mothers == 2 and num_n2_mothers == 0 and num_c1_mothers == 0: return 'dc_sl_sl'
    if num_sl_mothers == 0 and num_n2_mothers == 2 and num_c1_mothers == 0: return 'dc_n2_n2'
    if num_sl_mothers == 0 and num_n2_mothers == 0 and num_c1_mothers == 2: return 'dc_c1_c1'
    if num_sl_mothers == 1 and num_n2_mothers == 0 and num_c1_mothers == 1:
        if abs(mother_pdgid_list[0]) == 1000024: return 'dc_c1_sl'
        else:                                    return 'dc_sl_c1'
    if num_sl_mothers == 0 and num_n2_mothers == 1 and num_c1_mothers == 1:
        if abs(mother_pdgid_list[0]) == 1000024: return 'dc_c1_n2'
        else:                                    return 'dc_n2_c1'
    if num_sl_mothers == 1 and num_n2_mothers == 1 and num_c1_mothers == 0:
        if abs(mother_pdgid_list[0]) == 1000023: return 'dc_n2_sl'
        else:                                    return 'dc_sl_n2'

    return 'dc_none'

# ------------------------------------------------------------------------------
def getMll(el_list, mu_list):
    px = 0.
    py = 0.
    pz = 0.
    e = 0.

    for el in el_list:
        px += el.px
        py += el.py
        pz += el.pz
        e  += el.E
    for mu in mu_list:
        px += mu.px
        py += mu.py
        pz += mu.pz
        e  += mu.E

    m2 = (e*e - px*px - py*py - pz*pz)
    return math.copysign(math.sqrt(abs(m2)), m2)

# ------------------------------------------------------------------------------
def getPtll(el_list, mu_list):
    px = 0.
    py = 0.

    for el in el_list:
        px += el.px
        py += el.py
    for mu in mu_list:
        px += mu.px
        py += mu.py

    ptll2 = (px*px + py*py)
    return math.copysign(math.sqrt(abs(ptll2)), ptll2)

# ------------------------------------------------------------------------------
def getEmmaMt(el_list, mu_list):
    mll = getMll(el_list, mu_list)
    ptll = getPtll(el_list, mu_list)
    return math.sqrt(mll*mll + ptll*ptll)

# ------------------------------------------------------------------------------
def getParentPdgID(event, particle_index):
    original_pdgid = event.mc_pdgId.at(particle_index)
    current_index  = particle_index
    mother = None
    # prevent infinite loops!!!
    iteration = 0
    MAX_ITERATION = 100

    # Loop over truth particles looking for mother which is not the same particle
    # print 'searching for parent'
    while mother is None and iteration < MAX_ITERATION:
        iteration += 1

        # Loop over this particle's parents
        # print '  looping over parents (%s)' % event.mc_parent_index.at(current_index).size()
        for ll in xrange(event.mc_parent_index.at(current_index).size()):
            parent_index  = event.mc_parent_index.at(current_index).at(ll);
            parent_pdgid  = event.mc_pdgId.at(parent_index);
            parent_status = event.mc_status.at(parent_index);

            # if parent pdgid = original pdgid, the particle is the result of a scatter. take one step back...
            if parent_pdgid == original_pdgid:
                current_index = parent_index
                break
            # else we found the mother!!!
            else:
                mother = parent_pdgid

    # print 'parent pdgid: %s' % mother
    return mother

# ------------------------------------------------------------------------------
def getParentPdgIDFromBarcode(event, barcode):
    mc_index = None
    for index in xrange(event.mc_n):
        if event.mc_barcode.at(index) == barcode:
            mc_index = index
            # print 'found index matching this barcode -- index: %s  pdgid: %s' % (mc_index, event.mc_pdgId.at(mc_index))
            return getParentPdgID(event, mc_index)
    return None
