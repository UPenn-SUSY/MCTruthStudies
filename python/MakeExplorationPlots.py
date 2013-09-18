#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import ROOT
import rootlogon
import metaroot

# ==============================================================================
# helper definitions
flavor_channels = [ 'ee_os'
                  , 'ee_ss'
                  , 'mm_os'
                  , 'mm_ss'
                  , 'em_os'
                  , 'em_ss'
                  , 'me_os'
                  , 'me_ss'
                  , 'eee'
                  , 'eem'
                  , 'emm'
                  , 'mmm'
                  , 'multi'
                  , 'none'
                  ]
max_num_leptons = 4
max_num_jets = 2

el_prefix  = 'el_'
mu_prefix  = 'mu_staco_'
jet_prefix = 'jet_AntiKt4TruthJets_'

# ------------------------------------------------------------------------------
def doObjectSelection( event
                     , lep_pt_cut
                     , lep_eta_cut
                     , jet_pt_cut
                     , jet_eta_cut
                     , verbose = False
                     ):
    if verbose:
        print '----------------------------------------'
        print 'doing object selection for event: %s' % event.EventNumber

    # get baseline objects
    baseline = getBaselineObjects( event
                                 , lep_pt_cut
                                 , jet_pt_cut
                                 , verbose
                                 )

    overlap_removed = doOverlapRemoval(baseline, verbose)

    signal = getSignalObjects( event
                             , overlap_removed
                             , lep_pt_cut
                             , lep_eta_cut
                             , jet_pt_cut
                             , jet_eta_cut
                             , verbose
                             )

    return signal

# ------------------------------------------------------------------------------
def getBaselineObjects( event
                      , lep_pt_cut
                      , jet_pt_cut
                      , verbose = False
                      ):
    baseline_el = { 'num':0
                  , 'index':[]
                  , 'pt':[]
                  , 'eta':[]
                  , 'phi':[]
                  , 'charge':[]
                  }

    baseline_mu = { 'num':0
                  , 'index':[]
                  , 'pt':[]
                  , 'eta':[]
                  , 'phi':[]
                  , 'charge':[]
                  }

    baseline_jet = { 'num':0
                   , 'index':[]
                   , 'pt':[]
                   , 'eta':[]
                   , 'phi':[]
                   }

    # Get baseline electrons
    el_index_order = getPtSortedIndices(event.el_n, event.el_pt)
    for el_index in el_index_order:
        el_pt = event.el_pt.at(el_index)

        if el_pt < lep_pt_cut:
            if verbose:
                print '  electron %d failed pt cut (%f < %f)' % (el_index, el_pt, lep_pt_cut)
            continue
        baseline_el['num'] += 1
        baseline_el['index'].append(el_index)
        baseline_el['pt'].append(     el_pt)
        baseline_el['eta'].append(    event.el_eta.at(el_index))
        baseline_el['phi'].append(    event.el_phi.at(el_index))
        baseline_el['charge'].append( event.el_charge.at(el_index))

    # Get baseline muons
    mu_index_order = getPtSortedIndices(event.mu_staco_n, event.mu_staco_pt)
    for mu_index in mu_index_order:
        mu_pt = event.mu_staco_pt.at(mu_index)

        if mu_pt < lep_pt_cut:
            if verbose:
                print '  muon %d failed pt cut (%f < %f)' % (mu_index, mu_pt, lep_pt_cut)
            continue
        baseline_mu['num'] += 1
        baseline_mu['index'].append(mu_index)
        baseline_mu['pt'].append(     mu_pt)
        baseline_mu['eta'].append(    event.mu_staco_eta.at(mu_index))
        baseline_mu['phi'].append(    event.mu_staco_phi.at(mu_index))
        baseline_mu['charge'].append( event.mu_staco_charge.at(mu_index))

    # get baseline jets
    jet_index_order = getPtSortedIndices(event.jet_AntiKt4TruthJets_n, event.jet_AntiKt4TruthJets_pt)
    for jet_index in jet_index_order:
        jet_pt     = event.jet_AntiKt4TruthJets_pt.at(jet_index)

        if jet_pt < jet_pt_cut:
            if verbose:
                print '  jet %d failed pt cut (%f < %f)' % (jet_index, jet_pt, jet_pt_cut)
            continue
        baseline_jet['num'] += 1
        baseline_jet['index'].append(jet_index)
        baseline_jet['pt'].append(   jet_pt)
        baseline_jet['eta'].append(  event.jet_AntiKt4TruthJets_eta.at(jet_index))
        baseline_jet['phi'].append(  event.jet_AntiKt4TruthJets_phi.at(jet_index))

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
        print '    el:  %d - %s' % (overlap_removed_el['num'] , overlap_removed_el['index'])
        print '    mu:  %d - %s' % (overlap_removed_mu['num'] , overlap_removed_mu['index'])
        print '    jet: %d - %s' % (overlap_removed_jet['num'], overlap_removed_jet['index'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # do e-e overlap removal
    el_to_remove_ee = []
    dr_cut_ee = 0.1
    for i_ee in xrange(overlap_removed_el['num']):
        eta_i_ee = overlap_removed_el['eta'][i_ee]
        phi_i_ee = overlap_removed_el['eta'][i_ee]
        for j_ee in xrange(i_ee+1, overlap_removed_el['num']):
            eta_j_ee = overlap_removed_el['eta'][j_ee]
            phi_j_ee = overlap_removed_el['eta'][j_ee]
            if isOverlap(dr_cut_ee, eta_i_ee, phi_i_ee, eta_j_ee, phi_j_ee):
                pt_i_ee = overlap_removed_el['pt'][i_ee]
                pt_j_ee = overlap_removed_el['pt'][j_ee]
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
    for jet_it_ej in xrange(overlap_removed_jet['num']):
        eta_jet_ej = overlap_removed_jet['eta'][jet_it_ej]
        phi_jet_ej = overlap_removed_jet['phi'][jet_it_ej]
        for el_it_ej in xrange(overlap_removed_el['num']):
            eta_el_ej = overlap_removed_el['eta'][el_it_ej]
            phi_el_ej = overlap_removed_el['phi'][el_it_ej]
            if isOverlap(dr_cut_ej, eta_jet_ej, phi_jet_ej, eta_el_ej, phi_el_ej):
                if verbose:
                    pt_jet_ej = overlap_removed_jet['pt'][jet_it_ej]
                    pt_el_ej  = overlap_removed_el['pt'][el_it_ej]

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
    for el_it_je in xrange(overlap_removed_el['num']):
        eta_el_je = overlap_removed_el['eta'][el_it_je]
        phi_el_je = overlap_removed_el['phi'][el_it_je]
        for jet_it_je in xrange(overlap_removed_jet['num']):
            eta_jet_je = overlap_removed_jet['eta'][jet_it_je]
            phi_jet_je = overlap_removed_jet['phi'][jet_it_je]
            if isOverlap(dr_cut_je, eta_el_je, phi_el_je, eta_jet_je, phi_jet_je):
                if verbose:
                    pt_el_je  = overlap_removed_el['pt'][el_it_je]
                    pt_jet_je = overlap_removed_jet['pt'][jet_it_je]

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
    for mu_it_jm in xrange(overlap_removed_mu['num']):
        eta_mu_jm = overlap_removed_mu['eta'][mu_it_jm]
        phi_mu_jm = overlap_removed_mu['phi'][mu_it_jm]
        for jet_it_jm in xrange(overlap_removed_jet['num']):
            eta_jet_jm = overlap_removed_jet['eta'][jet_it_jm]
            phi_jet_jm = overlap_removed_jet['phi'][jet_it_jm]
            if isOverlap(dr_cut_jm, eta_mu_jm, phi_mu_jm, eta_jet_jm, phi_jet_jm):
                if verbose:
                    pt_mu_jm  = overlap_removed_mu['pt'][mu_it_jm]
                    pt_jet_jm = overlap_removed_jet['pt'][jet_it_jm]

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
    for el_it_em in xrange(overlap_removed_el['num']):
        eta_el_em = overlap_removed_el['eta'][el_it_em]
        phi_el_em = overlap_removed_el['phi'][el_it_em]
        for mu_it_em in xrange(overlap_removed_mu['num']):
            eta_mu_em = overlap_removed_mu['eta'][mu_it_em]
            phi_mu_em = overlap_removed_mu['phi'][mu_it_em]
            if isOverlap(dr_cut_em, eta_el_em, phi_el_em, eta_mu_em, phi_mu_em):
                if verbose:
                    pt_el_em = overlap_removed_el['pt'][el_it_em]
                    pt_mu_em  = overlap_removed_mu['pt'][mu_it_em]

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
    for i_mm in xrange(overlap_removed_mu['num']):
        eta_i_mm = overlap_removed_mu['eta'][i_mm]
        phi_i_mm = overlap_removed_mu['phi'][i_mm]
        for j_mm in xrange(i_mm+1, overlap_removed_mu['num']):
            eta_j_mm = overlap_removed_mu['eta'][j_mm]
            phi_j_mm = overlap_removed_mu['phi'][j_mm]
            if isOverlap(dr_cut_mm, eta_i_mm, phi_i_mm, eta_j_mm, phi_j_mm):
                if verbose:
                    pt_i_mm = overlap_removed_mu['pt'][i_mm]
                    pt_j_mm = overlap_removed_mu['pt'][j_mm]

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
        print '    el:  %d - %s' % (overlap_removed_el['num'] , overlap_removed_el['index'])
        print '    mu:  %d - %s' % (overlap_removed_mu['num'] , overlap_removed_mu['index'])
        print '    jet: %d - %s' % (overlap_removed_jet['num'], overlap_removed_jet['index'])
    return { 'el':overlap_removed_el
           , 'mu':overlap_removed_mu
           , 'jet':overlap_removed_jet
           }

# ------------------------------------------------------------------------------
def removeElements(particle_lists, to_remove):
    to_remove = list(set(to_remove))
    to_remove.sort(reverse=True)
    for tr in to_remove:
        particle_lists['num'] -= 1
        for pl in particle_lists:
            if isinstance(particle_lists[pl], list):
                del particle_lists[pl][tr]
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
def getSignalObjects( event
                    , baseline_objects
                    , lep_pt_cut
                    , lep_eta_cut
                    , jet_pt_cut
                    , jet_eta_cut
                    , verbose = False
                    ):
    signal_el  = baseline_objects['el']
    signal_mu  = baseline_objects['mu']
    signal_jet = baseline_objects['jet']

    # Get signal electrons
    to_remove_el = []
    for el_it in xrange(signal_el['num']):
        el_pt  = signal_el['pt'][el_it]
        el_eta = signal_el['eta'][el_it]

        if el_pt < lep_pt_cut:
            if verbose:
                print '  electron %d failed pt cut (%f < %f)' % (el_it, el_pt, lep_pt_cut)
            continue
        if abs(el_eta) < lep_eta_cut:
            if verbose:
                print '  electron %d failed pt cut (%f < %f)' % (el_it, el_pt, lep_pt_cut)
            continue
        to_remove_el.append(el_it)
    removeElements(signal_el, to_remove_el)

    # Get signal muons
    to_remove_mu = []
    for mu_it in xrange(signal_mu['num']):
        mu_pt  = signal_mu['pt'][mu_it]
        mu_eta = signal_mu['eta'][mu_it]

        if mu_pt < lep_pt_cut:
            if verbose:
                print '  muon %d failed pt cut (%f < %f)' % (mu_it, mu_pt, lep_pt_cut)
            continue
        if abs(mu_eta) < lep_eta_cut:
            if verbose:
                print '  muon %d failed eta cut (|%f| < %f)' % (mu_it, mu_eta, lep_eta_cut)
            continue
        to_remove_mu.append(mu_it)
    removeElements(signal_mu, to_remove_mu)

    # Get signal muons
    to_remove_jet = []
    for jet_it in xrange(signal_jet['num']):
        jet_pt  = signal_jet['pt'][jet_it]
        jet_eta = signal_jet['eta'][jet_it]

        if jet_pt < lep_pt_cut:
            if verbose:
                print '  jet %d failed pt cut (%f < %f)' % (jet_it, jet_pt, lep_pt_cut)
            continue
        if abs(jet_eta) < lep_eta_cut:
            if verbose:
                print '  jet %d failed eta cut (|%f| < %f)' % (jet_eta, jet_eta, lep_eta_cut)
            continue
        to_remove_jet.append(jet_it)
    removeElements(signal_jet, to_remove_jet)

    return {'el':signal_el, 'mu':signal_mu, 'jet':signal_jet}

# ------------------------------------------------------------------------------
def getFlavorChannel(signal_objects, verbose = False):
    if verbose:
        print 'getting flavor channel'

    num_el  = signal_objects['el']['num']
    num_mu  = signal_objects['mu']['num']
    num_lep = num_el+num_mu

    # 2-lepton events
    if num_lep == 2:
        # get charge product
        charge_product = 1
        for el_it in xrange(num_el):
            charge_product *= signal_objects['el']['charge'][el_it]
        for mu_it in xrange(num_mu):
            charge_product *= signal_objects['mu']['charge'][mu_it]

        if num_el == 2 and charge_product < 0: return 'ee_os'
        if num_el == 1 and charge_product < 0:
            if signal_objects['el']['pt'][0] >= signal_objects['mu']['pt'][0]:
                return 'em_os'
            return 'me_os'
        if num_el == 0 and charge_product < 0: return 'mm_os'
        if num_el == 2 and charge_product > 0: return 'ee_ss'
        if num_el == 1 and charge_product > 0:
            if signal_objects['el']['pt'][0] >= signal_objects['mu']['pt'][0]:
                return 'em_ss'
            return 'me_ss'
        if num_el == 0 and charge_product > 0: return 'mm_ss'

        print 'Oh no! Di-lepton event did not fall into any channel!!!'
        assert False

    # 3-lepton events
    if num_lep == 3:
        if num_el == 3: return 'eee'
        if num_el == 2: return 'eem'
        if num_el == 1: return 'emm'
        if num_el == 0: return 'mmm'

        print 'Oh no! Tri-lepton event did not fall into any channel!!!'
        assert False

    # multi-lepton events
    if num_lep >= 4:
        return 'multi'

    return 'none'

# ------------------------------------------------------------------------------
class hFlavorChannels(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_flavor_channel'
                , title = 'flavor channels'
                ):
        self.hist = ROOT.TH1F( name
                             , title
                             , len(flavor_channels)
                             , -0.5
                             , len(flavor_channels)-0.5
                             )
        # for fc in flavor_channels:
        for i, fc in enumerate(flavor_channels):
            # self.hist.GetXaxis().SetBinLabel(flavor_channels[fc]+1, fc)
            self.hist.GetXaxis().SetBinLabel(i+1, fc)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        # bin_num = flavor_channels[flavor_channel]
        bin_num = flavor_channels.index(flavor_channel)
        self.hist.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        self.hist.Write()
        self.hist.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        self.hist.Write()

# ------------------------------------------------------------------------------
class hPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_pt'
                , title = 'pt'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist_pt   = {}
        self.hist_diff = {}
        self.hist_2d   = {}
        for fc in flavor_channels:
            self.hist_pt[fc] = []
            for num_lep in xrange(max_num_leptons):
                self.hist_pt[fc].append( ROOT.TH1F( '%s_%s_%d' % (fc,name, num_lep)
                                                  , '%s - %s -- lepton %d; p_{T}^{%d} [GeV]' % (title, fc, num_lep, num_lep)
                                                  , num_bins, x_min, x_max
                                                  )
                                       )

            self.hist_diff[fc] = ROOT.TH1F( '%s_%s_diff' % (fc, name)
                                          , '%s - %s -- diff;p_{T}^{0} - p_{T}^{1} [GeV]' % (title, fc)
                                          , num_bins
                                          , -x_max
                                          , x_max
                                          )
            self.hist_2d[fc] = ROOT.TH2F( '%s_%s_2d' % (fc, name)
                                        , '%s - %s -- diff;p_{T}^{0};p_{T}^{1} [GeV]' % (title, fc)
                                        , num_bins
                                        , x_min
                                        , x_max
                                        , num_bins
                                        , x_min
                                        , x_max
                                        )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        lep_pt_list = []
        for el_pt in signal_objects['el']['pt']:
            lep_pt_list.append(el_pt/1000)
        for mu_pt in signal_objects['mu']['pt']:
            lep_pt_list.append(mu_pt/1000)

        for lep_it, lep_pt in enumerate(lep_pt_list):
            if lep_it == max_num_leptons: break
            self.hist_pt[flavor_channel][lep_it].Fill(lep_pt)
        if len(lep_pt_list) >= 2:
            self.hist_diff[flavor_channel].Fill(lep_pt_list[0]-lep_pt_list[1])
            self.hist_2d[flavor_channel].Fill(lep_pt_list[0],lep_pt_list[1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in flavor_channels:
            for num_lep in xrange(max_num_leptons):
                self.hist_pt[fc][num_lep].Write()
            self.hist_diff[fc].Write()
            self.hist_2d[fc].Write()

# ------------------------------------------------------------------------------
class hEta(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_eta'
                , title = 'eta'
                ):
        num_bins = 10
        x_min = -3.
        x_max = +3.

        self.hist_eta = {}
        for fc in flavor_channels:
            self.hist_eta[fc] = []
            for num_lep in xrange(max_num_leptons):
                self.hist_eta[fc].append( ROOT.TH1F( '%s_%s_%d' % (fc, name, num_lep)
                                                   , '%s - %s -- lepton %d; #eta^{%d}' % (title, fc, num_lep, num_lep)
                                                   , num_bins, x_min, x_max
                                                   )
                                        )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        lep_eta_list = []

        for el_eta in signal_objects['el']['eta']:
            lep_eta_list.append(el_eta)
        for mu_eta in signal_objects['mu']['eta']:
            lep_eta_list.append(mu_eta)
            
        for lep_it, lep_eta in enumerate(lep_eta_list):
            if lep_it == max_num_leptons: break
            self.hist_eta[flavor_channel][lep_it].Fill(lep_eta)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in flavor_channels:
            for num_lep in xrange(max_num_leptons):
                self.hist_eta[fc][num_lep].Write()

# ------------------------------------------------------------------------------
class hNumJet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_num_jet'
                , title = 'num jet'
                ):
        num_bins = 10
        x_min = -0.5
        x_max = num_bins - 0.5

        self.hist   = {}
        for fc in flavor_channels:
            self.hist[fc] = ROOT.TH1F( '%s_%s' % (fc, name)
                                     , '%s - %s; Jet Multiplicity' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        self.hist[flavor_channel].Fill(signal_objects['jet']['num'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in flavor_channels:
            self.hist[fc].Write()

# ------------------------------------------------------------------------------
class hJetPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_jet_pt'
                , title = 'jet pt'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist = {}
        for fc in flavor_channels:
            self.hist[fc] = []
            for num_jet in xrange(max_num_jets):
                self.hist[fc].append( ROOT.TH1F( '%s_%s_%d' % (fc, name, num_jet)
                                               , '%s - %s -- jet %d; p_{T}^{%d} [GeV]' % (title, fc, num_jet, num_jet)
                                               , num_bins, x_min, x_max
                                               )
                                    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        jet_pt_list = []
        for jet_pt in signal_objects['jet']['pt']:
            jet_pt_list.append(jet_pt/1000)

        # jet_pt_list.sort(reverse=True)
        for jet_it, jet_pt in enumerate(jet_pt_list):
            if jet_it == max_num_jets: break
            self.hist[flavor_channel][jet_it].Fill(jet_pt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in flavor_channels:
            for jet_itr in xrange(max_num_jets):
                self.hist[fc][jet_itr].Write()

# ------------------------------------------------------------------------------
class hMet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_met'
                , title = 'met'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist_int   = {}
        self.hist_int_w_mu = {}
        self.hist_noint = {}
        self.hist_diff = {}
        self.hist_diff_w_mu = {}
        for fc in flavor_channels:
            self.hist_int[fc] = ROOT.TH1F( '%s_%s_int' % (fc, name)
                                         , '%s int - %s; E_{T}^{miss, int} [GeV]' % (title, fc)
                                         , num_bins, x_min, x_max
                                         )
            self.hist_int_w_mu[fc] = ROOT.TH1F( '%s_%s_int_w_mu' % (fc, name)
                                              , '%s int - %s; E_{T}^{miss, int+#mu} [GeV]' % (title, fc)
                                              , num_bins, x_min, x_max
                                              )
            self.hist_noint[fc] = ROOT.TH1F( '%s_%s_noint' % (fc, name)
                                           , '%s no int - %s; E_{T}^{miss,no int} [GeV]' % (title, fc)
                                           , num_bins, x_min, x_max
                                           )
            self.hist_diff[fc] = ROOT.TH1F( '%s_%s_int_noint_diff' % (fc, name)
                                           , '%s diff - %s; E_{T}^{miss,no int} - E_{T}^{miss,int} [GeV]' % (title, fc)
                                           , num_bins, -x_max, x_max
                                           )
            self.hist_diff_w_mu[fc] = ROOT.TH1F( '%s_%s_int_noint_diff_w_mu' % (fc, name)
                                               , '%s diff - %s; E_{T}^{miss,no int} - E_{T}^{miss,int+#mu} [GeV]' % (title, fc)
                                               , num_bins, -x_max, x_max
                                               )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        # interacting type met
        met_etx_int = event.MET_Truth_Int_etx
        met_ety_int = event.MET_Truth_Int_ety
        met_int = math.sqrt( met_etx_int*met_etx_int
                           + met_ety_int*met_ety_int
                           )/1000.
        self.hist_int[flavor_channel].Fill(met_int)

        # interacting type met + muons pT
        met_etx_int_w_mu = met_etx_int
        met_ety_int_w_mu = met_ety_int
        for mu_index in signal_objects['mu']['index']:
            met_etx_int_w_mu -= event.mu_staco_px.at(mu_index)
            met_ety_int_w_mu -= event.mu_staco_py.at(mu_index)
        met_int_w_mu = math.sqrt( met_etx_int_w_mu*met_etx_int_w_mu
                                + met_ety_int_w_mu*met_ety_int_w_mu
                                )/1000.
        self.hist_int_w_mu[flavor_channel].Fill(met_int_w_mu)

        # non interacting type met
        met_etx_noint = event.MET_Truth_NonInt_etx
        met_ety_noint = event.MET_Truth_NonInt_ety
        met_noint = math.sqrt( met_etx_noint*met_etx_noint
                             + met_ety_noint*met_ety_noint
                             )/1000.
        self.hist_noint[flavor_channel].Fill(met_noint)

        # difference between interacting and non-interacting type mets
        self.hist_diff[flavor_channel].Fill(met_noint - met_int)
        self.hist_diff_w_mu[flavor_channel].Fill(met_noint - met_int_w_mu)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in flavor_channels:
            self.hist_int[fc].Write()
            self.hist_int_w_mu[fc].Write()
            self.hist_noint[fc].Write()
            self.hist_diff[fc].Write()
            self.hist_diff_w_mu[fc].Write()

# ------------------------------------------------------------------------------
def getInFile(in_file):
    f = None
    t = None
    try:
        if os.path.exists(in_file):
            f = ROOT.TFile.Open(in_file)
            t = f.Get('truth')
        if f is None:
            print 'Please provide valid file name'
            raise
        if t is None:
            print 'No tree named "truth" in input file'
            raise
    except:
        print 'in except'
        return None
    return {'file':f, 'tree':t}

# ------------------------------------------------------------------------------
def defineHists():
    hists = {}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['flavor_channel'] = ROOT.TH1I( 'flavor_channel'
                                       , 'flavor channels'
                                       , len(flavor_channels)
                                       , -0.5
                                       , len(flavor_channels)-0.5
                                       )

# ------------------------------------------------------------------------------
def plotTruth(tree):
    hists = {}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['channels'] = hFlavorChannels(name = 'channels')
    hists['pt']       = hPt(name             = 'pt'      )
    hists['eta']      = hEta(name            = 'eta'     )
    hists['num_jet']  = hNumJet(name         = 'num_jet' )
    hists['jet_pt']   = hJetPt(name          = 'jet_pt'  )
    hists['met']      = hMet(name            = 'met'     )
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    total_num_events = tree.GetEntries()
    for i, event in enumerate(tree):
        # print '======================================================'
        if i % 100 == 0:
            print 'Event %d of %d' % (i, total_num_events)
        num_el = event.el_n
        num_mu = event.mu_staco_n

        signal_objects = doObjectSelection( event
                                            , lep_pt_cut  = 10.e3
                                            , lep_eta_cut = 2.4
                                            , jet_pt_cut  = 20.e3
                                            , jet_eta_cut = 2.7
                                            # , verbose = True
                                            )
        flavor_channel = getFlavorChannel(signal_objects)

        for h in hists:
            hists[h].fill(flavor_channel, signal_objects, event)

    return hists

# ------------------------------------------------------------------------------
def readInputs():
    if len(sys.argv) != 3:
        print 'Incorrect number of inputs'
        print 'usage: '
        print '  %s <input file name> <output file name>' % sys.argv[0]
        return None

    in_file = sys.argv[1]
    out_file = sys.argv[2]

    return {'in':in_file, 'out':out_file}

# ------------------------------------------------------------------------------
def main():
    configs = readInputs()
    if configs is None: return

    input = getInFile(configs['in'])
    hists = plotTruth(input['tree'])

    out_file = ROOT.TFile(configs['out'], 'RECREATE')
    for h in hists:
        hists[h].writeToFile(out_file)
    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    main()
