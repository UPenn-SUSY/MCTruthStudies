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
baseline_el_pt_cut  = 10.e3
baseline_mu_pt_cut  = 10.e3
baseline_jet_pt_cut = 20.e3

baseline_el_eta_cut  = 2.4
baseline_mu_eta_cut  = 2.4
baseline_jet_eta_cut = 2.4


signal_el_pt_cut  = 10.e3
signal_mu_pt_cut  = 10.e3
signal_jet_pt_cut = 20.e3

signal_el_eta_cut  = 2.4
signal_mu_eta_cut  = 2.4
signal_jet_eta_cut = 2.4


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
                  # , 'theta':[]
                  , 'charge':[]
                  , 'px':[]
                  , 'py':[]
                  , 'pz':[]
                  , 'E':[]
                  }

    baseline_mu = { 'num':0
                  , 'index':[]
                  , 'pt':[]
                  , 'eta':[]
                  , 'phi':[]
                  # , 'theta':[]
                  , 'charge':[]
                  , 'px':[]
                  , 'py':[]
                  , 'pz':[]
                  , 'E':[]
                  }

    baseline_jet = { 'num':0
                   , 'index':[]
                   , 'pt':[]
                   , 'eta':[]
                   , 'phi':[]
                   , 'theta':[]
                   , 'px':[]
                   , 'py':[]
                   , 'pz':[]
                   , 'E':[]
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
        # baseline_el['theta'].append(    event.el_theta.at(el_index))
        baseline_el['charge'].append( event.el_charge.at(el_index))
        baseline_el['px'].append(    event.el_px.at(el_index))
        baseline_el['py'].append(    event.el_py.at(el_index))
        baseline_el['pz'].append(    event.el_pz.at(el_index))
        baseline_el['E'].append(     event.el_E.at(el_index))

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
        # baseline_mu['theta'].append(    event.mu_staco_theta.at(mu_index))
        baseline_mu['charge'].append( event.mu_staco_charge.at(mu_index))
        baseline_mu['px'].append(     event.mu_staco_px.at(mu_index))
        baseline_mu['py'].append(     event.mu_staco_py.at(mu_index))
        baseline_mu['pz'].append(     event.mu_staco_pz.at(mu_index))
        baseline_mu['E'].append(      event.mu_staco_E.at(mu_index))

    # get baseline jets
    jet_index_order = getPtSortedIndices(event.jet_AntiKt4TruthJets_n, event.jet_AntiKt4TruthJets_pt)
    for jet_index in jet_index_order:
        jet_pt     = event.jet_AntiKt4TruthJets_pt.at(jet_index)

        if jet_pt < jet_pt_cut:
            if verbose:
                print '  jet %d failed pt cut (%f < %f)' % (jet_index, jet_pt, jet_pt_cut)
            continue

        jet_phi = event.jet_AntiKt4TruthJets_phi.at(jet_index)
        jet_eta = event.jet_AntiKt4TruthJets_eta.at(jet_index)
        jet_theta = math.copysign(2*math.atan(math.exp(-abs(jet_eta))), jet_eta)

        baseline_jet['num'] += 1
        baseline_jet['index'].append(jet_index)
        baseline_jet['pt'].append(   jet_pt)
        baseline_jet['eta'].append(  jet_eta)
        baseline_jet['phi'].append(  jet_phi)
        baseline_jet['theta'].append(jet_theta)
        baseline_jet['px'].append(   jet_pt*math.cos(jet_phi))
        baseline_jet['py'].append(   jet_pt*math.sin(jet_phi))
        baseline_jet['pz'].append(   jet_pt*math.sin(jet_theta))
        baseline_jet['E'].append(    event.jet_AntiKt4TruthJets_E.at(jet_index))

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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get signal jets
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

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # get mll
    mll = getMll(signal_el, signal_mu)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # get met and met-related variables
    signal_met = {}

    met_etx_int = event.MET_Truth_Int_etx
    met_ety_int = event.MET_Truth_Int_ety
    met_int = math.sqrt( met_etx_int*met_etx_int
                        + met_ety_int*met_ety_int
                        )/1000.
    metrel_int = getMetRel( met_etx_int
                          , met_ety_int
                          , signal_el
                          , signal_mu
                          , signal_jet
                          , event
                          )/1000.

    met_etx_noint = event.MET_Truth_NonInt_etx
    met_ety_noint = event.MET_Truth_NonInt_ety
    met_noint = math.sqrt( met_etx_noint*met_etx_noint
                        + met_ety_noint*met_ety_noint
                        )/1000.
    metrel_noint = getMetRel( met_etx_noint
                            , met_ety_noint
                            , signal_el
                            , signal_mu
                            , signal_jet
                            , event
                            )/1000.

    signal_met = { 'int':met_int    , 'rel_int':metrel_int
                 , 'noint':met_noint, 'rel_noint':metrel_noint
                 }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return { 'el':signal_el
           , 'mu':signal_mu
           , 'jet':signal_jet
           , 'mll':mll, 'met':signal_met
           }

# ------------------------------------------------------------------------------
def getMetRel( met_etx, met_ety , signal_el, signal_mu, signal_jet, event):
    met_phi = math.atan2(met_ety, met_etx)

    min_dphi = 9999999999

    for el_index in signal_el['index']:
        el_phi = event.el_phi.at(el_index)
        dphi = abs(met_phi - el_phi)
        while dphi > 3.14159: dphi -= 3.14159
        if dphi < min_dphi: min_dphi = dphi

    for mu_index in signal_mu['index']:
        mu_phi = event.mu_staco_phi.at(mu_index)
        dphi = abs(met_phi - mu_phi)
        while dphi > 3.14159: dphi -= 3.14159
        if dphi < min_dphi: min_dphi = dphi

    for jet_index in signal_jet['index']:
        mu_phi = event.jet_AntiKt4TruthJets_phi.at(jet_index)
        dphi = abs(met_phi - mu_phi)
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
def getMll(el_list, mu_list):
    px = 0.
    py = 0.
    pz = 0.
    e = 0.

    for el_it in xrange(len(el_list['index'])):
        px += el_list['px'][el_it]
        py += el_list['py'][el_it]
        pz += el_list['pz'][el_it]
        e  += el_list['E' ][el_it]
    for mu_it in xrange(len(mu_list['index'])):
        px += mu_list['px'][mu_it]
        py += mu_list['py'][mu_it]
        pz += mu_list['pz'][mu_it]
        e  += mu_list['E' ][mu_it]

    m2 = (e*e - px*px - py*py - pz*pz)
    return math.copysign(math.sqrt(abs(m2)), m2)
