#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import ROOT
import rootlogon
import metaroot

import ObjectDefs as object_defs

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
        this_el = object_defs.Electron(event, el_index)
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
        this_mu = object_defs.Muon(event, mu_index)
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
        this_jet = object_defs.Jet(event, jet_index)
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

