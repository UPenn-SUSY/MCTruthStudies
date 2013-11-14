#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import ROOT
import rootlogon
import metaroot

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
def getTruthMbl(el_list, mu_list, jet_list):
    mbl_list = []

    return mbl_list

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
def getMbl(el_list, mu_list, b_jet_list):
    mbl_list = []

    for el in el_list:
        dr_min = None
        closest_b_jet = None
        for bjl in b_jet_list:
            dr = getDeltaR(el, bjl)
            if dr_min is None or dr < dr_min:
                dr_min = dr
                closest_b_jet = bjl
        if closest_b_jet is not None:
            mbl_list.append(getMll([el], [closest_b_jet]))
    for mu in mu_list:
        dr_min = None
        closest_b_jet = None
        for bjl in b_jet_list:
            dr = getDeltaR(mu, bjl)
            if dr_min is None or dr < dr_min:
                dr_min = dr
                closest_b_jet = bjl
        if closest_b_jet is not None:
            mbl_list.append(getMll([mu], [closest_b_jet]))

    return mbl_list

# ------------------------------------------------------------------------------
def getDeltaR(obj1, obj2):
    deta = abs(obj1.eta - obj2.eta)
    dphi = abs(obj1.phi - obj2.phi)
    while dphi > 3.14159:
        dphi -= 3.14159
    return math.sqrt(deta*deta + dphi*dphi)

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

