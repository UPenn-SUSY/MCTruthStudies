#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import array

from ctypes import cdll
from ctypes import c_double
from ctypes import c_void_p

import ROOT

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
    for index in reversed(xrange(event.mc_n)):
        if event.mc_barcode.at(index) == barcode:
            mc_index = index
            return getParentPdgID(event, mc_index)
    return None

# ------------------------------------------------------------------------------
def isBJet( my_jet
          , mc_pdg_id
          , mc_pt
          , mc_eta
          , mc_phi
          ):
    pt_cut_value = 10.e3
    delta_r_cut = 0.30

    my_jet_eta = my_jet.eta
    my_jet_phi = my_jet.phi

    for mc_it in xrange(len(mc_pdg_id)):
        if mc_pt.at(mc_it) < pt_cut_value or abs(mc_pdg_id.at(mc_it)) != 5:
            continue
        delta_eta = abs(abs(my_jet_eta) - abs(mc_eta.at(mc_it)))
        delta_phi = abs(abs(my_jet_phi) - abs(mc_phi.at(mc_it)))
        while delta_phi > 3.14159:
            delta_phi -= 3.14159
        delta_r = math.sqrt( delta_eta*delta_eta + delta_phi*delta_phi)

        if delta_r < delta_r_cut:
            return True
    return False
