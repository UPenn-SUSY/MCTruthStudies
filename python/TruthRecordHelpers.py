#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import array

from ctypes import cdll
from ctypes import c_double
from ctypes import c_void_p

import ROOT

# this_script_loc = os.path.realpath(__file__)
# lib_loc = '%s/../lib/libTruthRecordHelpers.so' % os.path.dirname(this_script_loc)
# print 'this_script_loc: %s' % this_script_loc
# print 'lib_loc: %s' % lib_loc
# truth_record_helpers_lib = cdll.LoadLibrary(lib_loc)

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
    # return truth_helpers.getParentPdgIdFromBarcode( barcode
    #                                               , event.mc_barcode
    #                                               , event.mc_pdgId
    #                                               , event.mc_parent_index
    #                                               , verbose = False
    #                                               )

    mc_index = None
    for index in reversed(xrange(event.mc_n)):
        if event.mc_barcode.at(index) == barcode:
            mc_index = index
            return getParentPdgID(event, mc_index)
    return None


# # ------------------------------------------------------------------------------
# def getParentPdgIdFromBarcode( barcode
#                              , mc_barcode
#                              , mc_pdg_id
#                              , mc_parent_index
#                              , verbose = False
#                              ):
#     return truth_record_helpers_lib.getParentPdgIdFromBarcode( barcode
#                                                              , mc_barcode
#                                                              , mc_pdg_id
#                                                              , mc_parent_index
#                                                              , verbose
#                                                              )
# 
# def getParentPdgId( mc_index
#                   , mc_mc_pdg_id
#                   , mc_parent_index
#                   , verbose
#                   ):
#     return truth_record_helpers_lib.getParentPdgId( mc_index
#                                                   , mc_mc_pdg_id
#                                                   , mc_parent_index
#                                                   , verbose
#                                                   )
