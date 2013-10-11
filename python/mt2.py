#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import array

from ctypes import cdll
from ctypes import c_double
from ctypes import c_void_p

this_script_loc = os.path.realpath(__file__)
lib_loc = '%s/../lib/libmt2.so' % os.path.dirname(this_script_loc)
print 'this_script_loc: %s' % this_script_loc
print 'lib_loc: %s' % lib_loc
mt2_lib = cdll.LoadLibrary(lib_loc)

# ------------------------------------------------------------------------------
# def getMT2(p0, p1, met):
def getMT2(signal_el, signal_mu, met_x, met_y, minv = 0., verbose = False):
    signal_leptons = []

    # add signal electrons
    for el_it in xrange(signal_el['num']):
        signal_leptons.append( { 'm':0.511
                               , 'px':signal_el['px'][el_it]
                               , 'py':signal_el['py'][el_it]
                               }
                             )
    # add signal muons
    for mu_it in xrange(signal_mu['num']):
        signal_leptons.append( { 'm':105.7
                               , 'px':signal_mu['px'][mu_it]
                               , 'py':signal_mu['py'][mu_it]
                               }
                             )

    # is number signal leptons != 2, return default value
    if len(signal_leptons) is not 2: return -999

    # get mT2
    mt2_lib.getMt2.restype = c_double
    # print 'getting mt2 value'
    mt2_value = mt2_lib.getMt2( c_double(signal_leptons[0]['m'])
                              , c_double(signal_leptons[0]['px'])
                              , c_double(signal_leptons[0]['py'])
                              , c_double(signal_leptons[1]['m'])
                              , c_double(signal_leptons[1]['px'])
                              , c_double(signal_leptons[1]['py'])
                              , c_double(met_x), c_double(met_y)
                              , c_double(minv)
                              , verbose
                              )
    if verbose:
        print 'got mt2 value: %s' % mt2_value

    return mt2_value
