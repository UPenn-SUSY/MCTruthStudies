#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import array

from ctypes import cdll
from ctypes import c_double
from ctypes import c_void_p

# mt2_lib = cdll.LoadLibrary('../mt2/obj/mt2_bisect.o')
# mt2_lib = cdll.LoadLibrary('MCTruthStudies/lib/libmt2.so')
mt2_lib = cdll.LoadLibrary('lib/libmt2.so')

class mt2_bisect(object):
    def __init__(self):
        # tmp = c_void_p(mt2_lib.get_new_mt2_bisect())
        tmp = mt2_lib.get_mt2_bisect()
        print tmp
        self.obj = tmp



        # self.obj = mt2_lib.get_new_mt2_bisect()
        print 'self.obj: %s' % self.obj
        # self.print_mt2()

    def mt2_bisect(self):
        mt2_lib.fun_mt2_bisect(self.obj)

    def mt2_massless(self):
        mt2_lib.fun_mt2_massless(self.obj)

    def set_momenta(self, pa0, pb0, pmiss0):
        mt2_lib.fun_set_momenta( self.obj
                               , double *pa0
                               , double *pb0
                               , double* pmiss0
                               )

    def set_mn(self, mn):
        mt2_lib.fun_set_mn(self.obj, mn)

    def get_mt2(self):
        return 'dummy!!!'
        return mt2_lib.fun_get_mt2(self.obj)

    def print_mt2(self):
        mt2_lib.fun_print(self.obj)

    def nevt():
        return mt2_lib.fun_nevt(self.obj)

    def set_momenta(self, pa0, pb0, pmiss0):
        mt2_lib.fun_set_momenta( self.obj,
                c_double(pa0['m']), c_double(pa0['px']), c_double(pa0['py']),
                c_double(pb0['m']), c_double(pb0['px']), c_double(pb0['py']),
                c_double(pmiss0['px']), c_double(pmiss0['py']))

# ------------------------------------------------------------------------------
# def getMT2(p0, p1, met):
def getMT2(signal_el, signal_mu, met_x, met_y):
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
    mt2_value = mt2_lib.getMt2( c_double(signal_leptons[0]['m'])  , c_double(signal_leptons[0]['px']) , c_double(signal_leptons[0]['py'])
                              , c_double(signal_leptons[1]['m'])  , c_double(signal_leptons[1]['px']) , c_double(signal_leptons[1]['py'])
                              , c_double(met_x), c_double(met_y)
                              )

    return mt2_value
