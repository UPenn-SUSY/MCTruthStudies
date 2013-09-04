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
def plotTruth(tree):
    for event in tree:
        num_el = event.el_n
        num_mu = event.mu_staco_n
        if num_el + num_mu <= 1: continue
        print 'Event: %d -- num_lep: %d' % (event.EventNumber, num_el+num_mu)
        for el_it in xrange(num_el):
            px = event.el_px.at(el_it)
            py = event.el_py.at(el_it)
            pt = math.sqrt(px*px + py*py)
            print '    el: %d  -  px: %d  -  py: %d  -  pt: %d' % (el_it, px, py, pt)
        for mu_it in xrange(num_mu):
            px = event.mu_staco_px.at(mu_it)
            py = event.mu_staco_py.at(mu_it)
            pt = math.sqrt(px*px + py*py)
            print '    mu: %d  -  px: %d  -  py: %d  -  pt: %d' % (mu_it, px, py, pt)


# ------------------------------------------------------------------------------
def main():
    inputs = getInFile('/afs/cern.ch/user/b/bjackson/work/public/MCGenTesting/144885/my.truth.ntup.root')
    if inputs is None: return
    print 'got in file'
    plotTruth(inputs['tree'])

# ==============================================================================
if __name__ == '__main__':
    main()
