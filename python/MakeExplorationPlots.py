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
# flavor_channels = { 'ee':0
#                   , 'mm':1
#                   , 'em':2
#                   , 'eee':3
#                   , 'eem':4
#                   , 'emm':5
#                   , 'mmm':6
#                   , 'multi':7
#                   , 'none':8
#                   }
flavor_channels = { 'ee_os':0
                  , 'ee_ss':1
                  , 'mm_os':2
                  , 'mm_ss':3
                  , 'em_os':4
                  , 'em_ss':5
                  , 'eee':6
                  , 'eem':7
                  , 'emm':8
                  , 'mmm':9
                  , '4l':10
                  , 'none':11
                  }

# ------------------------------------------------------------------------------
class hFlavorChannels(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__(self, name = 'h_flavor_chanel', title = 'flavor channels'):
        self.hist = ROOT.TH1F( name
                             , title
                             , len(flavor_channels)
                             , -0.5
                             , len(flavor_channels)-0.5
                             )
        for fc in flavor_channels:
            self.hist.GetXaxis().SetBinLabel(flavor_channels[fc]+1, fc)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        num_el = event.el_n
        num_mu = event.mu_staco_n

        bin = flavor_channels['none']

        if num_el+num_mu == 2:
            if num_el == 2:
                if event.el_charge.at(0)*event.el_charge.at(1) == -1:
                    bin = flavor_channels['ee_os']
                elif event.el_charge.at(0)*event.el_charge.at(1) == 1:
                    bin = flavor_channels['ee_ss']
                else:
                    print 'WARNING: ee charge is not OS or SS'
                    print '    q0: %d  --  q1: %d' % (event.el_charge.at(0), event.el_charge.at(1))
            elif num_mu == 2:
                if event.mu_staco_charge.at(0)*event.mu_staco_charge.at(1) == -1:
                    bin = flavor_channels['mm_os']
                elif event.mu_staco_charge.at(0)*event.mu_staco_charge.at(1) == 1:
                    bin = flavor_channels['mm_ss']
                else:
                    print 'WARNING: mm charge is not OS or SS'
                    print '    q0: %d  --  q1: %d' % (event.mu_staco_charge.at(0), event.mu_staco_charge.at(1))
            else:
                if event.el_charge.at(0)*event.mu_staco_charge.at(0) == -1:
                    bin = flavor_channels['em_os']
                elif event.el_charge.at(0)*event.mu_staco_charge.at(0) == 1:
                    bin = flavor_channels['em_ss']
                else:
                    print 'WARNING: em charge is not OS or SS'
                    print '    q0: %d  --  q1: %d' % (event.el_charge.at(0), event.mu_staco_charge.at(0))
        elif num_el+num_mu == 3:
            if num_el == 3:
                bin = flavor_channels['eee']
            elif num_el == 2:
                bin = flavor_channels['eem']
            elif num_el == 1:
                bin = flavor_channels['emm']
            else:
                bin = flavor_channels['mmm']
        elif num_el+num_mu == 4:
            bin = flavor_channels['none']

        self.hist.Fill(bin)

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
    hists['channels'] = hFlavorChannels()

    for event in tree:
        num_el = event.el_n
        num_mu = event.mu_staco_n
        if num_el + num_mu <= 1: continue
        print 'Event: %d -- num_lep: %d' % (event.EventNumber, num_el+num_mu)
        for el_it in xrange(num_el):
            px = event.el_px.at(el_it)
            py = event.el_py.at(el_it)
            pt = math.sqrt(px*px + py*py)
            print '    el: %d  -  charge: %d  -  px: %d  -  py: %d  -  pt: %d' % (el_it, event.el_charge.at(el_it), px, py, pt)
        for mu_it in xrange(num_mu):
            px = event.mu_staco_px.at(mu_it)
            py = event.mu_staco_py.at(mu_it)
            pt = math.sqrt(px*px + py*py)
            print '    mu: %d  -  charge: %d  -  px: %d  -  py: %d  -  pt: %d' % (mu_it, event.mu_staco_charge.at(mu_it), px, py, pt)

        hists['channels'].fill(event)

    return hists


# ------------------------------------------------------------------------------
def main():
    inputs = getInFile('/afs/cern.ch/user/b/bjackson/work/public/MCGenTesting/144885/my.truth.ntup.root')
    if inputs is None: return
    print 'got in file'
    hists = plotTruth(inputs['tree'])

    hists['channels'].hist.Draw()
    return hists

# ==============================================================================
if __name__ == '__main__':
    # main()
    hists = main()
