#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import ROOT
import rootlogon
import metaroot

import DiLeptonCutflow as cutflow
import TruthHists

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
                                       , len(cutflow.flavor_channels)
                                       , -0.5
                                       , len(cutflow.flavor_channels)-0.5
                                       )

# ------------------------------------------------------------------------------
def plotTruth(tree):
    hists = {}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['channels'] = TruthHists.hFlavorChannels()
    hists['pt']       = TruthHists.hPt()
    hists['eta']      = TruthHists.hEta()
    hists['num_jet']  = TruthHists.hNumJet()
    hists['jet_pt']   = TruthHists.hJetPt()
    hists['met']      = TruthHists.hMet()
    hists['mll']      = TruthHists.hMll()
    hists['mt2']      = TruthHists.hMt2()
    hists['ptll']     = TruthHists.hPtll()
    hists['emma_mt']  = TruthHists.hEmmaMt()
    hists['sr_ss']    = TruthHists.hSRSS()
    hists['sr_os']    = TruthHists.hSROS()
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    total_num_events = tree.GetEntries()
    for i, event in enumerate(tree):
        # print '======================================================'
        if i % 100 == 0:
            print 'Event %d of %d' % (i, total_num_events)
        num_el = event.el_n
        num_mu = event.mu_staco_n

        signal_objects = cutflow.doObjectSelection( event
                                          , lep_pt_cut  = 10.e3
                                          , lep_eta_cut = 2.4
                                          , jet_pt_cut  = 20.e3
                                          , jet_eta_cut = 2.7
                                          # , verbose = True
                                          )
        flavor_channel = cutflow.getFlavorChannel(signal_objects)

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
    print 'getting configs'
    configs = readInputs()
    if configs is None: return

    print 'getting input file'
    input = getInFile(configs['in'])
    print 'about to make plots!!!'
    hists = plotTruth(input['tree'])

    print 'creating output file'
    out_file = ROOT.TFile(configs['out'], 'RECREATE')
    print 'writing hists to file'
    for h in hists:
        hists[h].writeToFile(out_file)
    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    main()
