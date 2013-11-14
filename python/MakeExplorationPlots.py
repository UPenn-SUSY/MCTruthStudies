#!/usr/bin/env python
# ==============================================================================
# = usage ./MakeExplorationPlots.py <INPUT_TRUTH_NTUPLE> <OUT_FILE>
# ==============================================================================


import sys
import os.path
import optparse
import time
import math

import ROOT
import rootlogon
import metaroot

import BMinusLCutflow as cutflow
import TruthHists

# ------------------------------------------------------------------------------
def getInFile(in_file):
    f = None
    t = None
    try:
        print 'Looking for input file: %s' % in_file
        if os.path.exists(in_file) or 'root://' in in_file:
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
def addHistograms( hist_dict
                 , decay_category = 'dc_all'
                 , flavor_channel = 'fc_all'
                 ):
    selection_tag = '%s__%s' % (decay_category, flavor_channel)
    hist_dict[selection_tag] = {}
    # hist_dict[selection_tag]['channels']       = TruthHists.hFlavorChannels(selection_tag = selection_tag)
    # hist_dict[selection_tag]['decay_category'] = TruthHists.hDecayCategory( selection_tag = selection_tag)
    hist_dict[selection_tag]['num_lepton']     = TruthHists.hNumLepton(     selection_tag = selection_tag)
    hist_dict[selection_tag]['lepton_pt']      = TruthHists.hLeptonPt(      selection_tag = selection_tag)
    hist_dict[selection_tag]['lepton_eta']     = TruthHists.hLeptonEta(     selection_tag = selection_tag)
    hist_dict[selection_tag]['num_jet']        = TruthHists.hNumJet(        selection_tag = selection_tag)
    hist_dict[selection_tag]['jet_pt']         = TruthHists.hJetPt(         selection_tag = selection_tag)
    hist_dict[selection_tag]['met']            = TruthHists.hMet(           selection_tag = selection_tag)
    hist_dict[selection_tag]['mll']            = TruthHists.hMll(           selection_tag = selection_tag)
    hist_dict[selection_tag]['mbl']            = TruthHists.hMbl(           selection_tag = selection_tag)
    # hist_dict[selection_tag]['mt2']            = TruthHists.hMt2(           selection_tag = selection_tag)
    # hist_dict[selection_tag]['ptll']           = TruthHists.hPtll(          selection_tag = selection_tag)
    # hist_dict[selection_tag]['emma_mt']        = TruthHists.hEmmaMt(        selection_tag = selection_tag)
    # hist_dict[selection_tag]['sr_ss']          = TruthHists.hSRSS(          selection_tag = selection_tag)
    # hist_dict[selection_tag]['sr_os']          = TruthHists.hSROS(          selection_tag = selection_tag)
    # hist_dict[selection_tag]['pt_by_mother']   = TruthHists.hPtByMother(    selection_tag = selection_tag)
    # hist_dict[selection_tag]['eta_by_mother']  = TruthHists.hEtaByMother(   selection_tag = selection_tag)

# ------------------------------------------------------------------------------
def plotTruth(tree):
    hists = {}
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for dc in cutflow.decay_categories:
        for fc in cutflow.flavor_channels:
            addHistograms(hists, dc, fc)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    total_num_events = tree.GetEntries()
    for i, event in enumerate(tree):
        # print '======================================================'
        if i % 100 == 0:
            print 'Event %d of %d' % (i, total_num_events)
        if i > 500: break
        num_el = event.el_n
        num_mu = event.mu_staco_n

        ewk_cutflow = cutflow.BMinusLCutFlow(event)
        if not ewk_cutflow.valid_cutflow:
            continue

        for dc in cutflow.decay_categories:
            if dc != ewk_cutflow.decay_category and dc != 'dc_all': continue
            for fc in cutflow.flavor_channels:
                if fc != ewk_cutflow.flavor_channel and fc != 'fc_all': continue
                h_selection = '%s__%s' % (dc, fc)
                if not h_selection in hists: continue
                for h in hists[h_selection]:
                    hists[h_selection][h].fill(ewk_cutflow)

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
    for h_selection in hists:
        for h in hists[h_selection]:
            hists[h_selection][h].writeToFile(out_file)
    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    main()
