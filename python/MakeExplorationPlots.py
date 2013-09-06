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
                  , 'multi':10
                  , 'none':11
                  }

# ------------------------------------------------------------------------------
def getFlavorChannel(event, pt_cut, eta_cut):
    num_el = 0
    good_el_indices = []

    num_mu = 0
    good_mu_indices = []

    charge_product = 1

    for el_it in xrange(event.el_n):
        if event.el_pt.at(el_it) < pt_cut:
            print 'electron %d failed pt cut (%f)' % (el_it, event.el_pt.at(el_it))
            continue
        if abs(event.el_eta.at(el_it)) > eta_cut:
            print 'electron %d failed eta cut (%f)' % (el_it, event.el_eta.at(el_it))
            continue
        num_el += 1
        good_el_indices.append(el_it)
        charge_product *= event.el_charge.at(el_it)

    for mu_it in xrange(event.mu_staco_n):
        if event.mu_staco_pt.at(mu_it) < pt_cut:
            print 'muon %d failed pt cut (%f)' % (mu_it, event.mu_staco_pt.at(mu_it))
            continue
        if abs(event.mu_staco_eta.at(mu_it)) > eta_cut:
            print 'muon %d failed eta cut (%f)' % (mu_it, event.mu_staco_eta.at(mu_it))
            continue
        num_mu += 1
        good_mu_indices.append(mu_it)
        charge_product *= event.mu_staco_charge.at(mu_it)

    flavor_channel = 'none'
    if   num_el == 2 and num_mu == 0 and charge_product < 0:
        flavor_channel = 'ee_os'
    elif num_el == 2 and num_mu == 0 and charge_product > 0:
        flavor_channel = 'ee_ss'
    elif num_el == 0 and num_mu == 2 and charge_product < 0:
        flavor_channel = 'mm_os'
    elif num_el == 0 and num_mu == 2 and charge_product > 0:
        flavor_channel = 'mm_ss'
    elif num_el == 1 and num_mu == 1 and charge_product < 0:
        flavor_channel = 'em_os'
    elif num_el == 1 and num_mu == 1 and charge_product > 0:
        flavor_channel = 'em_ss'
    elif num_el == 3 and num_mu == 0:
        flavor_channel = 'eee'
    elif num_el == 2 and num_mu == 1:
        flavor_channel = 'eem'
    elif num_el == 1 and num_mu == 2:
        flavor_channel = 'emm'
    elif num_el == 0 and num_mu == 3:
        flavor_channel = 'mmm'
    elif num_el + num_mu > 3:
        flavor_channel = 'multi'

    if flavor_channel == 'none' and num_el+num_mu > 0:
        print 'WARNING: no flavor channel, but leptons'
        print '    num el: %d  --  num mu: %d' % (num_el, num_mu)

    return { 'channel':flavor_channel
           , 'el':good_el_indices
           , 'mu':good_mu_indices
           }

# ------------------------------------------------------------------------------
class hFlavorChannels(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_flavor_chanel'
                , title = 'flavor channels'
                , pt_cut = 10000
                , eta_cut = 2.4
                ):
        self.pt_cut  = pt_cut
        self.eta_cut = eta_cut

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

        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         )
        bin_num = flavor_channels[flavor_channel['channel']]
        self.hist.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        self.hist.Write()

# ------------------------------------------------------------------------------
class hPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_pt'
                , title = 'pt'
                , pt_cut = 10000
                , eta_cut = 2.4
                ):
        self.pt_cut  = pt_cut
        self.eta_cut = eta_cut

        num_bins = 20
        x_min = 0
        x_max = 200

        self.hist_0 = {}
        self.hist_1 = {}
        for fc in flavor_channels:
            self.hist_0[fc] = ROOT.TH1F( '%s_%s_0' % (fc, name)
                                       , '%s - %s -- leading;p_{T}^{lead} [GeV]' % (title, fc)
                                       , num_bins
                                       , x_min
                                       , x_max
                                       )
            self.hist_1[fc] = ROOT.TH1F( '%s_%s_1' % (fc, name)
                                       , '%s - %s -- subleading;p_{T}^{sublead} [GeV]' % (title, fc)
                                       , num_bins
                                       , x_min
                                       , x_max
                                       )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        pt_0 = -999
        pt_1 = -999

        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         )
        print flavor_channel
        if flavor_channel['channel'] == 'ee_os' or flavor_channel['channel'] == 'ee_ss':
            print 'found ee channel - getting pt values'
            pt_0 = event.el_pt.at(flavor_channel['el'][0])/1000
            pt_1 = event.el_pt.at(flavor_channel['el'][1])/1000
        elif flavor_channel['channel'] == 'mm_os' or flavor_channel['channel'] == 'mm_ss':
            print 'found mm channel - getting pt values'
            pt_0 = event.mu_staco_pt.at(flavor_channel['mu'][0])/1000
            pt_1 = event.mu_staco_pt.at(flavor_channel['mu'][1])/1000
        elif flavor_channel['channel'] == 'ee_os' or flavor_channel['channel'] == 'ee_ss':
            print 'found em channel - getting pt values'
            pt_0 = event.el_pt.at(      flavor_channel['el'][0])/1000
            pt_1 = event.mu_staco_pt.at(flavor_channel['mu'][0])/1000
        else:
            return

        self.hist_0[flavor_channel['channel']].Fill(pt_0)
        self.hist_1[flavor_channel['channel']].Fill(pt_1)
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in flavor_channels:
            self.hist_0[fc].Write()
            self.hist_1[fc].Write()

# ------------------------------------------------------------------------------
class hEta(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_eta'
                , title = 'eta'
                , pt_cut = 10000
                , eta_cut = 2.4
                ):
        self.pt_cut  = pt_cut
        self.eta_cut = eta_cut

        num_bins = 25
        x_min = -2.5
        x_max = 2.5

        self.hist_0 = {}
        self.hist_1 = {}
        for fc in flavor_channels:
            self.hist_0[fc] = ROOT.TH1F( '%s_%s_0' % (fc, name)
                                       , '%s - %s -- leading;#eta^{lead}' % (title, fc)
                                       , num_bins
                                       , x_min
                                       , x_max
                                       )
            self.hist_1[fc] = ROOT.TH1F( '%s_%s_1' % (fc, name)
                                       , '%s - %s -- subleading;#eta^{sublead}' % (title, fc)
                                       , num_bins
                                       , x_min
                                       , x_max
                                       )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        eta_0 = -999
        eta_1 = -999

        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         )
        print flavor_channel
        if flavor_channel['channel'] == 'ee_os' or flavor_channel['channel'] == 'ee_ss':
            print 'found ee channel - getting pt values'
            eta_0 = event.el_eta.at(flavor_channel['el'][0])/1000
            eta_1 = event.el_eta.at(flavor_channel['el'][1])/1000
        elif flavor_channel['channel'] == 'mm_os' or flavor_channel['channel'] == 'mm_ss':
            print 'found mm channel - getting pt values'
            eta_0 = event.mu_staco_eta.at(flavor_channel['mu'][0])/1000
            eta_1 = event.mu_staco_eta.at(flavor_channel['mu'][1])/1000
        elif flavor_channel['channel'] == 'ee_os' or flavor_channel['channel'] == 'ee_ss':
            print 'found em channel - getting pt values'
            eta_0 = event.el_eta.at(      flavor_channel['el'][0])/1000
            eta_1 = event.mu_staco_eta.at(flavor_channel['mu'][0])/1000
        else:
            return

        self.hist_0[flavor_channel['channel']].Fill(eta_0)
        self.hist_1[flavor_channel['channel']].Fill(eta_1)
        
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in flavor_channels:
            self.hist_0[fc].Write()
            self.hist_1[fc].Write()

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
    hists['channels_8']  = hFlavorChannels( name = 'channels_8'
                                          , pt_cut = 8000
                                          )
    hists['channels_10'] = hFlavorChannels( name = 'channels_10'
                                          , pt_cut = 10000
                                          )
    hists['pt_8']  = hPt(name = 'pt_8' , pt_cut = 8000 )
    hists['pt_10'] = hPt(name = 'pt_10', pt_cut = 10000)
    hists['eta_8']  = hEta(name = 'eta_8' , pt_cut = 8000 )
    hists['eta_10'] = hEta(name = 'eta_10', pt_cut = 10000)

    for event in tree:
        num_el = event.el_n
        num_mu = event.mu_staco_n
        if num_el + num_mu <= 1: continue

        for h in hists:
            hists[h].fill(event)
        # hists['channels_8'].fill(event)
        # hists['channels_10'].fill(event)
        # hists['pt_8'].fill(event)
        # hists['pt_10'].fill(event)

    return hists


# ------------------------------------------------------------------------------
def main():
    inputs = getInFile('/Users/bjackson/work/MCGenTesting/my.truth.ntup.root')
    # inputs = getInFile('/afs/cern.ch/user/b/bjackson/work/public/MCGenTesting/144885/my.truth.ntup.root')
    if inputs is None: return
    print 'got in file'
    hists = plotTruth(inputs['tree'])

    out_file = ROOT.TFile('out_hists.root', 'RECREATE')
    for h in hists:
        hists[h].writeToFile(out_file)
    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    main()
