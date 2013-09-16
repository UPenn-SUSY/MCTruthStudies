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
max_num_leptons = 4
max_num_jets = 2

el_prefix  = 'el_'
mu_prefix  = 'mu_staco_'
jet_prefix = 'jet_AntiKt4TruthJets_'

# ------------------------------------------------------------------------------
def getFlavorChannel(event, pt_cut, eta_cut, verbose = False):
    if verbose:
        print '-----------------------------------'
    num_el = 0
    good_el_indices = []
    good_el_pt = []

    num_mu = 0
    good_mu_indices = []
    good_mu_pt = []

    charge_product = 1

    for el_it in xrange(event.el_n):
        el_pt     = event.el_pt.at(el_it)
        el_eta    = event.el_eta.at(el_it)
        el_charge = event.el_charge.at(el_it)
        if verbose:
            print 'el[%d] -- pt: %f -- eta: %f -- charge: %d' % ( el_it
                                                                , el_pt
                                                                , el_eta
                                                                , el_charge
                                                                )

        if el_pt < pt_cut:
            if verbose:
                print '  electron %d failed pt cut (%f < %f)' % (el_it, el_pt, pt_cut)
            continue
        if abs(event.el_eta.at(el_it)) > eta_cut:
            if verbose:
                print '  electron %d failed eta cut (|%f| > %f)' % (el_it, el_eta, eta_cut)
            continue
        num_el += 1
        good_el_indices.append(el_it)
        good_el_pt.append(el_pt)
        charge_product *= el_charge

    for mu_it in xrange(event.mu_staco_n):
        mu_pt     = event.mu_staco_pt.at(mu_it)
        mu_eta    = event.mu_staco_eta.at(mu_it)
        mu_charge = event.mu_staco_charge.at(mu_it)
        if verbose:
            print 'mu[%d] -- pt: %f -- eta: %f -- charge: %d' % ( mu_it
                                                                , mu_pt
                                                                , mu_eta
                                                                , mu_charge
                                                                )

        if event.mu_staco_pt.at(mu_it) < pt_cut:
            if verbose:
                print '  muon %d failed pt cut (%f < %f)' % (mu_it, mu_pt, pt_cut)
            continue
        if abs(event.mu_staco_eta.at(mu_it)) > eta_cut:
            if verbose:
                print '  muon %d failed eta cut (|%f| > %f)' % (mu_it, mu_eta, eta_cut)
            continue
        num_mu += 1
        good_mu_indices.append(mu_it)
        good_mu_pt.append(mu_pt)
        charge_product *= mu_charge
        if verbose:
            print '  adding to num_mu: %d  --  %s' % (num_mu, good_mu_indices)

    if len(good_el_indices) > 0:
        tups = zip(good_el_pt, good_el_indices)
        tups.sort()
        tups.reverse()
        good_el_pt, good_el_indices = (list(t) for t in zip(*tups))
    if len(good_mu_indices) > 0:
        tups = zip(good_mu_pt, good_mu_indices)
        tups.sort()
        tups.reverse()
        good_mu_pt, good_mu_indices = (list(t) for t in zip(*tups))

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

    if verbose and flavor_channel == 'none' and num_el+num_mu > 0:
        print '  WARNING: no flavor channel, but leptons'
        print '      num el: %d  --  num mu: %d' % (num_el, num_mu)

    return { 'channel':flavor_channel
           , 'el':good_el_indices
           , 'mu':good_mu_indices
           }

# ------------------------------------------------------------------------------
def isOverlap(event, jet_eta, jet_phi, good_el_indices):
    for gei in good_el_indices:
        el_eta = event.el_eta.at(gei)
        el_phi = event.el_phi.at(gei)

        deta = jet_eta - el_eta
        dphi = abs(jet_phi - el_phi)
        while dphi > 3.14159:
            # print ' dphi: %f -- subtracting 3.14159' % dphi
            dphi -= 3.14159
        # print ' dphi: %f' % dphi
        dr = math.sqrt( deta*deta + dphi*dphi )

        return (dr < 0.2)

# ------------------------------------------------------------------------------
def getGoodJets( event
               , jet_pt_cut
               , jet_eta_cut
               , lep_pt_cut
               , lep_eta_cut
               , verbose = False
               ):
    good_jet_indices = []
    good_jet_pt = []

    flavor_channels = getFlavorChannel(event, lep_pt_cut, lep_eta_cut)

    for jet_it in xrange(event.jet_AntiKt4TruthJets_n):
        jet_pt  = event.jet_AntiKt4TruthJets_pt.at(jet_it)
        jet_eta = event.jet_AntiKt4TruthJets_eta.at(jet_it)
        jet_phi = event.jet_AntiKt4TruthJets_phi.at(jet_it)

        if jet_pt < jet_pt_cut:
            if verbose:
                print '  Jet %d failed pt cut (%f < %f)' % (jet_it, jet_pt, jet_pt_cut)
            continue
        if abs(jet_eta) > jet_eta_cut:
            if verbose:
                print '  Jet %d failed eta cut (|%f| > %f)' % (jet_it, jet_eta, jet_eta_cut)
            continue
        if isOverlap(event, jet_eta, jet_phi, flavor_channels['el']):
            if verbose:
                print '  Jet %d overlaps with electron' % jet_it
            continue

        good_jet_indices.append(jet_it)
        good_jet_pt.append(jet_pt)
        if verbose:
            print '  Adding to good jets: %d -- %s' % (len(good_jet_indices), good_jet_indices)

    if len(good_jet_indices) > 0:
        tups = zip(good_jet_pt, good_jet_indices)
        tups.sort()
        tups.reverse()
        good_jet_pt, good_jet_indices = (list(t) for t in zip(*tups))

    return good_jet_indices

# ------------------------------------------------------------------------------
def getInvisibles( event ):
    invisible_indices = {}

    for mc_it in xrange(event.mc_n):
        pass

# ------------------------------------------------------------------------------
class hFlavorChannels(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_flavor_channel'
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
                                         # , True
                                         )
        # print flavor_channel
        bin_num = flavor_channels[flavor_channel['channel']]
        self.hist.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        dir_name = 'pt_%s__eta_%s' % (self.pt_cut, self.eta_cut)
        if out_file.GetDirectory(dir_name) == None:
            out_file.mkdir(dir_name)
        out_file.cd(dir_name)
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

        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist_pt   = {}
        self.hist_diff = {}
        self.hist_2d   = {}
        for fc in flavor_channels:
            self.hist_pt[fc] = []
            for num_lep in xrange(max_num_leptons):
                self.hist_pt[fc].append( ROOT.TH1F( '%s_%s_%d' % (fc,name, num_lep)
                                                  , '%s - %s -- lepton %d; p_{T}^{%d} [GeV]' % (title, fc, num_lep, num_lep)
                                                  , num_bins, x_min, x_max
                                                  )
                                       )

            self.hist_diff[fc] = ROOT.TH1F( '%s_%s_diff' % (fc, name)
                                          , '%s - %s -- diff;p_{T}^{0} - p_{T}^{1} [GeV]' % (title, fc)
                                          , num_bins
                                          , -x_max
                                          , x_max
                                          )
            self.hist_2d[fc] = ROOT.TH2F( '%s_%s_2d' % (fc, name)
                                        , '%s - %s -- diff;p_{T}^{0};p_{T}^{1} [GeV]' % (title, fc)
                                        , num_bins
                                        , x_min
                                        , x_max
                                        , num_bins
                                        , x_min
                                        , x_max
                                        )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        pt = []

        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         # , True
                                         )
        for el_index in flavor_channel['el']:
            pt.append(event.el_pt.at(el_index)/1000.)
        for mu_index in flavor_channel['mu']:
            pt.append(event.mu_staco_pt.at(mu_index)/1000.)

        for num_lep, lep_pt in enumerate(pt):
            if num_lep == max_num_leptons: break
            self.hist_pt[flavor_channel['channel']][num_lep].Fill(lep_pt)
        if len(pt) >= 2:
            self.hist_diff[flavor_channel['channel']].Fill(pt[0]-pt[1])
            self.hist_2d[flavor_channel['channel']].Fill(pt[0],pt[1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        dir_name = 'pt_%s__eta_%s' % (self.pt_cut, self.eta_cut)
        if out_file.GetDirectory(dir_name) == None:
            out_file.mkdir(dir_name)
        out_file.cd(dir_name)
        for fc in flavor_channels:
            for num_lep in xrange(max_num_leptons):
                self.hist_pt[fc][num_lep].Write()
            self.hist_diff[fc].Write()
            self.hist_2d[fc].Write()

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

        num_bins = 10
        x_min = -3.
        x_max = +3.

        self.hist_eta = {}
        for fc in flavor_channels:
            self.hist_eta[fc] = []
            for num_lep in xrange(max_num_leptons):
                self.hist_eta[fc].append( ROOT.TH1F( '%s_%s_%d' % (fc, name, num_lep)
                                                   , '%s - %s -- lepton %d; #eta^{%d}' % (title, fc, num_lep, num_lep)
                                                   , num_bins, x_min, x_max
                                                   )
                                        )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        eta = []

        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         # , True
                                         )
        for el_index in flavor_channel['el']:
            eta.append(event.el_eta.at(el_index))
        for mu_index in flavor_channel['mu']:
            eta.append(event.mu_staco_eta.at(mu_index))

        for num_lep, lep_eta in enumerate(eta):
            if num_lep == max_num_leptons: break
            self.hist_eta[flavor_channel['channel']][num_lep].Fill(lep_eta)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        print 'hEta.writeToFile()'
        out_file.cd()
        dir_name = 'pt_%s__eta_%s' % (self.pt_cut, self.eta_cut)
        if out_file.GetDirectory(dir_name) == None:
            out_file.mkdir(dir_name)
        out_file.cd(dir_name)
        for fc in flavor_channels:
            for num_lep in xrange(max_num_leptons):
                self.hist_eta[fc][num_lep].Write()

# ------------------------------------------------------------------------------
class hNumJet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_num_jet'
                , title = 'num jet'
                , pt_cut = 10000
                , eta_cut = 2.4
                ):
        self.pt_cut  = pt_cut
        self.eta_cut = eta_cut

        num_bins = 10
        x_min = -0.5
        x_max = num_bins - 0.5

        self.hist   = {}
        for fc in flavor_channels:
            self.hist[fc] = ROOT.TH1F( '%s_%s' % (fc, name)
                                     , '%s - %s; Jet Multiplicity' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         # , True
                                         )

        num_jets = len(getGoodJets( event
                                  , 20000
                                  , 2.4
                                  , self.pt_cut
                                  , self.eta_cut
                                  # , True
                                  )
                      )
        self.hist[flavor_channel['channel']].Fill(num_jets)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        dir_name = 'pt_%s__eta_%s' % (self.pt_cut, self.eta_cut)
        if out_file.GetDirectory(dir_name) == None:
            out_file.mkdir(dir_name)
        out_file.cd(dir_name)
        for fc in flavor_channels:
            self.hist[fc].Write()

# ------------------------------------------------------------------------------
class hJetPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_jet_pt'
                , title = 'jet pt'
                , pt_cut = 10000
                , eta_cut = 2.4
                ):
        self.pt_cut  = pt_cut
        self.eta_cut = eta_cut

        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist = {}
        for fc in flavor_channels:
            self.hist[fc] = []
            for num_jet in xrange(max_num_jets):
                self.hist[fc].append( ROOT.TH1F( '%s_%s_%d' % (fc, name, num_jet)
                                               , '%s - %s -- jet %d; p_{T}^{%d} [GeV]' % (title, fc, num_jet, num_jet)
                                               , num_bins, x_min, x_max
                                               )
                                    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        pt = []

        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         # , True
                                         )
        jet_indicies = getGoodJets( event
                                  , 20000
                                  , 2.4
                                  , self.pt_cut
                                  , self.eta_cut
                                  )

        for jet_num, jet_index in enumerate(jet_indicies):
            if jet_num == max_num_jets: break
            self.hist[flavor_channel['channel']][jet_num].Fill(event.jet_AntiKt4TruthJets_pt.at(jet_index)/1000.)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        dir_name = 'pt_%s__eta_%s' % (self.pt_cut, self.eta_cut)
        if out_file.GetDirectory(dir_name) == None:
            out_file.mkdir(dir_name)
        out_file.cd(dir_name)
        for fc in flavor_channels:
            for jet_itr in xrange(max_num_jets):
                self.hist[fc][jet_itr].Write()

# ------------------------------------------------------------------------------
class hMet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , name = 'h_met'
                , title = 'met'
                , pt_cut = 10000
                , eta_cut = 2.4
                ):
        self.pt_cut  = pt_cut
        self.eta_cut = eta_cut

        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist   = {}
        for fc in flavor_channels:
            self.hist[fc] = ROOT.TH1F( '%s_%s' % (fc, name)
                                     , '%s - %s; E_{T}^{miss} [GeV]' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, event):
        flavor_channel = getFlavorChannel( event
                                         , self.pt_cut
                                         , self.eta_cut
                                         # , True
                                         )

        met_etx = event.MET_Truth_NonInt_etx
        met_ety = event.MET_Truth_NonInt_ety
        met = math.sqrt(met_etx*met_etx + met_ety*met_ety)/1000.
        self.hist[flavor_channel['channel']].Fill(met)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        dir_name = 'pt_%s__eta_%s' % (self.pt_cut, self.eta_cut)
        if out_file.GetDirectory(dir_name) == None:
            out_file.mkdir(dir_name)
        out_file.cd(dir_name)
        for fc in flavor_channels:
            self.hist[fc].Write()

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
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['channels']    = hFlavorChannels(name = 'channels'   , pt_cut = 0    , eta_cut = 3)
    hists['channels_8']  = hFlavorChannels(name = 'channels_8' , pt_cut = 8000 )
    hists['channels_10'] = hFlavorChannels(name = 'channels_10', pt_cut = 10000)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['pt']    = hPt(name = 'pt'   , pt_cut = 0    , eta_cut = 3)
    hists['pt_8']  = hPt(name = 'pt_8' , pt_cut = 8000 )
    hists['pt_10'] = hPt(name = 'pt_10', pt_cut = 10000)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['eta']    = hEta(name = 'eta'   , pt_cut = 0    , eta_cut = 3)
    hists['eta_8']  = hEta(name = 'eta_8' , pt_cut = 8000 )
    hists['eta_10'] = hEta(name = 'eta_10', pt_cut = 10000)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['num_jet']    = hNumJet(name = 'num_jet'   , pt_cut = 0    , eta_cut = 3)
    hists['num_jet_8']  = hNumJet(name = 'num_jet_8' , pt_cut = 8000 )
    hists['num_jet_10'] = hNumJet(name = 'num_jet_10', pt_cut = 10000)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['jet_pt']    = hJetPt(name = 'jet_pt'   , pt_cut = 0    , eta_cut = 3)
    hists['jet_pt_8']  = hJetPt(name = 'jet_pt_8' , pt_cut = 8000 )
    hists['jet_pt_10'] = hJetPt(name = 'jet_pt_10', pt_cut = 10000)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hists['met']    = hMet(name = 'met'   , pt_cut = 0    , eta_cut = 3)
    hists['met_8']  = hMet(name = 'met_8' , pt_cut = 8000 )
    hists['met_10'] = hMet(name = 'met_10', pt_cut = 10000)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    total_num_events = tree.GetEntries()
    for i, event in enumerate(tree):
        # print '======================================================'
        if i % 100 == 0:
            print 'Event %d of %d' % (i, total_num_events)
        num_el = event.el_n
        num_mu = event.mu_staco_n

        for h in hists:
            hists[h].fill(event)

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
    configs = readInputs()
    if configs is None: return

    input = getInFile(configs['in'])
    hists = plotTruth(input['tree'])

    out_file = ROOT.TFile(configs['out'], 'RECREATE')
    for h in hists:
        hists[h].writeToFile(out_file)
    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    main()
