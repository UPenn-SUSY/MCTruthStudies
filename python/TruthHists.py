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

# ------------------------------------------------------------------------------
class hFlavorChannels(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'flavor channels'
                ):
        self.hist = ROOT.TH1F( 'h__flavor_channel'
                             , title
                             , len(cutflow.flavor_channels)
                             , -0.5
                             , len(cutflow.flavor_channels)-0.5
                             )
        for i, fc in enumerate(cutflow.flavor_channels):
            self.hist.GetXaxis().SetBinLabel(i+1, fc)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        bin_num = cutflow.flavor_channels.index(flavor_channel)
        self.hist.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        self.hist.Write()
        self.hist.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        self.hist.Write()

# ------------------------------------------------------------------------------
class hPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'pt'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist_pt   = {}
        self.hist_diff = {}
        self.hist_2d   = {}
        for fc in cutflow.flavor_channels:
            self.hist_pt[fc] = []
            for num_lep in xrange(cutflow.max_num_leptons):
                self.hist_pt[fc].append( ROOT.TH1F( 'h__%s__lep_pt_%d' % (fc, num_lep)
                                                  , '%s - %s -- lepton %d; p_{T}^{%d} [GeV]' % (title, fc, num_lep, num_lep)
                                                  , num_bins, x_min, x_max
                                                  )
                                       )

            self.hist_diff[fc] = ROOT.TH1F( 'h__%s__lep_pt_diff' % fc
                                          , '%s - %s -- diff;p_{T}^{0} - p_{T}^{1} [GeV]' % (title, fc)
                                          , num_bins
                                          , -x_max
                                          , x_max
                                          )
            self.hist_2d[fc] = ROOT.TH2F( 'h__%s__lep_pt_2d' % fc
                                        , '%s - %s -- diff;p_{T}^{0};p_{T}^{1} [GeV]' % (title, fc)
                                        , num_bins
                                        , x_min
                                        , x_max
                                        , num_bins
                                        , x_min
                                        , x_max
                                        )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        lep_pt_list = []
        for el_pt in signal_objects['el']['pt']:
            lep_pt_list.append(el_pt/1000)
        for mu_pt in signal_objects['mu']['pt']:
            lep_pt_list.append(mu_pt/1000)

        for lep_it, lep_pt in enumerate(lep_pt_list):
            if lep_it == cutflow.max_num_leptons: break
            self.hist_pt[flavor_channel][lep_it].Fill(lep_pt)
        if len(lep_pt_list) >= 2:
            self.hist_diff[flavor_channel].Fill(lep_pt_list[0]-lep_pt_list[1])
            self.hist_2d[flavor_channel].Fill(lep_pt_list[0],lep_pt_list[1])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            for num_lep in xrange(cutflow.max_num_leptons):
                self.hist_pt[fc][num_lep].Write()
            self.hist_diff[fc].Write()
            self.hist_2d[fc].Write()

# ------------------------------------------------------------------------------
class hEta(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'eta'
                ):
        num_bins = 10
        x_min = -3.
        x_max = +3.

        self.hist_eta = {}
        for fc in cutflow.flavor_channels:
            self.hist_eta[fc] = []
            for num_lep in xrange(cutflow.max_num_leptons):
                self.hist_eta[fc].append( ROOT.TH1F( 'h__%s__lep_eta_%d' % (fc, num_lep)
                                                   , '%s - %s -- lepton %d; #eta^{%d}' % (title, fc, num_lep, num_lep)
                                                   , num_bins, x_min, x_max
                                                   )
                                        )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        lep_eta_list = []

        for el_eta in signal_objects['el']['eta']:
            lep_eta_list.append(el_eta)
        for mu_eta in signal_objects['mu']['eta']:
            lep_eta_list.append(mu_eta)

        for lep_it, lep_eta in enumerate(lep_eta_list):
            if lep_it == cutflow.max_num_leptons: break
            self.hist_eta[flavor_channel][lep_it].Fill(lep_eta)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            for num_lep in xrange(cutflow.max_num_leptons):
                self.hist_eta[fc][num_lep].Write()

# ------------------------------------------------------------------------------
class hNumJet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'num jet'
                ):
        num_bins = 10
        x_min = -0.5
        x_max = num_bins - 0.5

        self.hist   = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = ROOT.TH1F( 'h__%s__num_jets' % fc
                                     , '%s - %s; Jet Multiplicity' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        self.hist[flavor_channel].Fill(signal_objects['jet']['num'])

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            self.hist[fc].Write()

# ------------------------------------------------------------------------------
class hJetPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'jet pt'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = []
            for num_jet in xrange(cutflow.max_num_jets):
                self.hist[fc].append( ROOT.TH1F( 'h__%s__jet_pt_%d' % (fc, num_jet)
                                               , '%s - %s -- jet %d; p_{T}^{%d} [GeV]' % (title, fc, num_jet, num_jet)
                                               , num_bins, x_min, x_max
                                               )
                                    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        jet_pt_list = []
        for jet_pt in signal_objects['jet']['pt']:
            jet_pt_list.append(jet_pt/1000)

        # jet_pt_list.sort(reverse=True)
        for jet_it, jet_pt in enumerate(jet_pt_list):
            if jet_it == cutflow.max_num_jets: break
            self.hist[flavor_channel][jet_it].Fill(jet_pt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            for jet_itr in xrange(cutflow.max_num_jets):
                self.hist[fc][jet_itr].Write()

# ------------------------------------------------------------------------------
class hMet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'met'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.hist_met_int         = {}
        self.hist_metrel_int      = {}

        self.hist_met_int_w_mu    = {}
        self.hist_metrel_int_w_mu = {}

        self.hist_met_noint       = {}
        self.hist_metrel_noint    = {}

        self.hist_met_diff        = {}
        self.hist_met_diff_w_mu   = {}

        self.hist_metrel_diff        = {}
        self.hist_metrel_diff_w_mu   = {}

        for fc in cutflow.flavor_channels:
            self.hist_met_int[fc] = ROOT.TH1F( 'h__%s__met_int' % fc
                                             , '%s int - %s; E_{T}^{miss,int} [GeV]' % (title, fc)
                                             , num_bins, x_min, x_max
                                             )
            self.hist_metrel_int[fc] = ROOT.TH1F( 'h__%s__metrel_int' % fc
                                             , '%s int - %s; E_{T}^{miss,rel,int} [GeV]' % (title, fc)
                                             , num_bins, x_min, x_max
                                             )

            self.hist_met_int_w_mu[fc] = ROOT.TH1F( 'h__%s__met_int_w_mu' % fc
                                                  , '%s int - %s; E_{T}^{miss,int+#mu} [GeV]' % (title, fc)
                                                  , num_bins, x_min, x_max
                                                  )
            self.hist_metrel_int_w_mu[fc] = ROOT.TH1F( 'h__%s__metrel_int_w_mu' % fc
                                                  , '%s int - %s; E_{T}^{miss,rel,int+#mu} [GeV]' % (title, fc)
                                                  , num_bins, x_min, x_max
                                                  )

            self.hist_met_noint[fc] = ROOT.TH1F( 'h__%s__met_noint' % fc
                                               , '%s no int - %s; E_{T}^{miss,no int} [GeV]' % (title, fc)
                                               , num_bins, x_min, x_max
                                               )
            self.hist_metrel_noint[fc] = ROOT.TH1F( 'h__%s__metrel_noint' % fc
                                               , '%s no int - %s; E_{T}^{miss,rel,no int} [GeV]' % (title, fc)
                                               , num_bins, x_min, x_max
                                               )

            self.hist_met_diff[fc] = ROOT.TH1F( 'h__%s__met_int_noint_diff' % fc
                                              , '%s diff - %s; E_{T}^{miss,no int} - E_{T}^{miss,int} [GeV]' % (title, fc)
                                              , num_bins, -x_max, x_max
                                              )
            self.hist_met_diff_w_mu[fc] = ROOT.TH1F( 'h__%s__met_int_noint_diff_w_mu' % fc
                                                   , '%s diff - %s; E_{T}^{miss,no int} - E_{T}^{miss,int+#mu} [GeV]' % (title, fc)
                                                   , num_bins, -x_max, x_max
                                                   )

            self.hist_metrel_diff[fc] = ROOT.TH1F( 'h__%s__metrel_int_noint_diff' % fc
                                                 , '%s diff - %s; E_{T}^{miss,rel,no int} - E_{T}^{miss,rel,int} [GeV]' % (title, fc)
                                                 , num_bins, -x_max, x_max
                                                 )
            self.hist_metrel_diff_w_mu[fc] = ROOT.TH1F( 'h__%s__metrel_int_noint_diff_w_mu' % fc
                                                      , '%s diff - %s; E_{T}^{miss,rel,no int} - E_{T}^{miss,rel,int+#mu} [GeV]' % (title, fc)
                                                      , num_bins, -x_max, x_max
                                                      )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        # interacting type met
        met_int = signal_objects['met']['int']
        metrel_int = signal_objects['met']['rel_int']

        self.hist_met_int[flavor_channel].Fill(met_int)
        self.hist_metrel_int[flavor_channel].Fill(metrel_int)


        # non interacting type met
        met_noint = signal_objects['met']['noint']
        metrel_noint = signal_objects['met']['rel_noint']
        self.hist_met_noint[flavor_channel].Fill(met_noint)
        self.hist_metrel_noint[flavor_channel].Fill(metrel_noint)

        # difference between interacting and non-interacting type mets
        self.hist_met_diff[flavor_channel].Fill( met_noint
                                               - met_int
                                               )

        self.hist_metrel_diff[flavor_channel].Fill( metrel_noint
                                                  - metrel_int
                                                  )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            self.hist_met_int[fc].Write()
            self.hist_metrel_int[fc].Write()

            self.hist_met_int_w_mu[fc].Write()
            self.hist_metrel_int_w_mu[fc].Write()

            self.hist_met_noint[fc].Write()
            self.hist_metrel_noint[fc].Write()

            self.hist_met_diff[fc].Write()
            self.hist_met_diff_w_mu[fc].Write()

            self.hist_metrel_diff[fc].Write()
            self.hist_metrel_diff_w_mu[fc].Write()

# ------------------------------------------------------------------------------
class hMll(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'mll'
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.hist = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = ROOT.TH1F( 'h__%s__mll' % fc
                                     , '%s - %s -- mll; m_{ll} [GeV]' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        mll = signal_objects['mll']/1000
        self.hist[flavor_channel].Fill(mll)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            self.hist[fc].Write()

