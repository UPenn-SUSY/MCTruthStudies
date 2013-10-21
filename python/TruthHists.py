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
def writeToDir(hist, out_file, dir_name):
    out_file.cd()
    if dir_name is not '':
        if out_file.Get(dir_name) == None:
            out_file.mkdir(dir_name)
        out_file.cd(dir_name)
    hist.Write()

# ------------------------------------------------------------------------------
class hFlavorChannels(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'flavor channels'
                , dir_name = ''
                ):
        self.dir_name = dir_name

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
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'pt'
                , dir_name = ''
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.dir_name = dir_name

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
                writeToDir(self.hist_pt[fc][num_lep], out_file, self.dir_name)
            writeToDir(self.hist_diff[fc], out_file, self.dir_name)
            writeToDir(self.hist_2d[fc]  , out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hEta(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'eta'
                , dir_name = ''
                ):
        num_bins = 10
        x_min = -3.
        x_max = +3.

        self.dir_name = dir_name

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
                writeToDir(self.hist_eta[fc][num_lep], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hNumJet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'num jet'
                , dir_name = ''
                ):
        num_bins = 10
        x_min = -0.5
        x_max = num_bins - 0.5

        self.dir_name = dir_name

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
            writeToDir(self.hist[fc], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hJetPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'jet pt'
                , dir_name = ''
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.dir_name = dir_name

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
                writeToDir(self.hist[fc][jet_itr], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hMet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'met'
                , dir_name = ''
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.dir_name = dir_name

        self.hist_met_int         = {}
        self.hist_metrel_int      = {}

        self.hist_met_noint       = {}
        self.hist_metrel_noint    = {}

        self.hist_met_diff        = {}
        self.hist_metrel_diff        = {}

        for fc in cutflow.flavor_channels:
            self.hist_met_int[fc] = ROOT.TH1F( 'h__%s__met_int' % fc
                                             , '%s int - %s; E_{T}^{miss,int} [GeV]' % (title, fc)
                                             , num_bins, x_min, x_max
                                             )
            self.hist_metrel_int[fc] = ROOT.TH1F( 'h__%s__metrel_int' % fc
                                             , '%s int - %s; E_{T}^{miss,rel,int} [GeV]' % (title, fc)
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
            self.hist_metrel_diff[fc] = ROOT.TH1F( 'h__%s__metrel_int_noint_diff' % fc
                                                 , '%s diff - %s; E_{T}^{miss,rel,no int} - E_{T}^{miss,rel,int} [GeV]' % (title, fc)
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
            writeToDir(self.hist_met_int[fc], out_file, self.dir_name)
            writeToDir(self.hist_metrel_int[fc], out_file, self.dir_name)

            writeToDir(self.hist_met_noint[fc], out_file, self.dir_name)
            writeToDir(self.hist_metrel_noint[fc], out_file, self.dir_name)

            writeToDir(self.hist_met_diff[fc], out_file, self.dir_name)
            writeToDir(self.hist_metrel_diff[fc], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hMll(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'mll'
                , dir_name = ''
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name

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
            writeToDir(self.hist[fc], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hMt2(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'mt2'
                , dir_name = ''
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name

        self.hist = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = ROOT.TH1F( 'h__%s__mt2' % fc
                                     , '%s - %s -- mt2; m_{T2} [GeV]' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        mt2 = signal_objects['mt2']/1000
        self.hist[flavor_channel].Fill(mt2)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            writeToDir(self.hist[fc], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hPtll(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'ptll'
                , dir_name = ''
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name

        self.hist = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = ROOT.TH1F( 'h__%s__ptll' % fc
                                     , '%s - %s -- ptll; p_{T}^{ll} [GeV]' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        ptll = signal_objects['ptll']/1000
        self.hist[flavor_channel].Fill(ptll)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            writeToDir(self.hist[fc], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hEmmaMt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'emma_mt'
                , dir_name = ''
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name

        self.hist = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = ROOT.TH1F( 'h__%s__emma_mt' % fc
                                     , '%s - %s -- emma_mt; #sqrt{ (m_{ll})^{2} + (p_{T}^{ll})^{2} } [GeV]' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        emma_mt = signal_objects['emma_mt']/1000
        self.hist[flavor_channel].Fill(emma_mt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            writeToDir(self.hist[fc], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hSRSS(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'srss'
                , dir_name = ''
                ):
        num_bins = 7
        x_min = -0.5
        x_max = 6.5

        self.dir_name = dir_name

        self.hist = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = ROOT.TH1F( 'h__%s__srss' % fc
                                     , '%s - %s -- srss; SR channel' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        fill_bin = []

        if cutflow.isSRSS1(signal_objects): fill_bin.append(1)
        if cutflow.isSRSS2(signal_objects): fill_bin.append(2)
        if cutflow.isSRSS3(signal_objects): fill_bin.append(3)
        if cutflow.isSRSS4(signal_objects): fill_bin.append(4)
        if cutflow.isSRSS5(signal_objects): fill_bin.append(5)

        if len(fill_bin) == 0:
            self.hist[flavor_channel].Fill(0)
        else:
            for fb in fill_bin:
                self.hist[flavor_channel].Fill(fb)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            writeToDir(self.hist[fc], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hSROS(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'sros'
                , dir_name = ''
                ):
        num_bins = 7
        x_min = -0.5
        x_max = 6.5

        self.dir_name = dir_name

        self.hist = {}
        for fc in cutflow.flavor_channels:
            self.hist[fc] = ROOT.TH1F( 'h__%s__sros' % fc
                                     , '%s - %s -- sros; SR channel' % (title, fc)
                                     , num_bins, x_min, x_max
                                     )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, flavor_channel, signal_objects, event):
        fill_bin = []

        if cutflow.isSROSMT2a(signal_objects): fill_bin.append(1)
        if cutflow.isSROSMT2b(signal_objects): fill_bin.append(2)
        if cutflow.isSROSMT2c(signal_objects): fill_bin.append(3)

        if len(fill_bin) == 0:
            self.hist[flavor_channel].Fill(0)
        else:
            for fb in fill_bin:
                self.hist[flavor_channel].Fill(fb)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        out_file.cd()
        for fc in cutflow.flavor_channels:
            writeToDir(self.hist[fc], out_file, self.dir_name)

