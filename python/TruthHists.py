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
                , selection_tag = 'all'
                ):
        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = ROOT.TH1F( '%s__flavor_channel' % self.selection_tag
                             , title
                             , len(cutflow.flavor_channels)
                             , -0.5
                             , len(cutflow.flavor_channels)-0.5
                             )
        for i, fc in enumerate(cutflow.flavor_channels):
            self.hist.GetXaxis().SetBinLabel(i+1, fc)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        bin_num = cutflow.flavor_channels.index(ewk_cutflow.flavor_channel)
        self.hist.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hDecayCategory(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'decay_category'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist_all = ROOT.TH1F( '%s__decay_category__all' % self.selection_tag
                                 , title
                                 , len(cutflow.decay_categories)
                                 , -0.5
                                 , len(cutflow.decay_categories) - 0.5
                                 )
        for i, dc in enumerate(cutflow.decay_categories):
            self.hist_all.GetXaxis().SetBinLabel(i+1, dc)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        bin_num = cutflow.decay_categories.index(ewk_cutflow.decay_category)
        self.hist_all.Fill(bin_num)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist_all, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'pt'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist_pt = []
        for num_lep in xrange(cutflow.max_num_leptons):
            self.hist_pt.append( ROOT.TH1F( '%s__lep_pt_%d' % (self.selection_tag, num_lep)
                                          , '%s - %s -- lepton %d; p_{T}^{%d} [GeV]' % (title, self.selection_tag, num_lep, num_lep)
                                          , num_bins, x_min, x_max
                                          )
                               )

        self.hist_diff = ROOT.TH1F( '%s__lep_pt_diff' % self.selection_tag
                                  , '%s - %s -- diff;p_{T}^{0} - p_{T}^{1} [GeV]' % (title, self.selection_tag)
                                  , num_bins
                                  , -x_max
                                  , x_max
                                  )
        self.hist_2d = ROOT.TH2F( '%s__lep_pt_2d' % self.selection_tag
                                , '%s - %s -- diff;p_{T}^{0};p_{T}^{1} [GeV]' % (title, self.selection_tag)
                                , num_bins
                                , x_min
                                , x_max
                                , num_bins
                                , x_min
                                , x_max
                                )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        # get lepton pt values
        lep_pt_list = []
        for el_pt in [el.pt for el in ewk_cutflow.signal['el']]:
            lep_pt_list.append(el_pt/1000)
        for mu_pt in [mu.pt for mu in ewk_cutflow.signal['mu']]:
            lep_pt_list.append(mu_pt/1000)

        # fill histograms
        for lep_it, lep_pt in enumerate(lep_pt_list):
            if lep_it == cutflow.max_num_leptons: break
            self.hist_pt[lep_it].Fill(lep_pt)
        if len(lep_pt_list) >= 2:
            self.hist_diff.Fill( lep_pt_list[0]
                               - lep_pt_list[1]
                               )
            self.hist_2d.Fill( lep_pt_list[0]
                             , lep_pt_list[1]
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        for num_lep in xrange(cutflow.max_num_leptons):
            writeToDir(self.hist_pt[num_lep], out_file, self.dir_name)
        writeToDir(self.hist_diff, out_file, self.dir_name)
        writeToDir(self.hist_2d  , out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hEta(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'eta'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 10
        x_min = -3.
        x_max = +3.

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist_eta = []
        for num_lep in xrange(cutflow.max_num_leptons):
            self.hist_eta.append( ROOT.TH1F( '%s__lep_eta_%d' % (self.selection_tag, num_lep)
                                           , '%s - %s -- lepton %d; #eta^{%d}' % (title, self.selection_tag, num_lep, num_lep)
                                           , num_bins, x_min, x_max
                                           )
                                )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        lep_eta_list = []
        for el_eta in [el.eta for el in ewk_cutflow.signal['el']]:
            lep_eta_list.append(el_eta)
        for mu_eta in [mu.eta for mu in ewk_cutflow.signal['mu']]:
            lep_eta_list.append(mu_eta)

        for lep_it, lep_eta in enumerate(lep_eta_list):
            if lep_it == cutflow.max_num_leptons: break
            self.hist_eta[lep_it].Fill(lep_eta)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        for num_lep in xrange(cutflow.max_num_leptons):
            writeToDir(self.hist_eta[num_lep], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hPtByMother(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'pt_by_mother'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist_pt = {}
        for mother_type in cutflow.mother_type_list:
            self.hist_pt[mother_type] = []
            for num_lep in xrange(cutflow.max_num_leptons):
                self.hist_pt[mother_type].append( ROOT.TH1F( '%s__lep_pt_%d__mother_%s' % (self.selection_tag, num_lep, mother_type)
                                                            , '%s - %s -- lepton %d - mother: %s; pt^{%d}' % (title, self.selection_tag, num_lep, mother_type, num_lep)
                                                            , num_bins, x_min, x_max
                                                            )
                                                 )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        lep_pt_list = []
        lep_mother_list = []
        for el in ewk_cutflow.signal['el']:
            lep_pt_list.append(el.pt/1000.)
            lep_mother_list.append(el.parent_pdgid)
        for mu in ewk_cutflow.signal['mu']:
            lep_pt_list.append(mu.pt/1000.)
            lep_mother_list.append(mu.parent_pdgid)

        for lep_it in xrange(len(lep_pt_list)):
            if lep_it == cutflow.max_num_leptons: break

            lep_pt = lep_pt_list[lep_it]
            lep_mother = lep_mother_list[lep_it]
            if lep_mother is None: continue

            mother = None
            if   abs(lep_mother) >= 1000011 and abs(lep_mother) <= 1000016: mother = 'sl'
            elif abs(lep_mother) == 1000023: mother = 'n2'
            elif abs(lep_mother) == 1000024: mother = 'c1'
            else: mother = 'none'
            self.hist_pt[mother][lep_it].Fill(lep_pt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        for mother_type in cutflow.mother_type_list:
            for num_lep in xrange(cutflow.max_num_leptons):
                writeToDir(self.hist_pt[mother_type][num_lep], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hEtaByMother(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'eta_by_mother'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 10
        x_min = -3.
        x_max = +3.

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist_eta = {}
        for mother_type in cutflow.mother_type_list:
            self.hist_eta[mother_type] = []
            for num_lep in xrange(cutflow.max_num_leptons):
                self.hist_eta[mother_type].append( ROOT.TH1F( '%s__lep_eta_%d__mother_%s' % (self.selection_tag, num_lep, mother_type)
                                                            , '%s - %s -- lepton %d - mother: %s; #eta^{%d}' % (title, self.selection_tag, num_lep, mother_type, num_lep)
                                                            , num_bins, x_min, x_max
                                                            )
                                                 )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        lep_eta_list = []
        lep_mother_list = []
        for el in ewk_cutflow.signal['el']:
            lep_eta_list.append(el.eta)
            lep_mother_list.append(el.parent_pdgid)
        for mu in ewk_cutflow.signal['mu']:
            lep_eta_list.append(mu.eta)
            lep_mother_list.append(mu.parent_pdgid)

        # for lep_it, lep_eta, lep_mother in enumerate(lep_eta_list, lep_mother_list):
        for lep_it in xrange(len(lep_eta_list)):
            if lep_it == cutflow.max_num_leptons: break

            lep_eta = lep_eta_list[lep_it]
            lep_mother = lep_mother_list[lep_it]
            if lep_mother is None: continue

            mother = None
            if   abs(lep_mother) >= 1000011 and abs(lep_mother) <= 1000016: mother = 'sl'
            elif abs(lep_mother) == 1000023: mother = 'n2'
            elif abs(lep_mother) == 1000024: mother = 'c1'
            else: mother = 'none'
            self.hist_eta[mother][lep_it].Fill(lep_eta)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        for mother_type in cutflow.mother_type_list:
            for num_lep in xrange(cutflow.max_num_leptons):
                writeToDir(self.hist_eta[mother_type][num_lep], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hNumJet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'num jet'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 10
        x_min = -0.5
        x_max = num_bins - 0.5

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist   = []
        self.hist = ROOT.TH1F( '%s__num_jets' % self.selection_tag
                             , '%s - %s; Jet Multiplicity' % (title, self.selection_tag)
                             , num_bins, x_min, x_max
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        self.hist.Fill(len(ewk_cutflow.signal['jet']))

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hJetPt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'jet pt'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = []
        for num_jet in xrange(cutflow.max_num_jets):
            self.hist.append( ROOT.TH1F( '%s__jet_pt_%d' % (self.selection_tag, num_jet)
                                       , '%s - %s -- jet %d; p_{T}^{%d} [GeV]' % (title, self.selection_tag, num_jet, num_jet)
                                       , num_bins, x_min, x_max
                                       )
                            )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        jet_pt_list = []
        for jet_pt in [jet.pt for jet in ewk_cutflow.signal['jet']]:
            jet_pt_list.append(jet_pt/1000)

        # jet_pt_list.sort(reverse=True)
        for jet_it, jet_pt in enumerate(jet_pt_list):
            if jet_it == cutflow.max_num_jets: break
            self.hist[jet_it].Fill(jet_pt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        for jet_itr in xrange(cutflow.max_num_jets):
            writeToDir(self.hist[jet_itr], out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hMet(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'met'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 50
        x_min = 0
        x_max = 500

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist_met_int = ROOT.TH1F( '%s__met_int' % self.selection_tag
                                     , '%s int - %s; E_{T}^{miss,int} [GeV]' % (title, self.selection_tag)
                                     , num_bins, x_min, x_max
                                     )
        self.hist_metrel_int = ROOT.TH1F( '%s__metrel_int' % self.selection_tag
                                        , '%s int - %s; E_{T}^{miss,rel,int} [GeV]' % (title, self.selection_tag)
                                        , num_bins, x_min, x_max
                                        )

        self.hist_met_noint = ROOT.TH1F( '%s__met_noint' % self.selection_tag
                                       , '%s no int - %s; E_{T}^{miss,no int} [GeV]' % (title, self.selection_tag)
                                       , num_bins, x_min, x_max
                                       )
        self.hist_metrel_noint = ROOT.TH1F( '%s__metrel_noint' % self.selection_tag
                                          , '%s no int - %s; E_{T}^{miss,rel,no int} [GeV]' % (title, self.selection_tag)
                                          , num_bins, x_min, x_max
                                          )

        self.hist_met_diff = ROOT.TH1F( '%s__met_int_noint_diff' % self.selection_tag
                                      , '%s diff - %s; E_{T}^{miss,no int} - E_{T}^{miss,int} [GeV]' % (title, self.selection_tag)
                                      , num_bins, -x_max, x_max
                                      )
        self.hist_metrel_diff = ROOT.TH1F( '%s__metrel_int_noint_diff' % self.selection_tag
                                         , '%s diff - %s; E_{T}^{miss,rel,no int} - E_{T}^{miss,rel,int} [GeV]' % (title, self.selection_tag)
                                         , num_bins, -x_max, x_max
                                         )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        # interacting type met
        met_int = ewk_cutflow.met['int']
        metrel_int = ewk_cutflow.met['rel_int']

        self.hist_met_int.Fill(met_int)
        self.hist_metrel_int.Fill(metrel_int)

        # non interacting type met
        met_noint = ewk_cutflow.met['noint']
        metrel_noint = ewk_cutflow.met['rel_noint']

        self.hist_met_noint.Fill(met_noint)
        self.hist_metrel_noint.Fill(metrel_noint)

        # difference between interacting and non-interacting type mets
        self.hist_met_diff.Fill( met_noint
                               - met_int
                               )
        self.hist_metrel_diff.Fill( metrel_noint
                                  - metrel_int
                                  )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist_met_int, out_file, self.dir_name)
        writeToDir(self.hist_metrel_int, out_file, self.dir_name)

        writeToDir(self.hist_met_noint, out_file, self.dir_name)
        writeToDir(self.hist_metrel_noint, out_file, self.dir_name)

        writeToDir(self.hist_met_diff, out_file, self.dir_name)
        writeToDir(self.hist_metrel_diff, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hMll(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'mll'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = ROOT.TH1F( '%s__mll' % self.selection_tag
                             , '%s - %s -- mll; m_{ll} [GeV]' % (title, self.selection_tag)
                             , num_bins, x_min, x_max
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        mll = ewk_cutflow.mll/1000
        self.hist.Fill(mll)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hMt2(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'mt2'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = {}
        self.hist = ROOT.TH1F( '%s__mt2' % self.selection_tag
                             , '%s - %s -- mt2; m_{T2} [GeV]' % (title, self.selection_tag)
                             , num_bins, x_min, x_max
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        mt2 = ewk_cutflow.mt2/1000.
        self.hist.Fill(mt2)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hPtll(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'ptll'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = ROOT.TH1F( '%s__ptll' % self.selection_tag
                             , '%s - %s -- ptll; p_{T}^{ll} [GeV]' % (title, self.selection_tag)
                             , num_bins, x_min, x_max
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        ptll = ewk_cutflow.ptll/1000.
        self.hist.Fill(ptll)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hEmmaMt(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'emma_mt'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 40
        x_min = 0
        x_max = 400

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = ROOT.TH1F( '%s__emma_mt' % self.selection_tag
                             , '%s - %s -- emma_mt; #sqrt{ (m_{ll})^{2} + (p_{T}^{ll})^{2} } [GeV]' % (title, self.selection_tag)
                             , num_bins, x_min, x_max
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        emma_mt = ewk_cutflow.emma_mt/1000.
        self.hist.Fill(emma_mt)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hSRSS(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'srss'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 7
        x_min = -0.5
        x_max = 6.5

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = ROOT.TH1F( '%s__srss' % self.selection_tag
                             , '%s - %s -- srss; SR channel' % (title, self.selection_tag)
                             , num_bins, x_min, x_max
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        fill_bin = []

        if 'srss1' in ewk_cutflow.regions: fill_bin.append(1)
        if 'srss2' in ewk_cutflow.regions: fill_bin.append(2)
        if 'srss3' in ewk_cutflow.regions: fill_bin.append(3)
        if 'srss4' in ewk_cutflow.regions: fill_bin.append(4)
        if 'srss5' in ewk_cutflow.regions: fill_bin.append(5)

        for fb in fill_bin:
            self.hist.Fill(fb)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

# ------------------------------------------------------------------------------
class hSROS(object):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def __init__( self
                , title = 'sros'
                , dir_name = ''
                , selection_tag = 'all'
                ):
        num_bins = 7
        x_min = -0.5
        x_max = 6.5

        self.dir_name = dir_name
        self.selection_tag = selection_tag

        self.hist = {}
        self.hist = ROOT.TH1F( '%s__sros' % self.selection_tag
                             , '%s - %s -- sros; SR channel' % (title, self.selection_tag)
                             , num_bins, x_min, x_max
                             )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def fill(self, ewk_cutflow):
        fill_bin = []

        if 'srmt2a' in ewk_cutflow.regions: fill_bin.append(1)
        if 'srmt2b' in ewk_cutflow.regions: fill_bin.append(2)
        if 'srmt2c' in ewk_cutflow.regions: fill_bin.append(3)

        for fb in fill_bin:
            self.hist.Fill(fb)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def writeToFile(self, out_file):
        writeToDir(self.hist, out_file, self.dir_name)

