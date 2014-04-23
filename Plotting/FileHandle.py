#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import ROOT
#import rootlogon
# import metaroot

import yaml

# ==============================================================================
def readInputConfig(in_file_name):
    print 'finput config file: %s' % in_file_name
    config_file = open(in_file_name)
    config_dict = yaml.load(config_file)

    file_handles = []
    for cd in config_dict:
        file_handles.append( FileHandle( cd['label']
                                       , cd['file']
                                       , cd['dir']
                                       , cd['color']
                                       , cd['shape']
                                       , cd['line']
                                       , cd['xsec']
                                       )
                           )
    return file_handles

# ==============================================================================
# input files:
class FileHandle(object):
    def __init__(self, title, in_file, directory, color, shape, line, xsec):
        self.title = title
        self.color = color
        self.shape = shape
        self.line  = line
        self.xsec  = xsec

        self.in_file_name = in_file
        self.directory_name = directory

        self.in_file = ROOT.TFile(in_file)
        if directory == '':
            self.directory = self.in_file
        else:
            self.directory = self.in_file.GetDirectory(directory)

        keys = [lok.GetName() for lok in self.directory.GetListOfKeys()]
        print keys
        channel_name = [lok.GetName() for lok in self.directory.GetListOfKeys() if lok.GetName().startswith('fc_all__flavor_channel')]
        print 'channel_name: %s' % channel_name
        assert len(channel_name) == 1
        channel_name = channel_name[0]

        print channel_name
        print 'channels' in keys
        self.scale = self.directory.Get(channel_name).Integral()
        self.scale = 1/self.scale if self.scale != 0. else 1.

    def getListOfHists(self):
        hist_list = [obj.GetName() for obj in self.directory.GetListOfKeys()]
        return hist_list

    def getHist(self, hist_name, normalize = True, scale_to_xsec = 0.):
        this_hist_name = 'clone_%s' % hist_name
        if normalize:
            this_hist_name += '_norm'
        this_hist_name += '__%s' % self.title
        hist = self.directory.Get(hist_name).Clone(this_hist_name)
        hist.Sumw2()

        if not isinstance(hist, ROOT.TH2):
            moveOverflowToLastBin(hist)

        if normalize and hist.Integral() != 0.:
            hist.Scale(1./hist.Integral())
        if scale_to_xsec > 0.:
            hist.Scale(scale_to_xsec)
        return hist

# ------------------------------------------------------------------------------
def moveOverflowToLastBin(h, x_min=None, x_max=None):
    total_bins = h.GetNbinsX()
    total_entries = h.GetEntries()

    if x_min is None: x_min = h.GetXaxis().GetXmin()
    if x_max is None: x_max = h.GetXaxis().GetXmax()

    x_bins = []
    for i in xrange(total_bins):
        x_bins.append(h.GetBinLowEdge(i+1))

    # find total underflow
    underflow = h.GetBinContent(0)
    # if we want to truncate additional bins, loop over these bins
    min_bin = 1
    for i, x in enumerate(x_bins):
        if x < x_min:
            underflow += h.GetBinContent(i+1)
        else:
            min_bin = i+1
            break

    # move underflow to min_bin:
    for uf_bin in xrange(min_bin):
        h.SetBinContent(uf_bin, 0.)
    h.Fill(x_bins[min_bin-1], underflow)

    # find total overflow
    overflow = h.GetBinContent(total_bins+1)
    max_bin = 1
    for i, x in enumerate(x_bins):
        if x < x_max:
            max_bin = i+1
            continue
        overflow += h.GetBinContent(i+1)

    # move overflow to max_bin:
    for of_bin in xrange(max_bin+1, total_bins+2):
        h.SetBinContent(of_bin, 0.)
    h.Fill(x_bins[max_bin-1], overflow)

    h.SetEntries(total_entries)

