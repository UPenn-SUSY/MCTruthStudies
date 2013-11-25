#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import ROOT
import rootlogon
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
                                       )
                           )
    return file_handles

# ==============================================================================
# input files:
class FileHandle(object):
    def __init__(self, title, in_file, directory, color, shape, line):
        self.title = title
        self.color = color
        self.shape = shape
        self.line  = line

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

    def getHist(self, hist_name, normalize = True):
        this_hist_name = 'clone_%s' % hist_name
        if normalize:
            this_hist_name += '_norm'
        hist = self.directory.Get(hist_name).Clone(this_hist_name)
        if normalize and hist.Integral() != 0.:
            hist.Scale(1./hist.Integral())
        return hist

