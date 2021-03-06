#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import ROOT
import rootlogon
import metaroot

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
        # self.directory = None
        if directory == '':
            self.directory = self.in_file
        else:
            self.directory = self.in_file.GetDirectory(directory)

        keys = [lok.GetName() for lok in self.directory.GetListOfKeys()]
        print keys
        channel_name = [lok.GetName() for lok in self.directory.GetListOfKeys() if lok.GetName().startswith('dc_all__fc_all__flavor_channel')]
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

# ==============================================================================
class HistMerger(object):
    def __init__(self, file_handles, hist_name, normalize = True):
        self.name = hist_name
        self.file_names = [fh.title              for fh in file_handles]
        self.hists = [fh.getHist(hist_name, normalize) for fh in file_handles]
        self.setStyle(file_handles, self.hists)

        title = '%s;%s;%s' % ( self.hists[0].GetTitle()
                             , self.hists[0].GetXaxis().GetTitle()
                             , self.hists[0].GetYaxis().GetTitle()
                             )
        self.stack = ROOT.THStack('s_%s' % hist_name, title)

        leg_x1 = 0.80
        leg_x2 = 0.98
        leg_y1 = 0.98
        leg_y2 = leg_y1-(0.06*len(self.hists))

        big_leg_x1 = 0.05
        big_leg_x2 = 0.95
        big_leg_y1 = 0.98
        big_leg_y2 = leg_y1-(0.08*len(self.hists))

        self.legend = ROOT.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)
        self.big_legend = ROOT.TLegend(big_leg_x1, big_leg_y1, big_leg_x2, big_leg_y2)
        for i, h in enumerate(self.hists):
            self.stack.Add(h, 'HIST')
            self.legend.AddEntry(h, self.file_names[i], 'L')
            self.big_legend.AddEntry(h, self.file_names[i], 'L')

        canvas_name = hist_name.replace('h__', 'c__')
        if normalize:
            canvas_name = canvas_name + "_norm"
        self.canvas = ROOT.TCanvas(canvas_name)
        self.stack.Draw('nostack')
        self.legend.Draw()

        self.leg_canvas = ROOT.TCanvas('%s__leg' % canvas_name)
        self.big_legend.Draw()

    def setStyle(self, file_handles, hists):
        for fh,h in zip(file_handles, hists):
            h.SetMarkerStyle(fh.shape)
            h.SetMarkerColor(fh.color)
            h.SetLineStyle(fh.line)
            h.SetLineColor(fh.color)
            h.SetLineWidth(4)


# ------------------------------------------------------------------------------
def main():
    config_file_name = sys.argv[1]
    files = readInputConfig(config_file_name)

    out_file = ROOT.TFile.Open('truth_compare.canv.root', 'RECREATE')
    out_file.cd()
    list_of_hists = files[0].getListOfHists()
    num_hists = len(list_of_hists)
    for i, loh in enumerate(list_of_hists):
        print 'hist (%d of %d): %s' % (i, num_hists, loh)

        # uncomment if you want just a subset of plots
        # if 'srss' in loh: continue
        # if 'fc_ee_ss'  in loh: continue
        # if 'fc_em_ss'  in loh: continue
        # if 'fc_mm_ss'  in loh: continue
        # if 'fc_eee'    in loh: continue
        # if 'fc_eem'    in loh: continue
        # if 'fc_emm'    in loh: continue
        # if 'fc_mmm'    in loh: continue
        # if 'fc_multi'  in loh: continue
        # if 'fc_none'   in loh: continue
        # if 'dc_none'   in loh: continue
        # if 'dc_n2_n2'  in loh: continue
        # if 'dc_sl_n2'  in loh: continue
        # if 'dc_n2_sl'  in loh: continue
        # if 'dc_c1_n2'  in loh: continue
        # if 'dc_n2_c1'  in loh: continue
        # if 'lep_eta_2' in loh: continue
        # if 'lep_pt_2'  in loh: continue
        # if 'lep_eta_3' in loh: continue
        # if 'lep_pt_3'  in loh: continue
        if 'mother'    in loh and 'dc_all' not in loh: continue

        hm = HistMerger( files, loh, False )
        hm.canvas.Write()

        hm_norm = HistMerger( files, loh, True )
        hm_norm.canvas.Write()

        if 'channel' in loh:
            hm.leg_canvas.Write('c__leg')

        hm.canvas.Clear()
        hm_norm.canvas.Clear()
        hm.leg_canvas.Clear()
        for h in hm.hists:
            h.Delete()
        for h in hm_norm.hists:
            h.Delete()
    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    ROOT.gROOT.SetBatch()
    main()

