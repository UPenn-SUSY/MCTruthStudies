#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import ROOT
import rootlogon

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
