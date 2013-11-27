#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import ROOT
import rootlogon

# ==============================================================================
class HistMerger(object):
    def __init__(self, file_handles, hist_name, normalize = True, log = False):
        self.name = hist_name
        self.file_names = [fh.title              for fh in file_handles]
        self.hists = [fh.getHist(hist_name, normalize) for fh in file_handles]
        self.setStyle(file_handles, self.hists)
        self.normalize = normalize
        self.log = log

        self.title = '%s;%s;%s' % ( self.hists[0].GetTitle()
                                  , self.hists[0].GetXaxis().GetTitle()
                                  , self.hists[0].GetYaxis().GetTitle()
                                  )
        if isinstance(self.hists[0], ROOT.TH2):
            self.prep2DHistStack()
        elif isinstance(self.hists[0], ROOT.TH1):
            self.prep1DHistStack()

    # ------------------------------------------------------------------------------
    def prep1DHistStack(self):
        self.stack = ROOT.THStack('s_%s' % self.name, self.title)

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

        canvas_name = self.name.replace('h__', 'c__')
        if self.normalize:
            canvas_name = canvas_name + "__norm"
        if self.log:
            canvas_name = canvas_name + "__log_y"
        self.canvas = ROOT.TCanvas(canvas_name)
        if self.log:
            self.canvas.SetLogy()
        self.stack.Draw('nostack')
        self.legend.Draw()

        self.leg_canvas = ROOT.TCanvas('%s__leg' % canvas_name)
        self.big_legend.Draw()

    # ------------------------------------------------------------------------------
    def prep2DHistStack(self):
        canvas_name = self.name.replace('h__', 'c__')
        if self.normalize:
            canvas_name = canvas_name + "__norm"
        if self.log:
            canvas_name = canvas_name + "__log_z"
        self.canvas = ROOT.TCanvas(canvas_name)

        self.labels = []

        if self.log:
            self.canvas.SetLogz()

        if len(self.hists) == 1:
            self.hists[0].Draw('COLZ')

            this_label = self.getLabel(0)
            this_label.Draw()
            self.labels.append(this_label)
        else:
            if len(self.hists) == 2:
                self.canvas.Divide(1,2)
            elif len(self.hists) < 5:
                self.canvas.Divide(2,2)
            elif len(self.hists) < 7:
                self.canvas.Divide(2,3)
            elif len(self.hists) < 10:
                self.canvas.Divide(3,3)
            else:
                self.canvas.Divide(3,3)
                print 'Uh Oh! HistMerger cannot handle histogram lists of > 9 entries'

            for i, h in enumerate(self.hists):
                if i > 9: break
                self.canvas.cd(i+1)
                self.hists[i].Draw('COLZ')

                this_label = self.getLabel(i)
                this_label.Draw()
                self.labels.append(this_label)

        self.leg_canvas = None

    def setStyle(self, file_handles, hists):
        for fh,h in zip(file_handles, hists):
            h.SetMarkerStyle(fh.shape)
            h.SetMarkerColor(fh.color)
            h.SetLineStyle(fh.line)
            h.SetLineColor(fh.color)
            h.SetLineWidth(4)

    def getLabel(self, index):
        this_label = ROOT.TLatex(1,1, self.file_names[index])
        this_label.SetNDC()
        this_label.SetX(0.45)
        this_label.SetY(0.80)
        this_label.SetTextSize(0.05)
        return this_label
