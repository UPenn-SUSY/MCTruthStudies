#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import ROOT

# ==============================================================================
BUFFER = 1.15

# ==============================================================================
class HistMerger(object):
    def __init__(self, file_handles, hist_name, normalize = True, log = False, scale_to_xsec=False):
        self.name = hist_name
        self.normalize = normalize
        self.log = log
        self.scale_to_xsec = scale_to_xsec

        self.file_names = [fh.title for fh in file_handles]
        self.hists = [ fh.getHist( hist_name
                                 , normalize
                                 , 0 if not scale_to_xsec else fh.xsec
                                 )
                       for fh in file_handles
                     ]
        self.setStyle(file_handles, self.hists)

        self.title = '%s;%s;%s' % ( self.hists[0].GetTitle()
                                  , self.hists[0].GetXaxis().GetTitle()
                                  , self.hists[0].GetYaxis().GetTitle()
                                  )
        if isinstance(self.hists[0], ROOT.TH2):
            self.prep2DHistStack()
        elif isinstance(self.hists[0], ROOT.TH1):
            self.prep1DHistStack()
        elif isinstance(self.hists[0], ROOT.TGraph):
            self.prepGraphStack()

    # ------------------------------------------------------------------------------
    def prep1DHistStack(self):
        self.stack = ROOT.THStack('s_%s' % self.name, self.title)

        leg_x1 = 0.20
        leg_x2 = 0.43
        leg_y1 = 0.93
        leg_y2 = leg_y1-(0.05*len(self.hists))

        stat_box_x = 0.5
        stat_box_y = 0.90
        stat_box_text_size = 0.035

        big_leg_x1 = 0.05
        big_leg_x2 = 0.95
        big_leg_y1 = 0.98
        big_leg_y2 = leg_y1-(0.08*len(self.hists))

        # make legend
        self.legend = ROOT.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)
        self.big_legend = ROOT.TLegend(big_leg_x1, big_leg_y1, big_leg_x2, big_leg_y2)

        # set background color of legend
        self.legend.SetFillColor(0)
        self.big_legend.SetFillColor(0)

        self.legend.SetBorderSize(0)
        self.big_legend.SetBorderSize(0)

        # list of labels used as a stat box
        self.stat_box = []

        for i, h in enumerate(self.hists):
            self.stack.Add(h, 'HIST')
            self.legend.AddEntry(h, self.file_names[i], 'L')
            self.big_legend.AddEntry(h, self.file_names[i], 'L')

            num_entries = h.GetEntries()
            mean = getMeanSafe(h)
            rms = getRMSSafe(h)
            tmp_label = ROOT.TText( stat_box_x
                                  , stat_box_y
                                  , "Entries: %d  Mean: %.2f  RMS: %.2f" % ( num_entries
                                                                           , mean
                                                                           , rms
                                                                           )
                                  )
            tmp_label.SetTextSize(stat_box_text_size)
            tmp_label.SetTextColor(h.GetLineColor())
            tmp_label.SetNDC()
            tmp_label.SetX(stat_box_x)
            tmp_label.SetY(stat_box_y)
            self.stat_box.append(tmp_label)
            stat_box_y -= 0.05

        canvas_name = self.name.replace('h__', 'c__')
        if self.normalize:
            canvas_name = canvas_name + "__norm"
        if self.log:
            canvas_name = canvas_name + "__log_y"
        if self.scale_to_xsec > 0.:
            canvas_name = canvas_name + "__xsec"
        self.canvas = ROOT.TCanvas(canvas_name)
        if self.log:
            self.canvas.SetLogy()

        # add padding to top of plots
        self.stack.Draw('nostack')
        old_min = self.stack.GetHistogram().GetMinimum()
        old_max = self.stack.GetHistogram().GetMaximum()
        if self.log:
            new_max = math.pow( 10
                              , ( math.log(old_max, 10)
                                + (math.log(old_max, 10)-math.log(old_min, 10))*BUFFER/(1+BUFFER)
                                )
                              )
            self.stack.SetMaximum(new_max)
        else:
            new_max = (old_max + (old_max-old_min)*BUFFER/(1+BUFFER))
            self.stack.SetMaximum(new_max)

        # Draw plots and legend
        self.stack.Draw('nostack')
        if len(self.hists) <= 5:
            self.legend.Draw()

            for sb in self.stat_box:
                sb.Draw()

        self.leg_canvas = ROOT.TCanvas('%s__leg' % canvas_name)
        self.big_legend.Draw()

    # ------------------------------------------------------------------------------
    def prep2DHistStack(self):
        canvas_name = self.name.replace('h__', 'c__')
        if self.normalize:
            canvas_name = canvas_name + "__norm"
        if self.log:
            canvas_name = canvas_name + "__log_z"
        if self.scale_to_xsec > 0.:
            canvas_name = canvas_name + "__xsec"
        self.canvas = ROOT.TCanvas(canvas_name)

        self.labels = []

        if self.log:
            self.canvas.SetLogz()

        if len(self.hists) == 1:
            # add padding to top of plots
            self.hists[0].Draw('COLZ')
            old_min = self.hists[0].GetYaxis().GetXmin()
            old_max = self.hists[0].GetYaxis().GetXmax()
            new_max = (old_max + (old_max-old_min)*BUFFER/(1+BUFFER))
            self.hists[0].GetYaxis().SetRangeUser(old_min, new_max)

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
            elif len(self.hists) < 13:
                self.canvas.Divide(3,4)
            else:
                self.canvas.Divide(3,3)
                print 'Uh Oh! HistMerger cannot handle histogram lists of > 12 entries'

            for i, h in enumerate(self.hists):
                if i > 9: break
                self.canvas.cd(i+1)

                # add padding to top of plots
                self.hists[i].Draw('COLZ')
                old_min = self.hists[i].GetYaxis().GetXmin()
                old_max = self.hists[i].GetYaxis().GetXmax()
                new_max = (old_max + (old_max-old_min)*BUFFER/(1+BUFFER))
                self.hists[i].SetAxisRange(old_min, new_max, "Y")

                self.hists[i].Draw('COLZ')

                this_label = self.getLabel(i)
                this_label.Draw()
                self.labels.append(this_label)

        self.leg_canvas = None

    # ------------------------------------------------------------------------------
    def prepGraphStack(self):
        print 'prepGraphStack()'
        self.stack = []
        # # ROOT.THStack('s_%s' % self.name, self.title)

        leg_x1 = 0.20
        leg_x2 = 0.43
        leg_y1 = 0.93
        leg_y2 = leg_y1-(0.05*len(self.hists))

        # stat_box_x = 0.5
        # stat_box_y = 0.90
        # stat_box_text_size = 0.035

        big_leg_x1 = 0.05
        big_leg_x2 = 0.95
        big_leg_y1 = 0.98
        big_leg_y2 = leg_y1-(0.08*len(self.hists))

        # make legend
        self.legend = ROOT.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)
        self.big_legend = ROOT.TLegend(big_leg_x1, big_leg_y1, big_leg_x2, big_leg_y2)

        # set background color of legend
        self.legend.SetFillColor(0)
        self.big_legend.SetFillColor(0)

        self.legend.SetBorderSize(0)
        self.big_legend.SetBorderSize(0)

        # # list of labels used as a stat box
        # self.stat_box = []

        for i, h in enumerate(self.hists):
            # self.stack.Add(h, 'HIST')
            self.stack.append(h)
            # self.stack.Add(h, 'HIST')
            self.legend.AddEntry(h, self.file_names[i], 'L')
            self.big_legend.AddEntry(h, self.file_names[i], 'L')

        canvas_name = self.name.replace('h__', 'c__')
        if self.normalize:
            canvas_name = canvas_name + "__norm"
        if self.log:
            canvas_name = canvas_name + "__log_z"
        if self.scale_to_xsec > 0.:
            canvas_name = canvas_name + "__xsec"
        self.canvas = ROOT.TCanvas(canvas_name)
        if self.log:
            self.canvas.SetLogy()

        # # add padding to top of plots
        # self.stack.Draw('nostack')
        # old_min = self.stack.GetHistogram().GetMinimum()
        # old_max = self.stack.GetHistogram().GetMaximum()
        # if self.log:
        #     new_max = math.pow( 10
        #                       , ( math.log(old_max, 10)
        #                         + (math.log(old_max, 10)-math.log(old_min, 10))*BUFFER/(1+BUFFER)
        #                         )
        #                       )
        #     self.stack.SetMaximum(new_max)
        # else:
        #     new_max = (old_max + (old_max-old_min)*BUFFER/(1+BUFFER))
        #     self.stack.SetMaximum(new_max)

        # Draw plots and legend
        old_min = None
        old_max = None
        new_min = None
        new_max = None
        for i, ss in enumerate(self.stack):
            ss.Draw('%s' % 'apc' if i == 0 else 'pc')

            this_min = ss.GetHistogram().GetMinimum()
            this_max = ss.GetHistogram().GetMaximum()

            if old_min == None or this_min < old_min: old_min = this_min
            if old_max == None or this_max > old_max: old_max = this_max

        # add padding to top of plots
        if old_max is None: old_max = 1
        if self.log:
            if old_max is None or old_max < 0: old_max = 1
            if old_min is None or old_min < 0: old_min = 1e-6
            print 'old max: %s' % old_max
            print 'old min: %s' % old_min
            new_max = math.pow( 10
                              , ( math.log(old_max, 10)
                                + (math.log(old_max, 10)-math.log(old_min, 10))*BUFFER/(1+BUFFER)
                                )
                              )
            new_min = math.pow( 10
                              , ( math.log(old_min, 10)
                                - (math.log(old_max, 10)-math.log(old_min, 10))*(BUFFER/10)/(1+(BUFFER/10))
                                )
                              )
        else:
            print 'old max: %s' % old_max
            print 'old min: %s' % old_min
            if old_max is None: old_max = 1
            if old_min is None: old_min = 0
            print 'old max: %s' % old_max
            print 'old min: %s' % old_min
            new_max = (old_max + (old_max-old_min)*BUFFER/(1+BUFFER))
            new_min = (old_min - (old_max-old_min)*(BUFFER/10)/(1+(BUFFER/10)))
            if old_min > 0 and new_min < 0: new_min = 0

            print 'new max: %s' % new_max
            print 'new min: %s' % new_min

        for i, ss in enumerate(self.stack):
            if i == 0:
                self.frame = ROOT.TH1I( 'h_frame_%s' % self.name
                                      , self.title
                                      , 1
                                      , ss.GetHistogram().GetXaxis().GetXmin()
                                      , ss.GetHistogram().GetXaxis().GetXmax()
                                      )
                self.frame.GetYaxis().SetRangeUser(new_min, new_max)
                self.frame.Draw()
            ss.Draw('pc')

        # self.stack.Draw('nostack')
        if len(self.hists) <= 5:
            self.legend.Draw()

        self.leg_canvas = ROOT.TCanvas('%s__leg' % canvas_name)
        self.big_legend.Draw()


    # --------------------------------------------------------------------------
    def setStyle(self, file_handles, hists):
        for fh,h in zip(file_handles, hists):
            h.SetMarkerStyle(fh.shape)
            h.SetMarkerColor(fh.color)
            h.SetLineStyle(fh.line)
            h.SetLineColor(fh.color)
            h.SetLineWidth(4)

    # --------------------------------------------------------------------------
    def getLabel(self, index):
        this_label = ROOT.TLatex(1,1, self.file_names[index])
        this_label.SetNDC()
        # this_label.SetX(0.45)
        this_label.SetX(0.17)
        this_label.SetY(0.96)
        # this_label.SetTextSize(0.05)
        this_label.SetTextSize(0.045)
        return this_label

# ------------------------------------------------------------------------------
def getMeanSafe(h):
    integral = h.Integral()
    if integral == 0: return 0
    weighted_sum = 0
    for i in xrange(h.GetNbinsX()):
        bin_center = h.GetXaxis().GetBinCenter(i+1)
        weighted_sum += bin_center*h.GetBinContent(i+1)
    return weighted_sum/integral

# ------------------------------------------------------------------------------
def getRMSSafe(h):
    mean = getMeanSafe(h)
    integral = h.Integral()
    if integral == 0: return 0
    weighted_sum = 0
    for i in xrange(h.GetNbinsX()):
        bin_center = h.GetXaxis().GetBinCenter(i+1)
        weighted_sum += bin_center*bin_center*h.GetBinContent(i+1)
    # print 'weighted sum: %s' % weighted_sum
    # print 'integral: %s' % integral
    # print 'weighted sum/integral: %s' % (weighted_sum/integral)
    # print 'mean^2: %s' % (mean*mean)
    # print '(weighted_sum/integral - mean*mean): %s' % (weighted_sum/integral - mean*mean)
    if weighted_sum/integral < mean*mean: return 0
    return math.sqrt(weighted_sum/integral - mean*mean)
