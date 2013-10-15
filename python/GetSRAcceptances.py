#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import array

import ROOT
import rootlogon
import metaroot

import yaml

import DiLeptonCutflow as cutflow

# ==============================================================================
def readInputConfig(in_file_name):
    print 'finput config file: %s' % in_file_name
    config_file = open(in_file_name)
    config_dict = yaml.load(config_file)

    file_handles = []
    for cd in config_dict:
        print 'cd: %s' % cd
        file_handles.append( FileHandle( cd['label']
                                       , cd['files']
                                       , cd['color']
                                       , cd['shape']
                                       , cd['line']
                                       , cd['x_label']
                                       , cd['x_max']
                                       )
                           )
        print ''
    return file_handles

# ==============================================================================
# input files:
class FileHandle(object):
    def __init__(self, title, file_list, color, shape, line, x_label, x_max):
        self.title = title
        self.file_list = file_list
        self.color = color
        self.shape = shape
        self.line  = line
        self.x_label = x_label
        self.x_max   = x_max

# ------------------------------------------------------------------------------
def getSignalAcceptance(in_file_name, signal_region_function):
    f = ROOT.TFile.Open(in_file_name)
    tree = f.Get('truth')
    total_num_events = tree.GetEntries()
    num_sr = 0
    for i, event in enumerate(tree):
        if i % 100 == 0:
            print 'Event %d of %d' % (i, total_num_events)
        num_el = event.el_n
        num_mu = event.mu_staco_n

        signal_objects = cutflow.doObjectSelection( event
                                          , lep_pt_cut  = 10.e3
                                          , lep_eta_cut = 2.4
                                          , jet_pt_cut  = 20.e3
                                          , jet_eta_cut = 2.7
                                          # , verbose = True
                                          )
        flavor_channel = cutflow.getFlavorChannel(signal_objects)

        if signal_region_function(signal_objects):
            num_sr += 1
    return float(num_sr)/total_num_events

# ------------------------------------------------------------------------------
def getSignalAcceptanceCurve(file_handle, signal_region_function):
    print '-------------------------'
    x_value = []
    acc = []
    for fl in file_handle.file_list:
        print fl
        print fl['file']
        x_value.append(fl['x_value'])
        acc.append(getSignalAcceptance(fl['file'], signal_region_function))
    print x_value
    print acc

    acc_curve = ROOT.TGraph( len(x_value)
                           , array.array('d', x_value)
                           , array.array('d', acc)
                           )
    acc_curve.SetName('g__acc__%s' % file_handle.title)
    return acc_curve

# ------------------------------------------------------------------------------
def getSignalAcceptanceCanvas(file_handle_list, signal_region_function):
    # create legend
    leg_x1 = 0.80
    leg_x2 = 0.98
    leg_y1 = 0.98
    leg_y2 = leg_y1-(0.06*len(file_handle_list))

    # big_leg_x1 = 0.05
    # big_leg_x2 = 0.95
    # big_leg_y1 = 0.98
    # big_leg_y2 = leg_y1-(0.08*len(file_handle_list))

    leg     = ROOT.TLegend(leg_x1    , leg_y1    , leg_x2    , leg_y2    )
    # big_leg = ROOT.TLegend(big_leg_x1, big_leg_y1, big_leg_x2, big_leg_y2)

    # get acceptance values for each model, and add to TGraph
    acc_canv = ROOT.TCanvas('c__acc')
    acc_curve_list = []
    y_max = 0
    for fh in file_handle_list:
        label = fh.title
        acc_curve = getSignalAcceptanceCurve(fh, signal_region_function)
        acc_curve.SetLineStyle(fh.line)
        acc_curve.SetLineColor(fh.color)
        acc_curve.SetLineWidth(4)
        acc_curve.SetMarkerStyle(fh.shape)
        acc_curve.SetMarkerColor(fh.color)
        acc_curve_list.append(acc_curve)

        leg.AddEntry(acc_curve, fh.title, 'alp')
        # big_leg.AddEntry(acc_curve, fh.title, 'alp')

        # check maximum y value
        local_max = acc_curve.GetMaximum()
        y_max = max(y_max, local_max)

    # Draw TGraphs onto canvas
    drawn = False
    acc_canv.cd()
    for acl in acc_curve_list:
        if not drawn:
            acl.Draw('ALP')
            acl.GetHistogram().GetXaxis().SetTitle(fh.x_label)
            acl.GetHistogram().GetYaxis().SetTitle('fractional acceptance')

            # acc_curve.GetHistogram().GetYaxis().SetRangeUser(0.,1.1)
            # acl.GetHistogram().GetYaxis().SetRangeUser(0.,0.005)
            acc_curve.GetHistogram().GetYaxis().SetRangeUser(0.,1.2*y_max)

            acl.Draw('ALP')
            print acc_curve
            drawn = True
        else:
            acl.Draw('LPSAME')
    leg.Draw()

    return {'canv':acc_canv, 'curves':acc_curve_list, 'leg':leg}

# ------------------------------------------------------------------------------
def main():
    out_file = ROOT.TFile.Open('acceptances.canv.root', 'RECREATE')
    out_file.cd()

    config_file_name = sys.argv[1]
    file_handles = readInputConfig(config_file_name)

    # do sr-mt2a
    acc_curve_dict = getSignalAcceptanceCanvas(file_handles, cutflow.isSROSMT2a)
    out_file.cd()
    acc_curve_dict['canv'].Write('srmt2a')

    # do sr-mt2b
    acc_curve_dict = getSignalAcceptanceCanvas(file_handles, cutflow.isSROSMT2b)
    out_file.cd()
    acc_curve_dict['canv'].Write('srmt2b')

    # do sr-mt2c
    acc_curve_dict = getSignalAcceptanceCanvas(file_handles, cutflow.isSROSMT2c)
    out_file.cd()
    acc_curve_dict['canv'].Write('srmt2c')

    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    ROOT.gROOT.SetBatch()
    main()

