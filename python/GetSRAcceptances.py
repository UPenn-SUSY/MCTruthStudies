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
region_list = [ 'srss1'
              , 'srss2'
              , 'srss3'
              , 'srss4'
              , 'srss5'
              , 'srmt2a'
              , 'srmt2b'
              , 'srmt2c'
              ]

# ==============================================================================
def readInputConfig(in_file_name):
    print 'input config file: %s' % in_file_name
    config_file = open(in_file_name)
    config_dict = yaml.load(config_file)

    file_handles = []
    for cd in config_dict:
        file_handles.append( FileHandle( cd['label']
                                       , cd['files']
                                       , cd['color']
                                       , cd['shape']
                                       , cd['line']
                                       , cd['x_label']
                                       , cd['x_max']
                                       )
                           )
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
def getSignalAcceptance(in_file_name):
    print '----------------------------------------'
    print 'Getting signal acceptance for input file: %s' % in_file_name
    f = ROOT.TFile.Open(in_file_name)
    tree = f.Get('truth')
    total_num_events = float(tree.GetEntries())
    acc_in_regions = {r:0. for r in region_list}
    for i, event in enumerate(tree):
        if i % 100 == 0:
            print 'Event %d of %d' % (i, total_num_events)
        # if i > 500: break

        ewk_cutflow = cutflow.EwkCutFlow(event)
        if not ewk_cutflow.valid_cutflow: continue

        for region in region_list:
            if region in ewk_cutflow.regions:
                acc_in_regions[region] += 1
        # if 'srmt2a' in ewk_cutflow.regions:
        #     acc_in_regions['srmt2a'] += 1

    for air in acc_in_regions: acc_in_regions[air] /= total_num_events
    return acc_in_regions

# ------------------------------------------------------------------------------
def getSignalAcceptanceCurve(file_handle):
    x_value = []
    acc_in_regions = {r:[] for r in region_list}

    for fl in file_handle.file_list:
        x_value.append(fl['x_value'])
        acc = getSignalAcceptance(fl['file'])
        for r in region_list:
            acc_in_regions[r].append(acc[r])

    acc_curve = {r:ROOT.TGraph( len(x_value)
                               , array.array('d', x_value)
                               , array.array('d', acc_in_regions[r])
                               )
                for r in region_list
                }
    for r in region_list:
        acc_curve[r].SetName('g__acc__%s__%s' % (file_handle.title, r))

    max_values = {r:max(acc_in_regions[r]) for r in region_list}

    return {'curve':acc_curve, 'max':max_values}

# ------------------------------------------------------------------------------
# def getSignalAcceptanceCanvas(file_handle_list, signal_region_function):
def getSignalAcceptanceCanvas(file_handle_list):
    # create legend
    leg_x1 = 0.20
    leg_x2 = 0.38
    leg_y1 = 0.98
    leg_y2 = leg_y1-(0.06*len(file_handle_list))

    # big_leg_x1 = 0.05
    # big_leg_x2 = 0.95
    # big_leg_y1 = 0.98
    # big_leg_y2 = leg_y1-(0.08*len(file_handle_list))

    leg = ROOT.TLegend(leg_x1, leg_y1, leg_x2, leg_y2)
    leg.SetFillColor(0)
    # big_leg = ROOT.TLegend(big_leg_x1, big_leg_y1, big_leg_x2, big_leg_y2)
    # big_leg.SetFillColor(0)

    # get acceptance values for each model, and add to TGraph
    acc_canv = {r:ROOT.TCanvas('c__acc__%s' % r) for r in region_list}
    acc_curve_list = {r:[] for r in region_list}
    y_max = {r:0. for r in region_list}
    for fh in file_handle_list:
        label = fh.title
        # acc_curve = getSignalAcceptanceCurve(fh, signal_region_function)
        acc_curve = getSignalAcceptanceCurve(fh)
        for r in region_list:
            acc_curve['curve'][r].SetLineStyle(fh.line)
            acc_curve['curve'][r].SetLineColor(fh.color)
            acc_curve['curve'][r].SetLineWidth(4)
            acc_curve['curve'][r].SetMarkerStyle(fh.shape)
            acc_curve['curve'][r].SetMarkerColor(fh.color)
            acc_curve['curve'][r].SetMarkerSize(1.5)
            acc_curve_list[r].append(acc_curve['curve'][r])

        leg.AddEntry(acc_curve['curve'][region_list[0]], fh.title, 'alp')
        # big_leg.AddEntry(acc_curve, fh.title, 'alp')

        # check maximum y value
        local_max = {r:acc_curve['max'][r] for r in region_list}
        for r in region_list:
            y_max[r] = max(y_max[r], local_max[r])

    # Draw TGraphs onto canvas
    for r in region_list:
        drawn = False
        acc_canv[r].cd()
        for acl in acc_curve_list[r]:
            if not drawn:
                acl.Draw('ALP')
                acl.GetHistogram().GetXaxis().SetTitle(fh.x_label)
                acl.GetHistogram().GetYaxis().SetTitle('fractional acceptance')

                acl.GetHistogram().GetYaxis().SetRangeUser(0.,1.4*y_max[r])

                acl.Draw('ALP')
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

    # get acceptance curves and write to file
    acc_curve_dict = getSignalAcceptanceCanvas(file_handles)
    out_file.cd()
    for r in region_list:
        acc_curve_dict['canv'][r].Write(r)

    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    ROOT.gROOT.SetBatch()
    main()

