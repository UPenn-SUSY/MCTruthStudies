#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import ROOT
import rootlogon
# import metaroot

import yaml

import HistMerger
import FileHandle

# ------------------------------------------------------------------------------
def main():
    config_file_name = sys.argv[1]
    files = FileHandle.readInputConfig(config_file_name)
    print 'FILES: %s' % files

    out_file = ROOT.TFile.Open('truth_compare.canv.root', 'RECREATE')
    out_file.cd()
    list_of_hists = files[0].getListOfHists()
    num_hists = len(list_of_hists)
    for i, loh in enumerate(list_of_hists):
        if 'fc_all__' not in loh: continue
        for is_log in [True, False]:
            print 'hist (%d of %d): %s' % (i, num_hists, loh)

            hm = HistMerger.HistMerger( files, loh, False , is_log)
            hm.canvas.Write()

            hm_norm = HistMerger.HistMerger( files, loh, True, is_log)
            hm_norm.canvas.Write()

            if 'channel' in loh:
                if hm.leg_canvas is not None:
                    hm.leg_canvas.Write('c__leg')

            hm.canvas.Clear()
            hm_norm.canvas.Clear()
            if hm.leg_canvas is not None:
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

