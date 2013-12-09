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
            for norm in [True, False]:
                for xsec in [True, False]:
                    print 'hist (%d of %d): %s' % (i, num_hists, loh)

                    hm = HistMerger.HistMerger( files, loh, norm , is_log, xsec)
                    hm.canvas.Write()

                    if 'channel' in loh and is_log is False and norm is False and xsec is False:
                        if hm.leg_canvas is not None:
                            hm.leg_canvas.Write('c__leg')

                    hm.canvas.Clear()
                    if hm.leg_canvas is not None:
                        hm.leg_canvas.Clear()
                    for h in hm.hists:
                        h.Delete()
    out_file.Close()

# ==============================================================================
if __name__ == '__main__':
    ROOT.gROOT.SetBatch()
    main()

