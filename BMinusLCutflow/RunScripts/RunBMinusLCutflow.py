#!/usr/bin/env python

import sys
import os.path
import optparse
import time
import math

import glob

import ROOT

import os

# ==============================================================================
print 'loading libraries'
ROOT.gSystem.Load('${BASE_WORK_DIR}/lib/libProgressBar.so')
ROOT.gSystem.Load('${BASE_WORK_DIR}/lib/libmt2.so')
ROOT.gSystem.Load('${BASE_WORK_DIR}/lib/libTruthRecordHelpers.so')
ROOT.gSystem.Load('${BASE_WORK_DIR}/lib/libTruthNtupleLooper.so')
ROOT.gSystem.Load('${BASE_WORK_DIR}/lib/libHistogramHandlers.so')
ROOT.gSystem.Load('${BASE_WORK_DIR}/lib/libBMinusLCutflow.so')

print 'done loading libraries'

# ------------------------------------------------------------------------------
def getFileListFromDir(file_path):
    print 'getting files from dir: %s' % file_path
    file_list = glob.glob('%s/*' % file_path)
    return file_list

# ------------------------------------------------------------------------------
def getTChain(file_list, tree_name):
    t = ROOT.TChain(tree_name)
    for fl in file_list:
        print 'Adding file: %s' % fl
        t.AddFile(fl)
    return t

# ------------------------------------------------------------------------------
def runBMinusLCutflow( file_list
                     , tree_name = 'truth'
                     ):
    print "Adding files to TChain"
    t = getTChain(file_list, tree_name)

    # ==============================================================================
    print 'Creating BMinusL::Cutflow object'
    cf = ROOT.BMinusL.Cutflow(t)

    cf.Loop();

    bmlcf.writeToFile();

    # ==============================================================================
    print ''
    print ''

# ------------------------------------------------------------------------------
def main():
    input_file_or_dir = sys.argv[1]
    print input_file_or_dir

    file_list = getFileListFromDir(input_file_or_dir)

    runBMinusLCutflow(file_list)

    print ''
    print ''

# ==============================================================================
if __name__ == '__main__':
    main()
