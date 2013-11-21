#!/usr/bin/env python

import sys
import os.path
import optparse
import time

import ROOT

# ------------------------------------------------------------------------------
def pickEvent(tree, event_number):
    for event in tree:
        if event.EventNumber == event_number:
            print 'found event %d!' % event.EventNumber
            return event
    return None

# ------------------------------------------------------------------------------
def printTruthRecord(event, max_depth):
    out_file = file('event.gv', 'w')
    out_file.write('graph truth_record {\n')

    out_file.write(constructGraph(event, 0, 0, max_depth))

    # out_file.write('  1 [label="a"];\n')
    # out_file.write('  2 [label="b"];\n')
    # out_file.write('  3 [label="c"];\n')
    # out_file.write('  4 [label="d"];\n')

    # out_file.write('  1 -- 2;\n')
    # out_file.write('  1 -- 3;\n')
    # out_file.write('  3 -- 4;\n')

    out_file.write('}\n')
    out_file.close()

# ------------------------------------------------------------------------------
def constructGraph(event, current_index, current_depth, max_depth):
    label_block = ''
    labels = {}
    draw_labels = {}
    # this_graph_block = '\n'
    for mc_index in xrange(event.mc_barcode.size()):
        # print 'checking mcindex: %s' % mc_index
        current_barcode = event.mc_barcode.at(mc_index)
        current_pdgid   = event.mc_pdgId.at(mc_index)

        label_block += '  bc_%s [label="%s"];\n' % (current_barcode, current_pdgid)

    decay_block = ''
    # # max_itr = min(event.mc_barcode.size(), max_depth)
    max_itr = event.mc_barcode.size()
    for mc_index in xrange(max_itr):
        current_barcode = event.mc_barcode.at(mc_index)
        child_indices = event.mc_child_index.at(mc_index)
        parent_indices = event.mc_parent_index.at(mc_index)
        # print 'child indices: %s' % child_indices
        print 'num child indices: %s' % child_indices.size()
        print 'num parent indices: %s' % parent_indices.size()

        for ci in child_indices:
            print ci
            child_barcode = event.mc_barcode.at(ci)
            label_block += '  bc_%s -- bc_%s;\n' % (current_barcode, child_barcode)


    this_graph_block = '\n'
    this_graph_block += label_block
    this_graph_block += '\n'
    this_graph_block += decay_block

    return this_graph_block

# ------------------------------------------------------------------------------
def main():
    in_file_name = sys.argv[1]
    event_number = int(sys.argv[2])

    print 'getting file: %s' % in_file_name
    print 'printing truth record for event: %d' % event_number

    f = ROOT.TFile(in_file_name)
    t = f.Get('truth')
    event = pickEvent(t, event_number)

    printTruthRecord(event, 10)

if __name__ == '__main__':
    main()
