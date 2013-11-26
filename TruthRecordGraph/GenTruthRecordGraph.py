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
    out_file = file('event.dot', 'w')
    out_file.write('graph truth_record {\n')
    out_file.write('  size="30,30";\n')

    out_file.write(constructGraph(event, 0, 0, max_depth))

    out_file.write('}\n')
    out_file.close()

# ------------------------------------------------------------------------------
def getStyleString(current_pdgid):
    style_string = ''
    color = {'text':'black', 'fill':'white'}
    if abs(current_pdgid) >= 1e6:
        style_string = ', style=filled, color=lightblue'
    if abs(current_pdgid) == 5:
        style_string = ', style=filled, color=springgreen'
    if abs(current_pdgid) == 11 or abs(current_pdgid) == 13:
        style_string = ', style=filled, color=orange'
    if abs(current_pdgid) == 21:
        style_string = ', style=filled, color=snow3'
    if abs(current_pdgid) == 1 or abs(current_pdgid) == 2:
        style_string = ', style=filled, color=violet'

    return style_string

# ------------------------------------------------------------------------------
def constructGraph(event, current_index, current_depth, max_depth):
    # label_block = ''
    labels = {}
    draw_label = {}
    for mc_index in xrange(event.mc_barcode.size()):
        current_barcode     = event.mc_barcode.at(mc_index)
        current_pdgid       = event.mc_pdgId.at(mc_index)
        current_status_code = event.mc_status.at(mc_index)

        style_string = getStyleString(current_pdgid)

        draw_label[current_barcode] = False
        labels[current_barcode] = '  bc_%s [label="pdg: %s bc: %s sc: %s"%s];\n' % ( current_barcode
                                                                                   , current_pdgid
                                                                                   , current_barcode
                                                                                   , current_status_code
                                                                                   , style_string
                                                                                   )


    decay_block = ''
    # # max_itr = min(event.mc_barcode.size(), max_depth)
    max_itr = event.mc_barcode.size()
    for mc_index in xrange(max_itr):
        current_barcode = event.mc_barcode.at(mc_index)
        child_indices = event.mc_child_index.at(mc_index)
        parent_indices = event.mc_parent_index.at(mc_index)

        for ci in child_indices:
            child_barcode = event.mc_barcode.at(ci)
            decay_block += '  bc_%s -- bc_%s;\n' % (current_barcode, child_barcode)

            draw_label[current_barcode] = True
            draw_label[child_barcode] = True

    label_block = ''
    for l in labels:
        if draw_label[l]:
            label_block += labels[l]

    this_graph_block = '\n'
    this_graph_block += label_block
    this_graph_block += '\n'
    this_graph_block += decay_block
    this_graph_block += '\n'

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
