#!/usr/bin/env python

import sys
from os.path import isfile, basename, splitext, isdir, join
import os


def process_graph(graph_content, graph_name, graph_fpath):

    with open(graph_fpath, 'w') as f:
        for line in graph_content:
            f.write(line)

    num_components = int(graph_content[0].split(':')[1])
    mass = 0
    for comp in graph_content[1: num_components + 1]:
        mass += float(comp.split()[-1])



    num_bonds = int(graph_content[num_components + 1].split(':')[1])

    graph_name = graph_name.replace(' ', '_')
    lib_string = ' '.join([graph_fpath, graph_name, str(mass), str(num_components), 'Graphs'])
    return lib_string


def main():
    '''
    IN: (1) txt file with list of rule fragmented graph (optionally prefixed by Peptide names)
        (2) new DB root dir (<db_root>)
    OUT: <db_root>/graphs_dir  with .graph files per each compound; and <db_root>/library.info.graphs
    '''

    if len(sys.argv) < 3:
        print "Usage: this.py graphs.txt db_root"
        sys.exit(1)

    graphs_fpath = sys.argv[1]
    db_root = sys.argv[2]
    graphs_dir = join(db_root, "graphs_dir")
    lib_info = join(db_root, "library.info.graphs")

    if not isfile(graphs_fpath):
        print "Graphs list not found!", graphs_fpath
        sys.exit(1)

    if not isdir(graphs_dir):
        os.makedirs(graphs_dir)

    lib_info_content = []
    with open(graphs_fpath) as f:
        def __start_new_graph(cur_idx):
            cur_graph_fpath = join(graphs_dir, "compound_%s.graph" % str(cur_idx).zfill(5))
            cur_graph_name = "Compound_%s" % str(cur_idx).zfill(5)
            cur_graph_content = []
            return cur_graph_fpath, cur_graph_name, cur_graph_content

        cur_idx = 0
        cur_graph_fpath, cur_graph_name, cur_graph_content = __start_new_graph(cur_idx)
        next_graph_name = ""
        num_lines_to_read = 0
        for line in f:
            if num_lines_to_read <= 0:
                if line.startswith('number of components'):
                    if cur_graph_content:
                        lib_info_content.append(process_graph(cur_graph_content, cur_graph_name, cur_graph_fpath))
                    cur_idx += 1
                    cur_graph_fpath, cur_graph_name, cur_graph_content = __start_new_graph(cur_idx)
                    num_lines_to_read = int(line.split(':')[1]) + 1
                    if next_graph_name:
                        cur_graph_name = next_graph_name
                        next_graph_name = ""
                elif line.startswith('number of bonds'):
                    num_lines_to_read = int(line.split(':')[1]) + 1
                else:
                    next_graph_name = line.strip()
                    continue
            cur_graph_content.append(line)
            num_lines_to_read -= 1
        if cur_graph_content:
            lib_info_content.append(process_graph(cur_graph_content, cur_graph_name, cur_graph_fpath))

    with open(lib_info, 'w') as f:
        for line in lib_info_content:
            f.write(line + '\n')


if __name__ == '__main__':
    main()