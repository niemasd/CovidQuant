#!/usr/bin/env python3

# imports
from collections import deque
from datetime import datetime
from gzip import open as gopen
from os.path import isfile
from pysam import AlignmentFile
from sys import stderr, stdin
from zipfile import ZipFile

# useful constants
VERSION = '1.0.0'
RULE_DELIMS = ['==', '!=']
NUCS = {'A','C','G','T'}
ALIGNMENT_EXT_TO_QUAL = {
    'bam': 'rb',
    'cram': 'rc',
    'sam': 'r',
}

# error messages
ERROR_FILE_NOT_FOUND = "File not found"
ERROR_INVALID_ALIGNMENT_EXTENSION = "Invalid alignment file extension"
ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES = "Invalid PangoLEARN Decision Tree Rules file"

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# log printer
def print_log(s, end='\n'):
    print("[%s] %s" % (get_time(), s), file=stderr, end=end)

# helper tree node class
class Node:
    def __init__(self):
        self.children = dict() # children[(pos,delim,val)] = child Node
        self.read_count = 0

# helper tree class
class Tree:
    def __init__(self):
        self.root = Node()
        self.pos_val_to_nodes = dict() # map (position, nucleotide) tuple to nodes in the tree

    def traverse_edges(self):
        q = deque(); q.append(self.root)
        while len(q) != 0:
            n = q.popleft()
            for l in n.children:
                yield (n, n.children[l], l) # (u,v,l)
            q.extend(n.children.values())

    def add_pangolin_rule(self, pangolin_rule):
        curr_node = self.root
        lineage, decisions = [v.strip() for v in pangolin_rule.split('\t')]
        decisions = decisions.split(',')
        for decision in decisions:
            delim = None; delim_start = None
            for d in RULE_DELIMS:
                try:
                    delim_start = decision.index(d); delim = d; break
                except:
                    pass
            if delim is None:
                raise ValueError("Invalid decision: %s" % decision)
            pos = int(decision[:delim_start])
            val = decision[delim_start+len(delim):].strip().strip("'")
            label = (pos, delim, val)
            if label in curr_node.children:
                curr_node = curr_node.children[label]
            else:
                child = Node(); curr_node.children[label] = child
                if delim == '==':
                    pos_val_pairs = [(pos,val.upper())]
                elif delim == '!=':
                    pos_val_pairs = [(pos,nuc) for nuc in NUCS-{val.upper()}]
                else:
                    raise ValueError("Unsupported: %s" % delim)
                for pos_val in pos_val_pairs:
                    if pos_val in self.pos_val_to_nodes:
                        self.pos_val_to_nodes[pos_val].add(child)
                    else:
                        self.pos_val_to_nodes[pos_val] = {child}
                curr_node = child
        curr_node.lineage = lineage

# parse PangoLEARN decision tree rules txt
def parse_pangolearn_rules_txt(s, progress_percentage=10):
    tree = Tree(); lines = s.splitlines(); num_lines = len(lines); finished = progress_percentage
    print_log("Parsing %d lines from PangoLEARN rules file..." % num_lines)
    for i,line in enumerate(lines):
        if line.startswith('['):
            continue
        l = line.strip()
        if len(l) == 0:
            continue
        tree.add_pangolin_rule(l)
        percent_done = 100 * (i+1) / num_lines
        if percent_done >= finished:
            print_log("Finished parsing PangoLEARN rule line %d (%d%% done)" % (i, int(percent_done)))
            finished += progress_percentage
    print_log("Finished parsing PangoLEARN rules file")
    exit()
    return tree

# main execution
if __name__ == "__main__":
    # parse args
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_alignment', required=True, type=str, help="Input Alignment File (SAM/BAM)")
    parser.add_argument('-p', '--pangolearn_rules', required=True, type=str, help="Input PangoLEARN Decision Tree Rules") # https://github.com/cov-lineages/pangoLEARN/blob/master/pangoLEARN/data/decision_tree_rules.zip
    args = parser.parse_args()
    for fn in [args.input_alignment, args.pangolearn_rules]:
        if not isfile(fn):
            raise ValueError("%s: %s" % (ERROR_FILE_NOT_FOUND, fn))
    input_alignment_ext = args.input_alignment.lower().split('.')[-1].strip()
    if input_alignment_ext not in ALIGNMENT_EXT_TO_QUAL:
        raise ValueError("%s: %s" % (ERROR_INVALID_ALIGNMENT_EXTENSION, input_alignment_ext))
    print_log("CovidQuant v%s" % VERSION)
    print_log("Input Alignment: %s" % args.input_alignment)
    print_log("Input PangoLEARN Decision Tree Rules: %s" % args.pangolearn_rules)

    # load decision tree file
    print_log("Loading PangoLEARN decision tree rules...")
    if args.pangolearn_rules.lower().endswith('.zip'):
        rules_zip = ZipFile(args.pangolearn_rules)
        fns = [fn.strip() for fn in rules_zip.namelist() if fn.strip().lower().endswith('tree_rules.txt')]
        if len(fns) != 1:
            raise ValueError(ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES)
        tree = parse_pangolearn_rules_txt(rules_zip.open(fns[0]).read().decode())
    elif args.pangolearn_rules.lower().endswith('.txt.gz'):
        tree = parse_pangolearn_rules_txt(gopen(args.pangolearn_rules).read().decode())
    elif args.pangolearn_rules.lower().endswith('.txt'):
        tree = parse_pangolearn_rules_txt(open(args.pangolearn_rules).read())
    else:
        raise ValueError(ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES)

    # parse SAM/BAM and update tree
    sam = AlignmentFile(args.input_alignment, ALIGNMENT_EXT_TO_QUAL[input_alignment_ext])
    for alignment in sam.fetch(until_eof=True):
        if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
            continue
        for read_pos, ref_pos, val in alignment.get_aligned_pairs(with_seq=True):
            pos_val = (ref_pos, val)
            if pos_val in tree.pos_val_to_nodes:
                for node in tree.pos_val_to_nodes[pos_val]:
                    node.read_count += 1
