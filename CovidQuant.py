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
PROGRESS_PERCENTAGE = 10
PROGRESS_NUM_READS = 100000
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
    print("[%s] %s" % (get_time(), s), file=stderr, end=end); stderr.flush()

# helper tree node class
class Node:
    def __init__(self):
        self.children = dict() # children[(pos,delim,val)] = child Node
        self.read_count = 0

# helper tree class, adapted from TreeSwift
class Tree:
    def __init__(self):
        self.root = Node()
        self.pos_val_to_nodes = dict() # map (position, nucleotide) tuple to nodes in the tree

    def traverse_edges_levelorder(self):
        q = deque(); q.append(self.root)
        while len(q) != 0:
            n = q.popleft()
            for l in n.children:
                yield (n, n.children[l], l) # (u,v,l)
            q.extend(n.children.values())
    
    def traverse_nodes_postorder(self):
        s1 = deque(); s2 = deque(); s1.append(self.root)
        while len(s1) != 0:
            n = s1.pop(); s2.append(n); s1.extend(n.children.values())
        while len(s2) != 0:
            yield s2.pop()

    def traverse_nodes_preorder(self):
        s = deque(); s.append(self.root)
        while len(s) != 0:
            n = s.pop(); yield n; s.extend(n.children.values())

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

    def add_read_counts_from_sam(self, sam):
        for i,alignment in enumerate(sam.fetch(until_eof=True)):
            # progress and check that alignment is primary mapped
            i_plus_1 = i+1
            if i_plus_1 % PROGRESS_NUM_READS == 0:
                print_log("Parsing alignment %d..." % i_plus_1)
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue

            # valid alignment, so increment counts
            aligned_pairs = alignment.get_aligned_pairs(matches_only=True) # ignore indels (Pangolin seems to)
            seq = alignment.query_sequence
            for read_pos, ref_pos in aligned_pairs:
                if read_pos is None or ref_pos is None:
                    raise ValueError("Encountered None in SAM alignment")
                pos_val = (ref_pos, seq[read_pos])
                if pos_val in tree.pos_val_to_nodes:
                    nodes_to_inc = tree.pos_val_to_nodes[pos_val]
                    count_to_inc = 1. / len(nodes_to_inc)
                    for node in nodes_to_inc:
                        node.read_count += count_to_inc

    def quantify_lineages(self):
        # count num lineage nodes below each node (including the node itself)
        for node in self.traverse_nodes_postorder():
            node.num_lineage_nodes_below = sum(child.num_lineage_nodes_below for child in node.children.values())
            if hasattr(node, 'lineage'):
                node.num_lineage_nodes_below += 1

        # quantify lineage abundances
        abundance = dict(); num_nodes_per_lineage = dict()
        for node in self.traverse_nodes_preorder():
            # trickle down count at current node to children if applicable
            orig_read_count = node.read_count
            for child in node.children.values():
                count_to_trickle = orig_read_count * child.num_lineage_nodes_below / node.num_lineage_nodes_below
                node.read_count -= count_to_trickle; child.read_count += count_to_trickle
            # if this node is a lineage node, update abundance
            if hasattr(node, 'lineage'):
                if node.lineage in abundance:
                    abundance[node.lineage] += node.read_count
                    num_nodes_per_lineage[node.lineage] += 1
                else:
                    abundance[node.lineage] = node.read_count
                    num_nodes_per_lineage[node.lineage] = 1

        # normalize lineage abundances and return
        for lineage in abundance:
            abundance[lineage] /= num_nodes_per_lineage[lineage]
        tot = sum(abundance.values())
        for lineage in abundance:
            abundance[lineage] /= tot
        return abundance

# parse PangoLEARN decision tree rules txt
def parse_pangolearn_rules_txt(s):
    tree = Tree(); lines = s.splitlines(); num_lines = len(lines); finished = PROGRESS_PERCENTAGE
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
            print_log("Finished parsing PangoLEARN rule line %d (%d%% done)..." % (i, int(percent_done)))
            finished += PROGRESS_PERCENTAGE
    print_log("Finished parsing PangoLEARN rules file")
    return tree

# main execution
if __name__ == "__main__":
    # parse args
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_alignment', required=True, type=str, help="Input Alignment File (SAM/BAM)")
    parser.add_argument('-p', '--pangolearn_rules', required=True, type=str, help="Input PangoLEARN Decision Tree Rules") # https://github.com/cov-lineages/pangoLEARN/blob/master/pangoLEARN/data/decision_tree_rules.zip
    parser.add_argument('-o', '--output_abundances', required=False, type=str, default='stdout', help="Output Abundances (TSV)")
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
    print_log("Parsing alignment file...")
    sam = AlignmentFile(args.input_alignment, ALIGNMENT_EXT_TO_QUAL[input_alignment_ext])
    tree.add_read_counts_from_sam(sam)

    # quantify lineage abundances
    print_log("Quantifying lineage abundances...")
    abundances = tree.quantify_lineages()
    print_log("Finished quantifying lineage abundances")

    # output lineage abundances
    print_log("Outputting lineage abundances...")
    if args.output_abundances.lower() == 'stdout':
        from sys import stdout as out
    else:
        out = open(args.output_abundances, 'w')
    out.write("Lineage\tAbundance\n")
    for lineage in sorted(abundances.keys(), key=lambda x: abundances[x], reverse=True):
        out.write("%s\t%f\n" % (lineage, abundances[lineage]))
    out.close()
