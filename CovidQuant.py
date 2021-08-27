#!/usr/bin/env python3

# imports
from collections import deque
from datetime import datetime
from gzip import open as gopen
from json import load as jload
from os.path import isfile
from pysam import AlignmentFile
from sys import argv, setrecursionlimit, stderr, stdin
from urllib.request import urlopen
from zipfile import ZipFile

# useful constants
VERSION = '1.0.0'
PROGRESS_PERCENTAGE = 10
PROGRESS_NUM_READS = 100000
RULE_DELIMS = ['==', '!=']
NUCS = {'A','C','G','T','-'}
ALIGNMENT_EXT_TO_QUAL = {
    'bam': 'rb',
    'cram': 'rc',
    'sam': 'r',
}
INDENT = '  '
setrecursionlimit(10000)
CORRECT_NUC = { # https://github.com/cov-lineages/pangoLEARN/issues/13#issue-902259254
    '-': 'A',
    'A': 'C',
    'C': 'G',
    'G': 'T',
    'T': '-',
}

# error messages
ERROR_FILE_NOT_FOUND = "File not found"
ERROR_INVALID_ALIGNMENT_EXTENSION = "Invalid alignment file extension"
ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES = "Invalid PangoLEARN Decision Tree Rules file"

# convert a ViralMSA version string to a tuple of integers
def parse_version(s):
    return tuple(int(v) for v in s.split('.'))

# update CovidQuant to the newest version
def update_covidquant():
    tags = jload(urlopen(RELEASES_URL))
    newest = max(tags, key=lambda x: parse_version(x['name']))
    if parse_version(newest['name']) <= parse_version(VERSION):
        print("CovidQuant is already at the newest version (%s)" % VERSION); exit(0)
    old_version = VERSION; new_version = newest['name']
    url = 'https://raw.githubusercontent.com/niemasd/CovidQuant/%s/CovidQuant.py' % newest['commit']['sha']
    content = urlopen(url).read()
    try:
        with open(__file__, 'wb') as f:
            f.write(content)
    except PermissionError:
        print("ERROR: Received a permission error when updating CovidQuant. Perhaps try running as root?", file=stderr); exit(1)
    print("Successfully updated CovidQuant %s --> %s" % (old_version, new_version)); exit(0)

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# log printer
def print_log(s, end='\n'):
    print("[%s] %s" % (get_time(), s), file=stderr, end=end); stderr.flush()

# normalize dictionary such that values sum to 1
def normalize(d):
    out = dict(); tot = sum(d.values())
    for k in d:
        out[k] = d[k]/tot
    return out

# perform levelorder traversal on decision tree edges
def traverse_edges_levelorder(root):
    q = deque(); q.append(root)
    while len(q) != 0:
        n = q.popleft()
        for l in n:
            if l != '_LINEAGE':
                yield (n, n[l], l) # (u,v,l)
                q.append(n[l])

# perform postorder traversal on decision tree nodes
def traverse_nodes_postorder(root):
    s1 = deque(); s2 = deque(); s1.append(root)
    while len(s1) != 0:
        n = s1.pop(); s2.append(n); s1.extend(set(n.values())-{'_LINEAGE'})
    while len(s2) != 0:
        yield s2.pop()

# perform preorder traversal on decision tree nodes
def traverse_nodes_preorder(root):
    s = deque(); s.append(root)
    while len(s) != 0:
        n = s.pop(); yield n
        for l in n:
            s.append(n[l])

# draw the decision tree as lines of decisions
def print_decision_tree(root, num_indent=0):
    for l in root:
        if l == '_LINEAGE':
            print('%s- Lineage: %s' % (num_indent*INDENT, root[l]))
        else:
            pos,delim,symbol = l; print('%s- %s%s%s' % (num_indent*INDENT, pos, delim, symbol)); print_decision_tree(root[l], num_indent+1)

# print a GraphViz visualization of a decision tree
def graphviz(root):
    print("digraph G {")
    for u,v,l in traverse_edges_levelorder(root):
        u_label = id(u); v_label = id(v)
        if '_LINEAGE' in u:
            u_label = '%s (%s)' % (id(u),u['_LINEAGE'])
        if '_LINEAGE' in v:
            v_label = '%s (%s)' % (id(v),v['_LINEAGE'])
        print('"%s" -> "%s" [label="%s"];' % (u_label,v_label,l))
    print("}")

# parse PangoLEARN decision tree rules txt
def parse_pangolearn_rules_txt(s):
    root = dict(); lines = s.splitlines(); num_lines = len(lines); finished = PROGRESS_PERCENTAGE
    print_log("Parsing %d lines from PangoLEARN rules file..." % num_lines)
    for i,line in enumerate(lines):
        # skip lines starting with '[' and empty lines
        if line.startswith('['):
            continue
        l = line.strip()
        if len(l) == 0:
            continue

        # add this rule to the decision tree
        curr_node = root; lineage, decisions = [v.strip() for v in l.split('\t')]
        for decision in decisions.split(','):
            delim = None; delim_start = None
            for d in RULE_DELIMS:
                try:
                    delim_start = decision.index(d); delim = d; break
                except:
                    pass
            if delim is None:
                raise ValueError("Invalid decision: %s" % decision)
            pos = int(decision[:delim_start]); val = decision[delim_start+len(delim):].strip().strip("'")
            #label = (pos, delim, val)
            label = (pos, delim, CORRECT_NUC[val]) # TODO REPLACE WITH (pos,delim,val) ONCE PANGOLEARN FIXES THEIR TREE
            if label not in curr_node:
                curr_node[label] = dict()
            curr_node = curr_node[label]
        curr_node['_LINEAGE'] = lineage

        # output progress update
        percent_done = 100 * (i+1) / num_lines
        if percent_done >= finished:
            print_log("Finished parsing PangoLEARN rule line %d (%d%% done)..." % (i, int(percent_done))); finished += PROGRESS_PERCENTAGE
    print_log("Finished parsing PangoLEARN rules file")
    return root

# count symbols at each reference position from SAM/BAM/CRAM file
def read_counts_from_sam(sam):
    counts = dict() # counts[ref_pos][symbol] = count
    for i,alignment in enumerate(sam.fetch(until_eof=True)):
        # progress and check that alignment is primary mapped
        i_plus_1 = i+1
        if i_plus_1 % PROGRESS_NUM_READS == 0:
            print_log("Parsing alignment %d..." % i_plus_1)
        if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
            continue
        
        # valid alignment, so increment counts
        aligned_pairs = alignment.get_aligned_pairs(); seq = alignment.query_sequence
        for read_pos, ref_pos in aligned_pairs:
            # ignore insertions wrt reference (Pangolin seems to)
            if ref_pos is None:
                continue

            # get current read symbol
            if read_pos is None:
                symbol = '-'
            else:
                symbol = seq[read_pos].upper()

            # ignore N (Pangolin replaces it with ref symbol, but Pangolin uses a single consensus)
            if symbol == 'N':
                continue

            # increment count
            if ref_pos not in counts:
                counts[ref_pos] = dict()
            if symbol not in counts[ref_pos]:
                counts[ref_pos][symbol] = 0
            counts[ref_pos][symbol] += 1
    return counts

# best lineage (for sanity check)
def best_lineage(tree_root, counts):
    curr_node = tree_root
    while '_LINEAGE' not in curr_node:
        curr_counts = list()
        for label in curr_node:
            if label == '_LINEAGE':
                continue
            ref_pos, delim, symbol = label; label_count = 0
            if ref_pos in counts:
                if delim == '==' and symbol in counts[ref_pos]:
                    label_count = counts[ref_pos][symbol]
                elif delim == '!=':
                    for v in counts[ref_pos]:
                        if v != symbol:
                            label_count += counts[ref_pos][v]
            curr_counts.append((label_count,label))
        curr_node = curr_node[sorted(curr_counts, reverse=True)[0][1]]
    return curr_node['_LINEAGE']

# quantify lineages
def quantify_lineages(tree_root, counts):
    out = dict()
    def score(curr_node, curr_prob):
        if '_LINEAGE' in curr_node:
            if curr_node['_LINEAGE'] not in out:
                out[curr_node['_LINEAGE']] = 0
            out[curr_node['_LINEAGE']] += curr_prob; return
        for label in curr_node:
            ref_pos, delim, symbol = label
            if ref_pos in counts:
                label_count = 0
                if delim == '==' and symbol in counts[ref_pos]:
                    label_count = counts[ref_pos][symbol]
                elif delim == '!=':
                    for v in counts[ref_pos]:
                        if v != symbol:
                            label_count += counts[ref_pos][v]
                score(curr_node[label], curr_prob * label_count / sum(counts[ref_pos].values()))
    score(tree_root, 1.)
    return out

# main execution
if __name__ == "__main__":
    # check if user wants to update CovidQuant
    if '-u' in argv or '--update' in argv:
        update_covidquant()

    # parse args
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input_alignment', required=True, type=str, help="Input Alignment File (SAM/BAM)")
    parser.add_argument('-p', '--pangolearn_rules', required=True, type=str, help="Input PangoLEARN Decision Tree Rules") # https://github.com/cov-lineages/pangoLEARN/blob/master/pangoLEARN/data/decision_tree_rules.zip
    parser.add_argument('-o', '--output_abundances', required=False, type=str, default='stdout', help="Output Abundances (TSV)")
    parser.add_argument('--assign_top_lineage', action='store_true', help="Optional: Assign top lineage using decision tree (will only be included in log output)")
    parser.add_argument('--output_decision_coverage', required=False, type=str, default=None, help="Optional: Output coverage of each PangoLEARN genome position (TSV)")
    parser.add_argument('-u', '--update', action="store_true", help="Update CovidQuant (current version: %s)" % VERSION)
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
        tree_root = parse_pangolearn_rules_txt(rules_zip.open(fns[0]).read().decode())
    elif args.pangolearn_rules.lower().endswith('.txt.gz'):
        tree_root = parse_pangolearn_rules_txt(gopen(args.pangolearn_rules).read().decode())
    elif args.pangolearn_rules.lower().endswith('.txt'):
        tree_root = parse_pangolearn_rules_txt(open(args.pangolearn_rules).read())
    else:
        raise ValueError(ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES)

    # parse SAM/BAM to count symbols
    print_log("Parsing alignment file...")
    sam = AlignmentFile(args.input_alignment, ALIGNMENT_EXT_TO_QUAL[input_alignment_ext])
    counts = read_counts_from_sam(sam)

    # output coverage at each PangoLEARN genome position
    if args.output_decision_coverage is not None:
        print_log("Outputting coverage at each PangoLEARN decision genome position")
        out = open(args.output_decision_coverage, 'w')
        out.write("Position\tTotal Reads\t# A\t% A\t# C\t% C\t# G\t% G\t# T\t% T\t# -\t% -\n")
        for ref_pos in sorted({l[0] for u,v,l in traverse_edges_levelorder(tree_root)}):
            tmp = {'A':0, 'C':0, 'G':0, 'T':0, '-':0}
            for k in tmp:
                if k in counts[ref_pos]:
                    tmp[k] = counts[ref_pos][k]
            tot = sum(tmp.values())
            out.write("%d\t%d\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\n" % (ref_pos, tot, tmp['A'], tmp['A']/tot, tmp['C'], tmp['C']/tot, tmp['G'], tmp['G']/tot, tmp['T'], tmp['T']/tot, tmp['-'], tmp['-']/tot))
        out.close()

    # find the single best lineage using the simple decision tree algorithm
    if args.assign_top_lineage:
        print_log("Finding single best lineage...")
        print_log("Best lineage: %s" % best_lineage(tree_root,counts))

    # quantify lineage abundances
    print_log("Quantifying lineage abundances...")
    probabilities = quantify_lineages(tree_root, counts)
    abundances = normalize(probabilities)
    print_log("Finished quantifying lineage abundances")

    # output lineage abundances
    print_log("Outputting lineage abundances...")
    if args.output_abundances.lower() == 'stdout':
        from sys import stdout as out
    else:
        out = open(args.output_abundances, 'w')
    out.write("Lineage\tPath Probability\tAbundance\n")
    for lineage in sorted(abundances.keys(), key=lambda x: abundances[x], reverse=True):
        out.write("%s\t%f\t%f\n" % (lineage, probabilities[lineage], abundances[lineage]))
    out.close()
