#!/usr/bin/env python3

# imports
from gzip import open as gopen
from sys import stdin
from zipfile import ZipFile

# error messages
ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES = "Invalid PangoLEARN Decision Tree Rules"

# parse PangoLEARN decision tree rules txt
def parse_pangolearn_rules_txt(s):
    lines = [l.strip() for l in s.splitlines() if not l.startswith('[') and len(l.strip()) != 0]
    return None # TODO PARSE LINES

# main execution
if __name__ == "__main__":
    # parse args
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input PangoLEARN Decision Tree Rules")
    args = parser.parse_args()

    # load decision tree file
    if args.input.lower() == 'stdin':
        rules = parse_pangolearn_rules_txt(stdin.read())
    elif args.input.lower().endswith('.zip'):
        rules_zip = ZipFile(args.input)
        fns = [fn.strip() for fn in rules_zip.namelist() if fn.strip().lower().endswith('tree_rules.txt')]
        if len(fns) != 1:
            raise ValueError(ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES)
        rules = parse_pangolearn_rules_txt(rules_zip.open(fns[0]).read().decode())
    elif args.input.lower().endswith('.txt.gz'):
        rules = parse_pangolearn_rules_txt(gopen(args.input).read().decode())
    elif args.input.lower().endswith('.txt'):
        rules = parse_pangolearn_rules_txt(open(args.input).read())
    else:
        raise ValueError(ERROR_INVALID_PANGOLEARN_DECISION_TREE_RULES)
