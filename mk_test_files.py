#!/usr/bin/env python3
"""
Author : Ken Youens-Clark <kyclark@gmail.com>
Date   : 2024-12-10
Purpose: Recreate test files
"""

import argparse
import sys
from functools import partial
from subprocess import getstatusoutput
from typing import NamedTuple

SUFR = "./target/release/sufr"
SEQ1 = "data/inputs/1.fa"
SEQ2 = "data/inputs/2.fa"
SEQ3 = "data/inputs/3.fa"
UNIPROT = "data/inputs/uniprot.fa"
LONG = "data/inputs/long_dna_sequence.fa"
ABBA = "data/inputs/abba.fa"
OUTDIR = "./data/expected"

class Args(NamedTuple):
    """ Command-line arguments """
    build: bool
    verbose: bool


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Recreate test files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-b',
                        '--build',
                        help='Build release',
                        action='store_true')

    parser.add_argument('-v',
                        '--verbose',
                        help='Verbose',
                        action='store_true')

    args = parser.parse_args()

    return Args(build=args.build, verbose=args.verbose)


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()

    runner = partial(run, args.verbose)

    if args.build:
        runner("cargo build --release")

    runner(f'{SUFR} create --dna -o {OUTDIR}/1.sufr {SEQ1}')

    runner(f'{SUFR} create --dna -o {OUTDIR}/2.sufr {SEQ2}')

    runner(f'{SUFR} create --dna -o {OUTDIR}/3.sufr {SEQ3}')

    runner(f'{SUFR} create --dna --sequence-delimiter N '
           f'-o {OUTDIR}/2d.sufr {SEQ2}')

    runner(f'{SUFR} create -o {OUTDIR}/abba.sufr {ABBA}')

    runner(f'{SUFR} create --dna --allow-ambiguity -o {OUTDIR}/1n.sufr {SEQ1}')

    runner(f'{SUFR} create --dna --allow-ambiguity -o {OUTDIR}/2n.sufr {SEQ2}')

    runner(f'{SUFR} create --dna --ignore-softmask -o {OUTDIR}/2s.sufr {SEQ2}')

    runner(f'{SUFR} create --dna --allow-ambiguity --ignore-softmask '
           f'-o {OUTDIR}/2ns.sufr {SEQ2}')

    runner(f'{SUFR} create {LONG} --dna -o {OUTDIR}/long_dna_sequence.sufr')

    runner(f'{SUFR} create {LONG} --dna --allow-ambiguity '
           f'-o {OUTDIR}/long_dna_sequence_allow_ambiguity.sufr')

    runner(f'{SUFR} create -o {OUTDIR}/uniprot.sufr {UNIPROT}')

    runner(f'{SUFR} create -s 10111011 -o {OUTDIR}/uniprot-masked.sufr {UNIPROT}')


# --------------------------------------------------
def run(verbose: bool, cmd: str):
    """ Run command """

    if verbose:
        print(cmd)

    rv, out = getstatusoutput(cmd)
    if rv != 0:
        sys.exit(out)

# --------------------------------------------------
if __name__ == '__main__':
    main()
