#!/usr/bin/env python3

import re
import json
import argparse

parser = argparse.ArgumentParser(
    description="Reformat output from TRF",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "input", metavar="INPUT", type=str,
    help="Input dat or ngs file from TRF")
parser.add_argument(
    "--format", type=str, default="ngs",
    help="Format of input file, either `dat` or `ngs`")
parser.add_argument(
    "-o", "--out", type=str, default="test_out",
    help="Path to write output file")
parser.add_argument(
    "-f", "--outfmt", type=str, default="tsv",
    help="""
    Format for output file, either `tsv` or `gff`. Note that the start and end
    coordinates in both TSV and GFF files are 1-based end-inclusive, as they are
    in the original TRF output.
    """)

NGSFIELDS = ["start", "end", "per_len", "num_copies", "cons_len",
             "adj_pc_match", "adj_pc_indel", "aln_score",
             "pc_A", "pc_C", "pc_G", "pc_T", "entropy", "cons_seq",
             "seq_region", "flank_left", "flank_right"]

def parse_ngs(filename):
    out = []
    fh = open(filename, "r")
    current_record = None
    for line in fh:
        line = line.rstrip()
        headmatch = re.match(r"\@(.+)", line)
        if headmatch:
            current_record = headmatch.group(1)
            # print(f"Record {current_record}")
        else:
            linesplit = line.split(" ")
            if len(linesplit) != 17:
                raise Exception(f"Incorrect number of fields in line:\n{line}")
            else:
                linedict = {NGSFIELDS[i] : linesplit[i] for i in range(len(NGSFIELDS))}
                linedict['seqid'] = current_record
                out.append(linedict)
    fh.close()
    return(out)


def report_tsv(lines, filename):
    """Reformat parsed TRF output into TSV format

    Parameters
    ----------
    lines : list
        Output from parse_ngs()
    filename : str
        Path to filename to write output
    """
    report_fields = ["seqid", "start", "end", "per_len", "num_copies",
                     "cons_len", "adj_pc_match", "adj_pc_indel", "aln_score",
                     "pc_A", "pc_C", "pc_G", "pc_T", "entropy", "cons_seq"]
    # Skip reporting of extracted seq region and flanks, as these can be
    # extracted from the reference Fasta file
    fh = open(filename, "w")
    fh.write("\t".join(report_fields) + "\n")
    for linedict in lines:
        fh.write("\t".join([linedict[field] for field in report_fields]) + "\n")
    fh.close()


# main

args = parser.parse_args()

if args.input:
    lines = None
    if args.format == "ngs":
        lines = parse_ngs(args.input)
    if lines:
        if args.outfmt == "tsv":
            report_tsv(lines, args.out)
        elif args.outfmt == "gff":
            # placeholder
            print("not yet implemented")

