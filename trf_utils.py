#!/usr/bin/env python3

import re
import json
import argparse
import tempfile
# import logging
from subprocess import run
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Reformat output from TRF",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    "--input", type=str,
    help="Input dat or ngs file from TRF")
parser.add_argument(
    "--format", type=str, default="ngs",
    help="Format of input file, either `dat` or `ngs`")
parser.add_argument(
    "--min_length", type=int, default=100,
    help="Minimum length of a repeat region to report")
parser.add_argument(
    "--no_overlaps", action="store_true",
    help="""
    Report GFF without overlaps: if features overlap, pick highest scoring
    feature. In case of tied score, pick feature with the most repeat copies.
    """)
parser.add_argument(
    "-o", "--out", type=str, default="test_out",
    help="Path to write output file")
parser.add_argument(
    "-f", "--outfmt", type=str, default="gff",
    help="""
    Format for output file, either `tsv` or `gff`. Note that the start and end
    coordinates in both TSV and GFF files are 1-based end-inclusive, as they are
    in the original TRF output.
    """)

# logging.basicConfig(level=logging.DEBUG)

NGSFIELDS = ["start", "end", "per_len", "num_copies", "cons_len",
             "adj_pc_match", "adj_pc_indel", "aln_score",
             "pc_A", "pc_C", "pc_G", "pc_T", "entropy", "cons_seq",
             "seq_region", "flank_left", "flank_right"]


def parse_ngs(filename):
    out = defaultdict(list)
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
                if int(linedict['end']) - int(linedict['start']) + 1 > args.min_length:
                    out[current_record].append(linedict)
    fh.close()
    return(out)


def filter_overlaps(recs):
    """Filter TRF output records for overlapping features

    Drop features that are completely contained in another feature.  If two
    features overlap exactly (same start and end coordinates), pick the one
    with a higher alignment score. If both have the same alignment score, pick
    the one with more copies (i.e. shorter repeat period).  If two features
    overlap, keep both.

    Parameters
    ----------
    recs : dict
        Output from parse_args()

    Returns
    -------
    dict
        Same structure as input recs but with overlapping and duplicate
        features removed
    """
    out = {}
    for seqid in recs:
        # sort by start and end position
        recs[seqid] = sorted(recs[seqid], key=lambda x: int(x['start']))
        outstack = [] # store all entries
        for rec in recs[seqid]:
            if len(outstack) == 0:
                outstack.append(rec)
            else:
                if int(rec['start']) == int(outstack[-1]['start']):
                    if int(rec['end']) == int(outstack[-1]['end']):
                        # if coordinates overlap entirely
                        if float(rec['aln_score']) > float(outstack[-1]['aln_score']):
                            # pick higher scoring
                            discard = outstack.pop()
                            outstack.append(rec)
                        elif float(rec['aln_score']) == float(outstack[-1]['aln_score']):
                            # if scores identical
                            if float(rec['num_copies']) > float(outstack[-1]['num_copies']):
                                # pick repeat with more copies (i.e. shorter period)
                                discard = outstack.pop()
                                outstack.append(rec)
                    elif int(rec['end']) > int(outstack[-1]['end']):
                        # new record is entirely contained in old record
                        discard = outstack.pop()
                        outstack.append(rec)
                elif int(rec['start']) > int(outstack[-1]['start']):
                    if int(rec['start']) > int(outstack[-1]['end']):
                        # new record starts after old record, no overlap
                        # use > instead of >= because coordinates are end-inclusive
                        outstack.append(rec)
                    elif int(rec['end']) > int(outstack[-1]['end']):
                        # new record overlaps with old record
                        outstack.append(rec)
        out[seqid] = outstack
    return(out)


def report_tsv(recs, filename):
    """Reformat parsed TRF output into TSV format

    Parameters
    ----------
    recs : dict
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
    for seqid in recs:
        for linedict in recs[seqid]:
            fh.write("\t".join([linedict[field] for field in report_fields]) + "\n")
    fh.close()


def report_gff(recs, filename):
    """Reformat parsed TRF output to GFF3 format

    Parameters
    ----------
    recs : dict
        Output from parse_ngs()
    filename : str
        Path to filename to write output
    """
    gff_fields = ["seqid", "source", "type", "start", "end", "score", "strand",
                  "phase", "attributes"]
    attr_fields = ["per_len", "num_copies", "cons_len", "adj_pc_match",
                   "adj_pc_indel", "pc_A", "pc_C", "pc_G", "pc_T", "entropy",
                   "cons_seq"]
    fh = open(filename, "w")
    fh.write("##gff-version-3\n")
    for seqid in recs:
        lines = sorted(recs[seqid], key=lambda x: int(x["start"]))
        for linedict in lines:
            gff_out = {
                "seqid": linedict["seqid"],
                "source": "TRF",
                "type": "region",
                "start": linedict["start"],
                "end": linedict["end"],
                "score": linedict["aln_score"],
                "strand": ".",
                "phase": "."
                }
            gffid = "_".join([
                "TRF", linedict['seqid'], linedict['start'], linedict['end'],
                linedict['aln_score'], linedict['per_len']
                ])
            gff_out["attributes"] = "ID=" + gffid + ";" + ";".join([field+"="+str(linedict[field]) for field in attr_fields])
            fh.write("\t".join([gff_out[field] for field in gff_fields]) + "\n")
    fh.close()


if __name__ == "__main__":
    args = parser.parse_args()

    if args.input:
        recs = None
        if args.format == "ngs":
            recs = parse_ngs(args.input)
        elif args.format == "dat":
            print("not yet implemented")
        # Report output
        if recs:
            if args.no_overlaps:
                recs = filter_overlaps(recs)
            if args.outfmt == "tsv":
                report_tsv(recs, args.out)
            elif args.outfmt == "gff":
                report_gff(recs, args.out)

