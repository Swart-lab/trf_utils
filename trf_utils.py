#!/usr/bin/env python3

import re
import json
import argparse
import tempfile
import logging
from subprocess import run
from collections import defaultdict

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
    "--min_length", type=int, default=100,
    help="Minimum length of a repeat region to report")
parser.add_argument(
    "--no_overlaps", action="store_true",
    help="""
    Report GFF without overlaps: if features overlap, pick highest scoring
    feature. In case of tied score, pick feature with the most repeat copies.
    Requires `bedtools cluster` in path.
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

logging.basicConfig(level=logging.DEBUG)

NGSFIELDS = ["start", "end", "per_len", "num_copies", "cons_len",
             "adj_pc_match", "adj_pc_indel", "aln_score",
             "pc_A", "pc_C", "pc_G", "pc_T", "entropy", "cons_seq",
             "seq_region", "flank_left", "flank_right"]


def check_env():
    out = {}
    bedtools_vers = run(["bedtools", "--version"], capture_output=True)
    out['bedtools'] = bedtools_vers.stdout.decode().rstrip()
    return(out)


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


def report_gff_no_overlaps(recs, filename):
    """Pick best-scoring repeat region in case of overlaps

    Use bedtools cluster to generate overlap clusters
    For each cluster, pick repeat feature with highest alignment score
    If scores are tied, pick feature with the most repeat num_copies
    """
    # cluster overlapping features with bedtools cluster
    tf = tempfile.NamedTemporaryFile()
    report_gff(recs, tf.name)
    clusts = run(['bedtools','cluster','-i',tf.name],capture_output=True)
    tf.close()
    clusts = clusts.stdout.decode().split("\n")  # break on lines
    clusts = [line.split("\t") for line in clusts if line != ""]  # transform to list of lists
    # group by cluster number
    cd = defaultdict(list)
    for rec in clusts:
        cd[rec[9]].append(rec)
    # get highest scoring aln per cluster
    # TODO edge case where we have bookended features that are assigned to same cluster
    out = {}
    for clustid in cd:
        if len(cd[clustid]) == 1:
            logging.debug(
                f"cluster {str(clustid)} has only one member")
            out[clustid] = cd[clustid][0]
        else:
            topscore = [rec for rec in cd[clustid] if int(rec[5]) == max([int(r[5]) for r in cd[clustid]])]
            if len(topscore) == 1:
                logging.debug(
                    f"cluster {str(clustid)} has {str(len(cd[clustid]))} members, single top score")
                out[clustid] = topscore[0]
            else:
                ncmatches = []
                for rec in topscore:
                    ncmatch = re.search(r"num_copies=([^;]+);", rec[8])
                    ncmatch = float(ncmatch.group(1))
                    ncmatches.append((ncmatch, rec))
                # report record with highest num_copies
                logging.debug(
                    f"cluster {str(clustid)} has {str(len(cd[clustid]))} members, tied top score, picking highest num_copies")
                out[clustid] = sorted(ncmatches, key=lambda x: float(x[0]), reverse=True)[0][1]
    # write to file
    fh = open(filename, "w")
    fh.write("##gff-version-3\n")
    for clustid in out:
        fh.write("\t".join(out[clustid][0:9]) + "\n")  # strip cluster ID: last item in list
    fh.close()


if __name__ == "__main__":
    args = parser.parse_args()

    depvers = check_env()
    logging.debug(
        "Dependencies: " + " ".join([f"{tool}:{depvers[tool]}" for tool in depvers]))

    if args.input:
        recs = None
        if args.format == "ngs":
            recs = parse_ngs(args.input)
        elif args.format == "dat":
            print("not yet implemented")
        # Report output
        if recs:
            if args.outfmt == "tsv":
                report_tsv(recs, args.out)
            elif args.outfmt == "gff":
                if args.no_overlaps:
                    report_gff_no_overlaps(recs, args.out)
                else:
                    report_gff(recs, args.out)

