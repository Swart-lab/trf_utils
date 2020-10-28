Scripts to parse TRF output files
=================================

[Tandem Repeats Finder](https://github.com/Benson-Genomics-Lab/TRF) is a tool
for finding tandem repeats in DNA sequences. The default output is HTML but it
also produces reports in tabular format with the `-d` and `-ngs` options.

Here are scripts to parse the tabular output and reformat them into standard
formats (GFF, TSV) for downstream analyses.

Requires `bedtools cluster` to be in path.

For instructions run `python trf_utils.py --help`.
