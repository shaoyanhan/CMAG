#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import concurrent.futures


def calc_n50(nl_list):
    mid = sum(nl_list) * 0.5
    n50 = 0
    nl = 0
    for i in sorted(nl_list)[::-1]:
        nl += i
        if nl >= mid:
            return i


def calc_genome_stats(fasta_f):
    total_bases_num = 0
    ambiguous_bases_num = 0
    contigs_num = 0
    nl_list = []
    n50 = 0
    with SeqIO.parse(fasta_f, "fasta") as fah:
        for record in fah:
            contigs_num += 1
            total_bases_num += len(record.seq)
            nl_list.append(len(record.seq))
            ambiguous_bases_num += record.seq.count("N") + \
                record.seq.count("n")
    return (fasta_f, total_bases_num, ambiguous_bases_num, contigs_num, calc_n50(nl_list))


def drep_score(x, comw, conw, strw, n50w, sizew):
    '''
    compute drep score 
    '''
    score = (
        comw * x["completeness"]
        - conw * x["contamination"]
        + strw * x["contamination"] * (x["strain_heterogeneity"] / 100)
        + n50w * np.log10(x["N50"])
        + sizew * np.log10(x["size"])
    )
    return score


def galah_score(x):
    '''
    compute galah score 
    '''
    score = x["completeness"] - 5 * x["contamination"] - 5 * \
        x["contigs_num"] / 100 - 5 * x["ambiguous_bases_num"] / 100000
    return score


def set_score(df, comw=1, conw=5, strw=1, n50w=1, sizew=1, **kwargs):
    '''
    Args:
      df:  pandas dataframe, genome info,
        include header: 
          genome_path,
          completeness,
          contamination, 
          strain_heterogeneity,
          N50, size, contigs_num, ambiguous_bases_num

    Returns:
      df: dataframe,
          if set output, will save dataframe to output tsv file
    '''
    cancel = False
    for i in ["genome_path", "completeness", "contamination", "strain_heterogeneity", "N50", "size", "contigs_num", "ambiguous_bases_num"]:
        if i not in df.columns:
            print(f"{i} not in dataframe headers, please check genome info file")
            cancel = True
    if cancel:
        print("header error, exiting")
        sys.exit()

    df["drep_score"] = df.apply(lambda x: drep_score(
        x, comw, conw, strw, n50w, sizew), axis=0)
    df["galah_score"] = df.apply(lambda x: galah_score(x), axis=0)
    df = df.sort_values(["drep_score", "galah_score"], ascending=False)

    if ("output" in kwargs) and (kwars["output"] is not None):
        df.to_csv(kwargs["output"], sep='\t', index=False)
    return df


def set_stats(genome_info_file, threads=8):
    '''
    Args:
      genome_info_file: .tsv format, genome info file,
        include header:
          genome_path,
          completeness,
          contamination,
          strain_heterogeneity

    Returns:
      df: dataframe
    '''
    df = pd.read_csv(genome_info_file, sep='\t')

    calc_stats = False
    for col in ["N50", "size", "contigs_num", "ambiguous_bases_num"]:
        if not col in df.columns:
            calc_stats = True

    if calc_stats:
        stats_list = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            for stats_tuple in executor.map(calc_genome_stats, df["genome_path"].to_list()):
                stats_list.append(stats_tuple)
        df_stats = pd.DataFrame.from_records(stats_list, columns=[
                                             "genome_path", "size", "ambiguous_bases_num", "contigs_num", "N50"])
        df = df.merge(df_stats, how="inner")
    return df


def main():
    '''
    a drep/galah wrapper for picking representative genome 

    '''
    parser = argparse.ArgumentParser("pick representative genome")
    parser.add_argument("-comw", dest="comw",
                        default=1, type=float, help="completeness weight, default: 1")
    parser.add_argument("-conw", dest="conw",
                        default=5, type=float, help="contamination weight, default: 5")
    parser.add_argument("-strw", dest="strw",
                        default=1, fype=float, help="strain heterogeneity weight, default: 1")
    parser.add_argument("-n50w", dest="n50w",
                        default=1, fype=float, help="N50 weight, default: 1")
    parser.add_argument("-sizew", dest="sizew",
                        default=1, fype=float, help="genome size weight, default: 1")
    parser.add_argument("-t", dest="threads", default=8, type=int,
                        help="threads, used on calculate genomes stats")
    parser.add_argument("-gi", dest="gi",
                        required=True, help="genome info, tsv format")
    parser.add_argument("-o", dest="output", default=None,
                        help="output, default: None")
    args = parser.parse_args()

    df = set_stats(args.gi, args.threads)
    set_score(df, ags.comw, args.conw, args.strw,
              args.n50w, args.sizew, output=args.output)


if __name__ == "__main__":
    main()
