import argparse
from Rex import clonal_filter_graph as cfg, snippy_run as sr, snpdist_run as sdr
import pandas as pd
import os

def get_args():
    '''
    
    '''
    parser = argparse.ArgumentParser(description='Rex is a command-line tool that returns a list of non-clonal isolates based on user-defined criteria')

    cfilter = parser.add_argument_group('Clonal Filtering')
    
    cfilter.add_argument('--c-filter', '-sf', type=str, action='store_true', help='Run only the clonal filtering step, must provide molten snp-dist file with header from SNP-dists', default=False)
    cfilter.add_argument('--molten_input', '-mi', type=str, help='Path to the input file (CSV or TSV) containing a molten snp-dists file, with header, for clonal filtering step')
    cfilter.add_argument('--input', '-i', type=str, required=True, help='Path to the .tab or .tsv file containing names of isolats in the first column and the path to the isolate FASTA file or FASTQ files in the other column(s), with file names matching the isolateIDs in the metadata file')
    cfilter.add_argument('--metadata', '-m', type=str, required=True, help='Path to the file (CSV or TSV) containing metadata')
    cfilter.add_argument('--column', '-c', type=str, required=True, help='Name of the column in the metadata file corresponding to the isolateIDs, or isolate names')
    cfilter.add_argument('--output', '-o', type=str, required=True, help='Path to the output folder where results will be saved')

    snippy = parser.add_argument_group('SNIPPY run to generate core genome alignment')
    snippy.add_argument('--ram', '-r', type=int, help='RAM to be used for snippy run, in GB', default=8)
    snippy.add_argument('--threads', '-t', type=int, help='Number of CPUs to be used for snippy run', default=1)
    snippy.add_argument('--reference', '-R', help='Reference sequence to align isolates against, must be in genbank or gff format (.gbk or .gff)')

    snp_dist = parser.add_argument_group('Calculation of pairwise SNP distances using SNP-dists')
    snp_dist.add_argument('--snp_dists', '--sd', action='store_true', type=str, help='Argument for only running snp-dist and clonal filtering steps, must provide a core.aln file from snippy core-genome alignment', default=False)
    snp_dist.add_argument('--snpdist_input', '-sd', type=str, help='Path to the input file (CSV or TSV) containing core.aln file from snippy')
    

    return parser.parse_args()

def main():
    args = get_args()

    root = os.getcwd()

    if args.c-filter:
        cfg.ClonalFilterGraph(args.molten_input, args.metadata, args.column, args.output)

    if args.snp_dists:
        sdr.snpdist_runnation(args.snpdist_input, args.output)
        cfg.clonal_filter_graph(args.snpdist_input, args.metadata, args.column, args.output)

    else:
        sr.snippy_runnation(args.input, args.reference, args.threads, args.ram, root)
        sdr.snpdist_runnation(root)
        cfg.clonal_filter_graph(root, args.snpdist_input, args.metadata, args.column)
    
    

if __name__ == "__main__":
    main()