import argparse
from Rex import clonal_filter_graph as cfg, snippy_run as sr, snpdist_run as sdr
import pandas as pd


def get_args():
    '''
    
    '''
    parser = argparse.ArgumentParser(description='Rex is a command-line tool that returns a list of non-clonal isolates based on user-defined criteria')

    cfilter = parser.add_argument_group('Clonal Filtering')
    cfilter.add_argument('--reads', '-r', type=str, required=True, help='Path to the FASTQ files containing the raw reads for each isolate, with file names matching the isolateIDs in the metadata file')
    cfilter.add_argument('--assemblies', '-a', type=str, required=True, help='Path to the FASTA files containing the raw reads for each isolate, with file names matching the isolateIDs in the metadata file')
    cfilter.add_argument('--metadata', '-m', type=str, required=True, help='Path to the input file (CSV or TSV) containing metadata')
    cfilter.add_argument('--column', '-c', type=str, required=True, help='Name of the column in the metadata file corresponding to the isolateIDs, or isolate names')
    cfilter.add_argument('--output', '-o', type=str, required=True, help='Path to the output file where results will be saved')

    snippy = parser.add_argument_group('SNIPPY run to generate core genome alignment')
    snippy.add_argument('--fasta', '-fa', type=str, help='Path to the file (CSV or TSV) containing phenotype data for association analysis')
    snippy.add_argument_group.add_argument('--subtypes', '-s', help='path to file (CSV or TSV) containing subtypes to perform phenotype association analysis for, can be a single column (one subtyping method) or multiple columns (multiple subtyping methods)')

    snp_dist = parser.add_argument_group('Calculation of pairwise SNP distances using SNP-dists')
    snp_dist.add_argument('--snpdist_input', '-sd', type=str, help='Path to the input file (CSV or TSV) containing pairwise SNP distances for SNP distance analysis')


    return parser.parse_args()

def main():
    args = get_args()

    df = pd.read_csv(args.input, sep='\t' if args.input.endswith('.tsv') else ',')

    DI, DI_low, DI_high = di.DiscriminationIndex(df)

    if args.comparator is not None:
        if args.comparator not in df.columns:
            raise ValueError(f"Comparator column '{args.comparator}' not found in input file")
        comparator_col = df[args.comparator]

        AW_ab, AW_ba = aw.AdjustedWallace(df.iloc[:, 0], comparator_col)
    else:
        AW_ab, AW_ba = aw.AdjustedWallace(df.iloc[: 0], df.iloc[:, 1])


    results = pd.DataFrame({
        'DI': [DI],
        'CI_low': [DI_low],
        'CI_high': [DI_high],
        'AWC_AvsB' : [AW_ab],
        'AWC_BvsA' : [AW_ba]
    })

    if args.output.endswith('.tsv'):
        results.to_csv(args.output, sep='\t', index=False)
    else:
        results.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()