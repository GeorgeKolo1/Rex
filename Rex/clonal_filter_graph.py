import numpy as np
import pandas as pd
import random
import networkx as nx


def ClonalFilter(snp_dists, metadata, col, threshold=5):

    '''
    Function to filter out clonal isolates based on user-defined SNP distance threshold and metadata critera.
    
    Uses networkx connected components to create clusters, based on user-defined criteria, which consistent with the SNP address framework.

     Args:
        snp_dists (str): Path to a tab-separated file output by SNP-dists containing pairwise SNP distances between isolates.
        metadata (str): Path to a tab-separated or comma-seperated file containing metadata for the isolates, including columns for isolateID and one clustering criteria at minimum (e.g. host, year, etc...)
        threshold (int, optional): SNP distance threshold for defining clonal pairs. Defaults to 5.

    Returns:
        sample_keep (list): list of isolates to that are non-clonal, including one representative per clonal cluster plus all inherently non-clonal isolates

    Raises:

        '''
    
    random.seed(42)

    clonal_df = pd.read_csv(snp_dists, sep='\t')
    clonal_df = clonal_df.copy()
    clonal_df.columns = ['Isolate1','Isolate2','SNP_distance']
    clonal_unique = clonal_df['Isolate1'].unique()
    clonal_df['SNP_distance'] = clonal_df['SNP_distance'].astype(int)
    clonal_df = clonal_df[clonal_df['SNP_distance'] < threshold]

    if metadata.endswith('.csv'):
        metadata = pd.read_csv(metadata, sep=',')
    else:
        metadata = pd.read_csv(metadata, sep='\t')

    metadata = metadata[metadata[col].isin(clonal_unique)]
    meta_cols = [c for c in metadata.columns if c != col]
    print(f'Length of unique isolates in snp_dists file {len(clonal_unique)}, length of metadata after filtering {len(metadata)}')

    all_snp_isos = set(metadata[col]) - {'Reference'}
    isos_with_close_neighbour = (set(clonal_df['Isolate1'])
                                .union(set(clonal_df['Isolate2']))) - {'Reference'}
    truly_isolated = all_snp_isos - isos_with_close_neighbour

    print(f'Number of isolates (from metadata): {len(all_snp_isos)}')
    print(f'Isolates with at least one <{threshold} SNP neighbour: {len(isos_with_close_neighbour)}')
    print(f'Isolates with NO close neighbour anywhere: {len(truly_isolated)}')

    meta = metadata.set_index(col)
    clonal_df = clonal_df.merge(meta.add_suffix('_1'), on = 'Isolate1')
    clonal_df = clonal_df.merge(meta.add_suffix('_2'), on = 'Isolate2')

    #merge host metadata with clonal_df to get host information for both isolate1 and isolate2
    clonal_df = clonal_df.merge(host_lookup, left_on='Isolate1', right_on=col, how='left')
    clonal_df = clonal_df.rename(columns={'host' : 'host_1'}).drop(columns=col)
    
    clonal_df = clonal_df.merge(host_lookup, left_on='Isolate2', right_on=col, how='left')
    clonal_df = clonal_df.rename(columns={'host' : 'host_2'}).drop(columns=col)

    clonal_df = clonal_df.merge(year_lookup, left_on='Isolate1', right_on=col, how='left')
    clonal_df = clonal_df.rename(columns={'Year' : 'year_1'}).drop(columns=col)

    clonal_df = clonal_df.merge(year_lookup, left_on='Isolate2', right_on=col, how='left')
    clonal_df = clonal_df.rename(columns={'Year' : 'year_2'}).drop(columns=col)

    clonal_df = clonal_df.merge(DT_lookup, left_on='Isolate1', right_on=col, how='left')
    clonal_df = clonal_df.rename(columns={'PhageType' : 'dt_1'}).drop(columns=col)

    clonal_df = clonal_df.merge(DT_lookup, left_on='Isolate2', right_on=col, how='left')
    clonal_df = clonal_df.rename(columns={'PhageType' : 'dt_2'}).drop(columns=col)

    #Ensure clonal_df contains only isolates where host, year, and phage type match between the pair
    filtered = (
    (clonal_df['host_1'] == clonal_df['host_2']) &
    (clonal_df['year_1'] == clonal_df['year_2']) &
    (clonal_df['dt_1'] == clonal_df['dt_2']) &
    (clonal_df['Isolate1'] != clonal_df['Isolate2']) &
    (clonal_df['Isolate1'] != 'Reference') &
    (clonal_df['Isolate2'] != 'Reference') 
    )

    clonal_df = clonal_df[filtered]
    print(f'Number of clonal pairs after all filtering: {len(clonal_df)}')

    # Single linkage clustering via connected components
    sample_keep = []
    for (host, year, phagetype), group in clonal_df.groupby(['host_1', 'year_1', 'dt_1']):
        G = nx.Graph()
        for _, row in group.iterrows():
            G.add_edge(row['Isolate1'], row['Isolate2'])

        # Sort connected components for deterministic processing order
        for component in sorted(nx.connected_components(G), key=lambda c: sorted(c)):
            # One representative per connected component
            nodes = sorted(component)
            sample_keep.append(random.sample(nodes, 1)[0])

    sample_keep_deduped = list(set(sample_keep))
    clonal_isolates = set(clonal_df['Isolate1']).union(clonal_df['Isolate2'])
    singletons = all_snp_isos - clonal_isolates

    print(f'Representatives from clusters (raw): {len(sample_keep)}')
    print(f'Representatives from clusters (deduped): {len(sample_keep_deduped)}')
    print(f'Clonal isolates: {len(clonal_isolates)}')
    print(f'Singletons: {len(singletons)}')
    print(f'Overlap (reps also in singletons): {len(set(sample_keep_deduped) & singletons)}')

    sample_keep = list(set(sample_keep_deduped) | singletons)
    print(f'Final output: {len(sample_keep)}')

    return sample_keep
# EXAMPLE USAGE: sample_keep = ClonalFilter(clonal_df, df_meta2)