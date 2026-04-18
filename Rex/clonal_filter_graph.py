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
        col (int): The name of the column that contains isolateIDs, isolate names, or isolate references in the metadata file.

    Returns:
        sample_keep (list): list of isolates to that are non-clonal, including one representative per clonal cluster plus all inherently non-clonal isolates

    Raises:

        '''
    
    # Set random seed for reproducibility of representative selection from clusters
    random.seed(42)

    #Load and define snp-dist dataframe, filter to only include pairs below or equal to the user-defined threshold 
    clonal_df = pd.read_csv(snp_dists, sep='\t')
    clonal_df.columns = ['Isolate1','Isolate2','SNP_distance']
    clonal_unique = clonal_df['Isolate1'].unique()
    clonal_df['SNP_distance'] = clonal_df['SNP_distance'].astype(int)
    clonal_df = clonal_df[clonal_df['SNP_distance'] <= threshold]

    # Simple logic to parse file into dataframe depending on type of seperator
    if metadata.endswith('.csv'):
        metadata = pd.read_csv(metadata, sep=',')
    else:
        metadata = pd.read_csv(metadata, sep='\t')

    # Ensures that the returned isolates don't contain those that are in metadata but not in the snp-dists file, included this from experience..........
    metadata = metadata[metadata[col].isin(clonal_unique)]
    meta_cols = [c for c in metadata.columns if c != col]
    print(f'Length of unique isolates in snp_dists file {len(clonal_unique)}, length of metadata after filtering {len(metadata)}')

    all_snp_isos = set(metadata[col])
    isos_with_close_neighbour = (set(clonal_df['Isolate1'])
                                .union(set(clonal_df['Isolate2']))) - {'Reference'}
    truly_isolated = all_snp_isos - isos_with_close_neighbour

    print(f'Number of isolates (from metadata): {len(all_snp_isos)}')
    print(f'Isolates with at least one <{threshold} SNP neighbour: {len(isos_with_close_neighbour)}')
    print(f'Isolates with NO close neighbour anywhere: {len(truly_isolated)}')

    # Ensures the dynamic processing of metadata columns, joining on isolate 1 and isolate 2 to allow for filtering based on user-defined criteria (e.g. host, year, etc...)
    meta = metadata.set_index(col)
    clonal_df = clonal_df.join(meta.add_suffix('_1'), on = 'Isolate1')
    clonal_df = clonal_df.join(meta.add_suffix('_2'), on = 'Isolate2')

    #Mask to filter clonal_df such that only user-defined criteria are included in clonal_df
    mask = pd.Series(True, index=clonal_df.index)
    for base in meta_cols:
        mask &= (clonal_df[f'{base}_1'] == clonal_df[f'{base}_2'])

    mask &= (clonal_df['Isolate1'] != clonal_df['Isolate2'])
    mask &= (clonal_df['Isolate1'] != 'Reference')
    mask &= (clonal_df['Isolate2'] != 'Reference')

    clonal_df = clonal_df[mask]
    print(f'Number of clonal pairs after all filtering: {len(clonal_df)}')

    # Single linkage clustering via connected components using networkx
    sample_keep = []
    for _, group in clonal_df.groupby([f'{b}_1' for b in meta_cols]):
        G = nx.Graph()
        for _, row in group.iterrows():
            G.add_edge(row['Isolate1'], row['Isolate2'])

        # Sort connected components for deterministic processing order, also learned from experience that without this different isolates are returned each time...
        for component in sorted(nx.connected_components(G), key=lambda c: sorted(c)):
            # Random sample is taken from each "cluster" or connected component
            nodes = sorted(component)
            sample_keep.append(random.sample(nodes, 1)[0])

    sample_keep_deduped = list(set(sample_keep))
    clonal_isolates = set(clonal_df['Isolate1']).union(clonal_df['Isolate2'])
    singletons = all_snp_isos - clonal_isolates

    print(f'Representatives from clusters (raw): {len(sample_keep)}')
    print(f'Representatives from clusters (deduped): {len(sample_keep_deduped)}')
    print(f'Clonal isolates: {len(clonal_isolates)}')
    print(f'Singletons: {len(singletons)}')

    # There should be no overlap between sample_keep and singletons, as singletons should never be within clonal clusters, if they are something is broken....
    overlap = set(sample_keep_deduped) & singletons
    if overlap:
        raise RuntimeError(
            f"Something in core logical is broken... there are non-clonal isolates within the clonal clusters. Check the input files and the filtering criteria. "
            f"{len(overlap)} isolates appear as both representative and singleton: "
            f"{sorted(overlap)}"
    )

    sample_keep = list(set(sample_keep_deduped) | singletons)
    print(f'Final output: {len(sample_keep)}')

    return sample_keep