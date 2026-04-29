from Rex.clonal_filter_graph import ClonalFilter

def test_ClonalFilter():
    snp_dists = 'test/test_data/test_snp_dists.tsv'
    metadata = 'test/test_data/test_metadata.tsv'
    col = 'IsolateID'
    threshold = 5

    expected_output = ['IsolateA', 'IsolateC', 'IsolateE']
    output = ClonalFilter(snp_dists, metadata, col, threshold)

    assert set(output) == set(expected_output), f"Expected {expected_output}, but got {output}"