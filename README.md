# Rex
Rex is a command-line tool that identifies and clusters all isolates based on user-defined clonal criteria, returning a list of non-clonal isolates.

## Workflow
  1. Rex takes as input FASTA or FASTQ files, a metadata file, and a user set SNP threshold.
  2. Rex runs SNIPPY to produce a core-genome alignment
  3. Rex runs SNP-dists on the core.aln file from SNIPPY in order to calculate pairwise SNP distances
  4. Rex calculates clusters all isolates with SNP distances less than the user specified threshold and metadata, using networkX
  5. Am isolate from each user-defined clonal cluster is then sampled to act as a representative for each clonal cluster
  6. A final txt.file of non-clonal isolates is returned including the representatives isolates from non-clonal clusters as well as isolates that are inherently non-clonal.

## Quickstart

## Input


## Output

## Citation

## Why Rex?
The name originates from the fact that the CLI tool identifies user-defined clonal clusters. Star wars... Clones... Rex.... 

If you know, you know.
