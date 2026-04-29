def snipdist_runnation(listation, reference, cpus, ram, cwd):
    '''
    Function to run snippy on a list of isolates provided by the user

    Args:
        listation (str): Path to a tab-separated file containing names of isolates in the first column and the path to the isolate FASTA file or FASTQ files in the other column(s), with file names matching the isolateIDs in the metadata file
        reference (str): Path to the reference sequence to align isolates against, must be in genbank or gff format (.gbk or .gff)\
        cpus (int): Number of threads to be used for snippy run
        ram (int): RAM to be used for snippy run, in GB
    
    Returns:
        an output folder in working directory that contains the results of the snippy run
    Raises:


    '''

    # Construct the snippy command
    snippy_cmd = f"snippy-multi {listation} --ref {reference} --cpus {cpus} --ram {ram} --outdir {cwd}/snippy_output --force"

    # Run the snippy command
    try:
        subprocess.run(snippy_cmd, shell=True, check=True)
        print("Snippy run completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running snippy: {e}")