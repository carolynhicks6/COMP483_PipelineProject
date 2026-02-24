from Bio import SeqIO

# define paths of files
initial_files = snakemake.input.initial
filtered_files = snakemake.input.filtered
outfile = snakemake.output.results

with open(outfile, 'w') as out:
    # iterate through both lists at the same time
    for init_path, filt_path in zip(initial_files, filtered_files):

        # since snakemake passes full file path to python, extract sample name from path using string methods
        sample_name = init_path.split('/')[-1].split('_')[0]

        # counter to store numbers of paired reads
        initial_counter=0
        filtered_counter=0

        # count initial numbers of paired reads
        with open(init_path, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                initial_counter += 1
        
        # count filtered numbers of paired reads
        with open(filt_path, 'r') as handle:
            for record in SeqIO.parse(handle, "fastq"):
                filtered_counter += 1 
        
        # write to output file
        out.write(f"Sample {sample_name} had {initial_counter} read pairs before and {filtered_counter} read pairs after Bowtie2 filtering.\n")



        

