# define paths of input and output files
input_files=snakemake.input.tsv
output_file=snakemake.output.result

# make string of desired header indicating tab delimination
header="sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n "

# open the output file to write to it
with open(output_file, "w") as outfile:

    # # since snakemake passes full file path to python, extract sample name from path using string methods
    for path in input_files:
        sample_name=path.split('/')[-1].split('_')[0]

        # for each input path write the sample name and header to output file
        outfile.write(f"{sample_name}:\n")
        outfile.write(header)

        # open the input files
        with open(path, "r") as infile:
            # iterate through input and write the top 5 results to the output file
            for i, line in enumerate(infile):
                if i < 5:
                    outfile.write(line)
    
    # newline for formatting
    outfile.write("\n")
