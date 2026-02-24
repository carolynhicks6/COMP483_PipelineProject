from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# define file paths
infile=snakemake.input.fq
outfile=snakemake.output.longest_contig

# empty list to store contigs
contigs=[]

# parse through fasta and append sequences to contigs list
for record in SeqIO.parse(infile, "fasta"):
    contigs.append(record.seq)

# open outfile and write the longest contig to the file
with open (outfile, "w") as out:
    out.write(str(max(contigs, key=len)))

