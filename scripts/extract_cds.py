from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# empty lists to store sequences
records = []

# indicate input and output files in snakemake file
infile = snakemake.input.fasta
outfile = snakemake.output.index

for record in SeqIO.parse(infile, "fasta"): # parse input fasta
    description = record.description # this is the header for each sequence
    
    if "[protein_id=" not in description: # skip records w/o protein id
        continue

    # split at [protein=, extract text after it, split again at ]
    # only want to keep the protein id
    protein_id = description.split("[protein_id=")[1].split("]")[0]

    # writes new fasta file with just protein id
    records.append(SeqRecord(seq=record.seq, id=protein_id, description=""))

# use SeqIO to write transcriptome sequences to new output file
SeqIO.write(records, outfile, "fasta")





