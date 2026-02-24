samples=['SRR5660030', 'SRR5660033', 'SRR5660044', 'SRR5660045']
accession="GCF_000845245.1"
species="HCMV"
subfamily="Betaherpesvirinae"

# sort all expected output files
rule all:
    input:
        expand("paired_fastq/{sample}_1.fastq", sample=samples),
        expand("ref/{accession}_genomic.fna", accession=accession), 
        expand("ref/{accession}_cds.fna", accession=accession),
        expand("data/kallisto/{sample}/abundance.tsv", sample=samples),
        "results/sleuth/sleuth_results.tsv",
        expand("data/bowtie2/sam_files/{sample}.sam", sample=samples),
        expand("data/bowtie2/extracted_fastq/{sample}_R1.fastq",sample=samples),
        expand("data/bowtie2/extracted_fastq/{sample}_R2.fastq", sample=samples),
        "results/bowtie2/bowtie2_results.txt",
        expand("data/spades/{sample}/contigs.fasta", sample=samples),
        expand("data/spades/{sample}/scaffolds.fasta", sample=samples),
        expand("data/spades/{sample}_longest_contig", sample=samples),
        expand("results/blast/{sample}_blast.tsv", sample=samples),
        "results/blast/blast_results.tsv",
        "results/PipelineReport.txt"

# this rule downloads the genome and cds from NCBI for the indicated accession number
# this was accomplished by downloading and unzipping in a temporary directory
# then copying the unzipped files into the final directory 
# finally, the temporary directory was deleted
rule download_genome_and_cds:
    output:
        genome="ref/{accession}_genomic.fna",
        cds="ref/{accession}_cds.fna"
    shell:
        """
        mkdir -p temp/
        mkdir -p ref/

        datasets download genome accession {wildcards.accession} \
            --include genome,cds \
            --filename temp/ncbi_dataset.zip

        unzip -o temp/ncbi_dataset.zip -d temp
        
        find temp -name "*_genomic.fna" -exec cp {{}} {output.genome} \;
        find temp -name "cds_from_genomic.fna" -exec cp {{}} {output.cds} \;

        rm -rf temp/
        """

# the purpose of this rule is to build a transcriptome index for HCMV
# runs biopython and SeqIO to write index with just protein ids as header
# final result will be used as index for kallisto
rule make_index:
    input:
        fasta=f"ref/{accession}_cds.fna"
    output:
        index= "data/kallisto/index.txt"
    script:
        "scripts/extract_cds.py"

# the purpose of this rule is to make the kallisto index
# this allows kallisto to map reads to transcripts
rule kallisto_index:
    input:
        fasta="data/kallisto/index.txt"
    output:
        index="data/kallisto/index.idx"
    shell:
        "kallisto index -i {output.index} {input.fasta}"

# the purpose of this rule is to use kallisto to quantify transcripts per million for each sample
# the output is a directory for each sample which contains abundance files
rule quantify_tpm:
    input:
        index="data/kallisto/index.idx",
        fq1="paired_fastq/{sample}_1.fastq",
        fq2="paired_fastq/{sample}_2.fastq"
    output:
        out_dir=directory("data/kallisto/{sample}"),
        abundance="data/kallisto/{sample}/abundance.tsv"
    threads: 2
    shell:
        """
        kallisto quant -i {input.index} \
        -o {output.out_dir} -b 10 -t {threads} \
        {input.fq1} {input.fq2}
        """

# sleuth requires kallisto scripts to run - that's why it's input
# this tool is designed to compare the differential expression of RNA-seq data
# find statistical differences between the conditions and report to results file
rule run_sleuth:
    input:
        abundance = expand("data/kallisto/{sample}/abundance.tsv", sample=samples)
    output:
        "results/sleuth/sleuth_results.tsv"
    script:
        "scripts/sleuth.R"

# in order to run bowtie2, there needs to be an index of the reference genome
# this rule builds the index of 6 files
rule build_bowtie_index:
    input:
        genome=f"ref/{accession}_genomic.fna"
    output:
        "data/bowtie2/index/{species}.1.bt2",
        "data/bowtie2/index/{species}.2.bt2",
        "data/bowtie2/index/{species}.3.bt2",
        "data/bowtie2/index/{species}.4.bt2",
        "data/bowtie2/index/{species}.rev.1.bt2",
        "data/bowtie2/index/{species}.rev.2.bt2"
    shell:
        """
        mkdir -p data/bowtie2/index
        bowtie2-build {input.genome} data/bowtie2/index/{species}
        """

# this rule runs bowtie2
# need to wait until index is built before running bowtie2
# bowtie2 input is paired end fasta files
# output is sam files for each sample
# use bt2 prefix as parameter to take into account that there are 6 files for the bt2 index
rule run_bowtie2:
    input:
        idx = multiext(f"data/bowtie2/index/{species}", 
                       ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", 
                       ".rev.1.bt2", ".rev.2.bt2"),   
        fq1 = "paired_fastq/{sample}_1.fastq",
        fq2 = "paired_fastq/{sample}_2.fastq"
    output:
        "data/bowtie2/sam_files/{sample}.sam"
    params:
        bt2_prefix = f"data/bowtie2/index/{species}"
    shell:
        """
        mkdir -p data/bowtie2/sam_files
        bowtie2 -x {params.bt2_prefix} \
                -1 {input.fq1} \
                -2 {input.fq2} \
                -S {output}
        """

# in order to analyze results of bowtie2, the sam files need to be filtered and turned into paired end fastq files
# reads that didn't map to the genome are removed
rule filter_sam_extract_fastq:
    input:
        "data/bowtie2/sam_files/{sample}.sam"
    output:
        fq1="data/bowtie2/extracted_fastq/{sample}_R1.fastq",
        fq2="data/bowtie2/extracted_fastq/{sample}_R2.fastq"
    shell:
        """
        samtools view -b -F 4 {input} | \
        samtools sort -n | \
        samtools fastq -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -
        """ 

# given the filtered fastq files, the next step is to compare read pairs pre vs. post filtering
# this rule runs a python script to compare read pair numbers and write results to a file that summarizes all samples
rule bowtie2_results:
    input:
        initial=expand("paired_fastq/{sample}_1.fastq", sample=samples),
        filtered=expand("data/bowtie2/extracted_fastq/{sample}_R1.fastq", sample=samples)
    output:
        results="results/bowtie2/bowtie2_results.txt"
    script:
        "scripts/bowtie2_results.py"

# the next step is to run spades using bowtie2 output files
# this rule produces contigs and scaffolds fasta files in corresponding directory for each SRR
rule run_spades:
    input:
        fq1="data/bowtie2/extracted_fastq/{sample}_R1.fastq",
        fq2="data/bowtie2/extracted_fastq/{sample}_R2.fastq"
    output:
        contigs="data/spades/{sample}/contigs.fasta",
        scaffolds="data/spades/{sample}/scaffolds.fasta"
    params:
        outdir="data/spades/{sample}"
    shell:
        """
        mkdir -p data/spades
        spades.py -k 127 -t 2 --only-assembler \
        -1 {input.fq1} -2 {input.fq2} -o {params.outdir}
        """

# to find out which strains each assemly aligns to, the first thing to do is find the longest contig from spades assembly
# this rule does that using a python script, creating longest contig files for each SRR
rule spades_longest_contig:
    input:
        fq="data/spades/{sample}/contigs.fasta"
    output:
        longest_contig="data/spades/{sample}_longest_contig"
    script:
        "scripts/longest_contig.py"

# this rule uses shell commands to make a local ncbi database
# only include refseq genomes
# output is db directory with 7 different files
rule make_NCBI_database:
    output:
        db_files = multiext(f"db/{subfamily}", 
                       ".ndb", ".nhr", ".nin", ".not", 
                       ".nsq", ".ntf", ".nto")  
    shell:
        """
        mkdir -p db
        datasets download virus genome taxon {subfamily} --refseq --include genome --filename db/db.zip
        unzip -o db/db.zip -d db
        makeblastdb -in db/ncbi_dataset/data/*.fna -dbtype nucl -out db/{subfamily}
        """

# with the local database made, it is possible to run blast for spades output against this db
# the results are put into a tsv file 
rule blast:
    input:
        query="data/spades/{sample}_longest_contig",
        db_files = multiext(f"db/{subfamily}", 
                       ".ndb", ".nhr", ".nin", ".not", 
                       ".nsq", ".ntf", ".nto")
    output:
        results="results/blast/{sample}_blast.tsv"
    shell:
        """
        blastn -task megablast -query {input.query} -db db/{subfamily} \
        -out {output.results} -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle"
        """

# this rule takes all the tsv files generated from blast and organizes them into one tsv file
# this was accomplished using a python script that writes to desired output format
rule write_blast_results:
    input:
        tsv=expand("results/blast/{sample}_blast.tsv", sample=samples)
    output:
        result="results/blast/blast_results.tsv"
    script:
        "scripts/blast_results.py"

# this snakemake file wrote results for each of the tools used to separate files
# in order to tie everything together, this final rule takes the result files and appends them to the pipeline report
rule write_pipeline_report:
    input:
        cds="data/kallisto/index.txt",
        sleuth="results/sleuth/sleuth_results.tsv",
        bowtie2="results/bowtie2/bowtie2_results.txt",
        blast="results/blast/blast_results.tsv"
    output:
        report="results/PipelineReport.txt"
    shell:
        """
        echo "The HCMV genome ({accession}) has $(grep -c '^>' {input.cds}) CDS." >> {output.report}
        echo "" >> {output.report}
        paste <(cut -f1 {input.sleuth}) <(cut -f4 {input.sleuth}) <(cut -f2 {input.sleuth}) <(cut -f3 {input.sleuth})  >> {output.report}
        echo "" >> {output.report}
        cat {input.bowtie2} >> {output.report}
        echo "" >> {output.report}
        cat {input.blast} >> {output.report}
        """
