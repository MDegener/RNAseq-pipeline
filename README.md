# Lightweight preprocessing pipeline for RNAseq data

This lightweight pipeline is based on [nextflow](https://www.nextflow.io/) and was written for the preprocessing of RNAseq data. 

The following tools are implemented in the pipeline:

- Alignment of FASTQ files with [STAR](https://github.com/alexdobin/STAR)

- Sorting and indexing of BAM files with [samtools](http://www.htslib.org/) 

- Quality control with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://multiqc.info/) 

- Counting of exon-level and gene-level reads + calculation of gene TPM with [RNASeQC](https://github.com/getzlab/rnaseqc) 

- Estimation of percent spliced-in (PSI) with [MISO](https://miso.readthedocs.io/en/fastmiso/) (requires sorted/indexed BAM files)

## Getting started

1. Install [nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and [miniconda](https://docs.conda.io/en/latest/miniconda.html) 

2. Set up a conda environment with the required software: 

`conda create -p ./pipeline-env -c bioconda misopy samtools rna-seqc fastqc multiqc star`

3. Download genome annotation files to `RNAseq-pipeline/lib`:
```
# download comprehensive gene annotation from GENCODE v26
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz | gunzip -c > ./lib/gencode.v26.annotation.gtf

# download primary genome assembly from GENCODE v26
wget -O - ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh38.p10.genome.fa.gz | gunzip -c > ./lib/GRCh38.p10.genome.fa
```

4. Create a genome index for STAR, for example:

`STAR --runMode genomeGenerate --genomeDir ./lib/star_index_gencode_v26 --genomeFastaFiles ./lib/GRCh38.p10.genome.fa --sjdbGTFfile ./lib/gencode.v26.annotation.gtf  --sjdbOverhang 75`

5. Create MISO annotation files for hg38 with [rnaseqlib](https://rnaseqlib.readthedocs.io/en/clip/#creating-custom-gff-annotations-for-miso) (Note: may require you to manually install rnaseqlib)
```
# download GENCODE v26 gene predicition table (wgEncodeGencodeCompV26)
wget -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeCompV26.txt.gz | gunzip -c > ./lib/wgEncodeGencodeCompV26.txt

# delete all chromosomes with "_" (consider only reference chromosomes)
awk '$3 !~ /_/' ./lib/wgEncodeGencodeCompV26.txt > ./lib/ensGene.txt

# store gene prediction table in extra directory
mkdir ./lib/predGenesTable | mv ./lib/ensGene.txt ./lib/predGenesTable

# create annotation files based on gene prediction table (note: may require to manually install rnaseqlib to use this function)
python2 ./preprocessing/scripts/gff_make_annotation.py --flanking-rule commonshortest --genome-label hg38 ./lib/genePred/ ./lib

# rename directory
mv ./lib/commonshortest ./lib/miso_annotations_gencode_v26_refChrOnly

# index annotation file for skipped exons (= SE) with MISO
mkdir ./lib/miso_SE_indexed_hg38
index_gff --index ./lib/miso_annotations_gencode_v26_refChrOnly/SE.hg38.gff3 ./lib/miso_SE_indexed_hg38
```

## Running the pipeline

Exemplary command to run the pipeline in your terminal:
```
nextflow RNAseq-pipeline/main.nf -profile noggo -c RNAseq-pipeline/config/nextflow.config -w YOUR_WORK_DIRECTORY --flow="fastqc,star,sort,index,multiqc"` --inFiles=YOUR_INPUT_FILES --outDir=YOUR_OUTPUT_DIRECTORY --readLength=75 --paired=yes
```

By default, the pipeline expects FASTQ files. When you already have sorted and indexed BAM files and would like to use them as input, add the `bam` keyword to the flow parameter. Check the example below on how to correctly set the input directory parameter (i.e. `--inFiles`). 

In `config/nextflow.config`, you can define a profile for your processing environment. You can limit the maximum number of cores and memory to be used by the pipeline. See the [nextflow documentation on config profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles) for more infomation.

## Parameters

| Parameter     | Description                  | Example                                            |
|---------------|------------------------------|----------------------------------------------------|
| flow          |  determines which processes to run, flow keywords can be concatenated | `--flow="bam,miso,rnaseqc"` |
| readLength    |  average read length | `--readLength=75` |
| paired        |  paired-end RNAseq data? | `--paired="yes"` or `--paired="no"` | 
| libDir        |  directory that contains genome annotation files | `--libDir="RNAseq-pipeline/lib"` |
| inFiles       |  directory where input data is stored | FASTQ (single-end): `--inFiles="./*.fastq.gz"`<br> FASTQ (paired-end): `--inFiles="./fastq/*{1,2}.fastq.gz"` <br> BAM: `--inFiles="./bam/*{.bam,.bam.bai}"` |                     
| resDir        |  directory where output data will be stored  | `--resDir="./data"` | 

## Recommended flow definitions

1. To align your FASTQ files and sort/index the resulting BAM files: `--flow="fastqc,star,sort,index,multiqc"`

2. To only sort and index your BAM files: `--flow="bam,sort,index"`

3. To run MISO on sorted/indexed BAM files: `--flow="bam,miso"`

4. To run RNASeQC on BAM files: `--flow="bam,rnaseqc"`

5. Or 3. and 4. combined: `--flow="bam,miso,rnaseqc"`

Note: It is currently not possible to align FASTQ files, sort/index BAM files and run MISO or RNASeQC in one instance of the pipeline. Instead, you should first run the pipeline with `--flow="fastqc,star,sort,index,multiqc"` and subsequently with `--flow="bam,miso,rnaseqc"`.
