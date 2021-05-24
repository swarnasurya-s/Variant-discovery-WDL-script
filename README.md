# Variant Discovery Pipeline Using WDL-script
The variant discovery pipeline is used to find variants which includes SNPs and indels (insertion/deletions) from NGS (next-generation sequencing) reads. The pipeline includes the following:

Adapter trimming: Removal of contaminants such as adapters are done using Trimmomatic - Trimmomatic 0.39 http://www.usadellab.org/cms/?page=trimmomatic

Alignment: The sequence reads are aligned to the reference genome using BWA (Burrows-Wheeler Aligner) – bwa 0.7.17 http://bio-bwa.sourceforge.net/

BAM Processing: A series of intermediate steps to process and prepare the BAM file for variant calling. Picard Tools is used for BAM processing – Picard Tools 1.119 http://broadinstitute.github.io/picard/

Variant calling: Variants are called using the GATK Unified Genotyper – GATK 4.1.6.0 https://www.broadinstitute.org/gatk/

FastQC (v.0.11.5): A simple tool that provides an overview of raw data sequence - https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ Note in FastQC: if the program failed to export results, the module has an error. Therefore, you need to download the tool.
