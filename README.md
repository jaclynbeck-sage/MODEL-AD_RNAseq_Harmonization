# MODEL-AD RNAseq Harmonization Study

This repository harmonizes metadata and RNAseq data across multiple MODEL-AD studies, and performs differential expression analysis for each study.

General documentation about the project can be found on the [project's Synapse page](https://www.synapse.org/Synapse:syn25882591). Briefly, this pipeline works as follows:

1.  Metadata from each of the studies is harmonized so that key columns have the same names and information. Values are altered if necessary to be consistent across all studies.

2.  Raw fastq files from each study are aligned with the `nf-core/rnaseq` pipeline (STAR -\> RSEM) to get gene and transcript counts.

3.  The genotype of each sample is validated, where possible, with the `nf-core/sarek` pipeline (mpileup), using detection or non-detection of variants of altered genes. When a gene modification cannot be detected with RNAseq data, this step is ignored for mouse genes, and positive expression of human genes is used instead when the model contains human transgenes.

    1.  Detailed information on the development of the intervals.bed file used to define variant coordinates can be found here:

4.  The sex of each mouse is validated using expression of X- and Y-chromosome genes *Xist*, *Eif2s3y*, and *Ddx3y*.

5.  Differential expression is performed on the counts files, using only samples that passed QC.

Detailed documentation on running the pipeline end-to-end can be found here:
