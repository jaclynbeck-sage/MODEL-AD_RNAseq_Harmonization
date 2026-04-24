# Running the pipeline

This document is a developer resource with step-by-step instructions for running the RNASeq harmonization pipeline end-to-end.

> **Note:** The instructions contain Sage-specific internal links to Confluence and JIRA, and assumes access to the Seqera project and appropriate spaces on Synapse.
>
> -   **If you work at Sage** and do not have access to everything, these instructions should walk you through how to get access.
>
> -   **If you do not work at Sage** and want to reproduce the pipeline from the beginning, you will have to have your own Nextflow setup and possibly your own Synapse space to upload/download files if you are using the Nextflow steps.
>
>     -   If you only want to run the last few steps of the pipeline (R), you need a Synapse account to download any data in the [Model-AD RNAseq Harmonization Project](https://www.synapse.org/Synapse:syn25882591), and no special access is required.

**This pipeline trades off between R code and Seqera Platform several times.**

The general steps are:

1.  **R:** (1) Setup =\> (2) Metadata harmonization =\> (3) rnaseq samplesheet creation
2.  **Seqera:** (4) nf-core/rnaseq pipeline to align fastq files
3.  **R:** (5) sarek samplesheet creation
4.  **Seqera:** (6) nf-core/sarek pipeline to infer genotypes of each sample
5.  **R:** (7) Add provenance to counts files
6.  **R:** (8) Sample validation =\> (9) differential expression analysis

------------------------------------------------------------------------

## Set up Seqera Platform credentials

<https://sagebionetworks.jira.com/wiki/x/CICggg> is a great resource for getting set up with Seqera Platform. See that page for more details, but here is a brief summary of the steps you need to take to work on this pipeline:

1.  Check JumpCloud and make sure you have an item called `strides-ampad-workflows-towerviewer` in “All applications”. If you don’t, [create a JIRA ticket for IT](https://sagebionetworks.jira.com/jira/secure/CreateIssue.jspa?issuetype=3&pid=10083) to request membership.
2.  Our project is called `model-ad-rna-project`. If you do not have access to this project in Seqera Platform, [follow these instructions](https://sagebionetworks.jira.com/wiki/spaces/WF/pages/2191556616/Getting+Started+with+Nextflow+and+Seqera+Platform#Adding-Users-to-Existing-Seqera-Platform-Projects) to get yourself added, or ask Jaclyn to do it.
3.  Follow [these instructions](https://sagebionetworks.jira.com/wiki/spaces/WF/pages/2191556616/Getting+Started+with+Nextflow+and+Seqera+Platform#Deploying-to-Seqera-Platform-via-the-Web-UI) to set up credentials set up for Seqera Platform. You will need to:
    -   Set up a token/secret on Synapse that can be used for Seqera
    -   Create a profile in your computer's local AWS config file by following the instructions in [this Slack comment](https://sagebionetworks.slack.com/archives/C8SJHFCKT/p1661389316173499?thread_ts=1661387709.096529&cid=C8SJHFCKT).
    -   Get a SAML response token (covered in the above instructions)

------------------------------------------------------------------------

## Adding a new study to the pipeline

The pipeline depends on a config file that lists which studies to use. To add a new study, create a new row in `Model_AD_SynID_list.csv` with the study name and the Synapse IDs of the RNA seq assay metadata, biospecimen metadata, and individual metadata files. Make sure to specify a version for each Synapse ID.

The study name must match the value of “study” in the annotations of any of the metadata files, since it is used to query Synapse for files related to the study.

### Updating the intervals file

The intervals.bed file is used for genotype validation, and has the coordinates of where in a given gene the mutation/variant is expected to be. If the study has a new genotype that isn’t in the existing studies, let Jaclyn know so she can add coordinates for that variant to the intervals file. See also [Editing the intervals bed file](https://sagebionetworks.jira.com/wiki/spaces/MA/pages/4312563765/Sarek+Editing+the+intervals+bed+file?atlOrigin=eyJpIjoiZTMzZDU3ZGY3NjgwNDYxMjkyMjkzNTJhZmM2NjY4ZGMiLCJwIjoiYyJ9){.uri} for more detail on how this is done and what mutations are already covered in the current file.

------------------------------------------------------------------------

# Pipeline steps 1-3 (R)

The following steps are all in R. Each step will only use the studies listed in `Model_AD_SynID_list.csv`.

## Step 1 (Initial Synapse setup)

This step sets up the expected folder structure on Synapse for each study listed in `Model_AD_SynID_list.csv`, which makes it easier to go back and forth between Seqera and R results.

This step only needs to be run **once when you add a new data set** to the pipeline. All existing data sets are already set up.

Running the script extra times will not do anything to existing data on Synapse, so it is harmless to run it more than once.

> **File:** 01_SYN_CreateFolderStructure.R
>
> **Uploads:** None
>
> **Edits required if adding new study:** None

## Step 2 (Metadata harmonization)

This step pulls in all metadata files from all studies listed in `Model_AD_SynID_list.csv`, fixes some values in the older files, and merges all the data together. This script also prints out a warning if any metadata files have newer versions on Synapse than what is in the Synapse ID list.

**When adding a new study,** you must check the following columns for consistency:

-   `sex` should be “male” or “female”, lower-case

-   `ageDeath` should be numeric and a have reasonable numbers for age in months

-   `genotype` should only have values that exist in [rnaseq_genotype_label_map.csv](https://www.synapse.org/Synapse:syn69014401). This includes exactly matching capitalization and punctuation. Check with Jess Britton for the correct values to use for genotype if the new study hasn’t been added to that file yet.

-   All `specimenIDs` in the RNA seq assay metadata should exist in the biospecimen data, and all `individualIDs` associated with those `specimenIDs` should exist in the individual metadata.

The older Model-AD studies uploaded metadata before the template was standardized, and inconsistencies in those studies have already been accounted for in this script. Hopefully, newly-added studies should be more consistent and free of issues.

> **File:** 02_HarmonizeMetadata.R
>
> **Uploads:** One metadata per study to [Synapse](https://www.synapse.org/Synapse:syn61850200)
>
> **Edits required if adding a new study:**
>
> -   Add the study name to the script’s header comments.
>
> -   Examine the new metadata files for issues. Any fixes (fixing typos, renaming genotypes, etc) that aren’t already accounted for in the script should be added to the study-specific fixes section in the same style as the other studies with fixes, or to the genotype re-mapping section.
>
> -   **Comment thoroughly** any changes that are made, and why.
>
> -   If no fixes are required, add the study name to the comment at the end of the study-specific fixes section listing studies that need no corrections. This makes it clear that the study has been confirmed as problem-free rather than left unimplemented.

## Step 3 (nf-core/rnaseq sample sheet creation)

This script queries the AD Knowledge portal for all the fastq files related to each study, and creates sample sheets for the [nf-core/rnaseq](https://nf-co.re/rnaseq/3.14.0/) pipeline.

**When adding a new study,** you need to check the results of the `synTableQuery` for the new study and make sure that the number of fastq files is 2x the number of specimens, and that each file is annotated with a `specimenID` that exists in the metadata.

> **File:** 03_NF_CreateSampleSheets.R
>
> **Uploads:** Created sample sheets are uploaded to the [Sample Sheets](https://www.synapse.org/Synapse:syn62147112) folder on Synapse.
>
> **Edits required if adding new study:** If there are any issues with fastq file annotations or the query is not returning the right files, add any study-specific fixes in the same section as the other study fixes. Comment thoroughly any changes that are made, and why.

------------------------------------------------------------------------

# Pipeline Step 4 (Seqera Platform)

This is the step that aligns the fastq files to the reference genome, and is run on [Seqera Platform](https://tower.sagebionetworks.org/) under project `model-ad-rna-project`. This step uses the [nf-core/rnaseq](https://nf-co.re/rnaseq/3.14.0/) pipeline, version 3.14.0.

## Upload the sample sheet(s) to the AWS bucket

You will need to get a SAML response token, as described in [Set up tower credentials](https://sagebionetworks.jira.com/wiki/spaces/MA/pages/edit-v2/4295819317#Set-up-Tower-credentials). As a reminder: Go to JumpCloud, click the `strides-ampad-workflows-towerviewer` tile, and use the SAML response browser plug-in to get the token.

In a console, authenticate with:

```         
aws-saml --profile ampad
```

Paste your SAML response token when prompted. Upload the sample sheet:

```         
aws --profile ampad s3 cp <path_to_local_samplesheet> s3://model-ad-rna-project-tower-bucket/sample_sheets/
```

Check that it is there:

```         
aws --profile ampad s3 ls s3://model-ad-rna-project-tower-bucket/sample_sheets/
```

## Upload fastq files to AWS bucket

Launch the `synapse_stage` pipeline in Seqera Platform, with arguments:

-   `input`: path to the uploaded sample sheet on AWS.

-   `output`: desired output folder on AWS, which will be auto-created if it doesn't exist. The output folder should be in `s3://model-ad-rna-project-tower-scratch/work` so the fastqs get deleted after 30 days.

-   `entry`: synstage

## Run the RNASeq alignment pipeline

Launch the `rna_seq` pipeline in Seqera Tower, with arguments:

-   `input`: path to the *new* sample sheet created by `synapse_stage`. This sample sheet will be located in `<output folder from synapse_stage>/synstage_<name of synstage run on Seqera>`

-   `outdir`: path on AWS where the output of the pipeline should go. The data should go in `s3://model-ad-rna-project-tower-bucket/rnaseq_output/<study name>`.

-   You should not need to change any other settings. The pipeline has been set up to use the correct genome files and settings for Model-AD data.

Name the pipeline run `rnaseq_<study name>` or something else informative so we can quickly find previous runs for specific data sets. **Do not** use the auto-generated name.

The `rna_seq` pipeline is set up to use a spot instance, so it may fail and need to be resumed if the spot instance stops.

### Tag the successful run

When the run finally succeeds all the way, **tag that run** with “Harmonization” and “Final Data” so we know which one to refer back to. If the pipeline failed and needed to be resumed, tag only the final run that succeeded.

### Save the configuration of the run

After the `rna_seq` pipeline finishes, go to the run's page in Seqera Platform, click on the "Parameters" tab, and download the parameters as a JSON file. Name the file `rnaseq_<study_name>-params.json` and upload it to the [Configuration](https://www.synapse.org/Synapse:syn62147114) folder of the harmonization project on Synapse.

## Upload results to Synapse

You will need to upload 3 different sets of data: RSEM matrices, BAM files, and QC data. You will launch the `synapse_index` pipeline in Seqera Platform, once for each type of data with these arguments:

**RSEM matrices**

-   `filename_string`: rsem.merged

-   `parent_id`: the Synapse ID of the study's folder inside the [raw counts folder](https://www.synapse.org/Synapse:syn51132850). If a folder for the study doesn't already exist, make sure the study is in `Model_AD_SynID_list.csv` and run Step 1 to make the folder.

-   `s3_prefix`: \<outdir from rnaseq run\>/star_rsem/

-   `entry`: synindex

**BAM files**

-   `filename_string`: markdup.sorted.bam

-   `parent_id`: the Synapse ID of the study's folder inside the [BAM files folder](https://www.synapse.org/Synapse:syn63856101). If a folder for the study doesn't already exist, make sure the study is in `Model_AD_SynID_list.csv` and run Step 1 to make the folder.

-   `s3_prefix`: \<outdir from rna_seq\>/star_rsem/

-   `entry`: synindex

**QC files**

-   `parent_id`: the Synapse ID of the quality control folder for the study, which should be located inside the study's folder in the [Quality Control](https://www.synapse.org/Synapse:syn74490531) folder. If a folder doesn't exist for the study, make one with the study's name.

-   `s3_prefix`: \<outdir from rna_seq\>/multiqc/star_rsem/

-   `entry`: synindex

-   Delete the `filename_string` argument if it exists

------------------------------------------------------------------------

# Pipeline Step 5 (R)

### nf-core/sarek sample sheet creation

This step makes sample sheets of BAM files, for the [nf-core/sarek](https://nf-co.re/sarek/3.4.4/) pipeline to estimate genotypes.

> **File:** 05_NF_Create_BAM_Sample_Sheets.R
>
> **Uploads:** Created sample sheets are uploaded to the [Sample Sheets](https://www.synapse.org/Synapse:syn62147112) folder on Synapse.
>
> **Edits required if adding new study:** None

------------------------------------------------------------------------

# Pipeline Step 6 (Seqera Platform)

This is the step that estimates the genotype of each sample based on its aligned data, and is run on [Seqera Platform](https://tower.sagebionetworks.org/) under project `model-ad-rna-project`. This step uses the [nf-core/sarek](https://nf-co.re/sarek/3.4.4/) pipeline, version 3.4.4.

## Upload files to the AWS bucket

Upload the variant calling sample sheet(s) to AWS as in Step 4.

Then, launch the `synapse_stage` pipeline on Seqera Platform as in Step 4, pointing to the variant calling sample sheet.

Note: This uploads the BAM files from Synapse to the scratch bucket. This is a little inefficient, since the BAMs already exist in the main bucket, but it seemed easier to automate creation of the sample sheet and pipeline launch by staging them like this instead of trying to list all the BAM file paths via AWS queries.

## Run the Genotype calling pipeline

Launch the `genotype_calling` pipeline in Seqera Tower, with arguments:

-   `input`: path to the *new* sample sheet created by `synapse_stage`. This sample sheet will be located in `<output folder from synapse_stage>/synstage_<name of synstage run on Seqera>`

-   `outdir`: path on AWS where the output of the pipeline should go. The data should go in `s3://model-ad-rna-project-tower-bucket/genotype_calling/<study name>`.

-   `intervals`: s3://model-ad-rna-project-tower-bucket/reference_genomes/intervals_universal_genome.bed

    -   This is the default file for the pipeline and it already exists at that location. If you have uploaded a different intervals file elsewhere for testing, point to your file instead.

-   You should not need to change any other settings. The pipeline has been set up to use the correct settings for Model-AD data.

Name the pipeline run `genotype_<study name>` or something else informative. **Do not** use the auto-generated name.

The `genotype_calling` pipeline is set up to use a spot instance, so it may fail and need to be resumed if the spot instance stops.

As in Step 4, tag the final successful run with “Harmonization” and “Final Data”, and save the configuration of the run to the [Configuration](https://www.synapse.org/Synapse:syn62147114) folder. Name the configuration file `genotype_<study_name>-params.json`.

## Upload results to Synapse

Launch the `synapse_index` pipeline on Seqera Platform with arguments:

-   `filename_string`: bcftools.vcf.gz

-   `parent_id`: the Synapse ID of the study’s folder in the [Genotype Validation](https://www.synapse.org/Synapse:syn63913842) folder. If a folder for the study doesn't already exist, make sure the study is in `Model_AD_SynID_list.csv` and run Step 1 to make the folder.

-   `s3_prefix`: \<outdir from genotype_calling\>/variant_calling/bcftools/

-   `entry`: synindex

------------------------------------------------------------------------

# Pipeline Step 7 (R)

This step does some cleanup on Synapse files and merge all of the separate counts files into one big file. It also adds each study's name to its individual RSEM counts files, so each set of files has distinct names.

First, this step sets the provenance on the newly-created counts files to point to all the files that were used to generate it. Counts files are renamed in the process. Next, step gathers all the counts files for all the studies in `Model_AD_SynID_list.csv` and concatenates them into a single matrix. It does this for all 4 types of RSEM output.

**File:** 07_SYN_UpdateProvenance.R

**Uploads:** 4 files:

-   [Model-AD_all_studies.gene_counts.tsv](https://www.synapse.org/Synapse:syn62690577)

-   [Model-AD_all_studies.gene_tpm.tsv](https://www.synapse.org/Synapse:syn62690622)

-   [Model-AD_all_studies.transcript_counts.tsv](https://www.synapse.org/Synapse:syn62690714)

-   [Model-AD_all_studies.transcript_tpm.tsv](https://www.synapse.org/Synapse:syn62690807)

**Edits required if adding new study:** None

------------------------------------------------------------------------

# Pipeline Step 8 (R)

Documentation in progress.

------------------------------------------------------------------------

# Pipeline Step 9 (R)

Documentation in progress.
