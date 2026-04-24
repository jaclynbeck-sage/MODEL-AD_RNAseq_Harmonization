# Editing the intervals.bed file

## Table of Contents

1.  [Introduction]

2.  [Examples]

    1.  [Human gene mutation (London mutation of APP)](#human-gene-mutation-london-mutation-of-app)

    2.  [Mouse gene mutations]

3.  [Current mutations in the file]

    1.  [Human]

    2.  [Mouse]

4.  [Excluded mutations]

## Introduction

This pipeline does genotype verification where possible, using the `nf-core/sarek` workflow. The workflow requires a .bed file of intervals, where each interval is the genomic coordinates of the mutation to look for. Finding the coordinates for a mutation can be complicated and time consuming, so this document has a few examples done a few different ways to explain my process. **I show both a more complicated/involved way and an easier way for each example,** because not all mutations have the same information available.

We use a custom genome which was created by this repository: [Sage-Bionetworks/customReferenceMODEL-AD](https://github.com/Sage-Bionetworks/customReferenceMODEL-AD). All coordinates must be relative to the custom genome. For mouse genes, this is just the real genomic coordinates from the Ensembl v112 mouse genome. For human genes, there is some math involved:

-   The custom genome includes the sequences for the human genes APOE, APP, CLU, MAPT, and PSEN1. Each gene was extracted from the human genome and concatenated onto the mouse genome as its own chromosome.

-   Chromosome mapping:

    -   APOE =\> 20

    -   APP =\> 21

    -   MAPT =\> 22

    -   PSEN1 =\> 23

    -   CLU =\> 24

-   To get genomic coordinates for human genes, subtract the real coordinates of the start of the gene on the **human genome** from the coordinates of the mutation, and add 1.

    -   For example, APOE starts at coordinate 44905791 of chromosome 19 on the human genome. A mutation at coordinate 44905800 would be at 44905800 - 44905791 + 1 = 10 on chromosome 20 in our custom genome.

**The general process of finding the right coordinates is:**

1.  Go to the Jax page for the mouse model

2.  Try to find information about which specific amino acids or bases were changed

3.  Get a sequence to search for, either by looking at the Genotyping Protocols section or by figuring out where the amino acid / codon falls within the protein from UniProt/CCDS.

    1.  Alternatively, several mouse models have “[report cards](https://www.synapse.org/Synapse:syn25474060)” detailing exactly which part of the sequence was edited.

4.  Search for the sequence on the Ensembl page for the gene, looking at the exons of certain transcripts. **It’s important to use the [Ensembl archive for v112](https://may2024.archive.ensembl.org/index.html),** which corresponds to our custom genome.

    1.  **Pay attention to whether the gene is on the forward or reverse strand.** If it’s on the reverse strand, you may have to reverse-complement the sequence or count backward instead of forward depending on where you find your information.

5.  Get the coordinates of the target base(s) using the information provided for the matching exon. **The intervals file must be tab-separated, so do not edit it with software that replaces tabs with spaces.**

6.  After running the Sarek pipeline with your new intervals, look at the output **by hand** for a few known carriers and a few known WT mice and verify that the mutation is being recognized appropriately, shows the correct base change, and is at the right coordinates.

7.  Once you are sure that the intervals you put in are correct, update the [file on Synapse](https://www.synapse.org/Synapse:syn69046595) and **upload it to the Seqera tower bucket** under `s3://model-ad-rna-project-tower-bucket/reference_genomes/`.

> **Note:** Coordinates in the intervals file start at 0, while the coordinates for Ensembl and our genome start at 1. So all Ensembl coordinates need to have 1 subtracted from them for the intervals file.
>
> -   In the APOE example above, the start coordinate would be 10-1 = 9.
>
> The end coordinate in the intervals file is *not included* in the interval, so it must always be at least 1 base beyond the desired coordinates.
>
> Some part of the Sarek pipeline does not work if the start and end coordinates are only 1 base apart, so all 1-base mutations have an interval that is **2 bases wide**. I also add an extra base at the end of the other intervals just to be safe. This does not affect the outcome of the genotyping.
>
> -   In the APOE example above, the end coordinate would be 9+2 = 11.

------------------------------------------------------------------------

# Examples

## Human gene mutation (London mutation of APP) {#human-gene-mutation-london-mutation-of-app}

The “London” mutation in APP is one of the mutations in 5xFAD mice. The stock number for these mice (008730) is listed in the 5xFAD individual metadata (both UCI and Jax’s metadata), so we go to [jax.org](http://jax.org), search for that number, and end up at <https://www.jax.org/strain/008730>.

We are looking specifically for any information on how the mutations were introduced. This information might be under “Detailed Information”, “Development”, or the “Genetics” header, depending on the mouse model.

### Using the London mutation as an example

The London mutation is listed as V717I, which means that amino acid 717 of the human APP gene was changed from V → I. If you Google for “amino acid sequence table” and “amino acid abbreviations” you’ll find many tables that convert between amino acids and base triplets.

Looking at the codons for Val and Ile, we need to change the first G in the codon to an A to switch from V → I.

We have to convert this to genomic coordinates in the human genome (Ensembl release 112) and then to the coordinates of our custom genome (Ensembl release 112).

[The Ensembl release 112 archive page for APP](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000142192;r=21:25880535-26171128) shows that this gene is on the reverse strand, but all of our coordinates are relative to the forward strand, so there are some extra reversals that need to happen.

Possible places to find coordinates for this mutation:

-   The transcript table on the [gene page in Ensembl](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000142192;r=21:25880535-26171128) lists several transcripts. Find the transcript that seems to be the “main” one, in this case “APP-201”, and look at the links to CCDS and/or Uniprot.

-   The “Sequence” tab in the left-hand menu on the [gene page in Ensembl](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Sequence?db=core;g=ENSG00000142192;r=21:25880535-26171128) will give the full sequence of the gene with exons highlighted.

-   If you click on a transcript, the [“Exons” tab](https://may2024.archive.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000142192;r=21:25880535-26171128;t=ENST00000346798) under “Sequence” breaks the transcript into exons with genomic coordinates.

-   If the mutation is a well-known one (which this one is), you can just Google for the coordinates of the SNP. For example, [SNPedia](https://www.snpedia.com/index.php/Rs63750264) has a page for this mutation that provides coordinates. However, you need to make sure that the coordinates match what is in Ensembl version 112.

-   Look at the “Genotyping Protocols” section on the Jax page for the mouse model. There are typically links to the sequence used in Sanger sequencing or other methods, which gives you the sequence around the mutation. This section is less likely to have easily-readable information for human genes but is useful for mouse genes.

### Method 1 (the longer way)

In this case we will use a combination of the [UniProt page](https://www.uniprot.org/uniprotkb/P05067/entry#P05067-1) and the [CCDS page](https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS13576) for APP-201. On the UniProt page, there is a numbered amino acid sequence. If you look at amino acid #717, it is a V, so we are looking at the right variant. Copy the set of amino acids that are around 717 (`VIATVIVITL`) and search for that sequence on the CCDS page. (Alternatively, you could count to 717 along the amino acid sequence on the CCDS page). If you click on the V in the interactive amino acid sequence, it also highlights the codon in the sequence above it. So this particular V is codon “GTC” (GUC in RNA).

The sequence on this page alternates between black and blue text to indicate exons. Our V is in the second-to-last exon. The table above lists the real genomic coordinates of each exon, however those coordinates are relative to the forward strand but the sequence is on the reverse strand. This means that the second-to-last exon in the sequence is actually the second exon in the table above, **and** the start of the exon is at the end of the sequence. Counting backward from the end, the V is the 21st amino acid. The G in the codon is the only base that changes, so the location of the G in the sequence is 20\*3 + 2 = 62. Then, the real coordinate of this G is 25891722 + 62 = **25891784**. (Start coordinate of second exon + 62 bases).

Now we need to convert this to **our** genomic coordinates. When I made the custom genome, I appended each human gene as its own chromosome to the mouse genome, and each human gene starts at coordinate 1. We go back to the [Ensembl archive page](https://may2024.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000142192;r=21:25880535-26171128;t=ENST00000346798) for APP and find that it lists the genomic coordinates of the gene as 25,880,535-26,171,128. Make sure you are looking at the **gene** coordinates and not the transcript coordinates. These coordinates are relative to the forward strand, so we use 25,880,535 as our start. The coordinate of the G relative to our genome is 25891784 - 25880535 = **11249**. This number is already 0-based for the intervals file since we did not add 1 to it.

### Method 2 (shorter option with less counting)

Since this is a known SNP, we can take a shortcut. Find the codon on the CCDS page as above. Copy the sequence starting a few bases before the codon and ending a few bases after. I used `ATCGTCATC`. Go back to the Ensembl archive page for transcript APP-201 and click on the [“Exons” tab](https://may2024.archive.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000142192;r=21:25880535-26171128;t=ENST00000346798). Search for the sequence on the page. With all the markup it is a little hard to find the match, but you can see it in the second-to-last exon. Click directly on the G, and it will tell you that its genomic coordinate is **25891784]**. Then subtract 25880535 and as above.

### Add to the intervals file

Now that we have our coordinate of **11249**, we can add it to the `intervals_universal_genome.bed` file as a new line (tab-separated):

```         
21  11249 11251 APP-LonV717I
```

-   The “21” is the number of the chromosome **in our custom reference genome**, not the human chromosome (although in the case of APP, they are the same).

-   11249 is our start coordinate (0-based).

-   11251 is the end coordinate. The pipeline will look at the bases from [start, end), which means the start coordinate is included but the end coordinate is not.

-   “APP-LonV717I” is the label for the mutation. The label can be anything, but we have a lot of mutations in the file and more down the road, so I try to keep it as clear as possible which exact mutation corresponds to this line.

**Notes:**

-   Even though we only care about one base, the pipeline throws errors if the start and end coordinates are only 1 base apart (e.g. 11249 to 11250), so all of the mutations that are 1 base only actually have intervals that are **2 bases wide**.

-   These coordinates are 0-based coordinates. Since we subtracted the G coordinate from the start-of-gene coordinate, the final coordinate is already 0-based. However, if you were to use a given coordinate from the mouse genome, which doesn’t require subtracting anything, you would need to subtract 1 from the coordinate for the intervals file.

### Extra verification

Due to all of the coordinate systems, I do a final verification using R and the custom reference genome to make sure the base I want is at the location I think it is. Here is some basic R code to do this:

```         
BiocManager::install("rtracklayer")
library(synapser)

gtf_file <- synGet("syn62035250", downloadLocation = <convenient_filepath>)
fasta_file <- synGet("syn62035247", downloadLocation = <convenient_filepath>)

gtf <- rtracklayer::import(gtf_file$path, format = "gtf", feature.type = "gene")
       
subset(gtf, gene_name == "APP") # Will give you which chromosome the gene is on too
      
fasta <- rtracklayer::import(fasta_file$path, format = "fasta", type = "DNA")
    
app <- fasta[grepl("^21 dna:chromosome", names(fasta))][[1]]
    
# 0-based start coord is 11249, so 1-based coord is 11250. 
print(app[11250]) # Should be a "C" (complement of G on the reverse strand)
      
# Let's look at the whole codon, which should be from 11248-11250.
# The codon on the reverse strand, so we use reverseComplement().
print(reverseComplement(app[11248:11250]))
```

This prints out `GTC`, which matches what we expect. You can also print out a few more bases before and after to make sure the sequence matches what is in Ensembl.

------------------------------------------------------------------------

## Mouse gene mutations

Mouse genes are easier to add than human genes in some ways, mostly because we don’t have to convert between coordinate systems. Let’s use the **Abca7-V1613M** mutation as an example. No stock number was provided in the metadata of the UCI_ABCA7 study, but the [wiki page for the study](https://www.synapse.org/Synapse:syn27207345) has it. We end up at <https://www.jax.org/strain/035316> .

The [Ensembl v112 archive page for Abca7](https://may2024.archive.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000035722;r=10:79832328-79851406) says the gene is on the forward strand, so we don’t have to reverse complement anything.

### Method 1 (the longer way)

Find the transcript that seems to be the main one (*Abca7-203*), go to the UniProt page, and find the 1613th amino acid in the sequence. That AA is a Y, so **we are looking at the wrong transcript**. Instead we find that the correct transcript is *Abca7-202* (or *Abca7-201*, they are almost the same), and the 1613th AA on that UniProt page is a V. Copy the few AAs around the V (e.g. `IVVFI`) and search for it on the CCDS page for the *Abca7-202* transcript. Click on the correct V and find the highlighted codon. Copy the sequence a few bases around that codon (e.g. `GTGGTGTTC`).

Go back to the Ensembl page for *Abca7* and click on the *Abca7-202* transcript. Go to the “Exons” tab under “Sequence” in the left pane and search for the sequence. We find it in the middle of exon 35 but there’s no link to click on in the specific sequence to get the coordinates. Instead we count: starting from 79,845,985 (the start coordinate of the exon), we find the codon at 79846015 - 79846017. There are 2 G’s in our codon, so we have to look at an amino acid chart to figure out which G is being changed. Methionine (M) corresponds to the start codon AUG. So to change V → M, we change the first G of GUG to AUG. Our G of interest is at coordinate 79846015.

We can double-check we’re in the right spot by searching for the guide RNA provided by Jax (CTTGGTGGCAGTGTGCATAG). In some cases, the guide sequence doesn’t exactly match the real sequence so you have to search for parts of it, but in this case it does match, and is in exon 35 right before the codon we want.

### Method 2 (shorter option 1)

On the Jax page of the mouse model, scroll down to the “Genotyping Protocols” section and click on the link to Endpoint analysis. That page gives you the sequence around the mutation, in this case both the target mutation and two silent SNPs. Copy a short part of the sequence just after the (g/a) and search for it in the exon sequence of a transcript on the Ensembl archive page. In this case, any transcript will work. The sequence just appears in different exons depending on the transcript. Count where the target G is in the exon and add it to the start coordinate of the exon as above.

This method also gives us the location of two silent mutations, 1 and 4 bases before the G. We can include those silent mutations in the interval.

### Method 3 (shorter option 2)

This mouse model has a [development report card](https://www.synapse.org/Synapse:syn25474077) that gives the exact sequence that was edited. Search for a short part of the sequence that includes the mutation(s) on the Ensembl archive page as above, and find the coordinates.

### Add a line in the intervals file

```         
10  79846010  79846016  Abca7-V1613M
```

-   10 is the chromosome that *Abca7* is on

-   Since we are including the silent mutations, we start 4 bases before the G. The actual start coordinate would then be 79846011, which is 79846010 in the 0-based .bed format.

-   The 0-based coordinate of the target G is 79846015-1 = 79846014. Since the end coordinate in the file is excluded from the set of bases, we would add 1 to get 79846015 for the end coordinate. However, I added an extra base to be safe, given the issues I had with the pipeline for single-base mutations as mentioned above.

### Extra verification

As with the human gene above, we can verify that the coordinates we are using correspond to the bases we expect. Using the `fasta` object from the R code above:

```         
abca7 <- fasta[grepl("^10 dna:chromosome", names(fasta))][[1]]
print(abca7[79846011:79846015])
```

This should print out `AGTGG`, which covers our silent mutations and the target G.

------------------------------------------------------------------------

# Current mutations in the file

## Human

------------------------------------------------------------------------

### APP-LonV717I

The “London” familial AD mutation in human APP. This point mutation changes amino acid 717 of the APP-201 transcript from V → I (**G**TC → **A**TC).

**Location:** chr 21: 11250 (our genome), chr 21: 25891784 (Human genome)

**Genotype:** 5xFAD_carrier

**Jax:** <https://www.jax.org/strain/008730>

------------------------------------------------------------------------

### APP-FloI716V

The “Florida” familial AD mutation in human APP. Changes amino acid 716 of the APP-201 transcript from I → V (**A**TC → **G**TC).

**Location:** chr 21: 11253 (our genome), chr 21: 25891787 (Human genome)

**Genotype:** 5xFAD_carrier

**Jax:** <https://www.jax.org/strain/008730>

------------------------------------------------------------------------

### APP-SweM671L / APP-SweK670N

The “Swedish” familial AD mutation in human APP. This consists of two point mutations in adjacent amino acids of the APP-201 transcript: AA 671 changes from M → L (**A**TG → **C**TG) and AA 670 changes from K → N (AA**G** → AA**T**).

**Location:** chr 21: 17092 - 17093 (our genome), chr 21: 25897626 - 25897627 (Human genome)

**Genotype:** 5xFAD_carrier, 3xTg-AD_carrier

**Jax:** <https://www.jax.org/strain/008730>, <https://www.jax.org/strain/004807>

------------------------------------------------------------------------

### PSEN1-M146L

A familial AD mutation in human PSEN1. Changes AA 146 of the PSEN-201 transcript from M → L (**A**TG → **C**TG).

**Location:** chr 23: 37246 (our genome), chr14: 73173663 (Human genome)

**Genotype:** 5xFAD_carrier

**Jax:** <https://www.jax.org/strain/008730>

------------------------------------------------------------------------

### PSEN1-L286V

A familial AD mutation in human PSEN1. Changes AA 286 of the PSEN-201 transcript from L → V (**C**TC → **G**TC).

**Location:** chr 23: 61700 (our genome), chr14: 73198117 (Human genome)

**Genotype:** 5xFAD_carrier

**Jax:** <https://www.jax.org/strain/008730>

------------------------------------------------------------------------

## Mouse

------------------------------------------------------------------------

### Abca7-V1613M

A LOAD mutation in the mouse *Abca7* gene. Changes AA 1613 of the *Abca7-201* transcript from V → M (**G**TG → **A**TG). Introduction of the mutation also produced two silent mutations upstream. Note that the genotype below is labeled with the corresponding human mutation relative to the human gene (V1599M) rather than the mouse gene (V1613M).

**Location:** chr 10: 79846011 (silent), 79846014 (silent), 79846015 (target)

**Sequence in interval:** **A**GT**GG** → **C**GT**CA**

**Genotype:** Abca7-V1599M_homozygous

**Jax:** <https://www.jax.org/strain/035316>

**Model development:** <https://www.synapse.org/Synapse:syn25474077>

------------------------------------------------------------------------

### Abi3-S212F

A LOAD mutation in the mouse *Abi3* gene. Changes AA 212 of the *Abi3-201* transcript from S → F (T**CT** → T**TC**) on the reverse strand. It also introduces a silent mutation upstream. Note that the genotype below is labeled with the corresponding human mutation relative to the human gene (S209F) rather than the mouse gene (S212F).

**Location:** chr11: 95724846-95724847 (target), 95724849 (silent)

**Sequence in interval:** **AG**A**G** → **GA**A**T** (**C**T**CT** → AT**TC** on reverse strand)

**Genotype:** Abi3-S209F_homozygous

**Jax:** <https://www.jax.org/strain/035871>

**Model development:** <https://www.synapse.org/Synapse:syn25474072>

> This variant may not be reliably detectable, based on testing with the UCI Primary Screen. More testing will be needed if a full study using this genotype becomes available.

------------------------------------------------------------------------

### App-KI

A humanization of mouse *App* in which 3 point mutations are introduced into Exon 14 on the reverse strand. AA 676 of transcript *App-206* is changed from G → R (**G**GA → **C**GA), AA 681 is changed from F → Y (T**T**T → T**A**T), and AA 684 is changed from R → H (C**GC** → C**AT**).

**Location:** chr16:84762596, 84762580, 84762570 - 84762571

**Sequence at interval:** **GC**GGACTTCA**A**ATCCTGAATCATGTC**C** → **AT**GGACTTCA**T**ATCCTGAATCATGTC**G** (**G**GACATGATTCAGGAT**T**TGAAGTCC**GC** → **C**GACATGATTCAGGAT**A**TGAAGTCC**AT** on reverse strand)

**Genotypes:** hAbeta-KI_LoxP_homozygous, LOAD2

**Jax:** <https://www.jax.org/strain/030898>, <https://www.jax.org/strain/030670>

------------------------------------------------------------------------

### Bin1-K358R

A LOAD mutation in the mouse *Bin1* gene. Changes AA 537 of the *Bin1-201* transcript from K → R (A**AA** → A**GG**). It also introduces a silent mutation (C → G) one base upstream. Note that K358R is the notation for the human notation, not the mouse notation, which has several different names.

**Location:** chr18:32565427 (silent), 32565429-32565430 (target)

**Sequence in interval:** **C**A**AA** → **G**A**GG**

**Genotype:** Bin1-K358R_homozygous

**Jax:** <https://www.jax.org/strain/035872>

**Model development:** <https://www.synapse.org/Synapse:syn25474068>

> This interval has not been tested yet.

------------------------------------------------------------------------

### Il1rap-ex3-KO

This is a knockout of exon 3 of the mouse *Il1rap* gene. For detection we use the start coordinate of the forward [genotyping probe](https://www.jax.org/Protocol?stockNumber=032777&protocolID=41260) sequence through the end of exon 3. The coordinates on the page are relative to a different version of the mouse genome, but the correct coordinates for our genome are 26495665 - 26495744.

**Location:** chr16: 26495665 - 26495744

**Genotype:** LOAD2.IL1rap-Exon3KO_homozygous

**Jax:** <https://www.jax.org/strain/032777>

------------------------------------------------------------------------

### Il34\*Y213

A LOAD mutation in the mouse *Il34* gene. Changes AA 213 of the *Il34-201* transcript from Y → stop codon (T**AC** → T**GA**).

**Location:** chr8: 111468980 - 111468981

**Sequence at interval:** GT → TC (AC → GA on reverse strand)

**Genotype:** LOAD2.Il34-Y213_homozygous

**Jax:** <https://www.jax.org/strain/033541>

------------------------------------------------------------------------

### Picalm-H465R

A LOAD mutation in the mouse *Picalm* gene. Changes AA 465 of the *Picalm-206* transcript from H → R (CAT → AGA). It also introduces a silent mutation downstream, which is far enough away that we did not include it in the interval. Note that the genotype below is labeled with the corresponding human mutation relative to the human gene (H458R) rather than the mouse gene (H465R).

**Location:** chr7: 89831561 - 89831563

**Sequence in interval:** CAT → AGA

**Genotype:** Picalm-H458R_homozygous

**Jax:** <https://www.jax.org/strain/034037>

**Model development:** <https://www.synapse.org/Synapse:syn25474067>

> This variant may not be reliably detectable, based on testing with the UCI Primary Screen. More testing will be needed if a full study using this genotype becomes available.

------------------------------------------------------------------------

### Psen1-M146V

A familial AD mutation in the mouse *Psen1* gene. Changes AA 145 and 146 from IM → VV (**A**T**TA**TG → **G**T**GG**TG).

**Location:** chr12: 83761632, 83761634-83761635

**Sequence in interval:** **A**T**TA** → **G**T**GG**

**Genotype:** 3xTg-AD_carrier

**Jax:** <https://www.jax.org/strain/004807>

------------------------------------------------------------------------

### Ptprb\*D57N

A LOAD mutation in the mouse *Ptprb* gene. Changes AA 57 of the *Ptprb-203* transcript from D → N (**G**AT → **A**AT).

**Location:** chr10: 116113190

**Genotype:** LOAD2.Ptprb-D57N_KI/KI

**Jax:** <https://www.jax.org/strain/033867>

> This variant has not been tested yet.

------------------------------------------------------------------------

### Trem2-R47H

A LOAD mutation in the mouse *Trem2* gene that changes AA 47 from R → H (C**G**C → C**A**C). It also introduces 2 silent mutations downstream, which we leave out of the intervals file.

**Location:** chr17: 48655584

**Genotypes:** Trem2-R47H_homozygous, Trem2-R47H_CSS_homozygous, LOAD2

**Jax:** <https://www.jax.org/strain/028709>, <https://www.jax.org/strain/027918>, <https://www.jax.org/strain/030670>

------------------------------------------------------------------------

### Trem2-R47H-NSS

Identical to the Trem2-R47H mutation except that it changes two bases of AA 47 (C**GC** → C**AT**) instead of one, and introduces 9 silent mutations upstream, which we leave out of the intervals file.

**Location:** chr17: 48655584-48655585

**Sequence in interval:** GC → AT

**Genotype:** Trem2-R47H_NSS_homozygous

**Jax:** <https://www.jax.org/strain/034036>

**Model development:** <https://www.synapse.org/Synapse:syn25474065> (Trem2-R47H_NSS)

------------------------------------------------------------------------

## Excluded mutations

------------------------------------------------------------------------

### Adamts4 Enhancer KO

This is a deletion of \~2000 bases in a noncoding region of the genome near the mouse gene *Adamts4*. It cannot be detected from RNA seq data.

**Genotype:** LOAD2.Adamts4-KO_homozygous

**Jax:** <https://www.jax.org/strain/033868>

------------------------------------------------------------------------

### APOE4-KI

This is an insertion of a significant portion of the human APOE4 gene, and does not require an interval because it can be detected as expression of APOE4 in the RNA seq data.

**Genotype:** APOE4-KI_homozygous, LOAD2

**Jax:** <https://www.jax.org/strain/028709>, <https://www.jax.org/strain/030670>

------------------------------------------------------------------------

### Bin1 intron 1 SNP

This is a mutation in an intron of the mouse *Bin1* gene and cannot be detected from RNA seq data.

**Genotype:** LOAD2.Bin1-KI_KI/KI

**Jax:** <https://www.jax.org/strain/033869>

------------------------------------------------------------------------

### Cd2ap promoter SNP

This is a mutation in a noncoding region of the mouse *Cd2ap* gene and cannot be detected from RNA seq data.

**Genotype:** LOAD2.Cd2ap_KI/KI

**Jax:** <https://www.jax.org/strain/033873>

------------------------------------------------------------------------

### Clu-h2kbKI

A knock-in of \~2000 bases of the human CLU gene into the mouse *Clu* gene. This replaces exons 8 and 9 of the mouse gene with exons 8 and 9 of the human gene. The human sequence was further edited with SNP rs2279590 (A → G), which is in an intron and can’t be detected from RNA seq data. The knock-in can be detected as expression of human CLU in the RNA seq data.

**Genotype:** Clu-rs2279590_KI_homozygous

**Jax:** <https://www.jax.org/strain/037496>

------------------------------------------------------------------------

### **Epha1-exon1-SNP**

A LOAD mutation in the mouse *Epha1* gene. It is a silent mutation in AA 24 of the *Epha1-201* transcript (GC**G** → GC**T**, which are both Alanine).

**Location:** chr6: 42350073

**Genotype:** LOAD2.Epha1-KI_KI/KI

**Jax:** <https://www.jax.org/strain/033875>

> After testing, I found that this mutation is being detected in the Epha1-KI mice as a “deletion”, and is detected as a deletion in other genotypes as well, which indicates it’s not reliable for validation.

------------------------------------------------------------------------

### Ptk2b intron 5 SNP

This is a LOAD SNP in an intron of the mouse *Ptk2b* gene and is not detectable from RNA seq data.

**Genotype:** LOAD2.Ptk2b-intron5SNP_homozygous

**Jax:** <https://www.jax.org/strain/033913>

------------------------------------------------------------------------

### Scimp upstream SNP

This is a LOAD SNP in a noncoding region upstream of the mouse *Scimp* gene. Therefore it is not detectable from RNA seq reads.

**Genotype:** LOAD2.Scimp-upstreamSNP_homozygous

**Jax:** <https://www.jax.org/strain/032775>

------------------------------------------------------------------------

### Spi1\*rs1377416

This is a LOAD SNP in a noncoding region upstream of the mouse *Spi1* gene. Therefore it is not detectable from RNA seq reads.

**Genotype:** Spi1-rs1377416_homozygous

**Jax:** <https://www.jax.org/strain/035873>

**Model development:** <https://www.synapse.org/Synapse:syn25570580>

------------------------------------------------------------------------

### Tau-H2

This is likely a knock-in of human MAPT, H2 haplotype, but the [link to the Jax page](https://www.jax.org/strain/036601) on the LOAD2 wiki doesn’t work.

**Genotype:** LOAD2.hTau-H2_homozygous

------------------------------------------------------------------------

### tauP301L / MAPT

This is an insertion of human MAPT, and does not require an interval because it can be detected as expression of MAPT in the RNA seq data.

**Genotype:** 3xTg-AD_carrier

**Jax:** <https://www.jax.org/strain/004807>
