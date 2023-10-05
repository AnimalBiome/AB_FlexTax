# FlexTax: A flexible hybrid approach to assigning species-level taxonomy to full-length 16S rRNA gene sequences

## Description

This documentation describes AnimalBiome's process for assigning taxonomy to PacBio full-length 16S rRNA gene HiFi reads. We used both the sklearn and VSEARCH taxonomic classifiers in QIIME2.

We first discuss curation of the reference database, then the processing of sequences in QIIME2 with DADA2, then the implementation of the sklearn and VSEARCH classifiers, and lastly the creation of the final consensus taxonomy for our sequences.


## Process overview and output

![RD_SOP](https://github.com/AnimalBiome/AB_FlexTax/assets/67349833/ab40b5b8-ce6a-4fc3-b27f-f054558dff2f)

### Software used in this tutorial

*   [qiime2-2023.5](https://qiime2.org/)
    
*   [seqtk](https://github.com/lh3/seqtk) (if it doesn't come with qiime2)
    
*   [R](https://www.r-project.org/)
    
*   R package [tidyverse](https://www.tidyverse.org/)
    

## #1: Curate Reference Database

We have not tested this pipeline with uncurated reference databases, and results may vary depending on the level of curation applied. Uncurated reference databases should be used with caution.

We chose the [SILVA 138.1 NR99](https://www.arb-silva.de/no_cache/download/archive/release_138/Exports/) database and performed several curation steps to ensure that:

*   Reference sequences were full length or nearly full length versions of the 16S rRNA gene
    
*   Species-level taxa without formalized names or with uninformative names were removed
    
*   Non-prokaryotic taxa were removed
    

### Database curation process

1.1 Download the [Silva Reference](https://www.arb-silva.de/download/archive/)

**Length filtering**

1.2 Apply size filter to fasta file; we selected to retain sequences between 1300 and 1600 bp

```
seqtk comp Silva_nr99.fa | awk '{ if (($2 >= 1300) && ($2 <= 1600)) { print} }' | cut --fields 1 > SeqsToKeep.list";
seqtk subseq Silva_nr99.fa SeqsToKeep.list > Step1_Seqs.fa
```

**Uninformative taxonomy**

1.3 Manually review the filtered taxonomy file for uninformative designations, species without formalized names, etc. Add these to a text file ([exclude\_list.txt](https://github.com/AnimalBiome/16S_Pipeline_RD/blob/main/exclude_list.txt)). These were some of our names on our list:

```
s__metagenome
s__bacterium
s__uncultured
s__unclassified
s__unculturable
s__symbiont
s__soil
```

1.4 Exclude sequences from the taxonomy file that match the exclude\_list.txt

```
grep -f exclude_list.txt -v taxonomy_file > Step2_taxonomy.txt
```

**Non-prokaryotic taxonomy and mis-matched taxonomy**

1.5 Manually examine remaining unique taxonomy strings and exclude sequences with:

1.  Conflicting taxonomy between genus rank and species rank
    
    1.  ```
        d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Salmonella; s__Pantoea_agglomerans
        ```
        
2.  Eukaryotic species names
    
    1.  ```
        d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Enterococcaceae; g__Enterococcus; s__Heliconius_elevatus
        # this is an insect
        ```
        
3.  Uninformative species names (those not already caught by step 4)
    
    1.  ```
        d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Methyloligellaceae; g__uncultured; s__alpha_proteobacterium
        ```
        

1.6 Remove all of the discarded fasta sequences from all steps under 1.5 from the file `Step1_Seqs.fa` . Remove discarded entries from the Silva Taxonomy file.

_Overall, our manual curation reduced the size of the SILVA 138.1 NR99 database from 510,508 sequences to 116,173 sequences, and from 80,381 unique taxa to 17,746 unique taxa._

**Add any missing references or references of interest**

1.7 We added both COT (canine oral taxa \[Dewhirst et al., 2012\]) and FOT (feline oral taxa \[Dewhirst et al., 2015\]), _Peptoclostridium hiranonis_ and _Prevotella copri_ as “targeted” sequences to both the sequence and taxonomy files.

**Train QIIME2 sklearn classifier on this reference database**

1.8 Compile qiime archives

```
qiime tools import --type 'FeatureData[Taxonomy]' --input-path output_taxononmy_add.txt --output-path output_taxonomy.qza
qiime tools import --type 'FeatureData[Sequence]' --input-path output_seqs_add.fa --output-path output_seqs_addition.qza
```

1.9 Train model

```
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads output_seqs_addition.qza --i-reference-taxonomy output_taxonomy.qza --o-classifier output_naive_bayes_classifier.qza
```

## #2: Process PacBio sequences with QIIME2

Starting from a set of demultiplexed raw fastq files.

#### Step 2.1 - Trim sequences to only retain HiFi reads within appropriate length parameters

```
seqtk comp Sample1.hifi_read.fastq | awk '{ if (($2 >= 1300) && (\$2 <= 1600)) { print} }' | cut --fields 1 > Sample1.list";
seqtk subseq Sample1.hifi_read.fastq Sample1.list > Sample1_trimmed.hifi_reads.fq;
```

We chose to retain HiFi reads with lengths of 1300 to 1600 bp (i.e., those that are close to full-length 16S rRNA reads).

This step is recommended because we identified some concatenated amplicon fragments that were being incorrectly assigned to some barcodes. The succession of barcode sequences in line with the amplicon can be misinterpreted by the adaptor finder and therefore wrongly associate amplicon sequences to incorrect sampleIDs. In order to avoid this problem, we opted to only keep fragments within the target range of the 16S rRNA amplicon length.

#### Step 2.2 - Make a qiime manifest file

This manifest is a tab delimited file that needs to have a header with the following format:  
"sample-id"\[TAB\]"absolute-filepath"

#### Step 2.3 - Store fastq files as a qiime2 archive

This step is single threaded and can be slow. This will combine all the fastq files listed in your manifest file into a qiime2 archive and link each of them to a sample ID that will be referenced by qiime2 from now on.

```
qiime tools import \
       --type SampleData[SequencesWithQuality] \
       --input-path PacBioCCSmanifest.tsv \
       --output-path trimmed_reads.qza \
       --input-format SingleEndFastqManifestPhred33V2
```

#### Step 2.4 - Denoise reads using denoise-ccs

```
 qiime dada2 denoise-ccs \
   --i-demultiplexed-seqs trimmed_reads.qza \
   --p-trunc-len 0 \
   --p-max-ee 3 \
   --p-n-threads 45 \
   --p-n-reads-learn 100000 \
   --p-pooling-method 'pseudo' \
   --p-chimera-method 'pooled' \
   --p-front AGRGTTYGATYMTGGCTCAG \
   --p-adapter RGYTACCTTGTTACGACTT \
   --p-min-len 1300 \
   --p-max-len 1600 \
   --output-dir dada2_output \
   --verbose
```

This step will QC the fastq reads to orient them the same direction (PacBio does not output reads in a set direction) as well as identify the sequencing primers and remove them. No length filtering is applied in the adaptor trimming step, and only 3 errors are allowed. Any read not satisfying these conditions will be automatically discarded. We apply the same length filter as before however this is applied to the post QC read length (after adaptor trimming).  
The second part of this command will denoise the reads by evaluating the redundancies across the samples using different methods. We used pseudo pooling for this step. Refer to the qiime2 or dada2 documentation for details.  
The third part of this command refers to chimera detection where reads will be evaluated using the pooled method.  
Steps 2.3 and 2.4 and multithreaded and will use all 45 threads specified.

**Note:** Files with fewer than 10 HiFi reads should be removed prior to processing. Very low read counts may cause this step to fail if there are not enough nucleotides left in the array for QC evaluation, triggering an error.

**At end of this, you will end up with a table of unique ASV sequences and their counts in each sample.**

## #3: Assign Taxonomy to ASV sequences with the sklearn and VSEARCH qiime classifiers

#### Step 3.1 - Run sklearn classifier with default confidence threshold (0.7)

```
qiime feature-classifier classify-sklearn \
      --i-classifier output_naive_bayes_classifier.qza \
      --i-reads dada2_output/representative_sequences.qza \
      --o-classification sklearn_taxonomy.qza \
      --p-confidence 0.7 \
      --p-n-jobs 45
```

This step assigns a taxonomic label to each of the representative sequences detected by dada2 in the previous step. Each ASV is evaluated and assigned a label that will be as precise as possible while still maintaining a 0.7 confidence level of the results. If a 0.7 confidence level is not achieved, it will defer to the parent taxonomic rank (e.g. species label will be “Unassigned” but Genus will be assigned)

0.7 was picked in accordance to the author’s evaluation and descriptions in [https://doi.org/10.1186/s40168-018-0470-z](https://doi.org/10.1186/s40168-018-0470-z).

#### Step 3.2 - Run two rounds of VSEARCH classification

**Round 1 -** Very stringent on coverage and percent identity. We are using a top-hit approach here with lots of rejects to make sure we are evaluating enough candidate sequences. Because the thresholds are so stringent, we are not expecting to get many assignments out of this round of VSEARCH (i.e., many ASV sequences will be “Unassigned”). The label will be picked from the consensus of the 3 accepted results (51% majority).  
We can reduce processing time by targeting only the plus strand thanks to the denoising step that was run prior to this.

```
qiime feature-classifier classify-consensus-vsearch \
       --i-query dada2_output/representative_sequences.qza  \
       --i-reference-reads output_seqs_addition.qza \
       --i-reference-taxonomy output_taxonomy.qza  \
       --p-threads 45 --o-classification vsearch_R1.TAX.qza \
       --p-perc-identity 0.99 \
       --p-strand plus \
       --p-maxaccepts 3 \
       --p-query-cov 0.97 \
       --p-maxrejects 1000 \
       --p-top-hits-only \
       --output-dir Vsearch_outDir_Round1
```

**Round 2 -** Slightly more relaxed with larger consensus pool. This is a slightly more relaxed approach on the coverage. We are also accepting more results to build our consensus from. Because the thresholds are a bit more relaxed we are not needing to reject as many candidates, and keeping the maxrejects value low will improve run times.

```
qiime feature-classifier classify-consensus-vsearch \
       --i-query dada2_output/representative_sequences.qza  \
       --i-reference-reads output_seqs_addition.qza \
       --i-reference-taxonomy output_taxonomy.qza \
       --p-threads 45 \
       --o-classification vsearch_R2.TAX.qza \
       --p-perc-identity 0.99 \
       --p-strand plus \
       --p-maxaccepts 10 \
       --p-query-cov 0.9 \
       --p-maxrejects 250 \
       --p-top-hits-only \
       --output-dir Vsearch_outDir_Round2
```

#### Step 3.3 - Merge the two VSEARCH taxonomy tables

Merge the two VSEARCH taxonomy assignments together. The first data will take precedence over the next ones listed, so order matters here. We want to make sure our Round1 is being used before Round2’s labels.

```
unzip vsearch_R1.TAX.qza
unzip -d vsearch_round1_dir vsearch_R1.TAX.qza
grep -v 'Unassigned' vsearch_round1_dir/*/data/taxonomy.tsv > tax_R1_filteredTax.tsv
qiime tools import --type 'FeatureData[Taxonomy]' --input-path tax_R1_filteredTax.tsv --output-path tax_R1_filteredTax.qza
rm -r vsearch_round1_dir

qiime feature-table merge-taxa --i-data tax_R1_filteredTax.qza --i-data vsearch_R2.TAX.qza --o-merged-data vsearch_merged.qza
```

#### Step 3.4 - Prepare files for import into R

Next we want to prepare the directories so we can read the input files from R when creating the consensus taxonomy. Note that in the `biom_convert` command you will need to change the \[random\_string\] value to reflect the name of your unzipped `table.qza` object.

```
unzip -d vsearch_merged_dir vsearch_merged.qza
```

```
unzip -d sklearn_taxonomy_dir sklearn_taxonomy.qza
```

```
unzip -d biom_table_dir dada2_output/table.qza
biom convert --to-tsv -i biom_table_dir/[random_string]/data/feature-table.biom -o test.biom.tsv
```

## #4: Create Consensus Taxonomy

Resulting biom table was joined with the taxonomy tables from VSEARCH and sklearn by ASV identifier in R. The following code was used to create a consensus taxonomy from the sklearn and VSEARCH taxonomic assignments:

```
Taxon_MERGED_rolled_0.7conf = case_when(Taxon_VSEARCH == Taxon_SKLEARN ~ Taxon_VSEARCH,
                                        Taxon_VSEARCH != "Unassigned" & Taxon_VSEARCH != Taxon_SKLEARN & Consensus_numeric == 1 ~ Taxon_VSEARCH,
                                        Taxon_VSEARCH == "Unassigned"  ~ Taxon_SKLEARN,
                                        TRUE ~ Taxon_SKLEARN)
```

VSEARCH taxonomic assignments are highly accurate but because of the stringency, many are “Unassigned”. sklearn gives very accurate taxonomic assignment but sometimes not to species or genus level. We opted to use the sklearn taxonomy as a default, but replaced it with the VSEARCH assignment if the VSEARCH taxonomy was of high confidence (consensus=1). And that is what the above code is essentially doing.

## References

Bolyen, E., Rideout, J.R., Dillon, M.R. _et al._ Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. _Nat Biotechnol_ **37**, 852–857 (2019). [https://doi.org/10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9)

[Scikit-learn: Machine Learning in Python](https://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html), Pedregosa _et al._, JMLR 12, pp. 2825-2830, 2011.

Bokulich, N.A., Kaehler, B.D., Rideout, J.R. _et al._ Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. _Microbiome_ **6**, 90 (2018). [https://doi.org/10.1186/s40168-018-0470-z](https://doi.org/10.1186/s40168-018-0470-z)

Dewhirst, F. E., Klein, E. A., Thompson, E. C., Blanton, J. M., Chen, T., Milella, L., et al. (2012). The canine oral microbiome. _PLoS One_ 7:e36067.

Dewhirst, F. E., Klein, E. A., Bennett, M. L., Croft, J. M., Harris, S. J., and Marshall-Jones, Z. V. (2015). The feline oral microbiome: a provisional 16S rRNA gene based taxonomy with full-length reference sequences. _Vet. Microbiol._ 175, 294–303. doi: 10.1016/j.vetmic.2014.11.019

McKinney, W., & others. (2010). Data structures for statistical computing in python. In _Proceedings of the 9th Python in Science Conference_ (Vol. 445, pp. 51–56).

Bokulich, N.A., Kaehler, B.D., Rideout, J.R. _et al._ Optimizing taxonomic classification of marker-gene amplicon sequences with QIIME 2’s q2-feature-classifier plugin. _Microbiome_ **6**, 90 (2018). [https://doi.org/10.1186/s40168-018-0470-z](https://doi.org/10.1186/s40168-018-0470-z)

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. VSEARCH: a versatile open source tool for metagenomics. PeerJ. 2016 Oct 18;4:e2584. doi: 10.7717/peerj.2584. PMID: 27781170; PMCID: PMC5075697.

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods. 2016 Jul;13(7):581-3. doi: 10.1038/nmeth.3869. Epub 2016 May 23. PMID: 27214047; PMCID: PMC4927377.
