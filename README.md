# 16S_Pipeline_RD


## Description

AnimalBiome's taxonomic assignment methods for pacbio full legnth 16S Hifi reads.

### Environment preparation

* Install qiime2-2023.5
* Install seqtk (if it doesn't come with qiime2)
* Install R
* Install R packages tidyverse
* Install perl

### SOP

Starting from  a set of raw fastq files demultiplexed.



#### Step 0 - trim sequences to target hifi reads within length parameters for our amplicon
```
seqtk comp $file | awk '{ if ((\$2 >= 1300) && (\$2 <= 1600)) { print} }' | cut --fields 1 > reads_list/$2.list";
seqtk subseq $file reads_list/$2.list > reads_trimmed_1/demultiplex.$1."."_trim.fastq";
```

We use 1300 to 1600 bp

#### Step 1 - make a qiime manifest file (make_pacbio_manifest.pl)
This is a tab delimited file that needs to have special headers
"sample-id"[TAB]"absolute-filepath"

make_pacbio_manifest.pl loops through the files in the directory output from step0

#### Step 2 - Load fastq files into qiime2 archive format.
Single threaded can be slow
```
qiime tools import \
       --type SampleData[SequencesWithQuality] \
       --input-path PacBioCCSmanifest.tsv \
       --output-path reads_qza/trimmed_reads.qza \
       --input-format SingleEndFastqManifestPhred33V2
```
       
#### Step 3 - Denoise reads using denoise-ccs from qiime2.

```
 qiime dada2 denoise-ccs \
   --i-demultiplexed-seqs reads_qza/trimmed_reads.qza \
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

#### Step 4 - Assign Taxonomy
##### Step 4.1 - SKlearn

```
qiime feature-classifier classify-sklearn \
      --i-classifier step4_output_naive_bayes_classifier.qza \
      --i-reads dada2_output/representative_sequences.qza \
      --o-classification sklearn_taxonomy.qza \
      --p-confidence 0 \
      --p-n-jobs 45
```

##### Step 4.2 - 2 rounds of Vsearch

Round 1 - very stringent on coverage identity
``` 
qiime feature-classifier classify-consensus-vsearch \
       --i-query /home/qiime2/dada2_output/representative_sequences.qza  \
       --i-reference-reads /home/qiime2/step4_output_seqs_addition.qza \
       --i-reference-taxonomy /home/qiime2/step4_output_taxonomy.qza  \
       --p-threads 45 --o-classification vsearch_R1.TAX.2.qza \
       --p-perc-identity 0.99 \
       --p-strand plus \
       --p-maxaccepts 3 \
       --p-query-cov 0.97 \
       --p-maxrejects 1000 \
       --p-top-hits-only \
       --output-dir Vsearch_outDir_Round1.2
```

Round 2 - slightly more relaxed with larger consensus pool
``` 
qiime feature-classifier classify-consensus-vsearch \
       --i-query /home/qiime2/dada2_output/representative_sequences.qza  \
       --i-reference-reads /home/qiime2/step4_output_seqs_addition.qza \
       --i-reference-taxonomy /home/qiime2/step4_output_taxonomy.qza \
       --p-threads 45 \
       --o-classification vsearch_R2.TAX.2.qza \
       --p-perc-identity 0.99 \
       --p-strand plus \
       --p-maxaccepts 10 \
       --p-query-cov 0.9 \
       --p-maxrejects 250 \
       --p-top-hits-only \
       --output-dir Vsearch_outDir_Round2.2
```


##### Step 4.2.2 - merge and prepare the taxonomy and ASV table

* ```
  qiime feature-table merge-taxa --i-data vsearch_R1.TAX.3.qza --i-data vsearch_R2.TAX.3.qza --o-merged-data vsearch_merged.qza
  ```
* ```
  unzip -d vsearch_merged_dir vsearch_merged.qza
  ```
* ```
  unzip -d sklearn_taxonomy_dir sklearn_taxonomy.qza
  ```
* ```
  unzip -d biom_table_dir dada2_output/table.qza
  biom convert --to-tsv -i biom_table_dir/[random_string]/data/feature-table.biom -o test.biom
  ```
##### Step 5 - Roll back the Taxonomy

Hard coded paths to point the right files for now 

`Rscript Genus_Confidence_Rollback.Rmd`;
