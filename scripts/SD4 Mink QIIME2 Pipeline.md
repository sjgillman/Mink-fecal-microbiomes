# Bioinformatic pipeline for "Mink fecal microbiomes are influenced by sex, temperature and time post-defecation"
## QIIME2 verion: 2019.4
Samples are from Argonne National Laboratory
Pipeline adapted from qiime2 tutorial "Atacama soil microbiome" $ "Moving Pictures"
samples are EMP-Paired end multiplexed sequences with new primer set 
w/ barcodes read forward & no longer reversed in demux step

**Import**
```
qiime tools import \
--type EMPPairedEndSequences \
--input-path Reads \
--output-path mink-sequences.qza #you can name this whatever you want
##output aritfact: mink-sequences.qza
```

**Demultiplexing Sequences**
you will need metadata/mapping file
if you would like to only look at specific samples you can modify the metadata file
to only include those samples, the final demux.qza file will only contain those files
alternative for filter select samples after demultiplexing is provided further down
google sheets has an addin in called Keemei that can validate files
use this to ensure file is formatted correctly/ download as .tsv

```
qiime demux emp-paired \
--m-barcodes-file MinkMeta.tsv \
--m-barcodes-column BarcodeSequence \
--p-no-golay-error-correction \
--i-seqs mink-sequences.qza \
--o-per-sample-sequences minkdemux.qza \
--o-error-correction-details minkdemux-detail.qza

## make a summary visualization 

qiime demux summarize \
--i-data minkdemux.qza \ 
--o-visualization minkdemuxseq.qzv
```

**Denoising**
```
qiime dada2 denoise-paired \
--i-demultiplexed-seqs minkdemux.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 150 \
--p-trunc-len-r 150 \
--o-table minktable.qza \
--o-representative-sequences minkrep-seqs.qza \
--o-denoising-stats minkdenoising-stats.qza


## summary visulization table for determining sample depth for rarifying
qiime feature-table summarize \
--i-table minktable.qza \
--o-visualization minktable.qzv \
--m-sample-metadata-file MinkMeta.tsv


# make visualization artifacts of rep seq 
qiime feature-table tabulate-seqs \
--i-data minkrep-seqs.qza \
--o-visualization minkrep-seqs.qzv
#output visualization: rep-seq.qzv


#  view denoising stats
qiime metadata tabulate \
--m-input-file minkdenoising-stats.qza \
--o-visualization minkdenoising-stats.qzv
#output visualization: denoising-stats.qzv
```

**Generating a tree for phylogenetic diversity analyses**
```
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences minkrep-seqs.qza \
--o-alignment aligned-rep-mink-seqs.qza \
--o-masked-alignment masked-aligned-rep-mink-seqs.qza \
--o-tree unrooted-mink-tree.qza \
--o-rooted-tree rooted-mink-tree.qza
```

**Taxonomic Analysis sklearn**
import reference otus

```
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path SILVA_132_99_16S.fna \
--output-path SILVA_OTU.qza

## Import reference taxonomy file
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path taxonomy_7_levels.txt \
--output-path ref-taxonomy.qza



##  Extract reference reads
qiime feature-classifier extract-reads \
--i-sequences SILVA_OTU.qza \
--p-f-primer GTGCCAGCMGCCGCGGTAA \
--p-r-primer GGACTACHVGGGTWTCTAAT \
--p-trunc-len 150 \
--p-min-length 100 \
--p-max-length 400 \
--o-reads ref-seqs.qza

## Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy ref-taxonomy.qza \
--o-classifier mink-classifier.qza

## Test Classifier
qiime feature-classifier classify-sklearn \
--i-classifier mink-classifier.qza \
--i-reads minkrep-seqs.qza \
--o-classification mink-taxonomySILVA.qza

## fixing white spaces
qiime tools export \
--input-path mink-taxonomySILVA.qza \
--output-path taxonomy-with-spaces

qiime metadata tabulate \
--m-input-file taxonomy-with-spaces/taxonomy.tsv  \
--o-visualization taxonomy-as-metadata.qzv

qiime tools export \
--input-path taxonomy-as-metadata.qzv \
--output-path taxonomy-as-metadata

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-path taxonomy-as-metadata/metadata.tsv \
--output-path mink-taxonomy-without-spaces.qza


## create visualization
qiime metadata tabulate \
--m-input-file mink-taxonomy-without-spaces.qza \
--o-visualization mink-taxonomySILVA.qzv
```

**Filtering** 
```
#  filter out mitochondria and chloroplast
qiime taxa filter-table \
--i-table minktable.qza \
--i-taxonomy mink-taxonomy-without-spaces.qza \
--p-exclude mitochondria \
--o-filtered-table mink-table-filter.qza

qiime taxa filter-table \
--i-table mink-table-filter.qza \
--i-taxonomy mink-taxonomy-without-spaces.qza \
--p-exclude chloroplast \
--o-filtered-table clean-mink-table.qza

# get rid of unassigned
qiime taxa filter-table \
--i-table clean-mink-table.qza \
--i-taxonomy mink-taxonomy-without-spaces.qza \
--p-exclude Unassigned \
--o-filtered-table clean-mink-table-unassigned-rm.qza

# barplot
qiime taxa barplot \
--i-table clean-mink-table-unassigned-rm.qza \
--i-taxonomy mink-taxonomy-without-spaces.qza \
--m-metadata-file MinkMeta.tsv \
--o-visualization mink-taxa-bar-plotsSILVA-clean.qzv

#you can view any visual on qiime2view.com

#  visualize
qiime feature-table summarize \
--i-table clean-mink-table-unassigned-rm.qza \
--o-visualization clean-mink-table-unassigned-rm.qzv \
--m-sample-metadata-file MinkMeta.tsv

#  remove Bacteria only assigned
qiime taxa filter-table \
--i-table clean-mink-table-unassigned-rm.qza \
--i-taxonomy mink-taxonomy-without-spaces.qza \
--p-mode exact \
--p-exclude D_0__Bacteria \
--o-filtered-table clean-mink-table-unassigned_Unknown-rm.qza


#  remove Bacteria only assigned
qiime taxa filter-table \
--i-table clean-mink-table-unassigned_Unknown-rm.qza \
--i-taxonomy mink-taxonomy-without-spaces.qza \
--p-exclude D_0__Archaea \
--o-filtered-table clean-mink-table-unassigned_Unknown_Arch-rm.qza

#  barplot
qiime taxa barplot \
--i-table clean-mink-table-unassigned_Unknown_Arch-rm.qza \
--i-taxonomy mink-taxonomy-without-spaces.qza \
--m-metadata-file MinkMeta.tsv \
--o-visualization mink-taxa-bar-plotsSILVA-clean2.qzv
```




