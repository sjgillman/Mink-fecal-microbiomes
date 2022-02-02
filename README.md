# Mink (*Neovison vison*) fecal microbiomes are influenced by sex, temperature and time post-defecation 2022
#### Published Manuscript in Journal of Mammalogy; [DOI: gyab140](https://doi.org/10.1093/jmammal/gyab140)
Authors: Diana J.R. Lafferty, Sierra J. Gillman, Lane K. Jeakle, Brian J. Roell, Erin A. McKenney
All infromation can also be downloaded from [zenodo](https://zenodo.org/record/4560812#.YfnOc_XMJhE)

**Abstract**
Gut microbiomes encode myriad metabolic functions critical to mammalian ecology and evolution. While fresh fecal samples provide an efficient, noninvasive method of sampling gut microbiomes, collecting fresh feces from elusive species is logistically challenging. Nonfresh feces, however, may not accurately represent the gut microbiome of the host due to succession of gut microbial consortia postdefecation as well as colonization by microbes from the surrounding environment. Using American mink (Neovison vison) as a model species, we examined postdefecation microbial community succession to learn how ambient temperature and temporal sampling constraints influence the reliability of nonfresh feces to represent host gut microbiomes. To achieve our goal, we analyzed fresh mink feces (n = 5 females; n = 5 males) collected at the time of defecation from captive mink at a farm in the Upper Peninsula of Michigan and we subsequently subsampled each fecal specimen to investigate microbial community succession over five days, under both warm (21°C) and cold (–17°C to –1°C) temperature treatments. We found that both temperature and time influenced fecal microbiome composition; and we also detected significant sexual dimorphism in microbial community structures, with female mink microbiomes exhibiting significantly greater variation than males’ when exposed to the warm temperature treatment. Our results demonstrate that feces from unknown individuals can be a powerful tool for examining carnivore gut microbiomes, though rigorous study design is required because sex, ambient temperature, and time since defecation drive significant microbial variation and the sample size requirements necessary for detecting statistically significant differences between target populations is an important consideration for future ecologically meaningful research.

Directory structure | Description
--- | ---
minkr-gmb-JoM/
  README.md
  **data/** | **Description**
  *SD1 Mink Metadata.tsv* | the metadata file including sampleID, stable carbon and nitrogen values, sex, etc for each sample
  *SD3_Demultiplexed_Sequences.tar.gz* | demultiplexed EMP-paired end sequences demultiplexed on QIIME2
  *
  *taxonomy.qza* | QIIME2 artifact
  *rooted-tree.qza* | QIIME2 artifact created with MAFFT
  *phyloseq.rds* | phyloseq object created with MetaFile.tsv, OTU_table.qza, taxonomy.qza, and rooted-tree.qza and used for downstream analysis in R.
  *phyloseq_srs.rds* | normalized "phyloseq.rds" object.
  **scripts/** | **Description**
  *SD4 Mink QIIME2 Pipeline.md* | bioinformatic pipeline to prepare sequences for analysis in R
  * SD5 Mink Supplemental script.R* | code for analysis and figure creation in R.
  **images/**
  *gmb.png*

<p align="center">
<img src="images/gmb.png" width="500" />
  </p>


