### Mink fecal microbiomes are influenced by sex, temperature and time post-defecation
# Supplemental Analysis Script
#Diana J.R. Lafferty1, Sierra J. Gillman, Lane K. Jeakle, Brian J. Roell, Erin A. McKenney



library(microbiome) ## data analysis
library(qiime2R) # import data
library(phyloseq)
library(vegan) 
library(tidyverse)
library(ggplot2)
library(mctoolsr)
library(picante) ## faith's PD
library(Rmisc)## summary stats
library(picante)
library(cowplot)



setwd("~/Desktop/Projects/Mink/Longitudinal/Mink-R")
#### Import & create phyloseq dataframe with qiime2R and QIIME2 artifacts #####
## Following Tutorial: Integrating QIIME2 and R for data
# visualization and analysis using qiime2R by J. Bisanz
## you will need
# 1.) Metafile.tsv (Cold)
# 2.) taxonomy.qza
# 3.) table.qza (Cold)
# 4.) rooted.qza

metadata<-read_tsv("MinkMetaR.tsv")
SVs<-read_qza("clean-mink-table-unassigned_Unknown_Arch-rm.qza")
taxonomy<-read_qza("mink-taxonomy_renamed.qza")
taxtable<-taxonomy$data %>% as_tibble() %>%
  separate(Taxon, sep=";", c("Domain", "Phylum",
                             "Class", "Order",
                             "Family", "Genus",
                             "Species")) #convert the table into a tabular split version
tree<-read_qza("filter-rooted-mink-tree.qza")


## Create the phyloseq object
phy_objCom<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>%
              select(-Confidence) %>%
              column_to_rownames("Feature.ID") %>%
              as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))




##### Clean Taxonomy table #####
## Rename NAs to last known taxonomy rank
tax.clean <- data.frame(tax_table(phy_objCom))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
## import new taxonomy table
tax_table(phy_objCom) <- as.matrix(tax.clean)


## save phyloseq object
#saveRDS(phy_objCom, "physeqCom.rds")


set.seed(9242)
phyb.rar2 <- rarefy_even_depth(phy_objCom, sample.size = 10000)

summary(sample_sums(phyb.rar2)) ##checking to see they all have the same sequence depth
any(taxa_sums(phyb.rar2)== 0)


#### Community Composition ####

#subset to initial samples only
SubDat<-subset_samples(phyb.rar2, Temperature=="Initial")

##### Relative Abundance for whole community ####
## relative abundance
pseq.rel <- microbiome::transform(SubDat, "compositional")

## merge to phylum rank
phlyum <- tax_glom(pseq.rel, taxrank = "Phylum")
ntaxa(phlyum)

## melt
phylum_melt<- psmelt(phlyum)

## get summary statistics phyla
p_abund<-summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Phylum"))

## get summary stat for phyla for each sex

p_abundSex<-summarySE(phylum_melt, measurevar = "Abundance", groupvars =c("Phylum", "Sex"))


#export
#write_csv(p_abund, "PopLevelAbund.csv")
#write_csv(p_abundSex, "SexLevelAbund.csv")


##### population level diversity ####

##  pull metadata from physeq object
sam.meta <- meta(phyb.rar2)

## Add the rownames as a new colum for easy integration later.
sam.meta$sam_name <- rownames(sam.meta)

#### Non-phylogenetic diversities: Shannon ####
## calculated with microbiome package
## 
div_shan<- microbiome::alpha(phyb.rar2)
## can run index= "all" if you desire all alpha indices

## Add the rownames to diversity table
div_shan$sam_name <- rownames(div_shan)


#### Phylogenetic diversity: Faith's PD #####
#Phylogenetic diversity is calculated using the picante package.

## pull ASV table
phyb.rar.asvtab <- as.data.frame(phyb.rar2@otu_table)

## pull tree
phyb.rar.tree <- phyb.rar2@phy_tree

## We first need to check if the tree is rooted or not 

phyb.rar2@phy_tree
###rooted so we are good to go

## Getting the data ready
div_pd <- pd(t(phyb.rar.asvtab), phyb.rar.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the
#tree file comes from the first code  we used to read tree
#file (see making a phyloseq object section)

## Add the rownames to diversity table
div_pd$sam_name <- rownames(div_pd)

## STEP 4p. merge all of the alphas into one file
merged_table<-merge(div_pd,div_shan, by = "sam_name", all=T)
merged_table2<-merge(merged_table,sam.meta, by = "sam_name", all=T)
alpha_table <- merged_table2

#### summarize ###
ShanPop<-summarySE(alpha_table, measurevar = "diversity_shannon", groupvars =c("Temperature"))

PDPop<-summarySE(alpha_table, measurevar = "PD", groupvars =c("Temperature"))

################################## COLD Analysis ###############################################
## for cold vs warm analysis, samples were filtered in QIIME2 and re-rooted
# following was done in QIIME2:

qiime feature-table filter-seqs \
--i-data filtered-rep-seq.qza \
--i-table Cold-table.qza \
--o-filtered-data Cold-rep.qza


qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences Cold-rep.qza \
--o-alignment aligned-rep-mink-seqs.qza \
--o-masked-alignment masked-aligned-rep-Cold-seqs.qza \
--o-tree unrooted-Cold-tree.qza \
--o-rooted-tree Cold-tree.qza

## Data was then analyzed in R:

setwd("~/Desktop/Projects/Mink/Longitudinal/DiversityCold")

## import artifacts & metadata file
metadata<-read_tsv("ColdMeta-R.tsv")
SVs<-read_qza("Cold-table.qza")
taxonomy<-read_qza("mink-taxonomy_renamed.qza")
taxtable<-taxonomy$data %>% 
  as_tibble() %>% 
  separate(Taxon, sep=";", c("Domain", "Phylum",
                             "Class", "Order", 
                             "Family", "Genus",
                             "Species")) #convert the table into a tabular split version
tree<-read_qza("Cold-tree.qza")


## Create the phyloseq object
phy_objC<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% 
              select(-Confidence) %>% 
              column_to_rownames("Feature.ID") %>% 
              as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% 
                as.data.frame() %>% 
                column_to_rownames("SampleID")))


##### Clean Taxonomy table #####
## Rename NAs to last known taxonomy rank
tax.clean <- data.frame(tax_table(phy_objC))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
## import new taxonomy table
tax_table(phy_objC) <- as.matrix(tax.clean)


## save phyloseq object
#saveRDS(phy_objC, "physeqC.rds")

## if you ever want to pull back in
#phy_objC<- readRDS("physeqC.rds")

##  Equal sample sums
set.seed(9242) ## ensures rarifies the same each time script is run

summary(sample_sums(phy_objC)) ## helps determine depth for rarifying 

## rarefying: we already know our depth:
SubDatC <- rarefy_even_depth(phy_objC, sample.size = 10000)


####################################### Alpha Diversity Cold ######################################
##  pull metadata from physeq object
sam.metaC <- meta(SubDatC)

## Add the rownames as a new colum for easy integration later.
sam.metaC$sam_nameC <- rownames(sam.metaC)

#### Non-phylogenetic diversities: Shannon ####
## calculated with microbiome package

div_shanC<- microbiome::alpha(SubDatC)
## can run index= "all" if you desire all alpha indices

## Add the rownames to diversity table
div_shanC$sam_name <- rownames(div_shanC)

#### Phylogenetic diversity: Faith's PD #####
#Phylogenetic diversity is calculated using the picante package.

## pull ASV table
SubDatC.asvtab <- as.data.frame(SubDatC@otu_table)

## pull tree
SubDatC.tree <- SubDatC@phy_tree

## We first need to check if the tree is rooted or not 

SubDatC@phy_tree
###rooted so we are good to go

## Getting the data ready
div_pdC <- pd(t(SubDatC.asvtab), SubDatC.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the
#tree file comes from the first code  we used to read tree
#file (see making a phyloseq object section)

## Add the rownames to diversity table
div_pdC$sam_name <- rownames(div_pdC)

## STEP 4p. merge all of the alphas into one file
merged_tableC<-merge(div_pdC,div_shanC, by = "sam_name", all=T)
merged_table2C<-merge(merged_tableC,sam.metaC, by = "sam_name", all=T)
alpha_tableC <- merged_table2C

#write.csv(alpha_tableC, "cold_alphaC.csv")
#alpha_tableC<-read.csv("cold_alphaC.csv")

library(lmerTest)
library(car)
library(emmeans)
library(psycho)


ggplot(alpha_tableC,aes(x=PD))+geom_histogram()


pd_lme4C<-lmerTest::lmer(PD~Day*Sex+(1|Individual),data=alpha_tableC, REML=T)
anova(pd_lme4C)
##Type III Analysis of Variance Table with Satterthwaite's method
##Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
##Day     70.658 14.1316     5    40  0.9122 0.4829
##Sex     11.144 11.1438     1     8  0.7194 0.4210
##Day:Sex 37.742  7.5484     5    40  0.4873 0.7837

performance::r2(pd_lme4C)
#Conditional R2: 0.365
#Marginal R2: 0.102
summary(pd_lme4C)

resultsC <- analyze(pd_lme4C, CI = 95)

summary(resultsC) %>% 
  mutate(p = psycho::format_p(p))

PD_resC<-summary(resultsC) %>% 
  mutate(p = psycho::format_p(p))

write.csv(PD_resC, "PD_reC.csv")

print(resultsC)


## Assumption 1 - Linearity

## Graphically, plotting the model residuals (the difference 
#between the observed value and the model-estimated value) vs 
#the predictor 
alpha_tableC<-alpha_tableC
Plot.Model.F.Linearity<-plot(resid(pd_lme4C),alpha_tableC$PD)

## Assumption 2 Homogeneity of Variance
#Regression models assume that variance of the residuals
#is equal across groups. 

#extracts the residuals and places them in a new column in our original data table
alpha_tableC$lme10<- residuals(pd_lme4C) 
alpha_tableC$baslme10 <-abs(alpha_tableC$lme10) #creates a new column with the absolute value of the residuals
alpha_tableC$lme102 <- alpha_tableC$baslme10^2 #squares the absolute values of the residuals to provide the more robust estimate
pd_levenC <- lm(lme102 ~ Individual, data=alpha_tableC) #ANOVA of the squared residuals
anova(pd_levenC) #displays the results

##visually
plot(pd_lme4C) #creates a fitted vs residual plot

##Assumption 3: The residuals of the model are normally distributed.

#QQ plots 
library(lattice)
qqmath(pd_lme4C, id=0.05)
## overall looks good!!



ggplot(alpha_tableC,aes(x=diversity_shannon))+geom_histogram()


sh_lme4C<-lmerTest::lmer(diversity_shannon~Day*Sex+(1|Individual),data=alpha_tableC, REML=T)

anova(sh_lme4C)
#Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq  Mean Sq NumDF DenDF F value Pr(>F)
#Day     1.03461 0.206923     5    40  0.9466 0.4618
#Sex     0.14594 0.145944     1     8  0.6676 0.4375
#Day:Sex 0.11721 0.023442     5    40  0.1072 0.9901

performance::r2(sh_lme4C)
##Conditional R2: 0.812
#Marginal R2: 0.069
##
summary(sh_lme4C)




results <- analyze(sh_lme4C, CI = 95)

Shan_resC<-summary(resultsC) %>% 
  mutate(p = psycho::format_p(p))
write.csv(Shan_resC, "ShanreC.csv")
print(resultsC)



##Assumption 1 - Linearity

## Graphically, plotting the model residuals (the difference 
#between the observed value and the model-estimated value) vs 
#the predictor 

Plot.Model.F.Linearity<-plot(resid(sh_lme4C),alpha_tableC$diversity_shannon)

## Assumption 2 Homogeneity of Variance
#Regression models assume that variance of the residuals
#is equal across groups. 

#extracts the residuals and places them in a new column in our original data table
alpha_tableC$lme10<- residuals(sh_lme4C) 
alpha_tableC$baslme10 <-abs(alpha_tableC$lme10) #creates a new column with the absolute value of the residuals
alpha_tableC$lme102 <- alpha_tableC$baslme10^2 #squares the absolute values of the residuals to provide the more robust estimate
shan_levenC <- lm(lme102 ~ Individual, data=alpha_tableC) #ANOVA of the squared residuals
anova(pd_leven) #displays the results

##visually
plot(sh_lme4C) #creates a fitted vs residual plot

##Assumption 3: The residuals of the model are normally distributed.
qqmath(sh_lme4C, id=0.05)
## overall looks good!!


################################# Alpha Diversity Plotting For Cold ###############################
alphaC<-alpha_tableC
alphaC$Day<-factor(alphaC$Day, levels=c("day0", "day1", "day2", "day3","day4","day5"))


theme_USGS_box <- function(base_family = "Times New Roman", ...){
  theme_bw(base_family = base_family, ...) +
    theme(
      panel.grid = element_blank(),
      strip.background=element_rect(fill="white"),
      strip.text.x=element_text(size=12),
      plot.title = element_text(size = 8),
      axis.ticks.length = unit(-0.01, "in"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=12, color="black"),
      axis.text.y = element_text(size=11, color="black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")), 
      axis.text.x = element_text(angle = 30,hjust = .8,size=11, color="black",margin=unit(c(0.2,0.2,0.2,0.2), "cm")),
      axis.ticks.x = element_blank(),
      aspect.ratio = 1,
      legend.position="none"
    )}

p1<-ggplot(data = alphaC, aes(x = Day, y = PD, fill=Day)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6, color="black")+
  scale_fill_manual(values=c("grey27","deepskyblue", "blue","slateblue" ,"purple", "darkorchid4"))+
  theme_USGS_box()+scale_x_discrete(labels=c("day0"="Day 0", "day1"="Day 1", "day2"="Day 2", "day3"="Day 3", "day4"="Day 4", "day5"="Day 5"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,35))+
  ylab("Faith's PD")


p2<-ggplot(data = alphaC, aes(x = Day, y = diversity_shannon, fill=Day)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6, color="black")+
  scale_fill_manual(values=c("grey27","deepskyblue", "blue","slateblue" ,"purple", "darkorchid4"))+
  theme_USGS_box()+scale_x_discrete(labels=c("day0"="Day 0", "day1"="Day 1", "day2"="Day 2", "day3"="Day 3", "day4"="Day 4", "day5"="Day 5"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,4))+
  ylab("Shannon Diversity")

#tiff('Figure 1A_B.tiff', units="in", width=10, height=5, res=1200, compression = 'lzw')
plot_grid(p2, 
          p1,
          nrow = 2,rel_widths=c(1,1),rel_heights = c(1,1))
#dev.off()

#################################### Beta Diversity Cold #########################################

##permanova unweighted

uunifrac_dist = phyloseq::distance(SubDatC, method="unifrac", weighted=F)
adonis(uunifrac_dist ~ Temperature+Sex+Day, strata=sam.metaC$Individual,data = sam.metaC)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Temperature  1    0.1174 0.11738 0.81632  0.01372  0.112   
#Sex          1    0.3555 0.35551 2.47252 0.04156  0.003 **
#Day          4    0.4598 0.11496 0.79950 0.05376  0.009 **
#Residuals   53    7.6206 0.14379         0.89095          
#Total       59    8.5533                 1.0000        




## Test for hommogeneity
permutest(betadisper(uunifrac_dist, sam.meta$Day), strata=Individual)
#Response: Distances
##          Response: Distances
## Df  Sum Sq   Mean Sq     F N.Perm Pr(>F)
## Groups     5 0.03460 0.0069198 1.073    999  0.375
## Residuals 54 0.34823 0.0064488     

permutest(betadisper(uunifrac_dist, sam.meta$Sex), strata=Individual)

## Response: Distances
##             Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
## Groups     1 0.03733 0.037330 7.0448    999  0.008 **
##  Residuals 58 0.30734 0.005299                        

##weighted
wunifrac_dist = phyloseq::distance(SubDatC2, method="unifrac", weighted=T)
adonis(wunifrac_dist ~ Temperature+Day+Sex, strata=sam.meta$Individual,data = sam.meta)
#             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Temperature  1   0.00211 0.002113  0.1812 0.00317  0.232
# Day          4   0.00675 0.001687  0.1447 0.01013  0.439
# Sex          1   0.03943 0.039433  3.3813 0.05918  0.398
# Residuals   53   0.61808 0.011662         0.92753       
# Total       59   0.66638                  1.00000            


##### Figure 3A&B Beta: PCoA plot using the unweighted UniFrac as distance for Individuals with Sex ####
library(ape)
library(reshape2)

## pull ASV table
SubDatC.asvtab <- as.data.frame(SubDatC@otu_table)
## pull tree
SubDatC.tree <- SubDatC@phy_tree

# re-root tree
new_tre <- ape::multi2di(SubDatC.tree)
phy_tree(SubDatC)<-new_tre


## weighted Sex and Individual
wunifrac_dist = phyloseq::distance(SubDatC2, method="wunifrac")

wordination = ordinate(SubDatC2, method="PCoA", distance=wunifrac_dist)
p1<-phyloseq::plot_ordination(SubDatC2, uordination, color="Individual", shape="Sex", axes=c(1,2))+
  geom_point(size=3)+
  theme(legend.title=element_text(family="Times New Roman"),
        legend.text=element_text(family="Times New Roman"),
        legend.key = element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(),
        axis.title.y=element_text(size=12, family="Times New Roman"),
        axis.title.x=element_text(size=12, family="Times New Roman"),
        axis.text.y=element_text(size=10, family="Times New Roman", color="black"),
        axis.text.x=element_text(size=10, family="Times New Roman", color="black"))+
  scale_colour_manual(name = "Day",
                      labels = c("Day 0", "Day 1",
                                 "Day 2", "Day 3",
                                 "Day 4","Day 5"),
                      values = c("grey27","deepskyblue",
                                 "blue","slateblue" ,
                                 "purple", "darkorchid4"))




### same thing for sex and day unweighted
uunifrac_dist = phyloseq::distance(SubDatC2, method="unifrac")
uordination = ordinate(SubDatC2, method="PCoA", distance=uunifrac_dist)



p2<-phyloseq::plot_ordination(SubDatC2, uordination, color="Individual", shape="Sex", axes=c(1,2))+
  geom_point(size=3)+
  theme(legend.title=element_text(family="Times New Roman"),
        legend.text=element_text(family="Times New Roman"),
        legend.key = element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(),
        axis.title.y=element_text(size=12, family="Times New Roman"),
        axis.title.x=element_text(size=12, family="Times New Roman"),
        axis.text.y=element_text(size=10, family="Times New Roman", color="black"),
        axis.text.x=element_text(size=10, family="Times New Roman", color="black"))+
  scale_colour_manual(name = "Day",
                      labels = c("Day 0", "Day 1",
                                 "Day 2", "Day 3",
                                 "Day 4","Day 5"),
                      values = c("grey27","deepskyblue",
                                 "blue","slateblue" ,
                                 "purple", "darkorchid4"))


## group plots
library(extrafont)
#tiff('Figure 3A_B.tiff', units="in", width=10, height=5, res=1200, compression = 'lzw')
plot_grid(p1, 
          p2,
          nrow = 1,rel_widths=c(1,1.3),
          rel_heights = c(1,1),
          labels = c("A: Weighted", "B: Unweighted"),
          hjust=-.9, vjust=2.5)
#dev.off()


######################## Statisical Analysis for Beta Diversity Cold ############################

##permanova unweighted

uunifrac_dist <-UniFrac(SubDatC, weighted=F)
adonis(uunifrac_dist ~ Temperature+Sex+Day, strata=sam.metaC$Individual,data = sam.metaC)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Temperature  1    0.1174 0.11738 0.81632  0.01372  0.112   
#Sex          1    0.3555 0.35551 2.47252 0.04156  0.003 **
#Day          4    0.4598 0.11496 0.79950 0.05376  0.009 **
#Residuals   53    7.6206 0.14379         0.89095          
#Total       59    8.5533                 1.0000        




## Test for hommogeneity
permutest(betadisper(uunifrac_dist, sam.meta$Day), strata=Individual)
#Response: Distances
##          Response: Distances
## Df  Sum Sq   Mean Sq     F N.Perm Pr(>F)
## Groups     5 0.03460 0.0069198 1.073    999  0.375
## Residuals 54 0.34823 0.0064488     

permutest(betadisper(uunifrac_dist, sam.meta$Sex), strata=Individual)

## Response: Distances
##             Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
## Groups     1 0.03733 0.037330 7.0448    999  0.008 **
##  Residuals 58 0.30734 0.005299                        

##weighted
uunifrac_dist <-UniFrac(SubDatC, weighted=T)
adonis(wunifrac_dist ~ Temperature+Day+Sex, strata=sam.meta$Individual,data = sam.meta)
#             Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# Temperature  1   0.00211 0.002113  0.1812 0.00317  0.232
# Day          4   0.00675 0.001687  0.1447 0.01013  0.439
# Sex          1   0.03943 0.039433  3.3813 0.05918  0.398
# Residuals   53   0.61808 0.011662         0.92753       
# Total       59   0.66638                  1.00000            




####################################### WARM Analysis ###############################################

qiime feature-table filter-seqs \
--i-data filtered-rep-seq.qza \
--i-table Warm-table.qza \
--o-filtered-data Warm-rep.qza


qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences Warm-rep.qza \
--o-alignment aligned-rep-mink-seqs.qza \
--o-masked-alignment masked-aligned-rep-Warm-seqs.qza \
--o-tree unrooted-Warm-tree.qza \
--o-rooted-tree Warm-tree.qza

setwd("~/Desktop/Projects/Mink/Longitudinal/DiversityWarm")

## import artifacts & metadata file
metadata<-read_tsv("WarmMeta-R.tsv")
SVs<-read_qza("Warm-table.qza")
taxonomy<-read_qza("mink-taxonomy_renamed.qza")
taxtable<-taxonomy$data %>% as_tibble() %>% separate(Taxon, sep=";", c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) #convert the table into a tabular split version
tree<-read_qza("Warm-tree.qza")


## Create the phyloseq object
phy_objW<-phyloseq(
  otu_table(SVs$data, taxa_are_rows = T), 
  phy_tree(tree$data), 
  tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()), #moving the taxonomy to the way phyloseq wants it
  sample_data(metadata %>% as.data.frame() %>% column_to_rownames("SampleID")))


##### Clean Taxonomy table #####
## Rename NAs to last known group
tax.clean <- data.frame(tax_table(phy_objW))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
## import new taxonomy table
tax_table(phy_objW) <- as.matrix(tax.clean)


## save phyloseq object
#saveRDS(phy_objW, "physeqW.rds")

## if you ever want to pull back in
#phy_objW<- readRDS("physeqW.rds")

###############################################  Alpha Diversity Warm ##################################
##  Equal sample sums
set.seed(9242) ## ensures rarifies the same each time script is run

summary(sample_sums(phy_objW)) ## helps determine depth for rarifying 

SubDatW <- rarefy_even_depth(phy_objW, sample.size = 10000)

summary(sample_sums(SubDatW)) ##checking to see they all have the same sequence depth


##  pull metadata from physeq object
sam.metaW <- meta(SubDatW)

## Add the rownames as a new colum for easy integration later.
sam.metaW$sam_name <- rownames(sam.metaW)

#### Non-phylogenetic diversities: Shannon ####
## calculated with microbiome package
## 
div_shanW<- microbiome::alpha(SubDatW)
## can run index= "all" if you desire all alpha indices

## Add the rownames to diversity table
div_shanW$sam_name <- rownames(div_shanW)


#### Phylogenetic diversity: Faith's PD #####
#Phylogenetic diversity is calculated using the picante package.

## pull ASV table
phyb.rar.asvtab <- as.data.frame(SubDatW@otu_table)

## pull tree
phyb.rar.tree <- SubDatW@phy_tree


SubDatW@phy_tree
###rooted so we are good to go

## Getting the data ready
div_pdW <- pd(t(phyb.rar.asvtab), phyb.rar.tree,include.root=T) 
# t(ou_table) transposes the table for use in picante and the
#tree file comes from the first code  we used to read tree
#file (see making a phyloseq object section)

## Add the rownames to diversity table
div_pdW$sam_name <- rownames(div_pdW)

## STEP 4p. merge all of the alphas into one file
merged_tableW<-merge(div_pdW,div_shanW, by = "sam_name", all=T)
merged_table2W<-merge(merged_tableW,sam.metaW, by = "sam_name", all=T)
alpha_tableW <- merged_table2W



#write.csv(alpha_tableW, "warm_alpha.csv")

#alpha_tableW<-read.csv("warm_alpha.csv")


ggplot(alpha_tableW,aes(x=PD))+geom_histogram()


pd_lme4W<-lmerTest::lmer(PD~Day*Sex+(1|Individual),data=alpha_tableW, REML=T)

anova(pd_lme4W)
#Type III Analysis of Variance Table with Satterthwaite's method
#            Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
#Day     14.7769  2.9554     5    40  0.7118 0.6182
#Sex      1.0005  1.0005     1     8  0.2410 0.6367
#Day:Sex 26.8312  5.3662     5    40  1.2924 0.2864

#performance::r2(pd_lme4)
##Conditional R2: 0.601
##Marginal R2: 0.082

results <- analyze(pd_lme4W, CI = 95)

summary(pd_lme4W)
anova(pd_lme4W)
summary(resultsW) %>% 
  mutate(p = psycho::format_p(p))

PD_resW<-summary(resultsW) %>% 
  mutate(p = psycho::format_p(p))

write.csv(PD_resW, "PD_reW.csv")

print(resultsW)

emmeans(pd_lme4W, pairwise~Day)
##Assumption 1 - Linearity

## Graphically, plotting the model residuals (the difference 
#between the observed value and the model-estimated value) vs 
#the predictor 

Plot.Model.F.Linearity<-plot(resid(pd_lme4W),alpha_tableW$PD)

## Assumption 2 Homogeneity of Variance
#Regression models assume that variance of the residuals
#is equal across groups. 

#extracts the residuals and places them in a new column in our original data table
alpha_tableW$lme10<- residuals(pd_lme4W) 
alpha_tableW$baslme10 <-abs(alpha_tableW$lme10) #creates a new column with the absolute value of the residuals
alpha_tableW$lme102 <- alpha_tableW$baslme10^2 #squares the absolute values of the residuals to provide the more robust estimate
pd_levenW <- lm(lme102 ~ Individual, data=alpha_tableW) #ANOVA of the squared residuals
anova(pd_levenW) #displays the results

##visually
plot(pd_lme4W) #creates a fitted vs residual plot

##Assumption 3: The residuals of the model are normally distributed.

#QQ plots 
qqmath(pd_lme4W, id=0.05)
## overall looks good!!

## SHAN
ggplot(alpha_tableW,aes(x=diversity_shannon))+geom_histogram()
sh_lme4W<-lmerTest::lmer(diversity_shannon~Day*Sex+(1|Individual),data=alpha_tableW, REML=T)



anova(sh_lme4W)
##Type III Analysis of Variance Table with Satterthwaite's method
#          Sum Sq Mean Sq NumDF DenDF F value  Pr(>F)  
#Day     0.6434 0.12869     5    40  0.5278 0.75383  
#Sex     0.1567 0.15673     1     8  0.6428 0.44586  
#Day:Sex 3.9375 0.78750     5    40  3.2297 0.01528 *


emmeans(sh_lme4W, pairwise~Sex|Day)
#contrasts
#Day = day0:
#contrast      estimate    SE   df t.ratio p.value
#Female - Male    0.603 0.653 12.1  0.923  0.3740 

#Day = day1:
#contrast      estimate    SE   df t.ratio p.value
#Female - Male   -0.825 0.653 12.1 -1.264  0.2301 

#Day = day2:
#contrast      estimate    SE   df t.ratio p.value
#Female - Male   -0.590 0.653 12.1 -0.904  0.3836 

#Day = day3:
#contrast      estimate    SE   df t.ratio p.value
#Female - Male   -0.872 0.653 12.1 -1.336  0.2060 

#Day = day4:
#contrast      estimate    SE   df t.ratio p.value
#Female - Male   -0.796 0.653 12.1 -1.220  0.2458 

#Day = day5:
#contrast      estimate    SE   df t.ratio p.value
#Female - Male   -0.344 0.653 12.1 -0.526  0.6082 
#Degrees-of-freedom method: kenward-roger 


emmeans(sh_lme4W, pairwise~Day|Sex)

# Sex = Female:
# contrast    estimate    SE df t.ratio p.value
# day0 - day1   0.5606 0.312 40  1.795  0.4803 
# day0 - day2   0.2807 0.312 40  0.899  0.9445 
# day0 - day3   0.7017 0.312 40  2.247  0.2398 
# day0 - day4   0.5093 0.312 40  1.631  0.5838 
# day0 - day5   0.3141 0.312 40  1.006  0.9133 
# day1 - day2  -0.2799 0.312 40 -0.896  0.9452 
# day1 - day3   0.1411 0.312 40  0.452  0.9975 
# day1 - day4  -0.0513 0.312 40 -0.164  1.0000 
# day1 - day5  -0.2466 0.312 40 -0.790  0.9677 
# day2 - day3   0.4210 0.312 40  1.348  0.7566 
# day2 - day4   0.2286 0.312 40  0.732  0.9767 
# day2 - day5   0.0334 0.312 40  0.107  1.0000 
# day3 - day4  -0.1925 0.312 40 -0.616  0.9892 
# day3 - day5  -0.3877 0.312 40 -1.241  0.8139 
# day4 - day5  -0.1952 0.312 40 -0.625  0.9885 

# Sex = Male:
# contrast    estimate    SE df t.ratio p.value
# day0 - day1  -0.8668 0.312 40 -2.776  0.0826 
# day0 - day2  -0.9120 0.312 40 -2.920  0.0593 
# day0 - day3  -0.7731 0.312 40 -2.476  0.1559 
# day0 - day4  -0.8894 0.312 40 -2.848  0.0701 
# day0 - day5  -0.6321 0.312 40 -2.024  0.3474 
# day1 - day2  -0.0452 0.312 40 -0.145  1.0000 
# day1 - day3   0.0937 0.312 40  0.300  0.9996 
# day1 - day4  -0.0226 0.312 40 -0.072  1.0000 
# day1 - day5   0.2348 0.312 40  0.752  0.9739 
# day2 - day3   0.1389 0.312 40  0.445  0.9976 
# day2 - day4   0.0226 0.312 40  0.072  1.0000 
# day2 - day5   0.2800 0.312 40  0.896  0.9451 
# day3 - day4  -0.1163 0.312 40 -0.372  0.9990 
# day3 - day5   0.1411 0.312 40  0.452  0.9975 
# day4 - day5   0.2573 0.312 40  0.824  0.9613 






performance::r2(sh_lme4W)
#Conditional R2: 0.797
#Marginal R2: 0.112

resultsW <- analyze(sh_lme4W, CI = 95)

Shan_resW<-summary(resultsW) %>% 
  mutate(p = psycho::format_p(p))
write.csv(Shan_resW, "ShanreW.csv")
print(resultsW)

##Assumption 1 - Linearity

## Graphically, plotting the model residuals (the difference 
#between the observed value and the model-estimated value) vs 
#the predictor 

Plot.Model.F.Linearity<-plot(resid(sh_lme4W),alpha_table1W$diversity_shannon)

## Assumption 2 Homogeneity of Variance
#Regression models assume that variance of the residuals
#is equal across groups. 

#extracts the residuals and places them in a new column in our original data table
alpha_tableW$lme10<- residuals(sh_lme4W) 
alpha_tableW$baslme10 <-abs(alpha_table$lme10) #creates a new column with the absolute value of the residuals
alpha_tableW$lme102 <- alpha_tableW$baslme10^2 #squares the absolute values of the residuals to provide the more robust estimate
SH_levenW <- lm(lme102 ~ Individual, data=alpha_tableW) #ANOVA of the squared residuals
anova(SH_levenW) #displays the results

##visually
plot(sh_lme4W) #creates a fitted vs residual plot

##Assumption 3: The residuals of the model are normally distributed.
#QQ plots 
qqmath(sh_lme4W, id=0.05)

alpha2W<-alpha_tableW
alpha2W$Day<-factor(alpha2W$Day, levels=c("day0", "day1", "day2", "day3","day4","day5"))

################################# Alpha Diversity Plotting For WARM ###############################
P1<-ggplot(data = alpha2W, aes(x = Day, y = diversity_shannon, fill=Day)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6, color="black")+
  scale_fill_manual(values=c("grey27", "coral","indianred", "indianred4","firebrick2","red4"))+
  theme_USGS_box()+scale_x_discrete(labels=c("day0"="Day 0", "day1"="Day 1", "day2"="Day 2", "day3"="Day 3", "day4"="Day 4", "day5"="Day 5"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                     limits = c(0,4),sec.axis = dup_axis(name=NULL))+ylab("Shannon Diversity")

P2<-ggplot(data = alpha2W, aes(x = Day, y = PD, fill=Day)) +
  stat_boxplot(geom ='errorbar', width = 0.6) +
  geom_boxplot(width = 0.6, color="black")+
  scale_fill_manual(values=c("grey27", "coral","indianred", "indianred4","firebrick2","red4"))+
  theme_USGS_box()+scale_x_discrete(labels=c("day0"="Day 0", "day1"="Day 1", "day2"="Day 2", "day3"="Day 3", "day4"="Day 4", "day5"="Day 5"))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(0,18),sec.axis = dup_axis(name=NULL))+ylab("Faith's Phylogenetic Diversity")


#tiff('Figure 1C_D.tiff', units="in", width=10, height=5, res=1200, compression = 'lzw')
plot_grid(p2, 
          p1,
          nrow = 2,rel_widths=c(1,1),rel_heights = c(1,1))
#dev.off()


######################## Statisical Analysis for Beta Diversity Warm ############################

##permanova weighted
SubDatW.asvtab <- as.data.frame(SubDatW@otu_table)
## pull tree
SubDatW.tree <- SubDatW@phy_tree

# re-root tree
new_treW <- ape::multi2di(SubDatW.tree)
phy_tree(SubDatW)<-new_treW


wunifrac_dist <-UniFrac(SubDatW, weighted=T)
adonis(wunifrac_dist ~ Temperature+Sex+Day, strata=sam.metaW$Individual,data = sam.metaW)
#             Df S  umsOfSqs  MeanSqs  F.Model      R2      Pr(>F)    
#  Temperature  1   0.05041  0.050406  1.5109    0.02485  0.001 ***
#  Sex          1   0.18678  0.186780  5.5989     0.09207  0.011 *  
#  Day          4   0.02342  0.005855  0.1755     0.01154  0.729    
#Residuals   53   1.76810    0.033360   0.87154           
#Total       59   2.02870                  1.00000  
      

##hommogeneity
permutest(betadisper(wunifrac_dist, sam.metaW$Temperature), strata=Individual)
#Response: Distances
#Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.000105 0.00010465 0.0481    999  0.823
#Residuals 58 0.126151 0.00217501 

permutest(betadisper(wunifrac_dist, sam.metaW$Sex), strata=Individual)
#Response: Distances
#Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)   
#Groups     1 0.026117 0.0261173 10.086    999  0.002 **
#  Residuals 58 0.150181 0.0025893
##permanova unweighted

uunifrac_dist <-UniFrac(SubDatW, weighted=F)
adonis(uunifrac_dist ~ Temperature+Sex+Day, strata=sam.metaW$Individual,data = sam.metaW)
#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Temperature  1    0.1938 0.19383 1.23672 0.02088  0.014 *
#Sex          1    0.4112 0.41117 2.62350 0.04429  0.291  
#Day          4    0.3716 0.09290 0.59278 0.04003  0.895  
#Residuals   53    8.3064 0.15672         0.89480         
#Total       59    9.2830                 1.00000         




##hommogeneity
permutest(betadisper(uunifrac_dist, sam.metaW$Temperature), strata=Individual)
#Response: Distances
##          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     1 0.007963 0.0079627 2.1455    999  0.155
#Residuals 58 0.215260 0.0037114   


##### Figure 3C&D Beta: PCoA plot using the unweighted UniFrac as distance for Individuals with Sex ####

## pull ASV table
SubDatW.asvtab <- as.data.frame(SubDatW@otu_table)
## pull tree
SubDatW.tree <- SubDatW@phy_tree

# re-root tree
new_treW <- ape::multi2di(SubDatW.tree)
phy_tree(SubDatW)<-new_treW


## weighted Sex and Individual
wunifrac_dist = phyloseq::distance(SubDatW, method="wunifrac")

wordination = ordinate(SubDatW, method="PCoA", distance=wunifrac_dist)
p1<-phyloseq::plot_ordination(SubDatW, uordination, color="Individual", shape="Sex", axes=c(1,2))+
  geom_point(size=3)+
  theme(legend.title=element_text(family="Times New Roman"),
        legend.text=element_text(family="Times New Roman"),
        legend.key = element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(),
        axis.title.y=element_text(size=12, family="Times New Roman"),
        axis.title.x=element_text(size=12, family="Times New Roman"),
        axis.text.y=element_text(size=10, family="Times New Roman", color="black"),
        axis.text.x=element_text(size=10, family="Times New Roman", color="black"))+
  scale_colour_manual(name = "Day",
                      labels = c("Day 0", "Day 1",
                                 "Day 2", "Day 3",
                                 "Day 4","Day 5"),
                      values = c("grey27","deepskyblue",
                                 "blue","slateblue" ,
                                 "purple", "darkorchid4"))




### same thing for sex and day unweighted
uunifrac_dist = phyloseq::distance(SubDatW, method="unifrac")
uordination = ordinate(SubDatW, method="PCoA", distance=uunifrac_dist)



p2<-phyloseq::plot_ordination(SubDatW, uordination, color="Individual", shape="Sex", axes=c(1,2))+
  geom_point(size=3)+
  theme(legend.title=element_text(family="Times New Roman"),
        legend.text=element_text(family="Times New Roman"),
        legend.key = element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(),
        axis.title.y=element_text(size=12, family="Times New Roman"),
        axis.title.x=element_text(size=12, family="Times New Roman"),
        axis.text.y=element_text(size=10, family="Times New Roman", color="black"),
        axis.text.x=element_text(size=10, family="Times New Roman", color="black"))+
  scale_colour_manual(name = "Day",
                      labels = c("Day 0", "Day 1",
                                 "Day 2", "Day 3",
                                 "Day 4","Day 5"),
                      values = c("grey27","deepskyblue",
                                 "blue","slateblue" ,
                                 "purple", "darkorchid4"))


## group plots
#tiff('Figure 3C_D.tiff', units="in", width=10, height=5, res=1200, compression = 'lzw')
plot_grid(p1, 
          p2,
          nrow = 1,rel_widths=c(1,1.3),
          rel_heights = c(1,1),
          labels = c("A: Weighted", "B: Unweighted"),
          hjust=-.9, vjust=2.5)
#dev.off()





#### Figures  2 ####

###### plotting #####
sam.metaC <- meta(SubDatC)
sam.metaC$SampleID <- rownames(sam.metaC)
sam.newC<-sample_data(sam.metaC)


sample_data(SubDatC)<-sam.newC
phyb.rar.tree2 <- SubDatC@phy_tree
new_treC <- ape::multi2di(phyb.rar.tree2)
phyb.rar3<-SubDatC
phy_tree(phyb.rar3)<-new_treC

library(reshape2)
p = phyb.rar3
m = "wunifrac"
s = "SampleID"
d = "Day"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sam.metaC %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")
wu.sd$Type2<-factor(wu.sd$Type2,levels=c("day0", "day1","day2","day3","day4","day5"))
wu.sd$Type1<-factor(wu.sd$Type1,levels=c("day0", "day1","day2","day3","day4","day5"))

dfgtC<-wu.sd

dfgtC$Compared <- apply(dfgtC, 1, function(x) grepl(x['Type1'], x['Type2'],fixed=TRUE))
dfgtC$Compared <- as.character(dfgtC$Compared)
dfgtC$Compared[dfgtC$Compared == "TRUE"] <- "Within"
dfgtC$Compared[dfgtC$Compared == "FALSE"] <- "Between"


ColdWU<-summarySE(dfgtC, measurevar = "value", groupvars =c("Type1", "Compared"))
WithinCold<-subset(ColdWU, Compared=="Within")
WithinCold$Treatment<-"Cold"


##warm
sam.metaW <- meta(SubDatW)
sam.metaW$SampleID <- rownames(sam.metaW)
sam.newW<-sample_data(sam.metaW)
sample_data(SubDatW)<-sam.newW
phyb.rar.tree2 <- SubDatW@phy_tree
new_treW <- ape::multi2di(phyb.rar.tree2)
phyb.rar3<-SubDatW
phy_tree(phyb.rar3)<-new_treW


p = phyb.rar3
m = "wunifrac"
s = "SampleID"
d = "Day"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sam.metaW %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")
wu.sd$Type2<-factor(wu.sd$Type2,levels=c("day0", "day1","day2","day3","day4","day5"))
wu.sd$Type1<-factor(wu.sd$Type1,levels=c("day0", "day1","day2","day3","day4","day5"))

dfgt<-wu.sd

dfgt$Compared <- apply(dfgt, 1, function(x) grepl(x['Type1'], x['Type2'],fixed=TRUE))
dfgt$Compared <- as.character(dfgt$Compared)
dfgt$Compared[dfgt$Compared == "TRUE"] <- "Within"
dfgt$Compared[dfgt$Compared == "FALSE"] <- "Between"


WarmWU<-summarySE(dfgt, measurevar = "value", groupvars =c("Type1", "Compared"))
WithinWarm<-subset(WarmWU, Compared=="Within")
WithinWarm$Treatment<-"Warm"



Combined<-rbind(WithinWarm,WithinCold)

## plotting the within difference in hot and cold
p1<-ggplot(Combined, aes(x=Type1, y=value, group=Treatment, color=Treatment)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(0.05))+
  labs(x="", y = "")+
  scale_color_manual(values=c('cyan3','coral1'))+
  facet_wrap(~Treatment)+
  theme(legend.position="none",strip.background = element_rect(colour="black",fill="white"),
        legend.title=element_blank(),strip.text=element_text(size=12, family="Times New Roman"),
        legend.text=element_text(family="Times New Roman"),legend.key = element_blank(),
        panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(), axis.title.y=element_text(size=12, family="Times New Roman"),
        axis.title.x=element_text(size=12, family="Times New Roman"),
        axis.text.y=element_text(size=10, family="Times New Roman", color="black"),
        axis.text.x=element_text(size=10, family="Times New Roman", color="black", angle=45, vjust=0.5))+
        scale_x_discrete(labels=c("day0"="Day 0", "day1"="Day 1", "day2"="Day 2",
                           "day3"="Day 3", "day4"="Day 4", "day5"="Day 5"))+
  scale_y_continuous(breaks=c(0.05,0.1,.15,0.2,0.25,0.3, 0.35,0.4), limits=c(0.03,0.4))

p1


#### diferences in females vs males  from Warm Results ###
SubDatF<-subset_samples(SubDatW, Sex==c("Female"))
sam.metaF <- meta(SubDatF)
sam.metaF$SampleID <- rownames(sam.metaF)
sam.new<-sample_data(sam.metaF)
sample_data(SubDatF)<-sam.newF
phyb.rar.tree2 <- SubDatF@phy_tree
new_treF <- ape::multi2di(phyb.rar.tree2)
phyb.rar3<-SubDatF
phy_tree(phyb.rar3)<-new_treF


p = phyb.rar3
m = "wunifrac"
s = "SampleID"
d = "Day"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sam.metaF %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")
wu.sd$Type2<-factor(wu.sd$Type2,levels=c("day0", "day1","day2","day3","day4","day5"))
wu.sd$Type1<-factor(wu.sd$Type1,levels=c("day0", "day1","day2","day3","day4","day5"))

dfgtC<-wu.sd

dfgtC$Compared <- apply(dfgtC, 1, function(x) grepl(x['Type1'], x['Type2'],fixed=TRUE))
dfgtC$Compared <- as.character(dfgtC$Compared)
dfgtC$Compared[dfgtC$Compared == "TRUE"] <- "Within"
dfgtC$Compared[dfgtC$Compared == "FALSE"] <- "Between"


WarmFWU<-summarySE(dfgtC, measurevar = "value", groupvars =c("Type1", "Compared"))
WithinF<-subset(WarmFWU, Compared=="Within")
WithinF$Treatment<-"Female"


##Male
SubDatM<-subset_samples(SubDatW, Sex==c("Male"))
sam.metaM <- meta(SubDatM)
sam.metaM$SampleID <- rownames(sam.metaM)
sam.newM<-sample_data(sam.metaM)
sample_data(SubDatM)<-sam.newM
phyb.rar.tree2 <- SubDatM@phy_tree
new_treM <- ape::multi2di(phyb.rar.tree2)
phyb.rar3<-SubDatM
phy_tree(phyb.rar3)<-new_treM


p = phyb.rar3
m = "wunifrac"
s = "SampleID"
d = "Day"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sam.metaM %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")
wu.sd$Type2<-factor(wu.sd$Type2,levels=c("day0", "day1","day2","day3","day4","day5"))
wu.sd$Type1<-factor(wu.sd$Type1,levels=c("day0", "day1","day2","day3","day4","day5"))

dfgt<-wu.sd

dfgt$Compared <- apply(dfgt, 1, function(x) grepl(x['Type1'], x['Type2'],fixed=TRUE))
dfgt$Compared <- as.character(dfgt$Compared)
dfgt$Compared[dfgt$Compared == "TRUE"] <- "Within"
dfgt$Compared[dfgt$Compared == "FALSE"] <- "Between"


MWU<-summarySE(dfgt, measurevar = "value", groupvars =c("Type1", "Compared"))
WithinM<-subset(MWU, Compared=="Within")
WithinM$Treatment<-"Male"
Combined<-rbind(WithinM,WithinF)

## ggplot
p2<-ggplot(Combined, aes(x=Type1, y=value, group=Treatment, color=Treatment)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(0.05))+
  labs(x="", y = "")+
  scale_color_manual(values=c('gray45','gray18'))+ facet_wrap(~Treatment)+
  theme(legend.position="none",strip.background = element_rect(colour="black",fill="white"),
        strip.text=element_text(size=12, family="Times New Roman"),
        legend.title=element_blank(),legend.text=element_text(family="Times New Roman"),
        legend.key = element_blank(),panel.background=element_rect(color="black", fill="white"),
        panel.grid=element_blank(), axis.title.y=element_text(size=12, family="Times New Roman"),
        axis.title.x=element_text(size=12, family="Times New Roman"),
        axis.text.y=element_text(size=10, family="Times New Roman", color="black"),
        axis.text.x=element_text(size=10, family="Times New Roman", color="black", angle=45, vjust=.5))+
    scale_x_discrete(labels=c("day0"="Day 0", "day1"="Day 1", "day2"="Day 2",
                           "day3"="Day 3", "day4"="Day 4", "day5"="Day 5"))+
      scale_y_continuous(breaks=c(0.05,0.1,.15,0.2,0.25,0.3, 0.35,0.4,0.45), limits=c(0.03,0.45))
p2

#tiff('Figure 2.tiff', units="in", width=8, height=8, res=1200, compression = 'lzw')
plot_grid(p1, 
          p2,
          nrow = 2,rel_widths=c(1,1),rel_heights = c(1,1), hjust=-.9, vjust=2.5)
#dev.off()



#################################### POWER ######################################################

### POWER ANALYSIS MINK cold vs warm ####

library(pwr)
library(stats)
library(effectsize)

set.seed(1234)
setwd("~/Desktop/Projects/Mink/Longitudinal/POWER")

### T-Test for PD

## randomly subsample 10 from cold vs 10 warm
set.seed(1234)
#PDC<-read.csv("Cold_sub.csv")
#PDW<-read.csv("Warm_Sub.csv")

##sample

csamp<-PDC[sample(nrow(PDC),10),]
Wsamp<-PDW[sample(nrow(PDW),10),]



##remove shan diversity
keep <- c("PD")
Cpd2<-csamp[ , keep, drop=F]
Wpd2<-Wsamp[ , keep, drop=F]

##means
mc<-mean(Cpd2$PD)
mw<-mean(Wpd2$PD)

##pooled standard deviation
CSD<-sd(Cpd2$PD)
#2.86
WSD<-sd(Wpd2$PD)
#2.66


##combined
Cpd2$Type<-"Cold"
Wpd2$Type<-"Warm"

pdcom<-as.data.frame(rbind(Cpd2, Wpd2))


psd<-sd_pooled(PD~Type, data=pdcom)
##2.8
##check this is right

##cohen's d:
###what are ecologically meaningful differences?

(mc-mw)/psd
#.7

1.5/psd
#0.54
2.5/psd
#.90
3.5/psd
#1.26
4.5/psd
#1.63


ptab2<-cbind(NULL, NULL) 

for (i in c(0.54,.90,1.26,1.63)){
  pwrt2<-power.t.test(d=i,sd=psd,power=.8,sig.level=.05,type="two.sample",alternative="two.sided")
  ptab2<-as.data.frame(rbind(ptab2, cbind(pwrt2$d, pwrt2$n)))
}
ptab2

ptab2<- ptab2 %>% 
  mutate_if(is.numeric, round, digits = 0)
ptab2$Effect<-c(0.54,.90,1.26,1.63)
ptab2$Size<-ptab2$V2*2


### plot
labels<-c("n=826", "n=298", "n=154","n=92")
#tiff('PDPower.tiff', units="in", width=6.5, height=5, res=1200, compression = 'lzw')
p1<-ggplot(ptab2, aes(x=Effect, y=Size))+
  geom_line()+geom_point(size=2,shape=19 )+
  annotate(geom="text", x=1.3, y=360, label=("Faith's PD"), color="black", family="Times New Roman", size=3)+
  geom_text(label=labels, nudge_x =- 0.05, nudge_y = -6, family="Times New Roman", size=2.5)+
  ylab("Sample size")+xlab("Effect size\n(Difference in means)")+
  theme(text=element_text(family="Times New Roman"),legend.position="none",panel.grid=element_blank(),panel.background=element_rect(fill="white", color="black"),legend.title=element_blank(), legend.text=element_blank(),strip.text.x=element_text(color="black", family="Times New Roman", size=12),strip.background=element_rect(color="black", fill="white"),
        axis.text.y=element_text(size=10,family="Times New Roman", color="black"),axis.title.y=element_text(family="Times New Roman"),axis.text.x=element_text( hjust=.7,family="Times New Roman", color="black", size=9)) + 
  ggtitle(paste0(""))+scale_x_continuous( breaks = c(0.54,.90,1.26,1.63),  labels = c("54%\n(1.5 difference)","90%\n(2.5 difference)","126%\n(3.5 difference)", "163%\n(4.5 difference)"))+
  scale_y_continuous( breaks = c(100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850))
#dev.off()


### test
library(car)
library(exactRankTests)
library(RVAideMemoire)
library(plyr)

##test assumptions
byf.shapiro(PD~Type, data=pdcom)

##passed
leveneTest(PD~Type, data=pdcom)
##passed

t.test(Shan~Type, paired=F, var.equal=T,type="two.sample", data=shcom)
power.t.test(d=,sd=psds,power=.8,sig.level=.96,n=20,type="two.sample",alternative="two.sided")



#####Shannon
##remove shan diversity
keep <- c("Shan")
Csh2<-csamp[ , keep, drop=F]
Wsh2<-Wsamp[ , keep, drop=F]




## power 


##means
mcs<-mean(Csh2$Shan)
mws<-mean(Wsh2$Shan)
mcs-mws
#0.02
##pooled standard deviation
CSDs<-sd(Csh2$Shan)
#2.86
WSDs<-sd(Wsh2$Shan)
#2.66

##combined
Csh2$Type<-"Cold"
Wsh2$Type<-"Warm"

shcom<-as.data.frame(rbind(Csh2, Wsh2))

psds<-sd_pooled(Shan~Type, data=shcom)
##.99

0.02/.99
#
.5/psds
#0.5

1/psds
#1

1.5/psds
#1.52
2/psds
#2.0

ptabs<-cbind(NULL, NULL) 

for (i in c(.5, 1, 1.52, 2)){
  pwrts<-power.t.test(d=i,sd=psd,power=.8,sig.level=.05,type="two.sample",alternative="two.sided")
  ptabs<-as.data.frame(rbind(ptabs, cbind(pwrts$d, pwrts$n)))
}
ptabs

ptabs<- ptabs %>% 
  mutate_if(is.numeric, round, digits = 0)
ptabs$Effect<-c(.5, 1, 1.52, 2)
ptabs$Size<-ptabs$V2*2


### plot
labels<-c("n=964", "n=242", "n=106", "n=62")
#tiff('SHPOWER.tiff', units="in", width=6.5, height=5, res=1200, compression = 'lzw')
p2<-ggplot(ptabs, aes(x=Effect, y=Size))+
  geom_line()+geom_point(size=2,shape=19 )+
  annotate(geom="text", x=1.3, y=360, label=("Shannon Diversity"), color="black", family="Times New Roman", size=3)+
  geom_text(label=labels,  nudge_x =-.03, nudge_y =-30, family="Times New Roman", size=2.5)+
  ylab("")+xlab("Effect size\n(Difference in means)")+
  theme(text=element_text(family="Times New Roman"),legend.position="none",
        panel.grid=element_blank(),
        panel.background=element_rect(fill="white", color="black"),
        legend.title=element_blank(), legend.text=element_blank(),
        strip.text.x=element_text(color="black", family="Times New Roman", size=12),
        strip.background=element_rect(color="black", fill="white"),
        axis.text.y=element_text(size=10,family="Times New Roman", color="black"),
        axis.title.y=element_text(family="Times New Roman"),
        axis.text.x=element_text( family="Times New Roman", color="black", size=9)) + 
  ggtitle(paste0(""))+scale_x_continuous(breaks = c(0.5,1,1.52,2),labels = c("50%\n(0.5 difference)","100%\n(1 difference)","152%\n(1.5 difference)", "200%\n(2 difference)"))+
  scale_y_continuous(breaks = c(50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000))
#dev.off()



#tiff('Figure5.tiff', units="in", width=12, height=4, res=1200, compression = 'lzw')
plot_grid(p1,p2, 
          nrow = 1,rel_widths=c(1,1.1), labels = c("A", "B"))
#dev.off()



########################################### LEfSe #############################################
## get table ready for LEfSe
## conda activate qiime environment
## collapse to genus level
qiime taxa collapse \
--i-table clean_table.qza \
--o-collapsed-table collapse.table.qza \
--p-level 6 \
--i-taxonomy taxonomy_final.qza

## convert relative frequency table to bion & convert to txt
qiime feature-table relative-frequency \
--i-table collapse.table.qza \
--o-relative-frequency-table collapse.frequency.table.qza \
--output-dir exported/
  
  ## export to biom file
  qiime tools export \
--input-path collapse.frequency.table.qza \
--output-path exported/ 
  
  ## convert to txt/tsv with taxonomy
biom convert -i exported/feature-table.biom -o table.from_biom_w_taxonomy.txt --to-tsv --header-key taxonomy

## set up for lefse
#"|" instead of ";" subclass, main thing, Subject, no otu ID
## you've probably already done that to run it in the galaxy version
## before so you should know how to set it up
## FYI, after I was done I named my file: "MartenLtoG.txt" so you can follow along



#### LEfSe in python ####
## change environment
## to run LEfSe/Graphlan you will need to create a python 2 environment!
## activate python 2 environment
## you will need to install lefse and graphlan
## graphlan
conda install -c biobakery graphlan

## LEfSe
## https://github.com/biobakery/biobakery/wiki/lefse
conda install -c biobakery lefse


## run lefse, this will run just like if you were running it in the galaxy webpage
## you can even export the .res file from the webpage,
## but personally I feel this is much easier
format_input.py MartenLtoG.txt Sexm.in -c 1 -s 2 -u 3 -o 1000000

##Sexm.in is the outputfile you need to run the lefse.
run_lefse.py Sexm.in Sexm.res
###Sexm.res is the outputfile you need to run the grphlan.
## its the results of the differentially enriched taxa


# convert it with graphlan
## you can see what all of the different commands mean and
## how to change them in the top link
export2graphlan.py -i reduced.txt -o Sexm2.res -t tree.txt -a annot.txt --title "" --annotations 2,3 --external_annotations 4,5,6 --fname_row 0 --skip_rows 1,2 --ftop 200
# attach annotation to the tree
## the "annot.txt" is the textfile you will edit to make customized cladogram
## I like to run it as is first to note the things I want to change.

graphlan_annotate.py --annot annot.txt tree.txt outtree.txt
# generate the beautiful image
graphlan.py --dpi 300 --size 7.0 outtree.txt outimg.png --external_legends

###after initial lefse, I make a new clad table with only the enriched taxa and I do the same for the .res file!!!
## you can edit the ".res" file by opening it with a text editor.. maybe even excel!
## tinker around and you will get a feel
## makes the cladogram much less clustered!
## you can basically change whatever you want, node size, color of each group/each individual enriched taxa
## Rotation how many taxa you want labelled "A", "B","C", etc
## look through their webpage and you will see there is a sh*t ton you can do
##thats the reduced.txt & sexm2.res



