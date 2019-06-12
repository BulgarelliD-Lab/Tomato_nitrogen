#############################################################
#
# Ref to the ARTICLE
# 
# Federica Caradonia, Rodrigo Alegria Terrazas and Davide Bulgarelli
#
# federica.caradonia@unimore.it
# r.z.alegriaterrazas@dundee.ac.uk
# d.bulgarelli@dundee.ac.uk
#
# Revison June 2019
# 
# script to reproduce calculations and figures presented in the manuscript
# 
# Disclaimer: this is a temporary file, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#1st time installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#biocLite("DESeq2")
#biocLite("PMCMR")
#biocLite("VennDiagram")
#biocLite("UpSetR")
#biocLite("Tax4Fun")
#biocLite("qiimer")


#required packages 
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("PMCMR")
library("VennDiagram")
library("UpSetR")
library("plyr")
library("Tax4Fun")
#library("qiimer")
#library("biom")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set the working directory
#Davide cpu: remember from Federica's cpu is differnt
setwd("C:/Users/DB42008/Box Sync/Davide_lab_experiments/JH09/JH09_FC_R_analysis/")

#############################################################
#import the count matrix and the desing file
#############################################################


#OTU table generated using QIIME 1.9.0. 
#This file is designated Worksheet_ws2 in Supplementary Information
dat_info <- read.delim("JH09_FC_otu_table_nc2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the file 
dim(dat_info)
colnames(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Hv identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:86]))
OTU_97_reads

#total reads
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:86]))
OTU_97_reads_sum

#design file
#This file is designated Worksheet_ws1 in Supplementary Information
design <- read.delim("Map_JH09_FC_phyloseq.txt", sep = "\t", header=TRUE, row.names=1)
design

#remove chloroplast and mitochondria OTUs from the original dataset
#inspect the fist rows of the table to get insights into the column ConsensusLineage
dat_info[1:5, 87]

Chloroplast <- dat_info[grepl("Chloroplast", dat_info$ConsensusLineage), ]
dim(Chloroplast)

mitochondria <- dat_info[grepl("mitochondria", dat_info$ConsensusLineage), ]
dim(mitochondria)

#set a difference between the row names of the the three datasets: this information will be used to filter out Plant derived OTUs from the OTU table
noPlants <- setdiff(rownames(dat_info), c(rownames(Chloroplast), rownames(mitochondria)))

#inspect the results
length(rownames(dat_info))
length(noPlants)

#save the OTUids list generated at line 81. This will be used to usbset the OTU table and generate taxa-tables in QIIME
#write(noPlants, "JH09_FC_noPlant_OTUs_id.txt")
#The above mentioned file is designated Worksheet_ws3 in Supplementary Information

#generate a new OTU table which will be devoid of Chloroplast and Mitochondria OTUs
dat_info_noPlants <- dat_info[noPlants, ]

#create a new count matrix without OTUs assigned to Choloplast and Mitochondria
dat_count <- dat_info[, rownames(design)]
dat_count_noplants <- dat_info[noPlants, rownames(design)]
dim(dat_count_noplants)

#and a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_count_noplants), 87])
rownames(dat_tax_noPlants) <- rownames(dat_count_noplants)
#save the above file and in excel we will create a new tax table where each column represents a taxonomic rank
#write.table(dat_tax_noPlants, file="JH09_FC_dat_tax_noPlants.txt", sep="\t")

#check the effect of mitochondria/chloroplast depletion on the new OTU table
#with plant sequences
dim(dat_count)

#w/o plant-derived sequences
dim(dat_count_noplants)

#total number of reads w/o plant sequences
OTU_97_reads_noPlants <- colSums(dat_count_noplants)
OTU_97_reads_noPlants

#now sort per sample
sort(OTU_97_reads_noPlants)

#total number of reads
OTU_97_reads_noPlants_sum <- sum(OTU_97_reads_noPlants)
OTU_97_reads_noPlants_sum 

#Define the proportion of non-plant reads in the original dataset 
useful_reads <- (OTU_97_reads_noPlants_sum/OTU_97_reads_sum)*100
useful_reads

#create a dataset to visualise the proportion of reads per sample before and after removing OTUs assigned to chloroplast and mitochondria
OTU_97_reads_noPlants <- as.data.frame(OTU_97_reads_noPlants)
OTU_97_microbial_reads_proportion <- as.data.frame(colSums(dat_count_noplants)/colSums(dat_count))*100

#rename the columns in the generated datasets
colnames(OTU_97_reads_noPlants) <- c("reads")
colnames(OTU_97_microbial_reads_proportion) <- c("Microbial_OTUs_reads")

#combine these datasets with the design file
design_info <- cbind(design, OTU_97_reads_noPlants)
design_info_2 <- cbind(design_info, OTU_97_microbial_reads_proportion)

#calculate the max, min and mean number of reads for the dataset
mean(design_info_2$reads)
max(design_info_2$reads)
min(design_info_2$reads)

#############################################################
#Genererate the phyloseq object
#Data required: dat_count_noplants; design, JH09_FC_dat_tax_noPlants_ordered.txt, and 97_otus.tree.gz
#############################################################

#The OTU Table counts
JH09_FC_OTU <- otu_table(dat_count_noplants, taxa_are_rows=TRUE)

#The taxonomy information
#Note that the file JH02_JH03_dat_tax_noPlants_ordered.txt has been generated from the output of lines 96-100  
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
##The above mentioned file is designated Worksheet_ws4 in Supplementary Information
JH09_FC_taxa_ordered <- read.delim ("JH09_FC_dat_tax_noPlants_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
JH09_FC_taxa <- tax_table(as.matrix(JH09_FC_taxa_ordered))
dim(JH09_FC_taxa)

#The mapping file 
JH09_FC_map <- sample_data(design)

#The phylogenetic tree: the OTU table has been generated using a closed reference approach agains the greengenes 13_05 database use the corresponding phylogenetic tree 
JH09_FC_tree <- read_tree_greengenes("97_otus.tree.gz")

#check whether the tree is rooted
is.rooted(JH09_FC_tree)

#merge the files and create the phyloseq object
JH09_FC_data_phyloseq <- merge_phyloseq(JH09_FC_OTU, JH09_FC_taxa, JH09_FC_map,  JH09_FC_tree)
JH09_FC_data_phyloseq

#inspect the generated data
JH09_FC_data_phyloseq
sum(colSums(otu_table(JH09_FC_data_phyloseq)))
dim(dat_count_noplants)
sum(colSums(dat_count_noplants))

#remove the  samples with less than 1,000 HQ reads from the dataset
JH09_FC_data_phyloseq_2 <- subset_samples(JH09_FC_data_phyloseq, Microhabitat != "Out")
design_2 <- design[colnames(otu_table(JH09_FC_data_phyloseq_2)), ]
JH09_FC_data_phyloseq_2 

#H09_FC_data_phyloseq_2 <- subset_samples(JH09_FC_data_phyloseq_2, Treatments == "nothing")

#abundance filtering: remove OTUs tallyig less than 10 reads in at least 5% of the samples (=1 treatment)
JH09_FC_data_phyloseq_3 = filter_taxa(JH09_FC_data_phyloseq_2, function(x) sum(x > 10) > (0.05*length(x)), TRUE)
JH09_FC_data_phyloseq_3

#inspect the generated data
JH09_FC_data_phyloseq_2
sum(colSums(otu_table(JH09_FC_data_phyloseq_2)))
JH09_FC_data_phyloseq_3
sum(colSums(otu_table(JH09_FC_data_phyloseq_3)))

#proportion of retained OTUs
(length(rownames(otu_table(JH09_FC_data_phyloseq_3))))/(length(rownames(otu_table(JH09_FC_data_phyloseq_2)))) * 100
#proportion of retained reads
(sum(colSums(otu_table(JH09_FC_data_phyloseq_3))))/(sum(colSums(otu_table(JH09_FC_data_phyloseq_2)))) *100

#number of reads per sample
sort(colSums(otu_table(JH09_FC_data_phyloseq_3)))

#distribution
min(colSums(otu_table(JH09_FC_data_phyloseq_3)))
max(colSums(otu_table(JH09_FC_data_phyloseq_3)))
mean(colSums(otu_table(JH09_FC_data_phyloseq_3)))

#############################################################
#Figure 2: Alphadiversity calculations
#Data required: design_2; JH09_FC_data_phyloseq_rare_table_counts_2.txt
#############################################################

#rarefy the dataset
#JH09_FC_data_phyloseq_rare <- rarefy_even_depth(JH09_FC_data_phyloseq_3, rngseed=TRUE)

#extract and save the OTU table for reproducibility of the code
#JH09_FC_data_phyloseq_rare_table <- as.data.frame(otu_table(JH09_FC_data_phyloseq_rare))

#inspect the generated file
#class(JH09_FC_data_phyloseq_rare_table)
#dim(JH09_FC_data_phyloseq_rare_table)

#save the file for the reproducibility of the code
#write.table(JH09_FC_data_phyloseq_rare_table, file="JH09_FC_data_phyloseq_rare_table_counts.txt", sep="\t")
##The above mentioned file is designated Worksheet_ws5 in Supplementary Information

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_rare <- read.delim("JH09_FC_data_phyloseq_rare_table_counts_2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the generated file
dim(dat_count_rare)
colSums(dat_count_rare)

#generate a new phyloseq object wich will contain only the rarefied counts and the design file (only these two pieces of information are required for alphadiversity calculation)
JH09_FC_OTU_rare <- otu_table(dat_count_rare, taxa_are_rows=TRUE)
JH09_FC_map_rare <- sample_data(design_2)
JH09_FC_data_rare_phyloseq <- merge_phyloseq(JH09_FC_OTU_rare, JH09_FC_map_rare)

#Inspect the generated file
JH09_FC_data_rare_phyloseq
sample_sums(JH09_FC_data_rare_phyloseq)

#Index calculations
JH09_FC_alpha_rare <-  estimate_richness(JH09_FC_data_rare_phyloseq, measures = c("Observed", "Shannon", "Chao1"))
JH09_FC_alpha_rare

#generate a new dataframes for data visualisation

#Sample information

#Microhabitat
design_microhabitat <- as.data.frame(design_2[, 1])
rownames(design_microhabitat) <- rownames(design_2)
colnames(design_microhabitat) <- c("Microhabitat")

#Treatment
design_treatment <- as.data.frame(design_2[, 2])
rownames(design_treatment) <- rownames(design_2)
colnames(design_treatment) <- c("Treatment")

#data frame Genotype_Description
design_MT <- cbind(design_microhabitat, design_treatment)

#remove nursery sample
design_MT <- design_MT[colnames(dat_count_rare), ]

#Observed OTUs
JH09_FC_alpha_rare_Observed <- as.data.frame(JH09_FC_alpha_rare[ ,1])
rownames(JH09_FC_alpha_rare_Observed) <- rownames(JH09_FC_alpha_rare)
colnames(JH09_FC_alpha_rare_Observed) <- c("Observed")

#Combine the dataset sample description and Observed OTUs
JH09_FC_alpha_rare_Observed_MT <- cbind(design_MT, JH09_FC_alpha_rare_Observed)
JH09_FC_alpha_rare_Observed_MT <- as.data.frame(JH09_FC_alpha_rare_Observed_MT)
JH09_FC_alpha_rare_Observed_MT$Treatment
JH09_FC_alpha_rare_Observed_MT$Microhabitat

#Order the levels according to a defined order
#MTreatment
JH09_FC_alpha_rare_Observed_MT$Treatment <- ordered(JH09_FC_alpha_rare_Observed_MT$Treatment, levels=c("nursery", "nothing", "Liquiddigestate", "SRLiquiddigestate", "Pellet", "SCAMproduct", "Mineralfertilizer", "SRMineralfertilizer")) 

#plotting
#data visualisation: rhizosphere samples only
#remove vacant columns (bulk, OUT from the plot) Use droplevels to remove the empty levels from the list of levels
JH09_FC_alpha_rare_Observed_MT$Microhabitat <- droplevels(JH09_FC_alpha_rare_Observed_MT$Microhabitat)
boxplot(Observed ~ Treatment * Microhabitat, data=JH09_FC_alpha_rare_Observed_MT,  xlab = "Samples", ylab = "number of OTUs", main = "OTU richness", col = (c("gold", "brown", "brown", "brown", "brown", "brown", "brown","brown", "gold", "darkgreen", "darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen" )))

#remove nursery samples
JH09_FC_alpha_rare_Observed_MT_soil <- JH09_FC_alpha_rare_Observed_MT[JH09_FC_alpha_rare_Observed_MT$Treatment != "nursery",]
JH09_FC_alpha_rare_Observed_MT_soil
#assess the normality of distribution
shapiro.test(JH09_FC_alpha_rare_Observed_MT_soil$Observed)
#Microhabitat effect
Observed_OTUs_stats <- wilcox.test(Observed ~ Microhabitat, data = JH09_FC_alpha_rare_Observed_MT_soil,
                   exact = FALSE)
Observed_OTUs_stats 

#Chao1 OTUs
JH09_FC_alpha_rare_Chao1 <- as.data.frame(JH09_FC_alpha_rare[ ,2])
rownames(JH09_FC_alpha_rare_Chao1) <- rownames(JH09_FC_alpha_rare)
colnames(JH09_FC_alpha_rare_Chao1) <- c("Chao1")

#Combine the dataset sample description and Chao1 OTUs
JH09_FC_alpha_rare_Chao1_MT <- cbind(design_MT, JH09_FC_alpha_rare_Chao1)
JH09_FC_alpha_rare_Chao1_MT$Treatment

#Order the levels according to a defined order
#MTreatment
JH09_FC_alpha_rare_Chao1_MT$Treatment <- ordered(JH09_FC_alpha_rare_Chao1_MT$Treatment, levels=c("nursery", "nothing", "Liquiddigestate", "SRLiquiddigestate", "Pellet", "SCAMproduct", "Mineralfertilizer", "SRMineralfertilizer")) 

#plotting
#data visualisation: rhizosphere samples only
#remove vacant columns (bulk, OUT from the plot) Use droplevels to remove the empty levels from the list of levels
JH09_FC_alpha_rare_Chao1_MT$Microhabitat <- droplevels(JH09_FC_alpha_rare_Chao1_MT$Microhabitat)
boxplot(Chao1 ~ Treatment * Microhabitat, data=JH09_FC_alpha_rare_Chao1_MT,  xlab = "Samples", ylab = "number of OTUs", main = "OTU richness", col = (c("gold", "brown", "brown", "brown", "brown", "brown", "brown","brown", "gold", "darkgreen", "darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen" )))

#remove nursery samples
JH09_FC_alpha_rare_Chao1_MT_soil <- JH09_FC_alpha_rare_Chao1_MT[JH09_FC_alpha_rare_Chao1_MT$Treatment != "nursery",]
JH09_FC_alpha_rare_Chao1_MT_soil
#assess the normality of distribution
shapiro.test(JH09_FC_alpha_rare_Chao1_MT_soil$Chao1)
#Microhabitat effect
Chao1_OTUs_stats <- wilcox.test(Chao1 ~ Microhabitat, data = JH09_FC_alpha_rare_Chao1_MT_soil,
                                   exact = FALSE)
Chao1_OTUs_stats 

#Shannon OTUs
JH09_FC_alpha_rare_Shannon <- as.data.frame(JH09_FC_alpha_rare[ ,4])
rownames(JH09_FC_alpha_rare_Shannon) <- rownames(JH09_FC_alpha_rare)
colnames(JH09_FC_alpha_rare_Shannon) <- c("Shannon")

#Combine the dataset sample description and Shannon OTUs
JH09_FC_alpha_rare_Shannon_MT <- cbind(design_MT, JH09_FC_alpha_rare_Shannon)
JH09_FC_alpha_rare_Shannon_MT$Treatment

#Order the levels according to a defined order
#MTreatment
JH09_FC_alpha_rare_Shannon_MT$Treatment <- ordered(JH09_FC_alpha_rare_Shannon_MT$Treatment, levels=c("nursery", "nothing", "Liquiddigestate", "SRLiquiddigestate", "Pellet", "SCAMproduct", "Mineralfertilizer", "SRMineralfertilizer")) 

#plotting
#data visualisation: rhizosphere samples only
#remove vacant columns (bulk, OUT from the plot) Use droplevels to remove the empty levels from the list of levels
JH09_FC_alpha_rare_Shannon_MT$Microhabitat <- droplevels(JH09_FC_alpha_rare_Shannon_MT$Microhabitat)
boxplot(Shannon ~ Treatment * Microhabitat, data=JH09_FC_alpha_rare_Shannon_MT,  xlab = "Samples", ylab = "Shannon index", main = "Shannon", col = (c("gold", "brown", "brown", "brown", "brown", "brown", "brown","brown", "gold", "darkgreen", "darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen" )))

#remove nursery samples
JH09_FC_alpha_rare_Shannon_MT_soil <- JH09_FC_alpha_rare_Shannon_MT[JH09_FC_alpha_rare_Shannon_MT$Treatment != "nursery",]
JH09_FC_alpha_rare_Shannon_MT_soil
#assess the normality of distribution
shapiro.test(JH09_FC_alpha_rare_Shannon_MT_soil$Shannon)
#Microhabitat effect
Shannon_OTUs_stats <- wilcox.test(Shannon ~ Microhabitat, data = JH09_FC_alpha_rare_Shannon_MT_soil,
                                exact = FALSE)
Shannon_OTUs_stats 

#############################################################
#Figure 3: Betadiversity calculations
#Data required: design_2; JH02_JH03_RHM_data_phyloseq_3
#############################################################

#Transform the count in relative abundance cpm
JH09_FC_data_phyloseq_prop <- transform_sample_counts(JH09_FC_data_phyloseq_3,  function(x) 1e+03 * x/sum(x))

#PCoA bray distance
JH09_FC_data_phyloseq_prop_bray <- ordinate(JH09_FC_data_phyloseq_prop, "PCoA", "bray")
plot_ordination(JH09_FC_data_phyloseq, JH09_FC_data_phyloseq_prop_bray , color = "Treatments")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH09_FC_data_phyloseq_prop, JH09_FC_data_phyloseq_prop_bray , shape ="Microhabitat", color = "Treatments")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("green", "cyan", "brown", "gold", "darkorange", "purple", "darkgreen", "darkblue" ))
p + ggtitle("PCoA 16S data, Bray distance")

#Assess Microhabitat & Treatment effect (rhizosphere & root samples only)
#Subsetting: remove bulk samples
JH09_FC_data_phyloseq_prop_2 <- subset_samples(JH09_FC_data_phyloseq_prop, Microhabitat != "Bulk")
JH09_FC_data_phyloseq_prop_2
design_soil_grown <- design[colnames(otu_table(JH09_FC_data_phyloseq_prop_2)), ]

#BC distance adonis
BC <- phyloseq::distance(JH09_FC_data_phyloseq_prop_2, "bray")
adonis(BC ~ Treatments * Microhabitat, data= design_soil_grown, permutations = 5000)

#PCoA weighted unifrac distance
#info Unifrac: https://en.wikipedia.org/wiki/UniFrac
JH09_FC_data_phyloseq_prop_wunifrac <- ordinate(JH09_FC_data_phyloseq_prop, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH09_FC_data_phyloseq, JH09_FC_data_phyloseq_prop_wunifrac , color = "Treatments")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH09_FC_data_phyloseq_prop, JH09_FC_data_phyloseq_prop_wunifrac , shape ="Microhabitat", color = "Treatments")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = c("green", "cyan", "brown", "gold", "darkorange", "purple", "darkgreen", "darkblue" ))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#WU distance adonis
WU <- phyloseq::distance(JH09_FC_data_phyloseq_prop_2, "unifrac", weighted= TRUE)
adonis(WU ~ Treatments * Microhabitat, data= design_soil_grown, permutations = 5000)

#remove bulk soil sample and test for compartment effect
#WU distance adonis
WU <- phyloseq::distance(JH09_FC_data_phyloseq_prop_2, "unifrac", weighted= TRUE)
adonis(WU ~ Microhabitat *Timing, data= design_soil_grown, permutations = 5000)


############################################################################
#Hypothesis testing: create a new OTU table to be imported in qiime to test OTUs differentially regulated
#we use DEseq https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y
###########################################################################

#Microhabitat effect
#generate a new OTU table for calculation in QIIME
#extract OTU counts (non rarefied, yet filtered for abundance) from an abundance threshold phyloseq object
#the starting file is JH09_FC_data_phyloseq_3 and has to be pruned for nursery samples
JH09_FC_data_phyloseq_4 <- subset_samples(JH09_FC_data_phyloseq_3, Treatments != "nursery")
JH09_FC_data_phyloseq_4 <- subset_samples(JH09_FC_data_phyloseq_3, Treatments == "nothing")
#extract the count matrix
dat_count_threshold  <- as.data.frame(otu_table(JH09_FC_data_phyloseq_4))
#extract the cognate taxonomy information
dat_tax_threshold  <- as.data.frame(dat_tax_noPlants[rownames(dat_count_threshold  ), ])
colnames(dat_tax_threshold ) <- c("ConsensusLineage")
colnames(dat_tax_threshold)

#merge count and taxonomy files
dat_info_threshold <- cbind(dat_count_threshold , dat_tax_threshold )
dim(dat_info_threshold)
colnames(dat_info_threshold)
dat_info_threshold[1:5, ]

#generate a new design file for QIIME calculation
design_3 <- design[colnames(dat_count_threshold), ]
design_3

#save the files for the reproducibility of the code and the analysis in QIIME. Note: due to R formatting, remember to shift columns to the right and add the '#OTU ID' to the first column. Save the file as _2
#write.table(dat_info_threshold, file="JH09_FC_info_threshold.txt", sep = "\t")
#save mapping file for DESeq calculations. Note: due to R formatting, remember to shift columns to the right and add the '#SampleID' to the first column. Save the file as _2
#write.table(design_3, file="Map_JH09_FC_microhabitat.txt", sep = "\t") 

#Treatment effect: rhizosphere
#generate a new OTU table for calculation in QIIME
#extract OTU counts (non rarefied, yet filtered for abundance) from an abundance threshold phyloseq object
#the starting file is JH09_FC_data_phyloseq_4 (which does not contain nursery samples) and only rhizosphere samples retained
JH09_FC_data_phyloseq_5 <- subset_samples(JH09_FC_data_phyloseq_4, Microhabitat == "Rhizosphere")
JH09_FC_data_phyloseq_5
#extract the count matrix
dat_count_threshold_rhizosphere  <- as.data.frame(otu_table(JH09_FC_data_phyloseq_5))
#extract the cognate taxonomy information
dat_tax_threshold_rhizosphere  <- as.data.frame(dat_tax_noPlants[rownames(dat_count_threshold_rhizosphere  ), ])
colnames(dat_tax_threshold_rhizosphere) <- c("ConsensusLineage")
colnames(dat_tax_threshold_rhizosphere)

#merge count and taxonomy files
dat_info_threshold_rhizosphere <- cbind(dat_count_threshold_rhizosphere, dat_tax_threshold_rhizosphere)
dim(dat_info_threshold_rhizosphere)
colnames(dat_info_threshold_rhizosphere)
dat_info_threshold_rhizosphere[1:5, ]

#generate a new design file for QIIME calculation
design_4 <- design[colnames(dat_count_threshold_rhizosphere), ]
design_4

#save the files for the reproducibility of the code and the analysis in QIIME. Note: due to R formatting, remember to shift columns to the right and add the '#OTU ID' to the first column. Save the file as _2
#write.table(dat_info_threshold_rhizosphere, file="JH09_FC_info_threshold_rhizosphere.txt", sep = "\t")
#The above mentioned file is designated Worksheet_ws6 in Supplementary Information
#save mapping file for DESeq calculations. Note: due to R formatting, remember to shift columns to the right and add the '#SampleID' to the first column. Save the file as _2
#write.table(design_4, file="Map_JH09_FC_treatment_rhizosphere.txt", sep = "\t") 
#The above mentioned file is designated Worksheet_ws7 in Supplementary Information

#Treatment effect: roots
#generate a new OTU table for calculation in QIIME
#extract OTU counts (non rarefied, yet filtered for abundance) from an abundance threshold phyloseq object
#the starting file is JH09_FC_data_phyloseq_4 (which does not contain nursery samples) and only root samples retained
JH09_FC_data_phyloseq_6 <- subset_samples(JH09_FC_data_phyloseq_4, Microhabitat == "Roots")
JH09_FC_data_phyloseq_6
#extract the count matrix
dat_count_threshold_roots  <- as.data.frame(otu_table(JH09_FC_data_phyloseq_6))
#extract the cognate taxonomy information
dat_tax_threshold_roots  <- as.data.frame(dat_tax_noPlants[rownames(dat_count_threshold_roots  ), ])
colnames(dat_tax_threshold_roots) <- c("ConsensusLineage")
colnames(dat_tax_threshold_roots)

#merge count and taxonomy files
dat_info_threshold_roots <- cbind(dat_count_threshold_roots, dat_tax_threshold_roots)
dim(dat_info_threshold_roots)
colnames(dat_info_threshold_roots)
dat_info_threshold_roots[1:5, ]

#generate a new design file for QIIME calculation
design_5 <- design[colnames(dat_count_threshold_roots), ]
design_5

#save the files for the reproducibility of the code and the analysis in QIIME. Note: due to R formatting, remember to shift columns to the right and add the '#OTU ID' to the first column. Save the file as _2
#write.table(dat_info_threshold_roots, file="JH09_FC_info_threshold_roots.txt", sep = "\t")
#The above mentioned file is designated Worksheet_ws8 in Supplementary Information
#save mapping file for DESeq calculations. Note: due to R formatting, remember to shift columns to the right and add the '#SampleID' to the first column. Save the file as _2
#write.table(design_5, file="Map_JH09_FC_treatment_roots.txt", sep = "\t") 
#The above mentioned file is designated Worksheet_ws9 in Supplementary Information

#identify OTUs enriched in rhizosphere and roots responding to the treatments
#import the rhizosphere enriched table generated in QIIME
rhizosphere_enriched <- read.delim("JH09_FC_microhabitat_rhizosphere.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
rhizosphere_enriched[1:5, ]
#subset for rhizosphere enriched and adjusted p.value below 0.01
rhizosphere_eriched_001 <- rhizosphere_enriched[(rownames(rhizosphere_enriched)[which(rhizosphere_enriched$padj <0.01)]), ]
rhizosphere_eriched_001_2 <- rhizosphere_eriched_001[(rownames(rhizosphere_eriched_001)[which(rhizosphere_eriched_001$log2FoldChange < 0)]), ]
dim(rhizosphere_eriched_001)
dim(rhizosphere_eriched_001_2)

#import the treatment effect tables
#Pellet
pellet_enriched <- read.delim("JH09_FC_rhizosphere_Pellet.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
pellet_enriched[1:5, ]
#subset for pellet enriched and adjusted p.value below 0.01
pellet_eriched_001 <- pellet_enriched[(rownames(pellet_enriched)[which(pellet_enriched$padj <0.01)]), ]
pellet_eriched_001_2 <- pellet_eriched_001[(rownames(pellet_eriched_001)[which(pellet_eriched_001$log2FoldChange < 0)]), ]
dim(pellet_eriched_001)
dim(pellet_eriched_001_2)
pellet_eriched_001_2_rhizosphere = pellet_eriched_001_2
#define the impact of pellet on rhizosphere enriched bacteria
rhizosphere_enriched_pellet <- intersect(rownames(rhizosphere_eriched_001_2), rownames(pellet_eriched_001_2))
length(rhizosphere_enriched_pellet)
#inspect the taxonomies
pellet_enriched[rhizosphere_enriched_pellet, ]
#Liquiddigestate
Liquiddigestate_enriched <- read.delim("JH09_FC_rhizosphere_Liquiddigestate.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Liquiddigestate_enriched[1:5, ]
#subset for Liquiddigestate enriched and adjusted p.value below 0.01
Liquiddigestate_eriched_001 <- Liquiddigestate_enriched[(rownames(Liquiddigestate_enriched)[which(Liquiddigestate_enriched$padj <0.01)]), ]
Liquiddigestate_eriched_001_2 <- Liquiddigestate_eriched_001[(rownames(Liquiddigestate_eriched_001)[which(Liquiddigestate_eriched_001$log2FoldChange < 0)]), ]
dim(Liquiddigestate_eriched_001)
dim(Liquiddigestate_eriched_001_2)
#define the impact of Liquiddigestate on rhizosphere enriched bacteria
rhizosphere_enriched_Liquiddigestate <- intersect(rownames(rhizosphere_eriched_001_2), rownames(Liquiddigestate_eriched_001_2))
length(rhizosphere_enriched_Liquiddigestate)
#inspect the taxonomies
Liquiddigestate_enriched[rhizosphere_enriched_Liquiddigestate, ]
#SRLiquiddigestate
SRLiquiddigestate_enriched <- read.delim("JH09_FC_rhizosphere_SRLiquiddigestate.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
SRLiquiddigestate_enriched[1:5, ]
#subset for SRLiquiddigestate enriched and adjusted p.value below 0.01
SRLiquiddigestate_eriched_001 <- SRLiquiddigestate_enriched[(rownames(SRLiquiddigestate_enriched)[which(SRLiquiddigestate_enriched$padj <0.01)]), ]
SRLiquiddigestate_eriched_001_2 <- SRLiquiddigestate_eriched_001[(rownames(SRLiquiddigestate_eriched_001)[which(SRLiquiddigestate_eriched_001$log2FoldChange < 0)]), ]
dim(SRLiquiddigestate_eriched_001)
dim(SRLiquiddigestate_eriched_001_2)
#define the impact of SRLiquiddigestate on rhizosphere enriched bacteria
rhizosphere_enriched_SRLiquiddigestate <- intersect(rownames(rhizosphere_eriched_001_2), rownames(SRLiquiddigestate_eriched_001_2))
length(rhizosphere_enriched_SRLiquiddigestate)
#inspect the taxonomies
SRLiquiddigestate_enriched[rhizosphere_enriched_SRLiquiddigestate, ]
#SCAMproduct
SCAMproduct_enriched <- read.delim("JH09_FC_rhizosphere_SCAMproduct.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
SCAMproduct_enriched[1:5, ]
#subset for SCAMproduct enriched and adjusted p.value below 0.01
SCAMproduct_eriched_001 <- SCAMproduct_enriched[(rownames(SCAMproduct_enriched)[which(SCAMproduct_enriched$padj <0.01)]), ]
SCAMproduct_eriched_001_2 <- SCAMproduct_eriched_001[(rownames(SCAMproduct_eriched_001)[which(SCAMproduct_eriched_001$log2FoldChange < 0)]), ]
dim(SCAMproduct_eriched_001)
dim(SCAMproduct_eriched_001_2)
#define the impact of SCAMproduct on rhizosphere enriched bacteria
rhizosphere_enriched_SCAMproduct <- intersect(rownames(rhizosphere_eriched_001_2), rownames(SCAMproduct_eriched_001_2))
length(rhizosphere_enriched_SCAMproduct)
#inspect the taxonomies
SCAMproduct_enriched[rhizosphere_enriched_SCAMproduct, ]
#Mineralfertilizer
Mineralfertilizer_enriched <- read.delim("JH09_FC_rhizosphere_Mineralfertilizer.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Mineralfertilizer_enriched[1:5, ]
#subset for Mineralfertilizer enriched and adjusted p.value below 0.01
Mineralfertilizer_eriched_001 <- Mineralfertilizer_enriched[(rownames(Mineralfertilizer_enriched)[which(Mineralfertilizer_enriched$padj <0.01)]), ]
Mineralfertilizer_eriched_001_2 <- Mineralfertilizer_eriched_001[(rownames(Mineralfertilizer_eriched_001)[which(Mineralfertilizer_eriched_001$log2FoldChange < 0)]), ]
dim(Mineralfertilizer_eriched_001)
dim(Mineralfertilizer_eriched_001_2)
#define the impact of Mineralfertilizer on rhizosphere enriched bacteria
rhizosphere_enriched_Mineralfertilizer <- intersect(rownames(rhizosphere_eriched_001_2), rownames(Mineralfertilizer_eriched_001_2))
length(rhizosphere_enriched_Mineralfertilizer)
#inspect the taxonomies
Mineralfertilizer_enriched[rhizosphere_enriched_Mineralfertilizer, ]
#SRMineralfertilizer
SRMineralfertilizer_enriched <- read.delim("JH09_FC_rhizosphere_SRMineralfertilizer.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
SRMineralfertilizer_enriched[1:5, ]
#subset for SRMineralfertilizer enriched and adjusted p.value below 0.01
SRMineralfertilizer_eriched_001 <- SRMineralfertilizer_enriched[(rownames(SRMineralfertilizer_enriched)[which(SRMineralfertilizer_enriched$padj <0.01)]), ]
SRMineralfertilizer_eriched_001_2 <- SRMineralfertilizer_eriched_001[(rownames(SRMineralfertilizer_eriched_001)[which(SRMineralfertilizer_eriched_001$log2FoldChange < 0)]), ]
dim(SRMineralfertilizer_eriched_001)
dim(SRMineralfertilizer_eriched_001_2)
#define the impact of SRMineralfertilizer on rhizosphere enriched bacteria
rhizosphere_enriched_SRMineralfertilizer <- intersect(rownames(rhizosphere_eriched_001_2), rownames(SRMineralfertilizer_eriched_001_2))
length(rhizosphere_enriched_SRMineralfertilizer)
#inspect the taxonomies
SRMineralfertilizer_enriched[rhizosphere_enriched_SRMineralfertilizer, ]

#Prepare the data for UpSetR visualisation
#pellet
pellet_eriched_001_2_df <- as.data.frame(pellet_eriched_001_2[ ,1]) 
rownames(pellet_eriched_001_2_df) <- rownames(pellet_eriched_001_2)
colnames(pellet_eriched_001_2_df) <- c("counts_PE")
pellet_eriched_001_2_df[pellet_eriched_001_2_df > 1] <- 1
pellet_eriched_001_2_df
dim(pellet_eriched_001_2_df)
#Liquiddigestate
Liquiddigestate_eriched_001_2_df <- as.data.frame(Liquiddigestate_eriched_001_2[ ,1]) 
rownames(Liquiddigestate_eriched_001_2_df) <- rownames(Liquiddigestate_eriched_001_2)
colnames(Liquiddigestate_eriched_001_2_df) <- c("counts_LD")
Liquiddigestate_eriched_001_2_df[Liquiddigestate_eriched_001_2_df > 1] <- 1
Liquiddigestate_eriched_001_2_df
dim(Liquiddigestate_eriched_001_2_df)
#SRLiquiddigestate
SRLiquiddigestate_eriched_001_2_df <- as.data.frame(SRLiquiddigestate_eriched_001_2[ ,1]) 
rownames(SRLiquiddigestate_eriched_001_2_df) <- rownames(SRLiquiddigestate_eriched_001_2)
colnames(SRLiquiddigestate_eriched_001_2_df) <- c("counts_SRLD")
SRLiquiddigestate_eriched_001_2_df[SRLiquiddigestate_eriched_001_2_df > 1] <- 1
SRLiquiddigestate_eriched_001_2_df
dim(SRLiquiddigestate_eriched_001_2_df)
#SCAMproduct
SCAMproduct_eriched_001_2_df <- as.data.frame(SCAMproduct_eriched_001_2[ ,1]) 
rownames(SCAMproduct_eriched_001_2_df) <- rownames(SCAMproduct_eriched_001_2)
colnames(SCAMproduct_eriched_001_2_df) <- c("counts_SC")
SCAMproduct_eriched_001_2_df[SCAMproduct_eriched_001_2_df > 1] <- 1
SCAMproduct_eriched_001_2_df
dim(SCAMproduct_eriched_001_2_df)
#Mineralfertilizer
Mineralfertilizer_eriched_001_2_df <- as.data.frame(Mineralfertilizer_eriched_001_2[ ,1]) 
rownames(Mineralfertilizer_eriched_001_2_df) <- rownames(Mineralfertilizer_eriched_001_2)
colnames(Mineralfertilizer_eriched_001_2_df) <- c("counts_MF")
Mineralfertilizer_eriched_001_2_df[Mineralfertilizer_eriched_001_2_df > 1] <- 1
Mineralfertilizer_eriched_001_2_df
dim(Mineralfertilizer_eriched_001_2_df)
#SRMineralfertilizer
SRMineralfertilizer_eriched_001_2_df <- as.data.frame(SRMineralfertilizer_eriched_001_2[ ,1]) 
rownames(SRMineralfertilizer_eriched_001_2_df) <- rownames(SRMineralfertilizer_eriched_001_2)
colnames(SRMineralfertilizer_eriched_001_2_df) <- c("counts_SRMF")
SRMineralfertilizer_eriched_001_2_df[SRMineralfertilizer_eriched_001_2_df > 1] <- 1
SRMineralfertilizer_eriched_001_2_df
dim(SRMineralfertilizer_eriched_001_2_df)
#combine the datasets: note they have unequal values
#define a list of unique OTUs
OTU_list <- unique(c(c((c(rownames(pellet_eriched_001_2_df), rownames(Liquiddigestate_eriched_001_2_df))),
                        (c(rownames(SRLiquiddigestate_eriched_001_2_df), rownames(SCAMproduct_eriched_001_2_df)))), 
  (c(rownames(Mineralfertilizer_eriched_001_2_df), rownames(SRMineralfertilizer_eriched_001_2_df)))))
length(OTU_list)
#Pellet
PE_eriched_merging <- as.data.frame(pellet_eriched_001_2_df[OTU_list, ])
colnames(PE_eriched_merging) <- c("OTUs_PE")
row.names(PE_eriched_merging) <- as.vector(OTU_list)
#Liquid digestate
LD_eriched_merging <- as.data.frame(Liquiddigestate_eriched_001_2_df[OTU_list, ])
colnames(LD_eriched_merging) <- c("OTUs_LD")
row.names(LD_eriched_merging) <- as.vector(OTU_list)
#SRLiquid digestate
SRLD_eriched_merging <- as.data.frame(SRLiquiddigestate_eriched_001_2_df[OTU_list, ])
colnames(SRLD_eriched_merging) <- c("OTUs_SRLD")
row.names(SRLD_eriched_merging) <- as.vector(OTU_list)
#SCAM_product
SC_eriched_merging <- as.data.frame(SCAMproduct_eriched_001_2_df[OTU_list, ])
colnames(SC_eriched_merging) <- c("OTUs_SC")
row.names(SC_eriched_merging) <- as.vector(OTU_list)
#Mineralfertiliser
MF_eriched_merging <- as.data.frame(Mineralfertilizer_eriched_001_2_df[OTU_list, ])
colnames(MF_eriched_merging) <- c("OTUs_MF")
row.names(MF_eriched_merging) <- as.vector(OTU_list)
#SRMineralfertiliser
SRMF_eriched_merging <- as.data.frame(SRMineralfertilizer_eriched_001_2_df[OTU_list, ])
colnames(SRMF_eriched_merging) <- c("OTUs_SRMF")
row.names(SRMF_eriched_merging) <- as.vector(OTU_list)
#Merge the dataset
rhizosphere_treatment_OTUs <- cbind(PE_eriched_merging, LD_eriched_merging)
rhizosphere_treatment_OTUs <- cbind(rhizosphere_treatment_OTUs, SRLD_eriched_merging)
rhizosphere_treatment_OTUs <- cbind(rhizosphere_treatment_OTUs, SC_eriched_merging)
rhizosphere_treatment_OTUs <- cbind(rhizosphere_treatment_OTUs, MF_eriched_merging)
rhizosphere_treatment_OTUs <- cbind(rhizosphere_treatment_OTUs, SRMF_eriched_merging)
#set NA to 0
rhizosphere_treatment_OTUs[is.na(rhizosphere_treatment_OTUs)] <- 0
#visualisation
upset(rhizosphere_treatment_OTUs, sets = c("OTUs_PE", "OTUs_LD", "OTUs_SRLD", "OTUs_SC", "OTUs_MF", "OTUs_SRMF"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")


#identify OTUs enriched in roots and roots responding to the treatments
#import the roots enriched table generated in QIIME
roots_enriched <- read.delim("JH09_FC_microhabitat_roots.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
roots_enriched[1:5, ]
#subset for roots enriched and adjusted p.value below 0.05
roots_eriched_001 <- roots_enriched[(rownames(roots_enriched)[which(roots_enriched$padj <0.01)]), ]
roots_eriched_001_2 <- roots_eriched_001[(rownames(roots_eriched_001)[which(roots_eriched_001$log2FoldChange < 0)]), ]
dim(roots_eriched_001)
dim(roots_eriched_001_2)
#write.table(roots_eriched_001_2, file="JH09_FC_roots_enriched_001_2.txt", sep = "\t")
#define the the proportion of root enriched that are also rhizosphere enriched
plant_effect <- intersect(rownames(rhizosphere_eriched_001_2), rownames(roots_eriched_001_2))
length(plant_effect)
#identify the tanomy affiliation
plant_effect_taxa <- roots_eriched_001_2[plant_effect, ]
plant_effect_taxa
#write.table(plant_effect_taxa, file="JH09_FC_plant_eriched_001_2.txt", sep = "\t")
#define the root effect
root_enriched_only <- setdiff(rownames(roots_eriched_001_2), plant_effect)
root_enriched_only
root_enriched_taxa <- roots_eriched_001_2[root_enriched_only, ]
root_enriched_taxa
#write.table(root_enriched_taxa, file="JH09_FC_roots_enriched_only_001_2.txt", sep = "\t")

#create a phylogenetic tree for visualisation of rhizosphere and root enriched
JH09_FC_data_phyloseq_7 <- prune_taxa(unique(union(rownames(rhizosphere_eriched_001_2), rownames(roots_eriched_001_2))), JH09_FC_data_phyloseq_3)
JH09_FC_data_phyloseq_7
#phylogenetic tree
tree_FC = phy_tree(JH09_FC_data_phyloseq_7) 
#ape::write.tree(tree_FC, "JH09_FC_rhizo_root_enriched.tree")
#create the annotation files
#Rhizosphere enriched
tree_FC_taxonomic_annotation <- JH09_FC_taxa_ordered[ unique(union(rownames(rhizosphere_eriched_001_2), rownames(roots_eriched_001_2))), ]
#write.table(tree_FC_taxonomic_annotation, file="tree_FC_taxonomic_annotation.txt", sep = "\t")
#Microhabitat annotation
#Rhizosphere
rhizosphere_eriched_001_2_df <- as.data.frame(rhizosphere_eriched_001_2[ ,1]) 
rownames(rhizosphere_eriched_001_2_df) <- rownames(rhizosphere_eriched_001_2)
colnames(rhizosphere_eriched_001_2_df) <- c("counts_Rhizo")
rhizosphere_eriched_001_2_df[rhizosphere_eriched_001_2_df > 1] <- 1
rhizosphere_eriched_001_2_df
dim(rhizosphere_eriched_001_2_df)
#Roots
roots_eriched_001_2_df <- as.data.frame(roots_eriched_001_2[ ,1]) 
rownames(roots_eriched_001_2_df) <- rownames(roots_eriched_001_2)
colnames(roots_eriched_001_2_df) <- c("counts_Roots")
roots_eriched_001_2_df[roots_eriched_001_2_df > 1] <- 1
roots_eriched_001_2_df
dim(roots_eriched_001_2_df)

#combine the datasets: note they have unequal values
#define a list of unique OTUs
OTU_list_microhabitat <- unique(union(rownames(rhizosphere_eriched_001_2), rownames(roots_eriched_001_2)))
                       
#Rhizosphere
rhizo_eriched_merging <- as.data.frame(rhizosphere_eriched_001_2_df[OTU_list_microhabitat, ])
colnames(rhizo_eriched_merging) <- c("OTUs_rhizo")
row.names(rhizo_eriched_merging) <- as.vector(OTU_list_microhabitat)
#roots
roots_eriched_merging <- as.data.frame(roots_eriched_001_2_df[OTU_list_microhabitat, ])
colnames(roots_eriched_merging) <- c("OTUs_roots")
row.names(roots_eriched_merging) <- as.vector(OTU_list_microhabitat)

#Merge the dataset
microhabitat_OTUs <- cbind(rhizo_eriched_merging, roots_eriched_merging)
#set NA to 0
microhabitat_OTUs[is.na(microhabitat_OTUs)] <- 0
#write.table(microhabitat_OTUs, file="tree_FC_microhabitat_annotation.txt", sep = "\t")


#import the treatment effect tables: roots
#Pellet
pellet_enriched <- read.delim("JH09_FC_roots_Pellet.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
pellet_enriched[1:5, ]
#subset for pellet enriched and adjusted p.value below 0.05
pellet_eriched_001 <- pellet_enriched[(rownames(pellet_enriched)[which(pellet_enriched$padj <0.01)]), ]
pellet_eriched_001_2 <- pellet_eriched_001[(rownames(pellet_eriched_001)[which(pellet_eriched_001$log2FoldChange < 0)]), ]
dim(pellet_eriched_001)
dim(pellet_eriched_001_2)
pellet_eriched_001_2_roots = pellet_eriched_001_2
#Liquiddigestate
Liquiddigestate_enriched <- read.delim("JH09_FC_roots_Liquiddigestate.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Liquiddigestate_enriched[1:5, ]
#subset for Liquiddigestate enriched and adjusted p.value below 0.05
Liquiddigestate_eriched_001 <- Liquiddigestate_enriched[(rownames(Liquiddigestate_enriched)[which(Liquiddigestate_enriched$padj <0.01)]), ]
Liquiddigestate_eriched_001_2 <- Liquiddigestate_eriched_001[(rownames(Liquiddigestate_eriched_001)[which(Liquiddigestate_eriched_001$log2FoldChange < 0)]), ]
dim(Liquiddigestate_eriched_001)
dim(Liquiddigestate_eriched_001_2)
#SRLiquiddigestate
SRLiquiddigestate_enriched <- read.delim("JH09_FC_roots_SRLiquiddigestate.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
SRLiquiddigestate_enriched[1:5, ]
#subset for SRLiquiddigestate enriched and adjusted p.value below 0.05
SRLiquiddigestate_eriched_001 <- SRLiquiddigestate_enriched[(rownames(SRLiquiddigestate_enriched)[which(SRLiquiddigestate_enriched$padj <0.01)]), ]
SRLiquiddigestate_eriched_001_2 <- SRLiquiddigestate_eriched_001[(rownames(SRLiquiddigestate_eriched_001)[which(SRLiquiddigestate_eriched_001$log2FoldChange < 0)]), ]
dim(SRLiquiddigestate_eriched_001)
dim(SRLiquiddigestate_eriched_001_2)
#SCAMproduct
SCAMproduct_enriched <- read.delim("JH09_FC_roots_SCAMproduct.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
SCAMproduct_enriched[1:5, ]
#subset for SCAMproduct enriched and adjusted p.value below 0.05
SCAMproduct_eriched_001 <- SCAMproduct_enriched[(rownames(SCAMproduct_enriched)[which(SCAMproduct_enriched$padj <0.01)]), ]
SCAMproduct_eriched_001_2 <- SCAMproduct_eriched_001[(rownames(SCAMproduct_eriched_001)[which(SCAMproduct_eriched_001$log2FoldChange < 0)]), ]
dim(SCAMproduct_eriched_001)
dim(SCAMproduct_eriched_001_2)
#Mineralfertilizer
Mineralfertilizer_enriched <- read.delim("JH09_FC_roots_Mineralfertilizer.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
Mineralfertilizer_enriched[1:5, ]
#subset for Mineralfertilizer enriched and adjusted p.value below 0.05
Mineralfertilizer_eriched_001 <- Mineralfertilizer_enriched[(rownames(Mineralfertilizer_enriched)[which(Mineralfertilizer_enriched$padj <0.01)]), ]
Mineralfertilizer_eriched_001_2 <- Mineralfertilizer_eriched_001[(rownames(Mineralfertilizer_eriched_001)[which(Mineralfertilizer_eriched_001$log2FoldChange < 0)]), ]
dim(Mineralfertilizer_eriched_001)
dim(Mineralfertilizer_eriched_001_2)
#SRMineralfertilizer
SRMineralfertilizer_enriched <- read.delim("JH09_FC_roots_SRMineralfertilizer.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
SRMineralfertilizer_enriched[1:5, ]
#subset for SRMineralfertilizer enriched and adjusted p.value below 0.05
SRMineralfertilizer_eriched_001 <- SRMineralfertilizer_enriched[(rownames(SRMineralfertilizer_enriched)[which(SRMineralfertilizer_enriched$padj <0.01)]), ]
SRMineralfertilizer_eriched_001_2 <- SRMineralfertilizer_eriched_001[(rownames(SRMineralfertilizer_eriched_001)[which(SRMineralfertilizer_eriched_001$log2FoldChange < 0)]), ]
dim(SRMineralfertilizer_eriched_001)
dim(SRMineralfertilizer_eriched_001_2)

###################################################################################
#Data visualisation using UpSetR
###################################################################################

#Prepare the data for UpSetR visualisation
#pellet
pellet_eriched_001_2_df <- as.data.frame(pellet_eriched_001_2[ ,1]) 
rownames(pellet_eriched_001_2_df) <- rownames(pellet_eriched_001_2)
colnames(pellet_eriched_001_2_df) <- c("counts_PE")
pellet_eriched_001_2_df[pellet_eriched_001_2_df > 1] <- 1
pellet_eriched_001_2_df
dim(pellet_eriched_001_2_df)
#Liquiddigestate
Liquiddigestate_eriched_001_2_df <- as.data.frame(Liquiddigestate_eriched_001_2[ ,1]) 
rownames(Liquiddigestate_eriched_001_2_df) <- rownames(Liquiddigestate_eriched_001_2)
colnames(Liquiddigestate_eriched_001_2_df) <- c("counts_LD")
Liquiddigestate_eriched_001_2_df[Liquiddigestate_eriched_001_2_df > 1] <- 1
Liquiddigestate_eriched_001_2_df
dim(Liquiddigestate_eriched_001_2_df)
#SRLiquiddigestate
SRLiquiddigestate_eriched_001_2_df <- as.data.frame(SRLiquiddigestate_eriched_001_2[ ,1]) 
rownames(SRLiquiddigestate_eriched_001_2_df) <- rownames(SRLiquiddigestate_eriched_001_2)
colnames(SRLiquiddigestate_eriched_001_2_df) <- c("counts_SRLD")
SRLiquiddigestate_eriched_001_2_df[SRLiquiddigestate_eriched_001_2_df > 1] <- 1
SRLiquiddigestate_eriched_001_2_df
dim(SRLiquiddigestate_eriched_001_2_df)
#SCAMproduct
SCAMproduct_eriched_001_2_df <- as.data.frame(SCAMproduct_eriched_001_2[ ,1]) 
rownames(SCAMproduct_eriched_001_2_df) <- rownames(SCAMproduct_eriched_001_2)
colnames(SCAMproduct_eriched_001_2_df) <- c("counts_SC")
SCAMproduct_eriched_001_2_df[SCAMproduct_eriched_001_2_df > 1] <- 1
SCAMproduct_eriched_001_2_df
dim(SCAMproduct_eriched_001_2_df)
#Mineralfertilizer
Mineralfertilizer_eriched_001_2_df <- as.data.frame(Mineralfertilizer_eriched_001_2[ ,1]) 
rownames(Mineralfertilizer_eriched_001_2_df) <- rownames(Mineralfertilizer_eriched_001_2)
colnames(Mineralfertilizer_eriched_001_2_df) <- c("counts_MF")
Mineralfertilizer_eriched_001_2_df[Mineralfertilizer_eriched_001_2_df > 1] <- 1
Mineralfertilizer_eriched_001_2_df
dim(Mineralfertilizer_eriched_001_2_df)
#SRMineralfertilizer
SRMineralfertilizer_eriched_001_2_df <- as.data.frame(SRMineralfertilizer_eriched_001_2[ ,1]) 
rownames(SRMineralfertilizer_eriched_001_2_df) <- rownames(SRMineralfertilizer_eriched_001_2)
colnames(SRMineralfertilizer_eriched_001_2_df) <- c("counts_SRMF")
SRMineralfertilizer_eriched_001_2_df[SRMineralfertilizer_eriched_001_2_df > 1] <- 1
SRMineralfertilizer_eriched_001_2_df
dim(SRMineralfertilizer_eriched_001_2_df)
#combine the datasets: note they have unequal values
#define a list of unique OTUs
OTU_list <- unique(c(c((c(rownames(pellet_eriched_001_2_df), rownames(Liquiddigestate_eriched_001_2_df))),
                       (c(rownames(SRLiquiddigestate_eriched_001_2_df), rownames(SCAMproduct_eriched_001_2_df)))), 
                     (c(rownames(Mineralfertilizer_eriched_001_2_df), rownames(SRMineralfertilizer_eriched_001_2_df)))))
length(OTU_list)
#Pellet
PE_eriched_merging <- as.data.frame(pellet_eriched_001_2_df[OTU_list, ])
colnames(PE_eriched_merging) <- c("OTUs_PE")
row.names(PE_eriched_merging) <- as.vector(OTU_list)
#Liquid digestate
LD_eriched_merging <- as.data.frame(Liquiddigestate_eriched_001_2_df[OTU_list, ])
colnames(LD_eriched_merging) <- c("OTUs_LD")
row.names(LD_eriched_merging) <- as.vector(OTU_list)
#SRLiquid digestate
SRLD_eriched_merging <- as.data.frame(SRLiquiddigestate_eriched_001_2_df[OTU_list, ])
colnames(SRLD_eriched_merging) <- c("OTUs_SRLD")
row.names(SRLD_eriched_merging) <- as.vector(OTU_list)
#SCAM_product
SC_eriched_merging <- as.data.frame(SCAMproduct_eriched_001_2_df[OTU_list, ])
colnames(SC_eriched_merging) <- c("OTUs_SC")
row.names(SC_eriched_merging) <- as.vector(OTU_list)
#Mineralfertiliser
MF_eriched_merging <- as.data.frame(Mineralfertilizer_eriched_001_2_df[OTU_list, ])
colnames(MF_eriched_merging) <- c("OTUs_MF")
row.names(MF_eriched_merging) <- as.vector(OTU_list)
#SRMineralfertiliser
SRMF_eriched_merging <- as.data.frame(SRMineralfertilizer_eriched_001_2_df[OTU_list, ])
colnames(SRMF_eriched_merging) <- c("OTUs_SRMF")
row.names(SRMF_eriched_merging) <- as.vector(OTU_list)
#Merge the dataset
roots_treatment_OTUs <- cbind(PE_eriched_merging, LD_eriched_merging)
roots_treatment_OTUs <- cbind(roots_treatment_OTUs, SRLD_eriched_merging)
roots_treatment_OTUs <- cbind(roots_treatment_OTUs, SC_eriched_merging)
roots_treatment_OTUs <- cbind(roots_treatment_OTUs, MF_eriched_merging)
roots_treatment_OTUs <- cbind(roots_treatment_OTUs, SRMF_eriched_merging)
#set NA to 0
roots_treatment_OTUs[is.na(roots_treatment_OTUs)] <- 0
#visualisation
upset(roots_treatment_OTUs, sets = c("OTUs_PE", "OTUs_LD", "OTUs_SRLD", "OTUs_SC", "OTUs_MF", "OTUs_SRMF"), sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

#biocLite("Tax4Fun")
#biocLite("qiimer")
#library("Tax4Fun")
#library("qiimer")
#library("biom")

#########################################################################################################################################
#TAX4FUN
#PREDICTION OF FUNCTIONAL CAPABILITIES OF THE MICROBIOME
#data required: JH09_FC_otu_table_silva115_gg135filtered.txt, 
#silva database 115 folder
##########################################################################################################################################

OTU_TABLE_tax4fun <- importQIIMEData("JH09_FC_otu_table_silva115_gg135filtered.txt")
folderReferenceData <- ("C:\Users\ra42320\Box Sync2\Box Sync\Davide_lab_manuscripts\Federica_tomato_2018\tax4fun\Tax4Fun_input_files\SilvaSSURef_115_NR")
Tax4FunOutput <- Tax4Fun(OTU_TABLE_tax4fun, folderReferenceData, fctProfiling = F, refProfile = "UProC", shortReadMode = TRUE, normCopyNo = TRUE)
Tax4FunProfile<-Tax4FunOutput$Tax4FunProfile
Tax4FunProfile <- data.frame(t(Tax4FunOutput$Tax4FunProfile))
write.table(Tax4FunProfile, file="Tax4FunProfile_JH09_FC.txt", sep="\t")

########################################################################################################################################
