
setwd("/Users/denmanmoody/Documents/GitHub/mb2021_phyloseq")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

library("devtools")
library(phyloseq)
library(microbiome)

##### Error message that ggplot2 needed to be unloaded and reinstalled -- Fixed by:
remove.packages("ggplot2")
remove.packages("microbiome")

pseq<- readRDS("Denman_samples.rds")

pseq

## View sample data
summarize_phyloseq(pseq)

#create objects

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data

# check if any OTUs are not present in any samples (if any OTU slots are empty)
any(taxa_sums(pseq) == 0)

#source for removing chloro/mito/archaea; https://mibwurrepo.github.io/R_for_Microbial_Ecology/Microbiome_tutorial_V2.html#making-a-phyloseq-object

#remove chloroplast, mito, archaea
#the c__ o__ might not be necessary - left over from old code from Andy. Run the code and then check to see of chloro, mito, and arhaea get removed

pseq1 <- subset_taxa(pseq,Class!="c__Chloroplast")

pseq2 <- subset_taxa(pseq1,Order!="o__Mitochondria")

ps1 <- subset_taxa(pseq2,Kingdom!="Archaea")

#check to see if all chloroplasts, mitochondria, and archaea have been removed
View(ps1@tax_table)

#chloroplasts and mitochondria still in data under order level - remove

pseq3 <- subset_taxa(ps1,Order!="Chloroplast")

pseq4 <- subset_taxa(pseq3,Family!="Mitochondria")

View(pseq4@tax_table)

#Rename pseq4 to ps1 to use with the rest of the code
ps1 <- pseq4

#Making sure ranking are capitalized....??? --> Not super relevant for my project
rank_names(ps1)


###remove low prevalency reads --> Removing bacterial taxa that are in the tail
x1 = prune_taxa(taxa_sums(ps1) > 200, ps1) 
x2 = prune_taxa(taxa_sums(ps1) > 500, ps1) 
x3 = prune_taxa(taxa_sums(ps1) > 1000, ps1)

###Andy says to just prune to x1 (200) because this leaves us with 312 taxa and Andy usually tries to keep 300-500


#view ps1 - data before any pruning
plot(sort(taxa_sums(ps1), TRUE), type="h", ylim=c(0, 10000))

#view x1 - data after pruning off 200
plot(sort(taxa_sums(x1), TRUE), type="h", ylim=c(0, 10000))

#view x2 - data after pruning off 500
plot(sort(taxa_sums(x2), TRUE), type="h", ylim=c(0, 10000))

library(microbiome)

#looking at overview of the data after trimming  
summarize_phyloseq(ps1)
summarize_phyloseq(x1)
summarize_phyloseq(x2)

###Andy says to just prune to x1 (200) because this leaves us with 312 taxa and Andy usually tries to keep 300-500



##rarify data to make sequences an even read depth

#run code to check the minimum number of reads in data is now that it's been filtered and pruned
summarize_phyloseq(x1)

#min number of reads is 4759
#selecting read depth of 4750 = any samples with fewer than 4750 total reads will be removed

Rare <-rarefy_even_depth(x1, sample.size= 4750)

Rare
# 312 taxa = 312 ASVs

#check to see that all samples have been rarified to 4750 reads
summarize_phyloseq(Rare)

sample_depths <- sample_sums(Rare)
print(sample_depths)

## save filtered data
saveRDS(Rare, file="Denman_Samples_Rare_FINAL.rds")

##rename Rare to ps5 to use in future plots
ps5 <- Rare

sample_depths <- sample_sums(ps5)

print(sample_depths)

#############################################################




##################################

###calculate relative abundance at family, phylum, and genus level.

pseq_fam <- microbiome::aggregate_rare(ps5, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(ps5, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(ps5, level = "Genus", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")



##################################

#Perform plot ordination

#Set seed is code to keep axes on plots the same between different plots (plot ordination changes axes everytime you run it)
set.seed(4235421)

ord <- ordinate(ps5, "MDS", "bray")

#plot MDS/PcoA - can set "colour" and "shape" for any of your variables
#geompoint controls data point size on plot

plot_ordination(ps5, ord, color = "Sample.ID", shape = "Factor") + geom_point(size = 4)

plot_ordination(ps5, ord, color = "Factor", shape = "Sample.ID") + geom_point(size = 4) + scale_shape_binned()

##plot having issues because Family.1 is saved as an integer rather than a character --> Change column data type within metadata
ps5@sam_data$Family.1 <- as.character(ps5@sam_data$Family.1)


plot_ordination(ps5, ord, color = "Factor", shape = "Family.1") + geom_point(size = 4)

plot_ordination(ps5, ord, color = "Factor") + geom_point(size = 4)

plot_ordination(ps5, ord, color = "Factor") + geom_point(size = 4)

plot_ordination(ps5, ord, color = "Factor") + geom_point(size = 4)

plot_ordination(ps5, ord, color = "Family.1") + geom_point(size = 4)

plot_ordination(ps5, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)


p <- plot_ordination(ps5, ord, color = "Factor") + geom_point(size = 4)



