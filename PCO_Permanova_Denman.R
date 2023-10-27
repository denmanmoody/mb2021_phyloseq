
setwd("/Users/denmanmoody/Documents/GitHub/mb2021_phyloseq")

##Quick lib ---- 
#(if you have packages downloaded already)

library("devtools")
library(phyloseq)
library(microbiome)
library(vegan)

### To pull saved rarified samples
pseq<- readRDS("Denman_Samples_Rare_FINAL.rds")

## View sample data -- Rarified data should have 4750 reads for each ASV
summarize_phyloseq(pseq)


####################################################################

#Perform plot ordination

#Set seed is code to keep axes on plots the same between different plots (plot ordination changes axes everytime you run it)
set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")

##Change Family.1 from an integer to a character --> Change column data type within metadata
pseq@sam_data$Family.1 <- as.character(pseq@sam_data$Family.1)

#plot MDS/PcoA - can set "colour" and "shape" for any of your variables
#geompoint controls data point size on plot

plot_ordination(pseq, ord, color = "Family.1", shape = "Factor") + geom_point(size = 4)

plot_ordination(pseq, ord, color = "Factor") + geom_point(size = 4)


#not sure what this one is for...? Not sure what binned code is for
plot_ordination(pseq, ord, color = "Factor", shape = "Sample.ID") + geom_point(size = 4) + scale_shape_binned()

####################################################################


######### PERMANOVA ----
#source: https://microbiome.github.io/tutorials/PERMANOVA.html


#if any columns have missing values (NA), must remove

pseq_filtered <- subset_samples(pseq, !Family.1 %in% NA)

View(pseq_filtered@sam_data)

set.seed(423542)

Bray_dist<- phyloseq::distance(pseq_filtered, method = "bray", weighted = TRUE)

Sample_star <- data.frame(sample_data(pseq_filtered))

Test_bray <- adonis2(Bray_dist ~ Family.1*Factor, data = Sample_star)

Test_bray2 <- adonis2(Bray_dist ~ Factor*Family.1, data = Sample_star)

##################################################

