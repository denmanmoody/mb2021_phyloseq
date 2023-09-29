#Look at pseq.core

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)

#Load data ----


pseq <-  Marissa_MU42022_rare

#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#Day 1 only ----

pseq <- subset_samples(pseq.core, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

pseq_ <- subset_samples(pseq, !Sample.type %in% "Algae")


#Create pseq objects ----

pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")



#plot MDS/PcoA ----

set.seed(4235421)

ord <- ordinate(pseq.core, "MDS", "bray")

plot_ordination(pseq.core, ord, color = "Treatment", shape = "Age") + geom_point(size = 4)

#Day 1 only ----

pseq_core_filtered <- subset_samples(pseq.core, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

pseq_core_filtered <- subset_samples(pseq_core_filtered, !Sample.type %in% "Algae")

ord <- ordinate(pseq_core_filtered, "MDS", "bray")

plot_ordination(pseq_core_filtered, ord, color = "Treatment") + geom_point(size = 4)


#Extract core data ----

library(ggplot2)

file_path <- "C:/Users/maris/OneDrive/Documents/GitHub/pseq_core.csv"

# Write the otu_table to a CSV file
write.csv(pseq_core_filtered@otu_table, file = file_path)

if (file.exists(file_path)) {
  cat("OTU table CSV file has been created:", file_path, "\n")
} else {
  cat("Failed to create the CSV file.")
}

#Have to manually change sample names into treatments in excel sheet
#import excel sheet - have to make it with x, y setup

ggplot(pseq_core, aes(x=Treatment, y=Family, size = Abundance)) +
  geom_point(alpha=0.7)
