#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)
library(ggplot2)

#Load data ----
setwd("/Users/denmanmoody/Documents/GitHub/mb2021_phyloseq")

### To pull saved rarified samples
pseq<- readRDS("Denman_Samples_Rare_FINAL.rds")

## View sample data
summarize_phyloseq(pseq)

#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#create objects


pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

#convert to compositional data

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

##############################################################################

###PHYLUM LEVEL ----
###PHYLUM LEVEL - TOTAL ABUNDANCE

plot_bar(pseq, fill="Phylum")

#Sort the Phyla by abundance and pick the top 5

top5P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Phylum"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 5 Phyla

top5P = subset_taxa(pseq, Phylum %in% names(top5P.names))

#Plot of fouled vs non-fouled scallops with bar graph showing top 5 phyla in STACKED bar

plot_bar(top5P, x="Factor", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + labs(y = "Total Abundance")

###NOTE: To remove the grey block background from the plot above, use the + theme_classic() code:

plot_bar(top5P, x="Factor", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + labs(y = "Total Abundance") + theme_classic()


####################

#Plot of fouled vs non-fouled scallops with bar graph with top 5 phyla separated into their bar graphs

plot_bar(top5P, x="Factor", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + labs(y = "Total Abundance") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)) 



##############################################################################

###PHYLUM LEVEL - RELATIVE ABUNDANCE
# instead of using pseq (total abundance) plot with pseq.fam.rel (relative abundance at Family level)

#Sort the Phyla by RELATIVE abundance and pick the top 5
top5P.names = sort(tapply(taxa_sums(pseq.phy.rel), tax_table(pseq.phy.rel)[, "Phylum"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 5 Phyla
top5P = subset_taxa(pseq.phy.rel, Phylum %in% names(top5P.names))

#Plot of fouled vs non-fouled scallops with bar graph showing top 5 phyla in STACKED bar

plot_bar(top5P, x="Factor", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + labs(y = "Relative Abundance")


#Plot of fouled vs non-fouled scallops with bar graph with top 5 phyla separated into their bar graphs

plot_bar(top5P, x="Factor", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + labs(y = "Relative Abundance") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)) 




##############################################################################
















###FAMILY LEVEL ----
###FAMILY LEVEL - TOTAL ABUNDANCE - TOP 5 FAMILIES

#Sort the Families by total abundance and pick the top 5

top5P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Family"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 5 Families

top5P = subset_taxa(pseq, Family %in% names(top5P.names))

#Plot of fouled vs non-fouled scallops with bar graph showing top 5 families in STACKED bar

plot_bar(top5P, x="Factor", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Total Abundance", x = " ") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12))


#Plot of fouled vs non-fouled scallops with bar graph with top 5 families separated into their bar graphs

plot_bar(top5P, x="Factor", fill="Family", facet_grid = ~Family) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Total Abundance") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)) 



##############################################################################

###FAMILY LEVEL - TOTAL ABUNDANCE - TOP 10 FAMILIES

#Sort the Families by total abundance and pick the top 10

top10P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Family"], sum), TRUE)[1:10]

#Cut down the physeq data to only the top 10 Families

top10P = subset_taxa(pseq.fam.rel, Family %in% names(top10P.names))

#Plot of fouled vs non-fouled scallops with bar graph showing top 10 families in STACKED bar

plot_bar(top10P, x="Factor", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Total Abundance", x = " ") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12))


#Plot of fouled vs non-fouled scallops with bar graph with top 10 families separated into their bar graphs

plot_bar(top10P, x="Factor", fill="Family", facet_grid = ~Family) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Total Abundance") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)) 


##############################################################################


###FAMILY LEVEL - TOTAL ABUNDANCE - TOP 15 FAMILIES

#Sort the Families by total abundance and pick the top 15

top15P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Family"], sum), TRUE)[1:15]

#Cut down the physeq data to only the top 15 Families

top15P = subset_taxa(pseq.fam.rel, Family %in% names(top15P.names))

#Plot of fouled vs non-fouled scallops with bar graph showing top 15 families in STACKED bar

plot_bar(top15P, x="Factor", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Total Abundance", x = " ") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12))


#Plot of fouled vs non-fouled scallops with bar graph with top 10 families separated into their bar graphs

plot_bar(top15P, x="Factor", fill="Family", facet_grid = ~Family) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Total Abundance") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)) 


##############################################################################


###FAMILY LEVEL - RELATIVE ABUNDANCE - TOP 10 FAMILIES

#Sort the Families by relative abundance and pick the top 10

top10P.names = sort(tapply(taxa_sums(pseq.fam.rel), tax_table(pseq.fam.rel)[, "Family"], sum), TRUE)[1:10]

#Cut down the physeq data to only the top 10 Families

top10P = subset_taxa(pseq.fam.rel, Family %in% names(top10P.names))

#Plot of fouled vs non-fouled scallops with bar graph showing top 10 families in STACKED bar

plot_bar(top10P, x="Factor", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Relative Abundance", x = " ") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12))


#Plot of fouled vs non-fouled scallops with bar graph with top 10 families separated into their bar graphs

plot_bar(top10P, x="Factor", fill="Family", facet_grid = ~Family) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + labs(y = "Relative Abundance") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)) 

########################################################
#QUESTION FOR ANDY - When I pull the relative abundance plots I got "Other" and "Unknown" groups in my top 10 families that don't show up in the total abundance plots??? 
########################################################


##############################################################################

###GENUS LEVEL ----
###GENUS LEVEL - TOTAL ABUNDANCE - TOP 9 GENERA

#Sort genera by total abundance and pick the top 9
top9P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Genus"], sum), TRUE)[1:9]

#Cut down the physeq data to only the top 9 genera
top9P = subset_taxa(pseq, Genus %in% names(top9P.names))

#Plot of fouled vs non-fouled scallops with bar graph showing top 10 genera in STACKED bar

plot_bar(top9P, x="Factor", fill="Genus", facet_grid = ~Age) + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + labs(y = "Total Abundance", x = " ") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), strip.text = element_text(size = 12))

#Plot of fouled vs non-fouled scallops with bar graph with top 9 genera separated into their bar graphs

plot_bar(top9P, x="Factor", fill="Genus", facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + labs(y = "Total Abundance", x = " ") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), strip.text = element_text(size = 11))



###GENUS LEVEL - RELATIVE ABUNDANCE - TOP 9 GENERA

#Sort genera by relative abundance and pick the top 9
top9P.names = sort(tapply(taxa_sums(pseq.gen.rel), tax_table(pseq.gen.rel)[, "Genus"], sum), TRUE)[1:9]

#Cut down the physeq data to only the top 9 genera
top9P = subset_taxa(pseq.gen.rel, Genus %in% names(top9P.names))

#Plot of fouled vs non-fouled scallops with bar graph with top 9 genera separated into their bar graphs

plot_bar(top9P, x="Factor", fill="Genus", facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + labs(y = "Relative Abundance", x = " ") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 12), strip.text = element_text(size = 11))



#########################################################################

### PLOT A SINGLE GENUS - Example with Francisella (Total Abundance)

genus_to_plot <- "Francisella"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Total Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot


### PLOT A SINGLE GENUS - Example with Fusobacterium (Total Abundance)

family_to_plot <- "Fusobacteriaceae"
pseq_filtered <- subset_taxa(pseq, Family == family_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Family") +
  geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack") + labs(y = "Total Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot


### PLOT A SINGLE GENUS - Example with Fusobacterium (Relative Abundance)

family_to_plot <- "Fusobacteriaceae"
pseq_filtered <- subset_taxa(pseq.fam.rel, Family == family_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Family") +
  geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack") + labs(y = "Relative Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot



### PLOT ALL GENERA IN ONE FAMILY
# Example with Rhodobacteraceae (Family)

family_to_plot <- "Rhodobacteraceae"
pseq_filtered <- subset_taxa(pseq, Family == family_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot



# Example with Cyanobacteria (Phylum)

phyla_to_plot <- "Cyanobacteria"
pseq_filtered <- subset_taxa(pseq, Phylum == phyla_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Family") +
  geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot



#########################################################################

### PLOTTING HEAT MAPS

#Sort the OTUs by abundance and pick the top 20
top20OTU.names = names(sort(taxa_sums(pseq.fam.rel), TRUE)[1:20])

#Cut down the physeq.tree data to only the top 10 Phyla
top20OTU = prune_taxa(top20OTU.names, pseq.fam.rel)

top20OTU

plot_heatmap(top20OTU)

plot_heatmap(top20OTU, sample.label="Factor", sample.order="Factor")






