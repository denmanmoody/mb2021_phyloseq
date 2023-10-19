###top phylum/family/genera within microbiome data
##bar plots and violin plots

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)

#Load data ----
setwd("/Users/denmanmoody/Documents/GitHub/mb2021_phyloseq")

pseq<- readRDS("Denman_samples.rds")


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



#source for below: https://rstudio-pubs-static.s3.amazonaws.com/330760_8bba830836324bf6b100d4e76f49e3d2.html#other_visualizations 


#Bar plots ----

#PHYLUM LEVEL ----

##STARTING AT PHYLUM LEVEL - TOTAL ABUNDANCE

plot_bar(pseq, fill="Phylum")

#Sort the Phyla by abundance and pick the top 5

top5P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Phylum"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 10 Phyla

top5P = subset_taxa(pseq, Phylum %in% names(top5P.names))

#Plot

plot_bar(top5P, x="Age", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(top5P, x="Age", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

##code to change x-axis labels

plot <- plot_bar(top5P, x = "Age", fill = "Phylum") +
  geom_boxplot(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")

plot + scale_x_discrete(labels = c("Day 1", "Day 18", "Day 3", "Spat"))




#Plot

plot_bar(top5P, x="Treatment", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")


plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae"))

plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + 
  geom_boxplot(aes(color=Phylum, fill=Phylum), 
  stat="identity", 
  position="stack") 
+ theme_bw() 
+ theme(panel.grid = element_blank()) 
+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
+ theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)   
)



#Phylum relative abundance ----

###PHYLUM LEVEL - instead of using pseq (total abundance) plot with pseq.fam.rel (relative abundance at Family level)

top5P.names = sort(tapply(taxa_sums(pseq.phy.rel), tax_table(pseq.phy.rel)[, "Phylum"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 10 Phyla

top5P = subset_taxa(pseq.phy.rel, Phylum %in% names(top5P.names))


  labs(y = "Relative Abundance")

  plot_bar(top5P, x="Treatment", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
  
  
  plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12)   
  ) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Relative abundance", x = " ")


  
#FAMILY LEVEL ----

#Sort the Phyla by abundance and pick the top 5

top5P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Family"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 5 Family

top5P = subset_taxa(pseq, Family %in% names(top5P.names))

#Plot

plot_bar(top5P, x="Treatment", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)   
)plot_

#FAMILY LEVEL - plot by bar by treatment and plot relative abundance

top5P.names = sort(tapply(taxa_sums(pseq.fam.rel), tax_table(pseq.fam.rel)[, "Family"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 5 family

top5P = subset_taxa(pseq.fam.rel, Family %in% names(top5P.names))

plot_bar(top5P, x="Treatment", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Relative abundance", x = " ")



#GENUS LEVEL----
#- sort by treatment

top15P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Genus"], sum), TRUE)[1:15]

#Cut down the physeq data to only the top 15 Genus

top15P = subset_taxa(pseq, Genus %in% names(top15P.names))

#Plot


plot_bar(top15P, x="Treatment", fill="Genus", facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

#plot specific family

family_to_plot <- "Rhodobacteraceae"
pseq_filtered <- subset_taxa(pseq, Family == family_to_plot)

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ", title = "Rhodobacteraceae") +
  theme(plot.title = element_text(hjust = 0.5))


##PLOT SPECIFIC GENUS 

###Flavobacterium ----

# Filter data for only "Flavobacterium" genus
genus_to_plot <- "Flavobacterium"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

# Create the bar plot for "Flavobacterium" abundance between different "Tank_treatments"
plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack")

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")


###Celeribacter ----

genus_to_plot <- "Celeribacter"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

# Create the bar plot for "Flavobacterium" abundance between different "Tank_treatments"
plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack")

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

###Phaeobacter ----

genus_to_plot <- "Phaeobacter"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

# Create the bar plot for "Flavobacterium" abundance between different "Tank_treatments"
plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack")

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

###Vibrio ----

# Filter data for only "Vibrio" genus
genus_to_plot <- "Vibrio"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

#Denman's

plot <- plot_bar(top5P, x = "Factor", fill = "Family") +
  geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack") 

print(plot)


plot <- plot_bar(top5P, x = "Treatment", fill = "Family", facet_grid = ~Family) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme_bw()

plot + labs(y = "Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14), panel.grid = element_blank(),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank()) + theme_bw()

plot + scale_x_discrete(labels = c("Antibiotics", "Antibiotics + HT", "Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14), panel.grid = element_blank(),
                                                                                                                                                                                         panel.grid.major = element_blank(),
                                                                                                                                                                                         panel.grid.minor = element_blank()) + theme_bw()

plot + labs(y = "Relative Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), axis.text.y = element_text(size = 14), axis.title = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14), panel.grid = element_blank(),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank()) + theme_bw()

plot <- plot + scale_x_discrete(labels = c("Antibiotics", "Antibiotics + HT", "Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(plot)


###GENUS LEVEL - sort by treatment - only thing is that for Vibrio HS has outlier

top15P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Genus"], sum), TRUE)[1:15]

#Cut down the physeq data to only the top 15 Genus

top15P = subset_taxa(pseq, Genus %in% names(top15P.names))

#Plot

plot_bar(top15P, x="Treatment", fill="Genus", facet_grid=~Age) + geom_bar(aes(color=Genus, fill=Genus), stat="identity", facet_grid=~Age) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot <- plot + scale_x_discrete(labels = c("Antibiotics", "Antibiotics + HT", "Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(plot)


#denman's

plot <- plot_bar(top15P, x = "Factor", fill = "Genus", facet_grid = ~Genus)

plot <- plot + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(plot)


###################
  
  ##PLOT SPECIFIC GENUS 
  
  #START WITH FLAVOBACTERIUM
  
  # Filter data for only "Flavobacterium" genus - make sure name matches hot bacteria is spelled in silva database
  genus_to_plot <- "Flavobacterium"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus", facet_grid=~Age) +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot


##pseudoaltermonas

# Filter data for only "Pseudoalteromonas" genus
genus_to_plot <- "Pseudoalteromonas"
pseq_filtered <- subset_taxa(pseq_gen, Genus == genus_to_plot)

#pick colors

custom_color_palette <- ("Pseudoalteromonas" = "#1E90FF")

plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus", facet_grid=~Age) +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot


##rhodo

# Filter data for only Rhodo genus
# Code no longer works after OTU table was rarefied 

genus_to_plot <- "Phaeobacter"
pseq_filtered <- subset_taxa(pseq_gen, Genus == genus_to_plot)

plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus", facet_grid=~Age) +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot



##Vibrio

# Filter data for only "Vibrio" genus
genus_to_plot <- "Vibrio"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus", facet_grid=~Age) +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot


#Denman specific genus

#Francisella

genus_to_plot <- "Francisella"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot

#family level plot - Vibrio

family_to_plot <- "Vibrio"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot







##to plot all genera within a family
#using example of Rhodo

family_to_plot <- "Rhodobacteraceae"
pseq_filtered <- subset_taxa(pseq, Family == family_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot

#try with cyanobacteria (phylum level)

phyla_to_plot <- "Cyanobacteria"
pseq_filtered <- subset_taxa(pseq, Phylum == phyla_to_plot)

plot <- plot_bar(pseq_filtered, x = "Factor", fill = "Family") +
  geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack") + labs(y = "Abundance", x = "") + theme_bw(base_size = 17) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_text(size = 17), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 17), legend.text = element_text(size = 17), legend.title = element_text(size = 17), panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot



###NEXT SECTION - HEATMAPS (TOP ## OTU'S)

#Sort the OTUs by abundance and pick the top 20
top20OTU.names = names(sort(taxa_sums(pseq.fam.rel), TRUE)[1:20])

#Cut down the physeq.tree data to only the top 10 Phyla
top20OTU = prune_taxa(top20OTU.names, pseq.fam.rel)

top20OTU

plot_heatmap(top20OTU)

plot_heatmap(top20OTU, sample.label="Age", sample.order="Age")

plot_heatmap(top20OTU, sample.label="Age", sample.order="Age", taxa.label="Phylum", taxa.order="Family", low="white", high="purple", na.value="grey")

plot_heatmap(top20OTU, "NMDS", "bray", title="Bray-Curtis")


plot_heatmap(top20OTU, sample.label="Treatment", sample.order="Treatment")


plot_heatmap(top20OTU, "NMDS", "bray", title="Bray-Curtis", ample.label="Treatment", sample.order="Treatment")


##I want to look at heatmap of only larval samples

# Exclude specific data from top20OTU object

top20OTU = prune_taxa(top20OTU.names, pseq.fam.rel)

top20OTU

excluded_samples <- c("T10r1", "T10r2", "T10r3", "T11r1", "T11r3", "T12r1", "T12r2", "T12r3", "T13r1","T13r2","T13r3", "T14r1", "T14r2", "T15r1", "T15r2", "T16r1", "T16r2", "T16r3", "T1r1", "T1r2", "T1r3", "T2r1", "T2r3", "T3r1", "T3r2", "T3r3", "T4r1", "T4r2", "T9r1", "T9r3", "T11r2", "T14r3", "T2r2")

top20OTU <- prune_samples(!(sample_names(top20OTU) %in% excluded_samples), top20OTU)

plot_heatmap(top20OTU, sample.label="Tank_treatment", sample.order="Tank_treatment")

