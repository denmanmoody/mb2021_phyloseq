#plotting alpha diveristy 
#Sept 28th, 2023

#source: https://microbiome.github.io/tutorials/PlotDiversity.html

#Packages ----

library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)

#Get data ----

setwd("C:/Users/maris/OneDrive/Documents/GitHub/mb2021_phyloseq")

  
Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rare.rds")

pseq <- Marissa_MU42022_rare

#objects

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#Remove samples ----
#for treatment comparisons (below) need to remove "NA" values (algae)

pseq_filtered <- subset_samples(pseq, !Organism %in% "Algae")


#alpha diversity ----

ps1 <- prune_taxa(taxa_sums(pseq_filtered) > 0, pseq)

tab <- microbiome::alpha(ps1, index = "all")
kable(head(tab))

#below code require "Treatment" to be a factor rather than a character
#check if your variable is a character or factor

str(Metadata)

#Change chr to factor ----

Metadata$Treatment.fact <- as.factor(Metadata$Treatment)

Metadata$Age.fact <- as.factor(Metadata$Age)

str(Metadata)



#look for meta

ps1.meta <- meta(ps1)
kable(head(ps1.meta))

#combine meta and alpha diversity

ps1.meta$Shannon <- tab$diversity_shannon 
ps1.meta$InverseSimpson <- tab$diversity_inverse_simpson

#Treatment 

##### create a list of pairwise comparisons
Treatment1 <- levels(ps1.meta$Treatment)

Treatment <- levels(Metadata$Treatment.fact)
print(Treatment) 

Treatment.pairs <- combn(seq_along(Treatment), 2, simplify = FALSE, FUN = function(i)Treatment[i])
print(Treatment.pairs)

####Violin plots ----

p1 <- ggviolin(ps1.meta, x = "Treatment", y = "Shannon", add = "boxplot", 
               fill = "Treatment", 
               palette = c("#a6cee3", "#b2df8a", "#fdbf6f", "pink", "plum4"),
               ggtheme = theme_pubr(),
               font.label = 14)

p1 <- p1 + scale_x_discrete(labels = c("Control", "Probiotics", "Probiotics + HT", "High temperature (HT)", "Algae"))

p1 <- p1 + theme(legend.position = "none")

print(p1)

p1 <- p1 + stat_compare_means(comparisons = Treatment.pairs) 
print(p1)


#Age ----

Age <- levels(ps1.meta$Age)

Age <- levels(Metadata$Age.fact)
print(Age) 

Age.pairs <- combn(seq_along(Age), 2, simplify = FALSE, FUN = function(i)Age[i])
print(Age.pairs)

####Violin plots ----

p1 <- ggviolin(ps1.meta, x = "Age", y = "Shannon", add = "boxplot", 
fill = "Age", 
palette = c("#a6cee3", "#b2df8a", "#fdbf6f", "pink", "plum4"),
ggtheme = theme_pubr(),
font.label = 14)

p1 <- p1 + scale_x_discrete(limits = c("Day 01", "Day 03", "Day 06", "Day 15", "Spat"))


print(p1)

p1 <- p1 + stat_compare_means(comparisons = Age.pairs) 
print(p1)


#other alpha diversity indices ----

tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))

tab <- richness(pseq)
kable(head(tab))

tab <- dominance(pseq, index = "all")
kable(head(tab))

dominant(pseq)

tab <- rarity(pseq, index = "all")
kable(head(tab))

tab <- coverage(pseq, threshold = 0.5)
kable(head(tab))

tab <- core_abundance(pseq, detection = .1/100, prevalence = 90/100)

tab <- inequality(pseq)


tab <- evenness(pseq, "all")
kable(head(tab))


#Visualization ----

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Treatment"
                         )

p.shannon <- p.shannon + theme_classic() + 
  labs(x="", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Age"
)

p.shannon <- p.shannon + theme_classic() + 
  labs(x="", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Sample.type"
)

p.shannon <- p.shannon + theme_classic() + 
  labs(x="", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))
        
fill.colors = c(Control="cyan4",High.temperature = "deeppink4", Probiotics="forestgreen", Probiotics+HT= "slateblue2", NA="goldenrod1"))


palette = c("#a6cee3", "#b2df8a", "#fdbf6f", "pink", "plum4"),

p.shannon 


#Test significance (alpha diversity) ----

d <- meta(pseq)
d$diversity <- microbiome::diversity(pseq, "shannon")$shannon
# Split the values by group
spl <- split(d$diversity, d$Treatment)
# Kolmogorov-Smironv test
pv <- ks.test(spl$Control, spl$"Probiotics + HT", spl$"High temperature", spl$"Probiotics", spl$"NA")$p.value
# Adjust the p-value
padj <- p.adjust(pv)

#p value = 0.419

d <- meta(pseq)
d$diversity <- microbiome::diversity(pseq, "shannon")$shannon
# Split the values by group
spl <- split(d$diversity, d$Treatment)
# Kolmogorov-Smironv test
pv <- ks.test(spl$"Probiotics", spl$"High temperature")$p.value
# Adjust the p-value
padj <- p.adjust(pv)

#p value = 0.064 * most signif difference between treatments (alpha diversity)


# Split the values by group
spl <- split(d$diversity, d$Sample.type)
# Kolmogorov-Smironv test
pv <- ks.test(spl$"Algae", spl$"Larvae", spl$"Spat")$p.value
# Adjust the p-value
padj <- p.adjust(pv)

#pvalue = 0.01268 (algae, larvae, spat)

pv <- ks.test(spl$"Larvae", spl$"Spat")$p.value
# Adjust the p-value
padj <- p.adjust(pv)

#p value = 7.05e-14 (larvae and spat)


#Day 1 only ----
#Since high temp treatment is missing spat samples = diversity skewed lower - look at day 1 diversity between treatments

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Day 01")

p.shannon1 <- boxplot_alpha(pseq_Day1, 
                           index = "shannon",
                           x_var = "Treatment"
)

p.shannon1 <- p.shannon1 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 1") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon1 <- p.shannon1 + theme(legend.position = "none")

p.shannon1

#Day 3 only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Day 03")

p.shannon3 <- boxplot_alpha(pseq_Day1, 
                           index = "shannon",
                           x_var = "Treatment"
)

p.shannon3 <- p.shannon3 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 3") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon3 <- p.shannon3 + theme(legend.position = "none")

p.shannon3

#Day 6 only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Day 06")

p.shannon6 <- boxplot_alpha(pseq_Day1, 
                           index = "shannon",
                           x_var = "Treatment"
)

p.shannon6 <- p.shannon6 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 6") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon6 <- p.shannon6 + theme(legend.position = "none")

p.shannon6

#Day 15 only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Day 15")

p.shannon15 <- boxplot_alpha(pseq_Day1, 
                           index = "shannon",
                           x_var = "Treatment"
)

p.shannon15 <- p.shannon15 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 15") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon15 <- p.shannon15 + theme(legend.position = "none")

p.shannon15

#Spat only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Spat")

p.shannons <- boxplot_alpha(pseq_Day1, 
                           index = "shannon",
                           x_var = "Treatment"
)

p.shannons <- p.shannons + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Spat") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))


p.shannons <- p.shannons + theme(legend.position = "none")

p.shannons

#Align plots into grid ----

install.packages("gridExtra")
library(gridExtra)

grid.arrange(p.shannon1, p.shannon3, p.shannon6, p.shannon15, p.shannons, ncol = 2)
