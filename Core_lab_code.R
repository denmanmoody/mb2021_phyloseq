#https://microbiome.github.io/tutorials/Core.html
##########CORE
####Goals ----
#1) Subsetting samples
#2) Make stacked bar plots of taxonomy
#3) Calculate Core
#4) Makes stacked bar plots of core
#5) Make a table with the core bacteria (ask Andy)
#6) Search for at least two core in google scholar. 

#Load libraries ----
library(phyloseq)
library(microbiome)
library("devtools")

#Set data ----

Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rare.rds")

pseq <- Marissa_MU42022_rare

### Let's get into specific microbes 

#To calculate core, we need to subset the samples that we are interested in. For some groups there might be multiple groups.  Sub-setting data is really easy.

#Question 1: What samples are you determining the core for? All treatments? Specific treatments? One host? 

#####SUBSETTING DATA ----

New_object_name_SELECT<- subset_samples(pseq, Type == "Anenome")
New_object_name_REMOVE<-subset_samples(pseq, Type != "Anenome")

#Day 1 only ----
  
pseq <- subset_samples(pseq, !Age %in% c("Spat", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Spat only ----

pseq <- subset_samples(pseq, !Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"))

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#remove algae ----

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Larvae only with algae ----

pseq <- subset_samples(pseq, !Age %in% "Spat")

#Larvae only without algae ----

pseq <- subset_samples(pseq, !Age %in% "Spat")

pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#####Make a stacked bar plot of taxonomy

### Make data compositional ----
transform <- microbiome::transform

pseq <- transform(pseq, "compositional")
#Look at a specific level of taxonomy: must be over 1% relative abundance and found in half of samples (otherwise too messy)
pseq <- aggregate_rare(pseq, level = "phylum", detection = 1/100, prevalence = 50/100)
##Plot
Comp<- plot_composition(pseq, otu.sort = "abundance")
Comp
#Question 2: What does your plot show you? 

############Calculate Core###########››
#######Just to name Core
#Notice I'm not looking at all of the samples. So this is from different filtering than the above sample
#Make  compositional 
ps.rel <- microbiome::transform(New_object_name_SELECT, "compositional")
## For below, must be 1% rel abund or more and in 90% of samples. The object that is made tells you the names of the ASVs that fit that criteria (this can be helpful later for your scholar search)
core.taxa.standard_90 <- core_members(ps.rel, detection = 1/100, prevalence = 90/100)

#Question 3: What happens when you call up the last object that we made? 
#Question 4: What happens if you change th prevalence?  Try this from 50% to 100%. Let's draw a plot of this on the board that had # of ASVs on the Y and % prev on the X. 

#Make phyloseq object with core. Then we can make plots
Core1<- core(ps.rel, detection = 1/100, prevalence = 90/100)

#We can plot this:
Comp_core<-plot_composition(Core1, otu.sort = "abundance")

##### naming ASVs.  The psmelt function changes the data from a phyloseq format to regular format - this makes it easier to make plots 
Core2<- psmelt(Core1)

#### To better understand our data we can combine the ASV number (that links to a specific DNA sequence) to the  genus of that particular ASV. The following code does that with the mutate function. 
library(dplyr)
 
 Core3<- mutate(Core2, asv_gen= paste0(OTU, "-",genus))
##Plot core




#######Make plot with ggplot with better ASV names_ you can also use any of your meta data to collapse samples into groups. 
 head(Core3)
Best_plot<- ggplot(Core3, aes(x=Sample, y = Abundance, fill= asv_gen)) +geom_bar(stat="identity")

##Collapsed 
Best_plot_collapsed<- ggplot(Core3, aes(x=Type, y = Abundance, fill= asv_gen)) +geom_bar(stat="summary")
###Question 5: Not a question, but: Save a plot or multiple plots for your poster!
###Question 6: Search for at least two core in google scholar? Do those bacteria show up in  your search? If so where? It is worth spending a few minutes trying to figure out where the bacteria in your core are found. The environment? Other hosts? Same host? 