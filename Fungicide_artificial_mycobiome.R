#Code written by M. R. Aggerbeck. 
#mrag@envs.au.dk
#Creative commons, 2018. 

###----------------PACKAGES--------------
#Packages used in the following code. 
library("phyloseq"); packageVersion("phyloseq")
library("biomformat"); packageVersion("biomformat")
library("vegan"); packageVersion("vegan")
library("ggplot2"); packageVersion("ggplot2")
library("DESeq2"); packageVersion("DESeq2")
library("mvabund"); packageVersion("mvabund")
library("metacoder"); packageVersion("metacoder")
library("taxa"); packageVersion("taxa")
library("microbiomeSeq"); packageVersion("microbiomeSeq")
library("adespatial"); packageVersion("adespatial")
library("ggpubr"); packageVersion("ggpubr")
library("data.table"); packageVersion("data.table")
library("igraph"); packageVersion("igraph")
library("dplyr"); packageVersion("dplyr")
library("tidyr"); packageVersion("tidyr")
library("plotrix"); packageVersion("plotrix")
library("microbiome"); packageVersion("microbiome")


#This page of code is using Phyloseq as the primary package, and uses the others above as additional analyses. 
#Check out the Phyloseq tutorial and the various package manuals for more info on analyses. 

#Set working directory
setwd("D:/Marie/Documents/PhD/Alex_statistics/")
#Set user directory
uzdir <- "D:/Marie/Documents/PhD/Alex_statistics/"

###---------------------LOAD DATA-------------------------------------

#Biom file to import - preferably JSON
#Dataset 1
data_biom <- "feature-table_new_primer(1).biom"
#Dataset 2
#data_biom <- "dataset2.biom"

#Import data, removing any blank spaces in file path
biom_file <- paste(uzdir, data_biom, sep="")
#the actual act of importing the file into R
biom_otu_tax <- import_biom(biom_file)

#Import taxonomy table
tax1 <- as.matrix(read.table(file = 'dataset1_tax_name_un.csv', sep = '\t', header = TRUE))
#Use first row of table as rownames
rownames(tax1) <- tax1[,1]
tax1 <- tax1[,2:10]
head(tax1)

#Convert table to a Phyloseq taxonomy table
TAX = tax_table(tax1)
#Inspect
head(TAX)
#remove tax1 to declutter environment
rm(tax1)

#Combine OTU table and taxonomy table into a phyloseq object. 
physeq <- phyloseq(biom_otu_tax, TAX)
#Inspect
head(tax_table(physeq))
head(otu_table(physeq))

rank_names(tax_table(physeq))

taxa_names(physeq) <- tax_table(physeq)[,1]; head(tax_table(physeq)) 
rank_names(tax_table(physeq))
rownames(otu_table(physeq))[1:10]

#Import metadata file, removing blank spaces in path 
md_file <- paste(uzdir, "Metadata_ino.csv", sep="")   

#Has to be tsv. If csv, use code below. 
metadata <- import_qiime_sample_data(md_file)
metadata

#If the above doesn't produce a table, use this: 
metadata <- read.csv(md_file, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

#change names in metadata to fit proper nomenclature
temp <- as.matrix(metadata)
temp <- gsub("Mycobioma", "Skopobiota", temp)
temp <- gsub("Myco", "Skopo", temp)
temp <- gsub("Patogeno", "Pathogen", temp)
temp
metadata <- temp
temp <- colnames(metadata)
temp <- gsub("Mycobioma", "Skopobiota", temp)
temp <- gsub("Myco", "Skopo", temp)
temp <- gsub("Patogeno", "Pathogen", temp)
temp
colnames(metadata) = temp
metadata <- as.data.frame(metadata)

metadata <- sample_data(metadata)
rownames(metadata) = metadata$SampleID
head(metadata)

#Merge OTU data and metadata to single object - The Phyloseq Object(tm)
giovanni <- merge_phyloseq(physeq,metadata)

sample_names(giovanni) <- gsub("-", "_", sample_names(giovanni))
giovanni

#giovanni is now the s4 object - the phyloseq object - we will use throughout the rest of the code. 
#In theory, all other objects in the global environment could be removed now, but I tend to leave them, just in case I need them.

head(sample_data(giovanni))
head(tax_table(giovanni))
head(otu_table(giovanni))

#remove sample #46 in dataset 1
giovanni = subset_samples(giovanni, SampleID != "GF-T-ITS-46")
giovanni = prune_taxa(taxa_sums(giovanni) > 0, giovanni)

raw_data <- giovanni
gio_sp <- giovanni
rank_names(gio_sp)
tax_table(gio_sp) <- tax_table(gio_sp)[,-1]; head(tax_table(gio_sp))
rank_names(gio_sp)
gio_sp <- tax_glom(gio_sp, taxrank="X"); gio_sp
length(unique(taxa_names(gio_sp)))
taxa_names(gio_sp) <- tax_table(gio_sp)[,8]; taxa_names(gio_sp)[1:10]
tax_table(gio_sp) <- tax_table(gio_sp)[,-8]; head(tax_table(gio_sp))

tax_table(gio_sp) <-  gsub("s__unidentified", "", tax_table(gio_sp)); 
tax_table(gio_sp) <-  gsub("g__unidentified", "", tax_table(gio_sp)); 
tax_table(gio_sp) <-  gsub("f__unidentified", "", tax_table(gio_sp));
tax_table(gio_sp) <-  gsub("o__unidentified", "", tax_table(gio_sp));
tax_table(gio_sp) <-  gsub("c__unidentified", "", tax_table(gio_sp));
tax_table(gio_sp) <-  gsub("p__unidentified", "", tax_table(gio_sp));

sample_data(gio_sp)[1:5]
otu_table(gio_sp) [1:5]
tax_table(gio_sp)[1:5]

giovanni <- gio_sp

###------------------Filter out low-abundance reads-----------------

gio2 <- gio_sp
otu_table(gio2)[otu_table(gio2)<25 ] <- 0
gio2 = prune_taxa(taxa_sums(gio2) > 0, gio2)

###------------------SUBSETTING-------------------------------------

# #This for loop creates four subset from the fungicide treatments. 
# #We create a list of the names of the fungicides, and then use that to specify the subset command. 
# #Then we paste the name of the fungicide onto the subset and assign it as a separate phyloseq object.
# ID<- as.vector(unique(sample_data(giovanni)$Fungicide))
# for (i in 1:length(ID)){ 
#   temp <- subset_samples(giovanni, Fungicide==ID[i])
#   temp = prune_taxa(taxa_sums(temp) > 0, temp)
#   name <- paste(ID[i])
#   assign(name, temp)
# }
# 
# #Creates another set of subset, based on inoculation
# ID<- as.vector(unique(sample_data(giovanni)$Inoculation))
# for (i in 1:length(ID)){ 
#   temp <- subset_samples(giovanni, Inoculation==ID[i])
#   temp = prune_taxa(taxa_sums(temp) > 0, temp)
#   name <- paste(ID[i])
#   assign(name, temp)
# }

## Or you can make subsets of treatments. Specifies which metadata column, and which value in said column should be criteria. 
# No_ino <- subset_samples(giovanni, No_inoculation!="No")
# Control <- subset_samples(giovanni, Fungicide=="Control")
Cont_no_ino <- subset_samples(giovanni, No_inoculation=="No_inoculationC")

###----------------------Merge samples----------------


test = prune_taxa(taxa_sums(giovanni) > 0, giovanni)
test2 <- as.data.frame(unique(sample_data(test)[,4]))
groups_gio <- as.vector(as.character(test2$Treatment))
sample_data(test)$group_gio <- get_variable(test, "Treatment") %in% groups_gio

test2 <- as.data.frame(unique(sample_data(test)[,6]))
fungicide <- as.vector(as.character(test2$Fungicide))

test2 <- as.data.frame(unique(sample_data(test)[,13]))
inoculation <- as.vector(as.character(test2$Inoculation))
merged_sp = merge_samples(test, "Treatment")

merged_sp
sample_names(test)
sample_names(merged_sp)

sample_data(merged_sp)$Treatment = sample_names(merged_sp)
sample_data(merged_sp)$group_gio = sample_names(merged_sp) %in% groups_gio
merged_sp

test <- sample_data(merged_sp)$Fungicide
test <- gsub(1, "Blad", test); test <- gsub(2, "Control", test); test <- gsub(3, "CuS", test); test <- gsub(4, "Systemics", test)
sample_data(merged_sp)$Fungicide <- test

test <- sample_data(merged_sp)$Inoculation
test <- gsub(3, "Pathogen + Skopobiota", test); test <- gsub(4, "Skopobiota", test); test <- gsub(2, "Pathogen", test); test <- gsub(1, "Water", test)
sample_data(merged_sp)$Inoculation <- test
sample_data(merged_sp)

#gio_m_un <- merged_sp
#tax_table(gio_m_un) <-  gsub("s__unidentified", "", tax_table(gio_m_un)); 
#tax_table(gio_m_un) <-  gsub("g__unidentified", "", tax_table(gio_m_un)); 
#tax_table(gio_m_un) <-  gsub("f__unidentified", "", tax_table(gio_m_un));
#tax_table(gio_m_un) <-  gsub("o__unidentified", "", tax_table(gio_m_un));
#tax_table(gio_m_un) <-  gsub("c__unidentified", "", tax_table(gio_m_un));
#tax_table(gio_m_un) <-  gsub("p__unidentified", "", tax_table(gio_m_un));
#tax_table(gio_m_un)

#gio_un <- giovanni
#tax_table(gio_un) <-  gsub("s__unidentified", "", tax_table(gio_un)); 
#tax_table(gio_un) <-  gsub("g__unidentified", "", tax_table(gio_un)); 
#tax_table(gio_un) <-  gsub("f__unidentified", "", tax_table(gio_un));
#tax_table(gio_un) <-  gsub("o__unidentified", "", tax_table(gio_un));
#tax_table(gio_un) <-  gsub("c__unidentified", "", tax_table(gio_un));
#tax_table(gio_un) <-  gsub("p__unidentified", "", tax_table(gio_un));
#tax_table(gio_un)

###Merge unassigned taxa

sort(taxa_names(giovanni))

#IMPORTANT! Input manually any unwanted higher classes from the list produced by the sort() above
badtaxa = c("k__Fungi", "o__Agaricales", "o__Coniochaetales", "o__Pleosporales", "p__Ascomycota", "p__Basidiomycota",
            "c__Agaricomycetes", "c__Dothideomycetes", "c__Leotiomycetes", "c__Microbotryomycetes", "c__Sordariomycetes", 
            "c__Tremellomycetes")

gio_un_con <- merge_taxa(giovanni, badtaxa)
taxa_names(gio_un_con) <-  gsub("o__Agaricales", "Unassigned", taxa_names(gio_un_con)) 
head(tax_table(gio_un_con))
head(sample_data(gio_un_con))

sort(taxa_names(merged_sp))
badtaxa = c("k__Fungi", "o__Agaricales", "o__Coniochaetales", "o__Pleosporales", "p__Ascomycota", "p__Basidiomycota",
            "c__Agaricomycetes","c__Microbotryomycetes", "c__Sordariomycetes", "c__Tremellomycetes")

gio_m_un_con <- merge_taxa(merged_sp, badtaxa)
taxa_names(gio_m_un_con) <-  gsub("o__Agaricales", "Unassigned", taxa_names(gio_m_un_con)) 
head(tax_table(gio_m_un_con))
head(sample_data(gio_m_un_con))

sort(taxa_names(gio_m_un_con))

###----------Filter by relative abundance-----------

#Filters out anything below 0.01% of total abundance. 
minTotRelAbun = 0.1 
x = taxa_sums(giovanni)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
gio_filter1 = prune_taxa(keepTaxa, giovanni)

y = as.matrix(sort((x / sum(x))*100, decreasing = TRUE))

x = taxa_sums(merged_sp)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
gio_m_filter1 = prune_taxa(keepTaxa, merged_sp)


###-----------------Colour palettes-----------------------

# I will from now on use a colour scheme that

set_2 <- c("purple", "orange")
set_3 <- c("royalblue", "grey95","tomato")
set_ino <- c("cyan", "green","magenta", "gold")
set_fun <- c( "red", "violet", "brown","darkblue")


###----------------------Core biome----------------

write.table(otu_table(giovanni)[2], file='presence_pch.tsv', quote=FALSE, sep='\t')
write.table(y, file='taxalist_rel_ab.tsv', quote=FALSE, sep='\t')

core.taxa.standard <- core_members(Cont_no_ino, detection = 0, prevalence = 80/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_control_no_inoculation_80.tsv', quote=FALSE, sep='\t')

core.taxa.standard <- core_members(Cont_no_ino, detection = 0, prevalence = 50/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_control_no_inoculation_50.tsv', quote=FALSE, sep='\t')

# core.taxa.standard <- core_members(No_ino, detection = 0, prevalence = 80/100)
# core.taxa.standard
# write.table(core.taxa.standard, file='core_no_inoculation_80.tsv', quote=FALSE, sep='\t')
# 
# core.taxa.standard <- core_members(Control, detection = 0, prevalence = 80/100)
# core.taxa.standard
# write.table(core.taxa.standard, file='core_control_80.tsv', quote=FALSE, sep='\t')
# 
# core.taxa.standard <- core_members(Cont_no_ino, detection = 0, prevalence = 80/100)
# core.taxa.standard
# write.table(core.taxa.standard, file='core_control_no_inoculation_80.tsv', quote=FALSE, sep='\t')
# 
# core.taxa.standard <- core_members(No_ino, detection = 0, prevalence = 50/100)
# core.taxa.standard
# write.table(core.taxa.standard, file='core_no_inoculation_50.tsv', quote=FALSE, sep='\t')
# 
# core.taxa.standard <- core_members(Control, detection = 0, prevalence = 50/100)
# core.taxa.standard
# write.table(core.taxa.standard, file='core_control_50.tsv', quote=FALSE, sep='\t')
# 
# core.taxa.standard <- core_members(Cont_no_ino, detection = 0, prevalence = 50/100)
# core.taxa.standard
# write.table(core.taxa.standard, file='core_control_no_inoculation_50.tsv', quote=FALSE, sep='\t')

###----------------------PERMANOVA---------------------------

#Transform data to relative abundance ("compositional")
pseq.rel <- microbiome::transform(gio_filter1, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

p <- microbiome::plot_landscape(pseq.rel, method = "NMDS", distance = "jaccard", col = "Inoculation", size = 3)
print(p)

permanova <- vegan::adonis(t(otu)~Fungicide+Inoculation+Treatment,
                           data = meta, permutations=999, method = "jaccard")
write.table(as.data.frame(permanova$aov.tab), file='permanova_jaccard_all.tsv', quote=FALSE, sep='\t')

print(as.data.frame(permanova$aov.tab))

dist <- vegan::vegdist(t(otu))
temp <- anova(vegan::betadisper(dist, meta$Fungicide))
temp
write.table(temp, file='anova_fungicide.tsv', quote=FALSE, sep='\t')

#Under construction
coef <- coefficients(permanova)
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")


###-----------RICHNESS-----------

#This calculates richness for many difference diversity indices (check help file for specifics),
# but if calculated using microbiomeseq below, you get p-values on the figure as well. 
#Provides a combined plot (ggplot2) of three diversity incides and evenness (InvSimpson)
#the geom_boxplot adds median and std handles to the plot. 

#giovannirplot <- plot_richness(giovanni, x="Inoculation", color="Fungicide", measures=c("Shannon", "Simpson", "InvSimpson"), title = "Richness, giovanni dataset")
#giovannirplot + geom_boxplot(data = giovannirplot$data, aes(color = NULL), alpha = 0.1)

p_an <-plot_anova_diversity_pval(giovanni, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="Treatment",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_an_o
ggsave(p_an_o, file="richness_treatment.pdf")
write.table(p_an_pvalues, file='richness_treatment_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_o$data, file='richness_treatment_data.tsv', quote=FALSE, sep='\t')

#inoculation
p_an <-plot_anova_diversity_pval(giovanni, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="Inoculation",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p <- p_an_o+scale_colour_manual(values = set_ino )
p
ggsave(p, file="richness_inoculation.pdf")
write.table(p_an_pvalues, file='richness_inoculation_pval.tsv', quote=FALSE, sep='\t')
write.table(p$data, file='richness_inoculation_data.tsv', quote=FALSE, sep='\t')

#fungicide
p_an <-plot_anova_diversity_pval(giovanni, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="Fungicide",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p <- p_an_o+scale_colour_manual(values = set_fun )
p
ggsave(p, file="richness_fungicide.pdf")
write.table(p_an_pvalues, file='richness_fungicide_pval.tsv', quote=FALSE, sep='\t')
write.table(p$data, file='richness_fungicide_data.tsv', quote=FALSE, sep='\t')

#p1 <- plot_richness(Blad, x="Treatment", color="Inoculation", measures=c("Chao1", "Shannon", "Simpson", "InvSimpson"))
#p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)
#Blad_richness <- p1 + ggtitle("Richness, Blad")
#Blad_richness

###-----------Heatmaps------------

temp <- gio_m_un_con
tmp <- prune_taxa(names(sort(taxa_sums(temp),TRUE)[1:40]), temp)
p <- plot_heatmap(tmp, method="NMDS", distance="jaccard", sample.label="Treatment", 
                  sample.order = "Inoculation", title="40 Most Abundant Taxa") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank())
print(p)
ggsave(p, file="Heatmap_40_ino.pdf", width = 20, height = 25, units = "cm")

p <- plot_heatmap(tmp, method="NMDS", distance="jaccard", sample.label="Treatment", 
                  sample.order = "Fungicide", title="40 Most Abundant Taxa") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank())
print(p)
ggsave(p, file="Heatmap_40_fun.pdf", width = 20, height = 25, units = "cm")


###-----------Microbiomeseq------

#this package holds a lot of promise, but is still in beta, so if you use it, contact the author on how to reference it. 
#It implements a ton of statistical tools and returns nice figures, but is a little buggy. 

#Some of the commands in this packae requires a transposed matrix, so first, we'll transpose our phyloseq obejct:
# Extract abundance matrix from the phyloseq object
OTU2 <- otu_table(giovanni, taxa_are_rows = TRUE)
OTU2 <- t(OTU2)
#Merge new phyloseq object
TAX = as.matrix(tax_table(giovanni))
met <- as.data.frame(sample_data(giovanni))
pseq1 <- phyloseq(OTU2, TAX)
MBS_physeq<-merge_phyloseq(pseq1, met)

#OTU <- otu_table(MBS_Pch)
#TAX <- tax_table(MBS_Pch)
#met <- sample_data(MBS_Pch)

#write.table(OTU, file='otu_table.tsv', quote=FALSE, sep='\t')
#write.table(TAX, file='tax_table.tsv', quote=FALSE, sep='\t')
#write.table(met, file="sam_data.tsv", quote=FALSE, sep='\t')

#MBS_physeq is now the transposed version of giovanni. 

rm(OTU2, TAX, met, pseq1)

#subset transposed phyloseq object into inoculations
ID<- as.vector(unique(sample_data(MBS_physeq)$Inoculation))
for (i in 1:length(ID)){ 
  temp <- subset_samples(MBS_physeq, Inoculation==ID[i])
  temp = prune_taxa(taxa_sums(temp) > 0, temp)
  name <- paste("MBS_", ID[i], sep = "")
  assign(name, temp)
}
MBS_Pch

### RUN test_functions.R!!! 

#Pielou's evenness = Shannon evenness. 
p_an <-plot_anova_diversity_pval(MBS_physeq, method = c("simpson", "shannon", "evenness"),grouping_column ="Treatment",pValueCutoff=0.05)
p_an_o <- p_an[[1]]
p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
print(p_an_o)
#write.table(p_an_pvalues, file='richness_o_pval.tsv', quote=FALSE, sep='\t')
#write.table(p_an_o$data, file='richness_o_data.tsv', quote=FALSE, sep='\t')

p_an <-plot_anova_diversity_pval(MBS_Pch_Skopobiota, method = c("shannon", "evenness"),grouping_column ="Fungicide",pValueCutoff=0.05)
p_an_pm <- p_an[[1]]
p_an_pvalues <- p_an[[2]]
p_an_pm <- p_an_pm +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=set_fun)+
  ggtitle("Pathogen + Skopobiota")+
  theme(plot.title = element_text(hjust = 0.5))
write.table(p_an_pvalues, file='richness_pm_fun_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_pm$data, file='richness_pm_fun_data.tsv', quote=FALSE, sep='\t')

p1 <- plot_richness(MBS_Pch, x="Treatment", color="Fungicide", measures=c("Shannon", "InvSimpson"))+
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=set_fun)
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = Fungicide), alpha = 0.1)
#Pch_richness <- p1 + ggtitle("Richness, Pch")
#Pch_richness
p_an <-plot_anova_diversity_pval(MBS_Pch, method = c("shannon", "evenness"),grouping_column ="Fungicide",pValueCutoff=0.05)
p_an_p <- p_an[[1]]
p_an_pvalues <- p_an[[2]]

p_an_p <- p_an_p +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=set_fun)+
  ggtitle("Pathogen")+
  theme(plot.title = element_text(hjust = 0.5))
write.table(p_an_pvalues, file='richness_p_fun_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_p$data, file='richness_p_fun_data.tsv', quote=FALSE, sep='\t')
p_an_p

p_an <-plot_anova_diversity_pval(MBS_Skopobiota, method = c("shannon", "evenness"),grouping_column ="Fungicide",pValueCutoff=0.05)
p_an_m <- p_an[[1]]
p_an_pvalues <- p_an[[2]]
p_an_m <- p_an_m +  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=set_fun)+
  ggtitle("Skopobiota")+
  theme(plot.title = element_text(hjust = 0.5))
write.table(p_an_pvalues, file='richness_m_fun_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_m$data, file='richness_m_fun_data.tsv', quote=FALSE, sep='\t')
p_an_m

p_an <-plot_anova_diversity_pval(MBS_No_inoculation, method = c("shannon", "evenness"),grouping_column ="Fungicide",pValueCutoff=0.05)
p_an_n <- p_an[[1]]
p_an_pvalues <- p_an[[2]]
p_an_n <- p_an_n +  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=set_fun)+
  ggtitle("No Inoculation")+
  theme(plot.title = element_text(hjust = 0.5))
write.table(p_an_pvalues, file='richness_n_fun_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_n$data, file='richness_n_fun_data.tsv', quote=FALSE, sep='\t')
p_an_n

#arrange 4 plots in one

p <- ggarrange(p_an_pm,p_an_m,p_an_p,p_an_n,
          widths = 1, heights = 1.1,
          #labels = c("Pch+Skopobiota", "Skopobiota", "No inoculation", "Pch"),
          #label.x = 0.8, label.y = .90, hjust = 0, vjust = 0,
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom")
p
ggsave(p, file="Richness_ino_fun.pdf", width = 30, height = 20, units = "cm")

###--------------------------LCBD---------------------------

# Plot local contribution to beta diversity
physeq <- normalise_data(gio_m_un_con, norm.method = "relative")
#rownames(sample_data(g_sorted)) = sample_data(giovanni)[,13]

#NB: Using standard jaccard does not provide LCBD - use ab.jaccard or ab.sorensen to get abundance data for that.  
p <- plot_taxa(physeq,grouping_column="Inoculation",method="jaccard",number.taxa=20,filename='Relative_abundance_1_20_ino') 
p1 <- p + ggtitle("Relative abundance, 20 most common taxa (group: Inoculation)")
print(p1)
ggsave(p1, file="Relative_abundance_1_20_ino.pdf", width = 30, height = 20, units = "cm")

p <- plot_taxa(physeq,grouping_column="Fungicide",method="jaccard",number.taxa=20,filename='Relative_abundance_1_20_fun')
p1 <- p + ggtitle("Relative abundance, 20 most common taxa (group: Fungicide)") 
print(p1)
ggsave(p1, file="Relative_abundance_1_20_fun.pdf", width = 30, height = 20, units = "cm")

g_sorted_21 <- prune_taxa(names(sort(taxa_sums(physeq),TRUE)[21:50]), physeq)

p <- plot_taxa(g_sorted_21,grouping_column="Inoculation",method="jaccard",number.taxa=20,filename='Relative_abundance_21_40_ino')
p1 <- p + ggtitle("Relative abundance, 21-40 most common taxa (group: Inoculation)")
print(p1)
ggsave(p1, file="Relative_abundance_21_40_ino.pdf", width = 30, height = 20, units = "cm")

p <- plot_taxa(g_sorted_21,grouping_column="Fungicide",method="jaccard",number.taxa=20,filename='Relative_abundance_21_40_fun')
p1 <- p + ggtitle("Relative abundance, 21-40 most common taxa (group: Fungicide)") 
print(p1)
ggsave(p1, file="Relative_abundance_21_40_fun.pdf", width = 30, height = 20, units = "cm")

#---------------Ordination--------------

ord.res <- ordinate(giovanni, which_distance = "jaccard", method="NMDS", grouping_column="Treatment",pvalue.cutoff = 0.05)
p <- plot_ordination(giovanni, ord.res, type="samples", color="Inoculation", shape="Fungicide", title="Jaccard+NMDS ordination") +
  theme_bw()+
  scale_color_manual(values=c("cyan", "green","magenta", "gold", "red", "yellow", "black", "pink", "forestgreen", "purple")) +
  stat_ellipse(aes(group=Treatment),type = "norm")
print(p)
ggsave(p, file="Ordination_jaccard_nmds_treatment.pdf")
p1 <- p + facet_wrap(~Fungicide, 3)
print(p1)
ggsave(p, file="Ordination_jaccard_nmds_treatment_fw.pdf")

ord.res <- ordinate(giovanni, which_distance = "jaccard", method="NMDS", grouping_column="Treatment",pvalue.cutoff = 0.05)
p <- plot_ordination(giovanni, ord.res, type="samples", color="Inoculation", shape="Fungicide", title="Jaccard+NMDS ordination, ellipses = Inoculation") +
  theme_bw()+
  scale_color_manual(values=c("cyan", "green","magenta", "gold", "red", "yellow", "black", "pink", "forestgreen", "purple")) +
  stat_ellipse(aes(group=Inoculation),type = "norm")
print(p)
ggsave(p, file="Ordination_jaccard_nmds_inoculation.pdf")

ord.res <- ordinate(giovanni, which_distance = "jaccard", method="NMDS", grouping_column="Treatment",pvalue.cutoff = 0.05)
p <- plot_ordination(giovanni, ord.res, type="samples", color="Fungicide", shape="Inoculation", title="Jaccard+NMDS ordination, ellipses = Fungicide") +
  theme_bw()+
  scale_color_manual(values=c( "red", "violet", "brown","darkblue", "gold", "yellow", "black", "pink", "forestgreen", "purple")) +
  scale_shape_manual(values=c(2,3,15,10)) +
  stat_ellipse(aes(group=Fungicide),type = "norm")
print(p)
ggsave(p, file="Ordination_jaccard_nmds_fungicide.pdf")

dist1 <- as.matrix(phyloseq::distance(giovanni, "jaccard", type="sample"))
write.table(dist1, file='jaccard_dist_matrix.tsv', quote=FALSE, sep='\t')
dist1 <- as.matrix(phyloseq::distance(giovanni, "bray", type="sample"))
write.table(dist1, file='bray_dist_matrix.tsv', quote=FALSE, sep='\t')


#Mycobioma, Water, Pch, Pch_Mycobioma
#set_ino <- c("cyan", "green","magenta", "gold")
#Blad, Control, CuS, Systemics
#set_fun <- c( "red", "violet", "brown","darkblue")

#------------DeSeq2-------------

#use transposed phyloseq object

##conglomerate tips based on taxonomic level. Can ease computation in larger datasets. 
#temp <- taxa_level(MBS_physeq, "Family")

# deseq_sig_ino <- differential_abundance(MBS_physeq, grouping_column = "Treatment", 
#                                     pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "test_ino.csv")
# 
# 
# deseq_sig <- differential_abundance(MBS_No_inoculation, grouping_column = "Fungicide", 
#                                     pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "deseq_noino")
# deseq_sig <- differential_abundance(MBS_Pch, grouping_column = "Fungicide", 
#                                     pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "deseq_pch")
# deseq_sig <- differential_abundance(MBS_Pch_Skopobiota, grouping_column = "Fungicide", 
#                                     pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "deseq_pch_skopo")
# deseq_sig <- differential_abundance(MBS_Skopobiota, grouping_column = "Fungicide", 
#                                     pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "deseq_skopo")


OTU2 <- otu_table(gio_filter1, taxa_are_rows = TRUE)
OTU2 <- t(OTU2)
#Merge new phyloseq object
TAX = as.matrix(tax_table(gio2))
met <- as.data.frame(sample_data(gio2))
pseq1 <- phyloseq(OTU2, TAX)
MBS_physeq2<-merge_phyloseq(pseq1, met)
rm(OTU2, TAX, met, pseq1)

minTotRelAbun = 0.1 
x = taxa_sums(MBS_physeq2)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
gio3 = prune_taxa(keepTaxa, MBS_physeq2)
gio3 = MBS_physeq2

ID<- as.vector(unique(sample_data(gio3)$Inoculation))
for (i in 1:length(ID)){ 
  temp <- subset_samples(gio3, Inoculation==ID[i])
  temp = prune_taxa(taxa_sums(temp) > 0, temp)
  name <- paste("MBS_", ID[i], sep = "")
  cont <- subset_samples(temp, Fungicide=="Control")
  blad <- subset_samples(temp, Fungicide=="Blad")
  cus <- subset_samples(temp, Fungicide=="CuS")
  sys <- subset_samples(temp, Fungicide=="Systemics")
  tmp <- merge_phyloseq(cont,blad)
  name2 <- paste(name,"c-b")
  deseq_sig <- differential_abundance(tmp, grouping_column = "Fungicide", 
                                      pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = name2)
  tmp <- merge_phyloseq(cont,cus)
  name2 <- paste(name,"c-cu")
  deseq_sig <- differential_abundance(tmp, grouping_column = "Fungicide", 
                                      pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = name2)
  tmp <- merge_phyloseq(cont,sys)
  name2 <- paste(name,"c-s")
  deseq_sig <- differential_abundance(tmp, grouping_column = "Fungicide", 
                                      pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = name2)
  tmp <- merge_phyloseq(cus,blad)
  name2 <- paste(name,"cu_b")
  deseq_sig <- differential_abundance(tmp, grouping_column = "Fungicide", 
                                      pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = name2)
  tmp <- merge_phyloseq(cus,sys)
  name2 <- paste(name,"cu-s")
  deseq_sig <- differential_abundance(tmp, grouping_column = "Fungicide", 
                                      pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = name2)
  tmp <- merge_phyloseq(sys,blad)
  name2 <- paste(name,"s-b")
  deseq_sig <- differential_abundance(tmp, grouping_column = "Fungicide", 
                                      pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = name2)
}

# p <- plot_signif(deseq_sig_ino$plotdata, top.taxa = 100)
# print(p)
# 
# deseq_sig <- differential_abundance(giovanni, grouping_column = "Inoculation", output_norm = "log-relative", 
#                                     pvalue.threshold = 0.05, lfc.threshold = 0, filename = F)
# p <- plot_signif(deseq_sig$plotdata, top.taxa = 10)
# print(p)

###---------------Co-occurence networks--------------------

#co-occurence network are a nice way to illustrate which OTUs are likely to occur in the same samples, 
#and to illustrate patterns. 

tmp <- taxa_level(MBS_physeq, "Genus")

co_occrBlad<- co_occurence_network(MBS_physeq, grouping_column = "Fungicide", rhos = 0.35, method="cor", 
                                qval_threshold=0.05, select.condition = "Blad", scale.vertex.size=4, scale.edge.width=10, 
                                plotNetwork=T, plotBetweennessEeigenvalue=T)
co_occrCtrl <- co_occurence_network(tmp, grouping_column = "Fungicide", rhos = 0.35, method="cor", 
                                qval_threshold=0.05, select.condition = "Control", scale.vertex.size=4, scale.edge.width=10, 
                                plotNetwork=T, plotBetweennessEeigenvalue=T)
co_occrCuS <- co_occurence_network(tmp, grouping_column = "Fungicide", rhos = 0.35, method="cor", 
                                qval_threshold=0.05, select.condition = "CuS", scale.vertex.size=4, scale.edge.width=10, 
                                plotNetwork=T, plotBetweennessEeigenvalue=T)
co_occrSys <- co_occurence_network(tmp, grouping_column = "Fungicide", rhos = 0.35, method="cor", 
                                qval_threshold=0.05, select.condition = "Systemics", scale.vertex.size=4, scale.edge.width=10, 
                                plotNetwork=T, plotBetweennessEeigenvalue=T)



require(visNetwork)
g <- co_occrBlad$net$graph
h <- co_occrCtrl$net$graph
i <- co_occrCuS$net$graph
j <- co_occrSys$net$graph

data <- toVisNetworkData(g)
visNetwork(nodes = data$nodes, edges = data$edges, main="Giovanni dataset, Blad. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
data <- toVisNetworkData(h)
visNetwork(nodes = data$nodes, edges = data$edges, main="Giovanni dataset, Ctrl. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
data <- toVisNetworkData(i)
visNetwork(nodes = data$nodes, edges = data$edges, main="Giovanni dataset, CuS. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
data <- toVisNetworkData(j)
visNetwork(nodes = data$nodes, edges = data$edges, main="Giovanni dataset, Sys. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)


###----------DOESN'T WORK! - Phyloseq to Biom format---------

temp <- gio_filter1

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp) 

data = otu_table(temp) #OTU table
data
omd = tax_table(temp) #Taxonomy table
omd
smd = sample_data(gio_filter1)
colnames(smd)
smd <- smd[,1:13]
smdsmd <- smd[,-3]
smd$SampleID <- gsub("-", "_", smd$SampleID)
# Make a new biom object from component data
newbiom = make_biom(data, id = "fungicide_point1_rel_abun")

write_biom(newbiom, biom_file=file.path(uzdir, "fun_point1.biom"))
#write new metadata table
write.table(smd, file='metadata_lefse.tsv', quote=FALSE, sep='\t')
write.table(omd, file='taxonomy_lefse.tsv', quote=FALSE, sep='\t')

###-----------LEfSe---------------

library("yingtools2")
test <- lefse(gio_filter1, class="Fungicide", subclass = NA, subject = NA, anova.alpha = 0.05,
      wilcoxon.alpha = 0.05, lda.cutoff = 2,
      wilcoxon.within.subclass = FALSE, one.against.one = FALSE,
      mult.test.correction = 0, make.lefse.plots = FALSE,
      by_otus = FALSE, levels = rank_names(gio_filter1))

###-----------Test - create control column for heat tree------------

temp <- gio_filter1
sample_data(temp)$Control <- sample_data(temp)$No_inoculation
sample_data(temp)$Control <- gsub("No_inoculationB", "No", sample_data(temp)$Control); sample_data(temp)$Control 
sample_data(temp)$Control <- gsub("No_inoculationCS", "No", sample_data(temp)$Control); sample_data(temp)$Control 
sample_data(temp)$Control <- gsub("No_inoculationS", "No", sample_data(temp)$Control); sample_data(temp)$Control 
head(sample_data(temp))
head(tax_table(temp))

###-----------Metacoder----------

#Metacoder provides heat trees, which colour the branches of the tree that 
# has any significant change in abudance between two treatments. 

temp <- normalise_data(gio_filter1, norm.method = "relative")
temp <- tax_glom(temp, taxrank = "Genus")
#temp <- giovanni

#parse phyloseq object giovanni
obj <- parse_phyloseq(temp)


tissuegroup <- obj$data$sample_data$Treatment

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_treatment.pdf')


test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_treatment.tsv', quote=FALSE, sep='\t')

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_treatment_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_treatment_pval.tsv', quote=FALSE, sep='\t')

### Control

tissuegroup <- obj$data$sample_data$Control

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

# heat_tree_matrix(obj,
#                  data = "diff_table",
#                  node_size = n_obs,
#                  node_label = taxon_names,
#                  node_color = log2_median_ratio,
#                  node_color_trans = "linear",
#                  node_color_interval = c(-3, 3),
#                  edge_color_interval = c(-3, 3),
#                  layout = "davidson-harel",
#                  initial_layout = "reingold-tilford",
#                  node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
#                  #title = "Heat Tree Matrix",
#                  title_size = 0.05,
#                  node_size_axis_label = "Number of species",
#                  node_color_axis_label = "Log2 ratio median proportions",
#                  output_file = 'heattree_treatment.pdf')

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
          node_size_axis_label = "species count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
          title = "Red Nodes = Higher Abundance in Control, No Inoculation \n Blue Nodes = Higher Abundance in Treatment",
          title_size = 0.035,
          edge_label_size_range = c(1,5),
          output_file = 'heattree_control.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_treatment.tsv', quote=FALSE, sep='\t')

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
          node_size_axis_label = "species count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
          title = "Red Nodes = Significantly Higher Abundance in Control, No Inoculation \n Blue Nodes = Significantly Higher Abundance in Treatment",
          title_size = 0.035,
          edge_label_size_range = c(1,5),
          output_file = 'heattree_control_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_treatment_pval.tsv', quote=FALSE, sep='\t')


### Inoculaton

tissuegroup <- obj$data$sample_data$Inoculation

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_ino.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_ino.tsv', quote=FALSE, sep='\t')

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_ino_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_ino_pval.tsv', quote=FALSE, sep='\t')

###---------------HEAT TREES - FUNGICIDE SUBSETS---------------------

test <- gio_filter1

tmp <- normalise_data(test, norm.method = "relative")
tmp
#temp <- tax_glom(temp, taxrank = "Genus")
#temp <- giovanni

#Create subsets - if using for combined heat tree, DO NOT prune >0! 

ID<- as.vector(unique(sample_data(tmp)$Inoculation))
for (i in 1:length(ID)){ 
  temp <- subset_samples(tmp, Inoculation==ID[i])
  #temp = prune_taxa(taxa_sums(temp) > 0, temp)
  name <- paste("Filter1_relab_", ID[i], sep = "")
  assign(name, temp)
}

### No inoculation

obj <- parse_phyloseq(Filter1_relab_No_inoculation)
tissuegroup <- obj$data$sample_data$Fungicide

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_noino.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_noino.tsv', quote=FALSE, sep='\t')


#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_noino_pval.pdf')


test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_noino_pval.tsv', quote=FALSE, sep='\t')

### Pch

obj <- parse_phyloseq(Filter1_relab_Pch)
tissuegroup <- obj$data$sample_data$Fungicide

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_pch.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_pch.tsv', quote=FALSE, sep='\t')

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_pch_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_pch_pval.tsv', quote=FALSE, sep='\t')

### Pch + skopobiota

obj <- parse_phyloseq(Filter1_relab_Pch_Skopobiota)
tissuegroup <- obj$data$sample_data$Fungicide

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_pch_skopo.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_pch_skopo.tsv', quote=FALSE, sep='\t')

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_pch_skopo_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_pch_skopo_pval.tsv', quote=FALSE, sep='\t')

### Skopobiota

obj <- parse_phyloseq(Filter1_relab_Skopobiota)
tissuegroup <- obj$data$sample_data$Fungicide

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_skopo.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_skopo.tsv', quote=FALSE, sep='\t')


#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
                 #title = "Heat Tree Matrix",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_skopo_pval.pdf')


test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_skopo_pval.tsv', quote=FALSE, sep='\t')


# #Test p-values
# obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
# test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
# obj$data$diff_table <- test
# 
# heat_tree_matrix(obj,
#                  data = "diff_table",
#                  node_size = n_obs,
#                  node_label = taxon_names,
#                  node_color = log2_median_ratio,
#                  node_color_trans = "linear",
#                  node_color_interval = c(-3, 3),
#                  edge_color_interval = c(-3, 3),
#                  layout = "davidson-harel",
#                  initial_layout = "reingold-tilford",
#                  node_color_range = c("darkblue", "royalblue", "gray95", "red", "red4"),
#                  #title = "Heat Tree Matrix",
#                  title_size = 0.05,
#                  node_size_axis_label = "Number of species",
#                  node_color_axis_label = "Log2 ratio median proportions",
#                  output_file = 'heattree_fun_pval.pdf')
# 
# test <- obj$data$diff_table %>%
#   mutate(taxon_names = taxon_names(obj)[taxon_id], 
#          taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
#   arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
# write.table(test, file='heattree_fun_pval.tsv', quote=FALSE, sep='\t')

###



###------------Colour sign OTUs from LEfSe------

#Extract sign OTUs from LEfSe, set all other OTUs to log2 median ratio = 0 

#Import list of OTUs
#Get taxonomy data
#Combine tax data into string


#Tell metacoder to set all other otus to zero.
#Make list of OTUs that are not significant, set those to zero. 

#_______________object is "sign_otus"______________

###------------LEfSe------------

#This will make a subset and .tsv files of any OTUs found to be significantly more 
#abundant in either control or exclosure plots. 
#This is to create more precise graphics (the significant ones are flooded out by insignificant OTUs)

#Run a lefse in galaxy on the full dataset. 
#Get .tsv from lefse with all significant OTUs - format as needed in excel. 
#This bit of code gets the lefse .tsv, gets rownames as a list we can use to make a subset of the significant OTUs,
# and then writes .tsv files for further use - lefse and metacoder for graphics, for instance. 

sign <- read.table(file = 'Sign_lefse.tsv', sep = '\t', header = TRUE, row.names = 1)

sign_test <- rownames(sign)

sign_otus = prune_taxa(sign_test, giovanni)

#This is now a phyloseq object of only OTUs found significant in LEfSe
sign_otus

write.table(otu_table(sign_otus), file='LEfSe_sign_otus52.tsv', quote=FALSE, sep='\t')
write.table(tax_table(sign_otus), file='LEfSe_sign_tax52.tsv', quote=FALSE, sep='\t')


###-----------LEfSe-Metacoder--------------

#Uses the dataset of significant otus created from LEfSe. 

obj_bac_sign <- subset_taxa(sign_otus, Kingdom=="k:Bacteria")
obj_arc_sign <- subset_taxa(sign_otus, Kingdom=="k:Archaea")

##Archaea

obj <- parse_phyloseq(obj_arc_sign)

obj$data$otu_table <- calc_obs_props(obj,
                                     dataset = "otu_table",
                                     cols = obj$data$sam_data$sample_ids)

obj$data$otu_table[is.na(obj$data$otu_table)] <- 0

#test <- within(obj$data$otu_table, log2_median_ratio[log2_median_ratio == NaN] <- 0 )

obj$data$tax_abund <- calc_taxon_abund(obj, 
                                       dataset = "otu_table", 
                                       cols = obj$data$sam_data$sample_ids)

obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = obj$data$sam_data$Type)

obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund",
                                      cols = obj$data$tax_table$taxon_id,
                                      groups = obj$data$sam_data$Type)


obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value,
                                               method = "fdr")

hist(obj$data$diff_table$wilcox_p_value)

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          node_color_range = c("cyan", "gray", "tan"),
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          title = "Tan = Significantly higher abundance in control plots.", 
          node_label_size_range = c(0.04, 0.05), 
          title_size = 0.03)

##Bacteria

obj <- parse_phyloseq(obj_bac_sign)

obj$data$otu_table <- calc_obs_props(obj,
                                     dataset = "otu_table",
                                     cols = obj$data$sam_data$sample_ids)

obj$data$tax_abund <- calc_taxon_abund(obj, 
                                       dataset = "otu_table", 
                                       cols = obj$data$sam_data$sample_ids)

obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = obj$data$sam_data$Type)

obj$data$diff_table <- compare_groups(obj, dataset = "tax_abund",
                                      cols = obj$data$tax_table$taxon_id,
                                      groups = obj$data$sam_data$Type)

obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value,
                                               method = "fdr")

hist(obj$data$diff_table$wilcox_p_value)

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          node_color_range = c("cyan", "gray", "tan"),
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          title = "Tan = Significantly higher abundance in control plots.", 
          title_size = 0.03,
          node_label_size_range = c(0.01, 0.015),
          output_file = NULL)

###-----------Metacoder-Microbiomeseq------


# Extract abundance matrix from the phyloseq object
OTU2 <- otu_table(obj_bac_all, taxa_are_rows = TRUE)
OTU2 <- t(OTU2)
#Merge new phyloseq object
TAX = as.matrix(tax_table(giovanni))
met <- as.data.frame(metadata)
pseq1 <- phyloseq(OTU2, TAX)

bac_MSeq <-merge_phyloseq(pseq1, met)

#Runs an anova on the richness indices, and plots them with p-values and boxplots. 

p_an<-plot_anova_diversity(bac_MSeq, method = c("richness","simpson", "shannon", "evenness"),grouping_column ="Type",pValueCutoff=0.05)
print(p_an + labs(title = "Bacteria", subtitle = "* = p < 0.05, ** = p < 0.01") 
)

###----------------End of code-------------------
