#Code written by M. R. Aggerbeck. 
#mrag@envs.au.dk
#Creative commons, 2018. 

###----------------PACKAGES--------------
#Packages used in the following code. 
library("phyloseq"); packageVersion("phyloseq")
library("devtools"); packageVersion("devtools")
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
setwd("D:/Marie/Dropbox/Muskox_soil_microbial/WD/")
#Set user directory
uzdir <- "D:/Marie/Dropbox/Muskox_soil_microbial/WD/"

###---------------------LOAD DATA-------------------------------------

#Biom file to import - preferably JSON
#Dataset 
data_biom <- "otu_table_sintax_rar.biom"

#Import data, removing any blank spaces in file path
biom_file <- paste(uzdir, data_biom, sep="")
#the actual act of importing the file into R
biom_otu_tax <- import_biom(biom_file)
colnames(tax_table(biom_otu_tax)) = c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")
#inspect tax_table
head(tax_table(biom_otu_tax))
tax_table(biom_otu_tax) <- gsub(":", "_", tax_table(biom_otu_tax))
tax_table(biom_otu_tax) <- gsub("_sp._", "_sp_", tax_table(biom_otu_tax))
head(tax_table(biom_otu_tax))
head(otu_table(biom_otu_tax))


#Import metadata file, removing blank spaces in path 
md_file <- paste(uzdir, "mapping_file2.txt", sep="")   

#Has to be tsv. If csv, use code below. 
metadata <- import_qiime_sample_data(md_file)
metadata

###--------------Data preparation---------------

# #If the above doesn't produce a table, use this: 
# metadata <- read.csv(md_file, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")

#change names in metadata to fit proper nomenclature

temp <- as.matrix(metadata)
temp <- gsub("-", "_", temp)
temp
metadata <- temp
metadata <- as.data.frame(metadata)
metadata <- sample_data(metadata)

head(metadata)


#Merge OTU data and metadata to single object - The Phyloseq Object(tm)
physeq <- merge_phyloseq(biom_otu_tax,metadata)
sample_names(physeq) <- gsub("-", "_", sample_names(physeq))
sample_names(physeq)

#physeq is now the s4 object - the phyloseq object - we will use throughout the rest of the code. 
#In theory, all other objects in the global environment could be removed now, but I tend to leave them, just in case I need them.

head(sample_data(physeq))
head(tax_table(physeq))
head(otu_table(physeq))

#C4-5, 16-Tue_S16 is contaminated, proceed to remove from dataset. 
physeq = subset_samples(physeq, Sample_Name != "C4_5")
metadata = subset_samples(metadata, Sample_Name != "C4_5")
#/data loaded

raw_data <- physeq
#physeq <- physeqa
#physeqa <- physeq
##Reset to raw data in case of error
#physeq <- raw_data


rank_names(physeq)
sample_names(physeq) = sample_data(physeq)$Sample_Name

tax_table(physeq)[is.na(tax_table(physeq))] <- "unidentified"
tax_table(physeq)[1:10]

tax_table(physeq)[,7] <-  gsub("unidentified", "s_unidentified", tax_table(physeq)[,7]); 
tax_table(physeq)[,6] <-  gsub("unidentified", "g_unidentified", tax_table(physeq)[,6]); 
tax_table(physeq)[,5] <-  gsub("unidentified", "f_unidentified", tax_table(physeq)[,5]);
tax_table(physeq)[,4] <-  gsub("unidentified", "o_unidentified", tax_table(physeq)[,4]);
tax_table(physeq)[,3] <-  gsub("unidentified", "c_unidentified", tax_table(physeq)[,3]);
tax_table(physeq)[,2] <-  gsub("unidentified", "p_unidentified", tax_table(physeq)[,2]);

physeq <- tax_glom(physeq, taxrank="Species"); physeq
length(unique(taxa_names(physeq)))
length(unique(tax_table(physeq)[,7]))

tax_table(physeq) <-  gsub("s_unidentified", "", tax_table(physeq)); 
tax_table(physeq) <-  gsub("g_unidentified", "", tax_table(physeq)); 
tax_table(physeq) <-  gsub("f_unidentified", "", tax_table(physeq));
tax_table(physeq) <-  gsub("o_unidentified", "", tax_table(physeq));
tax_table(physeq) <-  gsub("c_unidentified", "", tax_table(physeq));
tax_table(physeq) <-  gsub("p_unidentified", "", tax_table(physeq));

met <- as.matrix(sample_data(physeq))

renam <- as.matrix(sample_data(physeq2)[,6])

renam <-  gsub("Control", "Grazed", renam)
renam <-  gsub("Exclosure", "Ungrazed", renam)

met[,6] <- renam
met <- as.data.frame(met)
met <- sample_data(met)

sample_data(MBS_physeq) <- met

###-----------------Make unique names for taxa----------------------

mynames = NULL

for (i in 1:length(taxa_names(physeq))){
  name <- makeTaxLabel(taxa_names(physeq)[i],physeq)
  mynames <- rbind(mynames, c(name))
  
}

#Find duplicates
n_occur <- data.frame(table(mynames))
n_occur[n_occur$Freq > 1,]


taxa_names(physeq) = make.unique(mynames, sep="_")
taxa_names(physeq)[1:10]

physeq
sample_data(physeq)[1:5]
otu_table(physeq) [1:5]
tax_table(physeq)[1:15]

###------------------Filter out low-abundance reads-----------------

physeq10 <- physeq
otu_table(physeq10)[otu_table(physeq10)<10 ] <- 0
physeq10 = prune_taxa(taxa_sums(physeq10) > 0, physeq10)
physeq = prune_taxa(taxa_sums(physeq) > 0, physeq)
physeq10
physeq

###------------------SUBSETTING-------------------------------------

# #This for loop creates four subset from the fungicide treatments (other dataset). 
# #We create a list of the names of the fungicides, and then use that to specify the subset command. 
# #Then we paste the name of the fungicide onto the subset and assign it as a separate phyloseq object.
# ID<- as.vector(unique(sample_data(physeq)$Fungicide))
# for (i in 1:length(ID)){ 
#   temp <- subset_samples(physeq, Fungicide==ID[i])
#   temp = prune_taxa(taxa_sums(temp) > 0, temp)
#   name <- paste(ID[i])
#   assign(name, temp)
# }
# 
# #Creates another set of subset, based on inoculation
# ID<- as.vector(unique(sample_data(physeq)$Inoculation))
# for (i in 1:length(ID)){ 
#   temp <- subset_samples(physeq, Inoculation==ID[i])
#   temp = prune_taxa(taxa_sums(temp) > 0, temp)
#   name <- paste(ID[i])
#   assign(name, temp)
# }

## Or you can make subsets of treatments. Specifies which metadata column, and which value in said column should be criteria. 
# No_ino <- subset_samples(physeq, No_inoculation!="No")
# Control <- subset_samples(physeq, Fungicide=="Control")
#Cont_no_ino <- subset_samples(physeq, No_inoculation=="No_inoculationC")

###----------------------Merge samples----------------


temp = prune_taxa(taxa_sums(physeq) > 0, physeq)
temp2 <- as.data.frame(unique(sample_data(temp)[,5]))
groups_phy <- as.vector(as.character(temp2$Sample))
sample_data(temp)$group_phy <- get_variable(temp, "Sample") %in% groups_phy

merged_phy = merge_samples(temp, "Sample")

merged_phy
sample_names(temp)
sample_names(merged_phy)

sample_data(merged_phy)
sample_data(merged_phy)$Sample_Name = sample_names(merged_phy)
merged_phy
physeq

tmp <- sample_data(merged_phy)

class(tmp)
write.table(tmp, file='merged_metadata.tsv', quote=FALSE, sep='\t')

#It's a hassle to change every column in R - much faster to just copy-paste block metadata in excel. 
#Do this, then return and import new file. 

md_file <- paste(uzdir, "merged_metadata_alt.tsv", sep="")   
meta <- import_qiime_sample_data(md_file)

otu <- otu_table(merged_phy)
tax <- tax_table(merged_phy)

temp <- merge_phyloseq(otu,tax,meta)
merged_phy <- temp
sample_data(merged_phy)

rm(otu,tax,temp,temp2,tmp)

###Merge unassigned taxa

sort(taxa_names(physeq))

#IMPORTANT! Input manually any unwanted higher classes from the list produced by the sort() above
badtaxa = c("k__Fungi", "o__Agaricales", "o__Coniochaetales", "o__Pleosporales", "p__Ascomycota", "p__Basidiomycota",
            "c__Agaricomycetes", "c__Dothideomycetes", "c__Leotiomycetes", "c__Microbotryomycetes", "c__Sordariomycetes", 
            "c__Tremellomycetes")

phy_un_con <- merge_taxa(physeq, badtaxa)
taxa_names(phy_un_con) <-  gsub("o__Agaricales", "Unassigned", taxa_names(phy_un_con)) 
head(tax_table(phy_un_con))
head(sample_data(phy_un_con))

sort(taxa_names(merged_phy))
badtaxa = c("k__Fungi", "o__Agaricales", "o__Coniochaetales", "o__Pleosporales", "p__Ascomycota", "p__Basidiomycota",
            "c__Agaricomycetes","c__Microbotryomycetes", "c__Sordariomycetes", "c__Tremellomycetes")

phy_m_un_con <- merge_taxa(merged_phy, badtaxa)
taxa_names(phy_m_un_con) <-  gsub("o__Agaricales", "Unassigned", taxa_names(phy_m_un_con)) 
head(tax_table(phy_m_un_con))
head(sample_data(phy_m_un_con))

sort(taxa_names(phy_m_un_con))

# ###----------Create unique names (?)---------
# 
# temp = prune_taxa(taxa_sums(physeq) > 0, physeq)
# temp2 <- as.data.frame(unique(sample_data(temp)[,11]))
# groups_phy <- as.vector(as.character(temp2$Sample))
# sample_data(temp)$group_phy <- get_variable(temp, "Sample") %in% groups_phy
# 
# merged_phy_2 = merge_samples(temp, "Treatment")
# 
# merged_phy_2
# sample_names(temp)
# sample_names(merged_phy_2)
# 
# sample_data(merged_phy_2)
# sample_data(merged_phy_2)$Sample_Name = sample_names(merged_phy_2)
# merged_phy_2
# physeq


###----------Filter by relative abundance-----------

physeq

#Filters out anything below 0.01% of total abundance. 
minTotRelAbun = 0.01 
x = taxa_sums(physeq)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
phy_filter1 = prune_taxa(keepTaxa, physeq)

y = as.matrix(sort((x / sum(x))*100, decreasing = TRUE))

x = taxa_sums(merged_phy)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
phy_m_filter1 = prune_taxa(keepTaxa, merged_phy)


###-----------------Colour palettes-----------------------

# I will from now on use a colour scheme that

set_2 <- c("royalblue","tomato")
set_3 <- c("royalblue", "grey95","tomato")
set_block <- c("cyan", "green","magenta", "gold")
set_fun <- c( "red", "violet", "brown","darkblue")


###----------------------Core biome----------------

#write.table(otu_table(physeq)[2], file='presence_pch.tsv', quote=FALSE, sep='\t')
write.table(y, file='taxalist_rel_ab.tsv', quote=FALSE, sep='\t')

core.taxa.standard <- core_members(physeq, detection = 0, prevalence = 75/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_90.tsv', quote=FALSE, sep='\t')

core = prune_taxa(core.taxa.standard, physeq)

excl_phy = subset_samples(physeq, Treatment=="Exclosure")
excl_phy = prune_taxa(taxa_sums(excl_phy) > 0, excl_phy)
excl_phy

ctrl_phy = subset_samples(physeq, Treatment=="Control")
ctrl_phy = prune_taxa(taxa_sums(ctrl_phy) > 0, ctrl_phy)
ctrl_phy

allTaxa = taxa_names(physeq)
#myTaxa <- allTaxa[!(allTaxa %in% core.taxa.standard)]
ctrlTaxa <- taxa_names(ctrl_phy)
exclTaxa <- taxa_names(excl_phy)
myTaxa <- allTaxa[!(allTaxa %in% ctrlTaxa)]
unique_excl = prune_taxa(myTaxa, physeq)
unique_excl = prune_taxa(taxa_sums(unique_excl) > 0, unique_excl)
unique_excl
ns_unique_excl = prune_taxa

myTaxa <- allTaxa[!(allTaxa %in% exclTaxa)]
unique_ctrl = prune_taxa(myTaxa, physeq)
unique_ctrl = prune_taxa(taxa_sums(unique_ctrl) > 0, unique_ctrl)
unique_ctrl

###----------------------Make file for LEfSe---------------------

test <- physeq
tax_table(test)[is.na(tax_table(test))] <- "unidentified"
test <- tax_glom(test, taxrank="Species"); test
otu_lefse <- as.data.frame(otu_table(test)) 
tax_lefse <-tax_table(test)
met_lefse <-sample_data(test)
met_lefse <- as.matrix(met_lefse)

names_lefse <- paste(tax_lefse[,1],tax_lefse[,2],tax_lefse[,3],tax_lefse[4],tax_lefse[,5],tax_lefse[,6],tax_lefse[,7],"x", sep="å")
names_lefse <-  gsub("s_", "", names_lefse); 
names_lefse <-  gsub("f_", "", names_lefse); 
names_lefse <-  gsub("g_", "", names_lefse); 
names_lefse <-  gsub("o_", "", names_lefse); 
names_lefse <-  gsub("c_", "", names_lefse); 
names_lefse <-  gsub("p_", "", names_lefse);
names_lefse <-  gsub("k_", "", names_lefse);
names_lefse <-  gsub("åååx", "", names_lefse);
names_lefse <-  gsub("ååx", "", names_lefse); 
names_lefse <-  gsub("åx", "", names_lefse); 
names_lefse <-  gsub("å", "|", names_lefse); 
head(names_lefse)
head(rownames(otu_lefse))
head(colnames(otu_lefse))
colnames_otu_lefse <- met_lefse[,8]
colnames(otu_lefse) = colnames_otu_lefse
length(names_lefse)
length(unique(names_lefse))

names_lefse_u <- make.unique(names_lefse)
head(names_lefse_u)
samdata_lefse <- rbind(met_lefse[,8],met_lefse[,6])
rownames(samdata_lefse) = cbind("sampleID", "Treatment")
otu_lefse <- as.matrix(otu_lefse)
rownames(otu_lefse) = names_lefse_u
temp <- rbind(samdata_lefse,otu_lefse)
head(otu_lefse)
head(temp)

write.table(temp, file='lefse_sintax_rar_taxglom_raw.tsv', quote=FALSE, sep='\t', col.names=FALSE)

###----------------------Make subset of significant lefse otus---------------------------------------

#load file from lefse
lefse_sign_raw  <- as.character(as.matrix(LEfSe_sintax_rarefied_glom_raw_SIGN[,1]))

test

otu <- otu_table(test)
tax <- tax_table(test)
met <- sample_data(test)

temp <- rownames(tax)
tmp <- cbind(temp,tax)
head(tmp)
rownames(tmp) = names_lefse_u
head(tmp)[1]

tmp <- tax_table(tmp)

temp <- merge_phyloseq(otu,tmp,met)
temp

taxa_names(temp) = names_lefse_u

head(sample_data(temp))
head(tax_table(temp))
head(otu_table(temp))

subset_lefse <- prune_taxa(lefse_sign_raw, test)
subset_lefse



###----------------------PERMANOVA---------------------------

#Transform data to relative abundance ("compositional")
pseq.rel <- microbiome::transform(core, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

p <- microbiome::plot_landscape(pseq.rel, method = "NMDS", distance = "jaccard", col = "Treatment", size = 3)
print(p)

permanova <- vegan::adonis(t(otu)~Treatment,
                           data = meta, permutations=999, method = "jaccard")
write.table(as.data.frame(permanova$aov.tab), file='permanova_jaccard_all.tsv', quote=FALSE, sep='\t')

print(as.data.frame(permanova$aov.tab))

dist <- vegan::vegdist(t(otu))
temp <- anova(vegan::betadisper(dist, meta$Treatment))
temp
write.table(temp, file='anova_treatment.tsv', quote=FALSE, sep='\t')

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

#physeqrplot <- plot_richness(physeq, x="Inoculation", color="Fungicide", measures=c("Shannon", "Simpson", "InvSimpson"), title = "Richness, physeq dataset")
#physeqrplot + geom_boxplot(data = physeqrplot$data, aes(color = NULL), alpha = 0.1)

p_an <-plot_anova_diversity_pval(physeq2, method = c("richness","shannon", "evenness"),
                                 grouping_column ="Treatment",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("* = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p <- p_an_o+scale_colour_manual(values = set_2, labels = c("Grazed   ", "Ungrazed   "))
p <- p + theme(legend.title = element_blank())+
  theme(legend.position="bottom")
p
ggsave(p, file="richness_treatment.pdf")
write.table(p_an_pvalues, file='richness_treatment_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_o$data, file='richness_treatment_data.tsv', quote=FALSE, sep='\t')

#Blocks
p_an <-plot_anova_diversity_pval(physeq, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="Sample",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_an_o
ggsave(p, file="richness_inoculation.pdf")
write.table(p_an_pvalues, file='richness_inoculation_pval.tsv', quote=FALSE, sep='\t')
write.table(p$data, file='richness_inoculation_data.tsv', quote=FALSE, sep='\t')

p1 <- plot_richness(physeq, x="Treatment", color="Treatment", measures=c("Chao1", "Shannon", "Simpson", "InvSimpson"))
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = NULL), alpha = 0.1)
p1
#Blad_richness <- p1 + ggtitle("Richness, Blad")
#Blad_richness

###-----------Heatmaps------------

temp <- phy_m_filter1
tmp <- temp
tmp <- prune_taxa(names(sort(taxa_sums(temp),TRUE)[1:100]), temp)
p <- plot_heatmap(unique_excl, method="NMDS", distance="jaccard", sample.label="Treatment", 
                  sample.order = "Treatment", title="Core Taxa") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank())
print(p)
ggsave(p, file="Heatmap_core.pdf", width = 20, height = 25, units = "cm")

###-----------Microbiomeseq------

#this package holds a lot of promise, but is still in beta, so if you use it, contact the author on how to reference it. 
#It implements a ton of statistical tools and returns nice figures, but is a little buggy. 

#Some of the commands in this packae requires a transposed matrix, so first, we'll transpose our phyloseq obejct:
# Extract abundance matrix from the phyloseq object
OTU2 <- otu_table(physeq, taxa_are_rows = TRUE)
OTU2 <- t(OTU2)
#Merge new phyloseq object
TAX = as.matrix(tax_table(physeq))
met <- as.data.frame(sample_data(physeq))
pseq1 <- phyloseq(OTU2, TAX)
MBS_physeq<-merge_phyloseq(pseq1, met)

#OTU <- otu_table(MBS_Pch)
#TAX <- tax_table(MBS_Pch)
#met <- sample_data(MBS_Pch)

#write.table(OTU, file='otu_table.tsv', quote=FALSE, sep='\t')
#write.table(TAX, file='tax_table.tsv', quote=FALSE, sep='\t')
#write.table(met, file="sam_data.tsv", quote=FALSE, sep='\t')

#MBS_physeq is now the transposed version of physeq. 

rm(OTU2, TAX, met, pseq1)

MBS_bac <- subset_taxa(MBS_physeq, Kingdom=="k_Bacteria")
MBS_arc <- subset_taxa(MBS_physeq, Kingdom=="k_Archaea")

#subset transposed phyloseq object into inoculations
ID<- as.vector(unique(sample_data(MBS_physeq)$Block))
for (i in 1:length(ID)){ 
  temp <- subset_samples(MBS_physeq, Inoculation==ID[i])
  temp = prune_taxa(taxa_sums(temp) > 0, temp)
  name <- paste("MBS_", ID[i], sep = "")
  assign(name, temp)
}
MBS_

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

p1 <- plot_richness(core, x="Treatment", color="Block", measures=c("Shannon", "InvSimpson"))+
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
p1 <- p1 + geom_boxplot(data = p1$data, aes(color = Treatment), alpha = 0.1)
Richness <- p1 + ggtitle("Richness, Pch")
Richness
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
temp <- normalise_data(physeq, norm.method = "relative")
temp <- tax_glom(temp, taxrank = "Family")
taxa_names(temp) = make.unique(tax_table(temp)[,5])
#rownames(sample_data(g_sorted)) = sample_data(physeq)[,13]
tmp <- prune_taxa(names(sort(taxa_sums(physeq),TRUE)[1:25]), physeq)
tmp <- normalise_data(physeq, norm.method = "relative")

#NB: Using standard jaccard does not provide LCBD - use ab.jaccard or ab.sorensen to get abundance data for that.  
p <- plot_taxa(temp,grouping_column="Type",method="ab.sorensen",number.taxa=21,filename='Relative_abundance_1_20') 
p1 <- p + ggtitle("Relative abundance, 20 most common taxa")
print(p)
ggsave(p, file="Relative_abundance_1_20_treatment_ab_all.pdf", width = 30, height = 20, units = "cm")
write.table(p$data, file='lcbd_data_all.tsv', quote=FALSE, sep='\t')

g_sorted_21 <- prune_taxa(names(sort(taxa_sums(temp),TRUE)[21:50]), physeq)

p <- plot_taxa(g_sorted_21,grouping_column="Treatment",method="jaccard",number.taxa=20,filename='Relative_abundance_21_40_ino')
p1 <- p + ggtitle("Relative abundance, 21-40 most common taxa")
print(p1)
ggsave(p1, file="Relative_abundance_21_40_treatment.pdf", width = 30, height = 20, units = "cm")

#---------------Ordination--------------

ord.res <- ordinate(physeq, which_distance = "bray", method="NMDS", grouping_column="Treatment",pvalue.cutoff = 0.05)
p <- plot_ordination(physeq, ord.res, type="samples", color="Type", shape="Block") +
  theme_bw()+
  scale_color_manual(values=c(set_2)) +
  stat_ellipse(aes(group=Treatment),type = "norm")
p <- p + theme(legend.title = element_blank())+
  theme(legend.position="bottom")
print(p)
ggsave(p, file="Ordination_bray_nmds_treatment.pdf")

meta_test <- as.data.frame(sample_data(physeq))

ord.res <- ordination(physeq, which_distance = "bray", method = "NMDS", grouping_column = "Treatment", pvalue.cutoff = 0.05)
ord.res.sol <- ord.res$solution
p <- plot_ordisurf(sol = ord.res.sol, meta_table = meta_test, env.variable = 'pH', grouping_column = "Treatment")
print(p)

p1 <- p + facet_wrap(~Block, 3)
print(p1)
ggsave(p, file="Ordination_jaccard_nmds_treatment_fw.pdf")

#------------DeSeq2-------------

#use transposed phyloseq object

##conglomerate tips based on taxonomic level. Can ease computation in larger datasets. 
#temp <- taxa_level(MBS_physeq, "Family")

deseq_sig <- differential_abundance(MBS_physeq2, grouping_column = "Treatment", 
                                   pvalue.threshold = 0.05, lfc.threshold = 0.1, filename = "Deseq2_bac_")
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


OTU2 <- otu_table(phy_filter1, taxa_are_rows = TRUE)
OTU2 <- t(OTU2)
#Merge new phyloseq object
TAX = as.matrix(tax_table(physeq10))
met <- as.data.frame(sample_data(physeq10))
pseq1 <- phyloseq(OTU2, TAX)
MBS_physeq2<-merge_phyloseq(pseq1, met)
rm(OTU2, TAX, met, pseq1)

minTotRelAbun = 0.1 
x = taxa_sums(MBS_physeq2)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
phy3 = prune_taxa(keepTaxa, MBS_physeq2)
phy3 = MBS_physeq2

ID<- as.vector(unique(sample_data(phy3)$Inoculation))
for (i in 1:length(ID)){ 
  temp <- subset_samples(phy3, Inoculation==ID[i])
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

p <- plot_signif(deseq_sig$plotdata, top.taxa = 100)
print(p)
# 
# deseq_sig <- differential_abundance(physeq, grouping_column = "Inoculation", output_norm = "log-relative", 
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
visNetwork(nodes = data$nodes, edges = data$edges, main="physeq dataset, Blad. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
data <- toVisNetworkData(h)
visNetwork(nodes = data$nodes, edges = data$edges, main="physeq dataset, Ctrl. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
data <- toVisNetworkData(i)
visNetwork(nodes = data$nodes, edges = data$edges, main="physeq dataset, CuS. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
data <- toVisNetworkData(j)
visNetwork(nodes = data$nodes, edges = data$edges, main="physeq dataset, Sys. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)

###-----------Metacoder----------

#Metacoder provides heat trees, which colour the branches of the tree that 
# has any significant change in abundance between two treatments. 

temp <- normalise_data(core, norm.method = "relative")
#remember to clear cache (i.e. delete obj before doingn this, otherwise it messes up sometimes. ) 
temp <- tax_glom(subset_lefse, taxrank = "Family")
#temp <- physeq
subset_lefse
temp1 <- subset_taxa(temp, Kingdom=="k_Bacteria")
temp1 = prune_taxa(taxa_sums(temp1) > 0, temp1)
temp2 <- subset_taxa(temp, Kingdom=="k_Archaea")
temp2 = prune_taxa(taxa_sums(temp2) > 0, temp2)

#parse phyloseq object physeq
obj <- parse_phyloseq(temp1)


tissuegroup <- obj$data$sample_data$Treatment

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

heat_tree(obj,
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_range = c("red4", "red", "gray95", "royalblue", "darkblue"),
                 #title = "Blue = Upregulated in Control \n Red = Upregulated in Exclosures",
                 #title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_treatment_signlefse_f.pdf')



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
          node_size = n_obs,
          node_label = taxon_names,
          node_color = log2_median_ratio,
          node_color_trans = "linear",
          node_color_interval = c(-3, 3),
          edge_color_interval = c(-3, 3),
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
          node_color_range = c("red4", "red", "gray95", "royalblue", "darkblue"),
          #title = "Heat Tree Matrix",
          title_size = 0.05,
          node_size_axis_label = "Number of species",
          node_color_axis_label = "Log2 ratio median proportions",
          output_file = 'heattree_treatment.pdf')

###--------------------heat tree matrices----------------------------

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

test <- phy_filter1

tmp <- normalise_data(test, norm.method = "relative")
tmp
#temp <- tax_glom(temp, taxrank = "Genus")
#temp <- physeq

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

sign_otus = prune_taxa(sign_test, physeq)

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
TAX = as.matrix(tax_table(physeq))
met <- as.data.frame(metadata)
pseq1 <- phyloseq(OTU2, TAX)

bac_MSeq <-merge_phyloseq(pseq1, met)

#Runs an anova on the richness indices, and plots them with p-values and boxplots. 

p_an<-plot_anova_diversity(bac_MSeq, method = c("richness","simpson", "shannon", "evenness"),grouping_column ="Type",pValueCutoff=0.05)
print(p_an + labs(title = "Bacteria", subtitle = "* = p < 0.05, ** = p < 0.01") 
)

###----------------End of code-------------------
