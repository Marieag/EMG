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
library("pairwiseAdonis"); packageVersion("pairwiseAdonis")

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
temp <- gsub("Mycobioma", "SK_ACEA1", temp)
temp <- gsub("Myco", "Sk", temp)
temp <- gsub("Patogeno", "Pathogen", temp)
temp
metadata <- temp
temp <- colnames(metadata)
temp <- gsub("Mycobioma", "SK_ACEA1", temp)
temp <- gsub("Myco", "Sk", temp)
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


###------------------Filter out low-abundance reads-----------------

gio2 <- gio_sp
otu_table(gio2)[otu_table(gio2)<25 ] <- 0
gio2 = prune_taxa(taxa_sums(gio2) > 0, gio2)

giovanni <- gio2

###----------------------Merge samples----------------


temp = prune_taxa(taxa_sums(giovanni) > 0, giovanni)
temp2 <- as.data.frame(unique(sample_data(temp)[,4]))
groups_gio <- as.vector(as.character(temp2$Treatment))
sample_data(temp)$group_gio <- get_variable(temp, "Treatment") %in% groups_gio

temp2 <- as.data.frame(unique(sample_data(temp)[,6]))
fungicide <- as.vector(as.character(temp2$Fungicide))

temp2 <- as.data.frame(unique(sample_data(temp)[,13]))
inoculation <- as.vector(as.character(temp2$Inoculation))
merged_sp = merge_samples(temp, "Treatment")

merged_sp
sample_names(temp)
sample_names(merged_sp)

sample_data(merged_sp)$Treatment = sample_names(merged_sp)
sample_data(merged_sp)$group_gio = sample_names(merged_sp) %in% groups_gio
merged_sp

temp <- sample_data(merged_sp)$Fungicide
temp <- gsub(1, "Blad", temp); temp <- gsub(2, "Control", temp); temp <- gsub(3, "CuS", temp); temp <- gsub(4, "Systemics", temp)
sample_data(merged_sp)$Fungicide <- temp

temp <- sample_data(merged_sp)$Inoculation
temp <- gsub(3, "Pathogen + SK_ACEA1", temp); temp <- gsub(4, "SK_ACEA1", temp); temp <- gsub(2, "Pathogen", temp); temp <- gsub(1, "Water", temp)
sample_data(merged_sp)$Inoculation <- temp
sample_data(merged_sp)

###Merge unassigned taxa

sort(taxa_names(giovanni))

#IMPORTANT! Manually input any unwanted higher classes from the list produced by the sort() above
badtaxa = c("k__Fungi", "o__Agaricales", "o__Coniochaetales", "o__Pleosporales", "p__Ascomycota", "p__Basidiomycota",
            "c__Agaricomycetes", "c__Dothideomycetes", "c__Leotiomycetes", "c__Microbotryomycetes", "c__Sordariomycetes", 
            "c__Tremellomycetes")

gio_un_con <- merge_taxa(giovanni, badtaxa)
taxa_names(gio_un_con) <-  gsub("o__Agaricales", "Unassigned", taxa_names(gio_un_con)) 
head(tax_table(gio_un_con))
head(sample_data(gio_un_con))

sort(taxa_names(merged_sp))
#IMPORTANT! Manually input any unwanted higher classes from the list produced by the sort() above
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

###----------Get mean abundance + std---------------

#Make list of means, std and medians for all EVS in dataset

tissuegroups<- as.vector(as.matrix(unique(sample_data(giovanni)[,4])))
mean_std_list<- matrix(nrow=dim(otu_table(giovanni))[1], "")

for(i in 1:length(tissuegroups)){
  temp <- subset_samples(giovanni, Treatment==tissuegroups[i])
  name <- paste("Means_", tissuegroups[i], sep="")
  name <- gsub("\\s*\\([^\\)]+\\)_","",name)
  name2 <- paste(tissuegroups[i])
  name2 <- gsub("\\s*\\([^\\)]+\\)_","",name2) 
  
  x <- cbind(
    rownames(otu_table(temp)),
    rowMeans2(otu_table(temp)),
    rowSds(otu_table(temp)),
    rowMedians(otu_table(temp)),
    ""
  )
  colnames(x) = c(name2, "Mean", "St.Dev.", "Median", "")
  
  assign(name, x) 
  mean_std_list <- cbind(mean_std_list, x)
}

head(mean_std_list)

write.table(mean_std_list, file='EVS_mean_std_med.tsv', quote=FALSE, sep='\t')
###


###-----------------Colour palettes-----------------------

# Colour schemes to match colours used during sampling etc. 

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

###----------------------PERMANOVA---------------------------

#Transform data to relative abundance ("compositional")
pseq.rel <- microbiome::transform(giovanni, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

p <- microbiome::plot_landscape(pseq.rel, method = "NMDS", distance = "jaccard", col = "Inoculation", size = 3)
print(p)

permanova <- vegan::adonis(t(otu)~Fungicide+Inoculation+Treatment, data = meta, permutations=999, method = "jaccard")
write.table(as.data.frame(permanova$aov.tab), file='permanova_jaccard_all.tsv', quote=FALSE, sep='\t')
print(as.data.frame(permanova$aov.tab))

#post-hoc test

permanova <- vegan::adonis(t(otu)~Fungicide+Inoculation+Treatment, data = meta, permutations=999, method = "jaccard")
permanova$aov.tab

post_hoc_permanova_fung <- pairwise.adonis(t(otu), meta$Fungicide, sim.function = "vegdist",sim.method = "jaccard", p.adjust.m = "bonferroni", reduce = NULL, perm = 999)
post_hoc_permanova_fung
write.table(as.data.frame(post_hoc_permanova_fung), file='permanova_jaccard_fung_post_hoc_bonferroni.tsv', quote=FALSE, sep='\t')


post_hoc_permanova_ino <- pairwise.adonis(t(otu), meta$Inoculation, sim.function = "vegdist", sim.method = "jaccard", p.adjust.m = "bonferroni", reduce = NULL, perm = 999)
post_hoc_permanova_ino
write.table(as.data.frame(post_hoc_permanova_ino), file='permanova_jaccard_ino_post_hoc_bonferroni.tsv', quote=FALSE, sep='\t')


post_hoc_permanova_treatment <- pairwise.adonis(t(otu), meta$Treatment, sim.function = "vegdist", sim.method = "jaccard", p.adjust.m = "bonferroni", reduce = NULL, perm = 999)
post_hoc_permanova_treatment
write.table(as.data.frame(post_hoc_permanova_treatment), file='permanova_jaccard_treatment_post_hoc_bonferroni.tsv', quote=FALSE, sep='\t')

###-----------Transpose phyloseq objects---------------

#Some of the commands in this packae requires a transposed matrix, so first, we'll transpose our phyloseq obejct:
# Extract abundance matrix from the phyloseq object
OTU2 <- otu_table(giovanni, taxa_are_rows = TRUE)
OTU2 <- t(OTU2)
#Merge new phyloseq object
TAX = as.matrix(tax_table(giovanni))
met <- as.data.frame(sample_data(giovanni))
pseq1 <- phyloseq(OTU2, TAX)
MBS_physeq<-merge_phyloseq(pseq1, met)


#MBS_physeq is now the transposed version of giovanni. 

rm(OTU2, TAX, met, pseq1)

###------------Subset transposed phyloseq objects-------------------

#subset transposed phyloseq object into inoculations
ID<- as.vector(unique(sample_data(MBS_physeq)$Inoculation))
for (i in 1:length(ID)){ 
  temp <- subset_samples(MBS_physeq, Inoculation==ID[i])
  temp = prune_taxa(taxa_sums(temp) > 0, temp)
  name <- paste("MBS_", ID[i], sep = "")
  assign(name, temp)
}

#subset transposed phyloseq object into fungicide
ID<- as.vector(unique(sample_data(MBS_physeq)$Fungicide))
for (i in 1:length(ID)){ 
  temp <- subset_samples(MBS_physeq, Fungicide==ID[i])
  temp = prune_taxa(taxa_sums(temp) > 0, temp)
  name <- paste("MBS_", ID[i], sep = "")
  assign(name, temp)
}

###-----------alpha diversity-----------

grouping_column ="Fungicide"
pValueCutoff=0.05


#Inoculation subsets: interchange "physeq =" with:
 physeq = MBS_No_inoculation
# physeq = MBS_SK_ACEA1
# physeq = MBS_Pch
# physeq = MBS_Pch_SK_ACEA1

temp1 <- estimate_richness(physeq, split = TRUE, measures = "Shannon")
measure <- rep("Shannon", dim(temp1)[1])
fung <- sample_data(physeq)[,6]
#sample, value, measure, fungicide
temp3 <- cbind(rownames(temp1),data.frame(temp1, row.names=NULL), measure, fung)
colnames(temp3) = c("sample", "value", "measure", "Fungicide"); colnames(temp3)

a1 <- aov(temp3$value ~ temp3$Fungicide)
posthoc <- TukeyHSD(x=a1, 'temp3$Fungicide', conf.level=0.95)
posthoc_shan <- as.data.frame(posthoc$`temp3$Fungicide`)
colnames(posthoc_shan) =c("diff","lwr","upr","padj"); posthoc_shan

#Save this - remember to change filename according to physeq object!
write.table(as.data.frame(posthoc_shan), file='shannon_posthoc_fungicide_noino.tsv', quote=FALSE, sep='\t')


posthoc_shan2 <- cbind(measure = paste0("Shannon"), rownames(posthoc_shan), data.frame(posthoc_shan, row.names=NULL))
colnames(posthoc_shan2)[2] = "name"; colnames(posthoc_shan2)
posthoc_shan2 <- posthoc_shan2 %>%
  separate(name, c("from", "to"), "-")
ph_shan<-posthoc_shan2[!(posthoc_shan2$padj>0.05),]; ph_shan


#evenness

temp2 <- evenness(physeq, index = 'pielou')
measure <- rep("Pielou's evenness", dim(temp2)[1])
fung <- sample_data(physeq)[,6]
#sample, value, measure, fungicide
temp4 <- cbind(rownames(temp1),data.frame(temp2, row.names=NULL), measure, fung)
colnames(temp4) = c("sample", "value", "measure", "Fungicide"); colnames(temp4)

a1 <- aov(temp4$value ~ temp4$Fungicide)
posthoc <- TukeyHSD(x=a1, 'temp4$Fungicide', conf.level=0.95)
posthoc_even <- as.data.frame(posthoc$`temp4$Fungicide`)
colnames(posthoc_even) =c("diff","lwr","upr","padj"); posthoc_even

#Save this - remember to change filename according to physeq object!
write.table(as.data.frame(posthoc_shan), file='evenness_posthoc_fungicide_noino.tsv', quote=FALSE, sep='\t')


posthoc_even2 <- cbind(measure = paste0("Pielou's evenness"), rownames(posthoc_even), data.frame(posthoc_even, row.names=NULL))
colnames(posthoc_even2)[2] = "name"; colnames(posthoc_even2)
posthoc_even2 <- posthoc_even2 %>%
  separate(name, c("from", "to"), "-")
ph_even<-posthoc_even2[!(posthoc_even2$padj>0.05),];ph_even

df_pw <- rbind(ph_shan, ph_even); df_pw

df <-  rbind(temp3, temp4)

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
p<-p+theme_bw()
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("")+xlab("")
p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("")

#This loop will generate the lines and signficances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,
                 aes(x,y),
                 data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),
                                         which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), 
                                   y = c((2-(i*0.1)),(2-(i*0.1))),
                                   #y = c(as.numeric(as.character(df_pw[i,"upr"])),as.numeric(as.character(df_pw[i,"upr"]))), 
                                   measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
                 color="black",lineend = "butt",
                 arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  
  p<-p+geom_text(inherit.aes=F,
                 aes(x=x,y=y,label=label),
                 data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+
                                      which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,
                                 y = (2-(i*0.1)),
                                 #y=as.numeric(as.character(df_pw[i,"upr"])),
                                 measure=as.character(df_pw[i,"measure"]),
                                 label=as.character(cut(as.numeric(as.character(df_pw[i,"padj"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}
p

#Polish off plot

#REMEMBER! Change to correct plot title.

p <- p +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values=set_fun)+
   ggtitle("No Inoculation")+
  # ggtitle("SK-ACEA1")+
  # ggtitle("Pathogen")+
  # ggtitle("Pathogen + SK-ACEA1")+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(ylim = c(0,2))
p

#REMEMBER! Save each plot as a different object: 

 p_n <- p; p_n #No inoculation
# p_m <- p; p_m #Mycobiome SK-ACEA1
# p_p <- p; p_p #Pathogen
# p_pm <- p; p_pm #Pathogen + mycobiome 

p <- ggarrange(p_m,p_n,p_pm,p_p,
               widths = 1, heights = 1.1,
               ncol = 2, nrow = 2,
               common.legend = TRUE, legend = "bottom")
p
ggsave(p, file="Richness_ino_fun.pdf", width = 30, height = 20, units = "cm")


#---


#Fungicide subsets: interchange "physeq =" with:

physeq = MBS_Control
# physeq = MBS_Blad
# physeq = MBS_CuS
# physeq = MBS_Systemics

grouping_column ="Inoculation"
pValueCutoff=0.05

temp1 <- estimate_richness(physeq, split = TRUE, measures = "Shannon")
measure <- rep("Shannon", dim(temp1)[1])
ino <- sample_data(physeq)[,13]
#sample, value, measure, Inoculation
temp3 <- cbind(rownames(temp1),data.frame(temp1, row.names=NULL), measure, ino)
colnames(temp3) = c("sample", "value", "measure", "Inoculation"); colnames(temp3)

a1 <- aov(temp3$value ~ temp3$Inoculation)
posthoc <- TukeyHSD(x=a1, 'temp3$Inoculation', conf.level=0.95)
posthoc_shan <- as.data.frame(posthoc$`temp3$Inoculation`)
colnames(posthoc_shan) =c("diff","lwr","upr","padj"); posthoc_shan

#Save this - remember to change filename according to physeq object!
write.table(as.data.frame(posthoc_shan), file='shannon_posthoc_inoculation_ctrl.tsv', quote=FALSE, sep='\t')


posthoc_shan2 <- cbind(measure = paste0("Shannon"), rownames(posthoc_shan), data.frame(posthoc_shan, row.names=NULL))
colnames(posthoc_shan2)[2] = "name"; colnames(posthoc_shan2)
posthoc_shan2 <- posthoc_shan2 %>%
  separate(name, c("from", "to"), "-")
ph_shan<-posthoc_shan2[!(posthoc_shan2$padj>0.05),]; ph_shan


#evenness

temp2 <- evenness(physeq, index = 'pielou')
measure <- rep("Pielou's evenness", dim(temp2)[1])
ino <- sample_data(physeq)[,13]
#sample, value, measure, Inoculation
temp4 <- cbind(rownames(temp1),data.frame(temp2, row.names=NULL), measure, ino)
colnames(temp4) = c("sample", "value", "measure", "Inoculation"); colnames(temp4)

a1 <- aov(temp4$value ~ temp4$Inoculation)
posthoc <- TukeyHSD(x=a1, 'temp4$Inoculation', conf.level=0.95)
posthoc_even <- as.data.frame(posthoc$`temp4$Inoculation`)
colnames(posthoc_even) =c("diff","lwr","upr","padj"); posthoc_even

#Save this - remember to change filename according to physeq object!
write.table(as.data.frame(posthoc_even), file='evenness_posthoc_inoculation_ctrl.tsv', quote=FALSE, sep='\t')


posthoc_even2 <- cbind(measure = paste0("Pielou's evenness"), rownames(posthoc_even), data.frame(posthoc_even, row.names=NULL))
colnames(posthoc_even2)[2] = "name"; colnames(posthoc_even2)
posthoc_even2 <- posthoc_even2 %>%
  separate(name, c("from", "to"), "-")
ph_even<-posthoc_even2[!(posthoc_even2$padj>0.05),];ph_even

df_pw <- rbind(ph_shan, ph_even); df_pw

df <-  rbind(temp3, temp4)

#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column),data=df)
p<-p+geom_boxplot()+geom_jitter(position = position_jitter(height = 0, width=0))
p<-p+theme_bw()
p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))
p<-p+facet_wrap(~measure,scales="free_y",nrow=1)+ylab("Observed Values")+xlab("")
p<-p+theme(strip.background = element_rect(fill = "white"))+xlab("")
#p

#This loop will generate the lines and signficances
for(i in 1:dim(df_pw)[1]){
  p<-p+geom_path(inherit.aes=F,
                 aes(x,y),
                 data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),
                                         which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), 
                                   y = c((2-(i*0.1)),(2-(i*0.1))),
                                   #y = c(as.numeric(as.character(df_pw[i,"upr"])),as.numeric(as.character(df_pw[i,"upr"]))), 
                                   measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), 
                 color="black",lineend = "butt",
                 arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
  
  p<-p+geom_text(inherit.aes=F,
                 aes(x=x,y=y,label=label),
                 data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+
                                      which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,
                                 y = (2-(i*0.1)),
                                 #y=as.numeric(as.character(df_pw[i,"upr"])),
                                 measure=as.character(df_pw[i,"measure"]),
                                 label=as.character(cut(as.numeric(as.character(df_pw[i,"padj"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))))
}
p

#Polish off plot

#REMEMBER! Change to correct plot title.

p <- p +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values=set_fun)+
  ggtitle("Control")+
  # ggtitle("Blad")+
  # ggtitle("CuS")+
  # ggtitle("Systemics")+
  theme(plot.title = element_text(hjust = 0.5))+
  coord_cartesian(ylim = c(0,2))
p

#REMEMBER! Save each plot as a different object: 

p_c <- p; p_c #Control
# p_b <- p; p_b #Blad
# p_cus <- p; p_cus #CuS
# p_s <- p; p_s #Systemics 

p <- ggarrange(p_c,p_b,p_cus,p_s,
               widths = 1, heights = 1.1,
               ncol = 2, nrow = 2,
               common.legend = TRUE, legend = "bottom")
p
ggsave(p, file="Richness_fun_ino.pdf", width = 30, height = 20, units = "cm")


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


###--------------------------LCBD---------------------------

tax_table(gio_m_un_con)[3,6] = "g__Epicoccum"
tax_table(gio_m_un_con)[3,7] = "s__Epicoccum_nigrum"
taxa_names(gio_m_un_con)[3] = "s__Epicoccum_nigrum"

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


#------------DeSeq2-------------

tax_table(MBS_physeq)[3,6] = "g__Epicoccum"
tax_table(MBS_physeq)[3,7] = "s__Epicoccum_nigrum"
taxa_names(MBS_physeq)[3] = "s__Epicoccum_nigrum"


OTU2 <- otu_table(gio_filter1, taxa_are_rows = TRUE)
OTU2 <- t(OTU2)
#Merge new phyloseq object
TAX = as.matrix(tax_table(gio_filter1))
met <- as.data.frame(sample_data(gio_filter1))
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

###---------------Co-occurence networks--------------------

#co-occurence network are a nice way to illustrate which OTUs are likely to occur in the same samples, 
#and to illustrate patterns. 

tmp <- taxa_level(MBS_physeq, "Genus")

tmp <- MBS_physeq
tmp


temp = gsub(pattern = "s__", replacement = "", x = taxa_names(tmp))
temp = as.character(make.unique(temp))
taxa_names(tmp) = temp
tmp

co_occrBlad<- co_occurence_network(tmp, grouping_column = "Fungicide", rhos = 0.35, method="cor", 
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

###-----------Create control column for heat tree------------

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

# Convert counts to proportions, then calculates per-taxon proportions and per-group occurence 
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
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
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_treatment_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_treatment_pval.tsv', quote=FALSE, sep='\t')

### Control

tissuegroup <- obj$data$sample_data$Control

# Convert counts to proportions, then calculates per-taxon proportions and per-group occurence 
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Note the heat_tree and not heat_tree_matrix here, as we're comparing two groups - control vs. everything else.
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
  mutate(taxon_names = taxon_names(obj)[taxon_id], taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
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

# Convert counts to proportions, then calculates per-taxon proportions and per-group occurence 
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
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

tax_table(gio_filter1)[3,6] = "g__Epicoccum"
tax_table(gio_filter1)[3,7] = "s__Epicoccum_nigrum"
taxa_names(gio_filter1)[3] = "s__Epicoccum_nigrum"

test <- gio_filter1

tmp <- normalise_data(test, norm.method = "relative")
tmp

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

# Convert counts to proportions, then calculates per-taxon proportions and per-group occurence 
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
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

# Convert counts to proportions, then calculates per-taxon proportions and per-group occurence 
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
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

### Pch + SK-ACEA1

obj <- parse_phyloseq(Filter1_relab_Pch_SK_ACEA1)
tissuegroup <- obj$data$sample_data$Fungicide

# Convert counts to proportions, then calculates per-taxon proportions and per-group occurence 
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
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
                 output_file = 'heattree_pch_Sk.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_pch_Sk.tsv', quote=FALSE, sep='\t')

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
                 output_file = 'heattree_pch_Sk_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_pch_Sk_pval.tsv', quote=FALSE, sep='\t')

### SK-ACEA1

obj <- parse_phyloseq(Filter1_relab_SK_ACEA1)
tissuegroup <- obj$data$sample_data$Fungicide

# Convert counts to proportions, then calculates per-taxon proportions and per-group occurence 
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)
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
                 output_file = 'heattree_Sk.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_Sk.tsv', quote=FALSE, sep='\t')


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
                 output_file = 'heattree_Sk_pval.pdf')


test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_Sk_pval.tsv', quote=FALSE, sep='\t')


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

###----------------End of code-------------------
