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
setwd("D:/Marie/Documents/PhD/Giovanni_statistics/")
#Set user directory
uzdir <- "D:/Marie/Documents/PhD/Giovanni_statistics/"

###---------------------LOAD DATA-------------------------------------

#Biom file to import - preferably JSON

data_biom <- "table_1_new.biom"

#Create file path, removing any blank spaces
biom_file <- paste(uzdir, data_biom, sep="")
#the actual act of importing the file into R
biom_otu_tax <- import_biom(biom_file)

#inspect file, prune colnmes as needed
head(biom_otu_tax)
#colnames(biom_otu_tax)[1:15] <- paste("C", colnames(biom_otu_tax)[1:15], sep = "")
#head(biom_otu_tax)

#Import taxonomy table - remember to prune metadata so that there are only one type of delimiter. 
tax1 <- as.matrix(read.table(file = 'taxonomy_ITS_genus_unassigned2.csv', sep = '\t', header = TRUE))

rownames(tax1) <- tax1[,1]; head(tax1)
rank_names(tax1)


#Check contents of new table

#otherwise, use "taxmat = as.matrix(#TABLE FILE#, rownames.force=TRUE)"

#Convert table to a Phyloseq taxonomy table
TAX = tax_table(tax1); head(TAX)
#remove tax1 to declutter environment
rm(tax1)

#Combine OTU table and taxonomy table into a phyloseq object. 
physeq <- phyloseq(biom_otu_tax, TAX); physeq
#Inspect
head(tax_table(physeq))
head(otu_table(physeq))
rank_names(tax_table(physeq))

taxa_names(physeq) <- tax_table(physeq)[,13]; head(tax_table(physeq)) 
rank_names(tax_table(physeq))

#Retain only columns with taxonomic data - pick correct columns, check file first!
#tax_table(physeq) <- tax_table(physeq)[,-10]; head(tax_table(physeq))
tax_table(physeq) <- tax_table(physeq)[,3:12]; head(tax_table(physeq)) 
rank_names(tax_table(physeq))

#Import metadata file, removing blank spaces in path
md_file <- paste(uzdir, "Metadata.tsv", sep="")   
metadata <- import_qiime_sample_data(md_file); head(metadata)
#metadata$X.SampleID <- gsub("-", "_", metadata$X.SampleID); head(metadata)

##If the above doesn't produce a table, use this: 
#metadata <- read.csv(md_file, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
#metadata <- sample_data(metadata2)
#rownames(metadata2) = metadata2$SampleID
##Use whichever one works. I haven't been able to figure out why the import_qiime_sample_data() 
## only works sometimes. Probably a unicode/tabulator thing.  

#Merge OTU data and metadata to single object - The Phyloseq Object(tm)
giovanni <- merge_phyloseq(physeq,metadata)

#sample_names(giovanni) <- gsub("-", "_", sample_names(giovanni)); sample_names(giovanni)

 

giovanni

raw_data <- giovanni

#giovanni is now the s4 object - the phyloseq object - we will use throughout the rest of the code. 
#In theory, all other objects in the global environment could be removed now, but I tend to leave them, just in case I need them.
#Use sample_names and taxa_names to change names of samples and taxa. 

#---------Inspect dataset-----------

#inspect contents of Phyloseq metaobjects
sample_data(giovanni)[1:5]
otu_table(giovanni) [1:5]
tax_table(giovanni)[1:5]

#(R practice)
#inspect number of observations in each sample.
colSums(otu_table(giovanni)); plot(colSums(otu_table(giovanni)))
#rowSums(otu_table(giovanni)); plot(rowSums(otu_table(giovanni)))
#inspect first ten values of first column of tax_table.
tax_table(giovanni) [,1][1:10]
#find number of unique values/names in rows of tax_table:
unique(sample_data(giovanni) [,2])
length(unique(tax_table(giovanni)[,10]))
length(unique(tax_table(giovanni)[,7]))

#/data inspection

###------------------CONGLOMERATING-------------------------------------

##IMPORTANT!!!
##Half of these commands (the ones removing and naming columns) are inputfile-specific!!!
##Inspect and adjust accordingly! 

#Conglomerating tips to species level
gio_sp <- giovanni
rank_names(gio_sp)
tax_table(gio_sp) <- tax_table(gio_sp)[,-7]; head(tax_table(gio_sp))
tax_table(gio_sp) <- tax_table(gio_sp)[,-7]; head(tax_table(gio_sp))
rank_names(gio_sp)
gio_sp <- tax_glom(gio_sp, taxrank="Name_sp"); gio_sp
taxa_names(gio_sp) <- tax_table(gio_sp)[,8]; taxa_names(gio_sp)[1:10]
tax_table(gio_sp) <- tax_table(gio_sp)[,-8]; head(tax_table(gio_sp))

sample_data(gio_sp)[1:5]
otu_table(gio_sp) [1:5]
tax_table(gio_sp)[1:5]

#Conglomerating tips to genus level
gio_g <- giovanni
rank_names(gio_g)
tax_table(gio_g) <- tax_table(gio_g)[,1:8]; head(tax_table(gio_g))
length(unique(tax_table(gio_g)[,7]))
taxa_names(gio_g) <- tax_table(gio_g)[,8]; taxa_names(gio_g)[1:10]
tax_table(gio_g) <- tax_table(gio_g)[,-8]; head(tax_table(gio_g))
gio_g <- tax_glom(gio_g, taxrank="Name_g"); gio_g
rank_names(gio_g)
taxa_names(gio_g) <- tax_table(gio_g)[,7]; taxa_names(gio_g)[1:10]
#write.table(tax_table(gio_g), file='gio_g_tax.tsv', quote=FALSE, sep='\t')
tax_table(gio_g) <- tax_table(gio_g)[,-7]; head(tax_table(gio_g))

gio_g_un <- gio_g
tax_table(gio_g_un) <-  gsub("s__unidentified", "", tax_table(gio_g_un)); 
tax_table(gio_g_un) <-  gsub("g__unidentified", "", tax_table(gio_g_un)); 
tax_table(gio_g_un) <-  gsub("f__unidentified", "", tax_table(gio_g_un));
tax_table(gio_g_un) <-  gsub("o__unidentified", "", tax_table(gio_g_un));
tax_table(gio_g_un) <-  gsub("c__unidentified", "", tax_table(gio_g_un));
tax_table(gio_g_un) <-  gsub("p__unidentified", "", tax_table(gio_g_un));
tax_table(gio_g_un)


sample_data(gio_g)[1:5]
otu_table(gio_g) [1:5]
tax_table(gio_g)[1:5]

gio_f <- tax_glom(gio_g, taxrank="Family"); gio_f

###--------------CHOOSING OBJECT AS PRIMARY PHYSEQ-------------------

#From now on, "giovanni" will be used as the nominator for the physeq being used.
gio_sp_1<- gio_sp
gio_sp_1 <-subset_samples(gio_sp, retain_discard_1=="retain")
gio_sp_1 = prune_taxa(taxa_sums(gio_sp_1) > 0, gio_sp_1); gio_sp_1
gio_sp_2 <- subset_samples(gio_sp, retain_discard_2=="retain")
gio_sp_2 = prune_taxa(taxa_sums(gio_sp_2) > 0, gio_sp_2); gio_sp_2
gio_sp_3 <- subset_samples(gio_sp, retain_discard_3=="retain")
gio_sp_3 = prune_taxa(taxa_sums(gio_sp_3) > 0, gio_sp_3); gio_sp_3

gio_g_1<- gio_g_un
gio_g_1 <-subset_samples(gio_g, retain_discard_1=="retain")
gio_g_1 = prune_taxa(taxa_sums(gio_g_1) > 0, gio_g_1); gio_g_1
gio_g_2 <- subset_samples(gio_g, retain_discard_2=="retain")
gio_g_2 = prune_taxa(taxa_sums(gio_g_2) > 0, gio_g_2); gio_g_2
gio_g_3 <- subset_samples(gio_g, retain_discard_3=="retain")
gio_g_3 = prune_taxa(taxa_sums(gio_g_3) > 0, gio_g_3); gio_g_3
 
gio_g3_c <- subset_samples(gio_g_3, trunk_canes_1=="canes")
gio_g3_c = prune_taxa(taxa_sums(gio_g3_c) > 0, gio_g3_c); gio_g3_c
gio_g3_a <- subset_samples(gio_g_3, trunk_canes_1=="trunk")
gio_g3_a = prune_taxa(taxa_sums(gio_g3_a) > 0, gio_g3_a); gio_g3_a


temp <- as.data.frame(sample_data(gio_sp_3)[,7])
temp$Tissue_3 <- ifelse(grepl("arm", temp$Description_3), "arm", "cane")
sample_data(gio_sp_3)[,8] <- temp$Tissue_3

gio_sp3_c <- subset_samples(gio_sp_3, trunk_canes_1=="canes")
gio_sp3_c = prune_taxa(taxa_sums(gio_sp3_c) > 0, gio_sp3_c); gio_sp3_c
gio_sp3_a <- subset_samples(gio_sp_3, trunk_canes_1=="trunk")
gio_sp3_a = prune_taxa(taxa_sums(gio_sp3_a) > 0, gio_sp3_a); gio_sp3_a

##Retained for reference purposes
#This for loop creates four subset from treatments ("Fungicide" from older dataset). 
#We create a list of the names of the fungicides, and then use that to specify the subset command. 
#Then we paste the name of the fungicide onto the subset and assign it as a separate phyloseq object.

#ID<- as.vector(unique(sample_data(giovanni)$Fungicide))
#for (i in 1:length(ID)){ 
#  temp <- subset_samples(giovanni, Fungicide==ID[i])
#  name <- paste(ID[i])
#  assign(name, temp)
#}

## Or you can make subsets of treatments. Specifies which metadata column, and which value in said column should be criteria. 
#Subset <- subset_samples(physeq, data_column=="value to subset - e.g. control")

###----------Filter by reads------------------------

temp <- gio_sp_1

otu_table(temp)[otu_table(temp)<25 ] <- 0
y = as.matrix(sort((x / sum(x))*100, decreasing = TRUE))
write.table(y, file='taxalist_rel_ab_tissue_pt1.tsv', quote=FALSE, sep='\t')

###----------Filter by relative abundance-----------

#Filters out anything below 0.01% of total abundance. 
minTotRelAbun = 0.1 
x = taxa_sums(gio_sp_1)
keepTaxa = rownames(as.data.frame(which(((x / sum(x))*100) > minTotRelAbun)))
gio_sp1_filter1 = prune_taxa(keepTaxa, gio_sp_1)

###----------Get mean abundance + std---------------

#Actual abundance

temp <- normalise_data(merged_sp_2, norm.method = "relative")
tmp <- boxplot(t(otu_table(merged_sp_2)))
temp1 <- t(otu_table(gio_sp_1))

names_list <- as.vector(colnames(temp1))
names_list <- names_list[!is.na(names_list)]

mean_std <- list()
mean_std_list<- matrix(nrow=0, ncol=2)
colnames(mean_std_list) = c("Mean", "St.Dev.")

for (i in 1:length(names_list)){

Mean <- mean(temp1[,i])
St_Dev <- sqrt(var(temp1[,i]))

mean_std_list_int <-cbind(Mean, St_Dev)
mean_std_list <- rbind(mean_std_list_int, mean_std_list) 
}

norm <- mean_std_list
norm

write.table(mean_std_list, file='mean_std_2.tsv', quote=FALSE, sep='\t')

#Abundance, 

boxplot(t(otu_table(merged_sp_3)))
temp1 <- t(otu_table(merged_sp_3))

names_list <- as.vector(colnames(temp1))
names_list <- names_list[!is.na(names_list)]

mean_std <- list()
mean_std_list<- matrix(nrow=0, ncol=2)


for (i in 1:length(names_list)){
  
  Mean <- mean(temp1[,i])
  St_Dev <- sqrt(var(temp1[,i]))
  
  mean_std_list_int <-cbind(Mean, St_Dev)
  mean_std_list <- rbind(mean_std_list_int, mean_std_list) 
}

colnames(mean_std_list) = c("Mean", "St.Dev.")
mean_std_list

write.table(mean_std_list, file='mean_std_3.tsv', quote=FALSE, sep='\t')

#EVS per sample

#Make list of means, std and medians for all EVS in dataset 1

tmp <- as.matrix(unique(sample_data(gio_sp_1)[,3]))
tissuegroups<- as.vector(tmp)
mean_std_list<- matrix(nrow=341, ncol=0)

for(i in 1:length(tissuegroups)){
  temp <- subset_samples(gio_sp_1, trunk_canes_1==tissuegroups[i])
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

write.table(mean_std_list, file='EVS_mean_std_med_1.tsv', quote=FALSE, sep='\t')

#Make list of means, std and medians for all EVS in dataset 2

tmp <- as.matrix(unique(sample_data(gio_sp_2)[,5]))
tissuegroups<- as.vector(tmp)
mean_std_list<- matrix(nrow=296, ncol=0)

for(i in 1:length(tissuegroups)){
  temp <- subset_samples(gio_sp_2, Tissue_2==tissuegroups[i])
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

write.table(mean_std_list, file='EVS_mean_std_med_2.tsv', quote=FALSE, sep='\t')

#Make list of means, std and medians for all EVS in dataset 3

tmp <- as.matrix(unique(sample_data(gio_sp_3)[,7]))
tissuegroups<- as.vector(tmp)
mean_std_list<- matrix(nrow=219, ncol=0)

for(i in 1:length(tissuegroups)){
  temp <- subset_samples(gio_sp_3, Description_3==tissuegroups[i])
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

write.table(mean_std_list, file='EVS_mean_std_med_3.tsv', quote=FALSE, sep='\t')
###


for (i in 1:length(names_list)){
  
  Mean <- mean(temp1[,i])
  St_Dev <- sqrt(var(temp1[,i]))
  
  mean_std_list_int <-cbind(Mean, St_Dev)
  mean_std_list <- rbind(mean_std_list_int, mean_std_list) 
}


###----------Colour palettes------------

# I will from now on use a colour scheme that
# red = symptomatic, 
# blue = asymptomatic, 
# purple = anything that combines the two (such as trunks in general), 
# orange/yellow for anything to do with the canes, 
# brown/burgundy for anything that combines trunks and canes.

#For trunks vs canes
set_2 <- c("purple", "orange")
#trunks vs canes, metacoder
set_3p <- c("purple", "maroon", "orange")
set_3w <- c("purple", "white", "orange")
set_3_rb <- c("royalblue", "purple","red")
#For leaf symptoms
set_6 <- c("royalblue", "steelblue2", "skyblue", "purple", "magenta", "red")
#For tissue groups
set_8 <- c("darkblue","blue","purple", "maroon", "darkorange2", "tomato3", "orange", "gold")

###------------------------Richness---------------------

#Pielou's evenness = Shannon evenness. 
p_an <-plot_anova_diversity_pval(merged_taxunassigned1, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="trunk_canes_1",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_an_o+scale_colour_manual(values = set_2)
p1 <- p_an_o+scale_colour_manual(values = set_2)
ggsave(p1, file="Richness_1sp.pdf", width = 30, height = 20, units = "cm")
write.table(p_an_pvalues, file='richness_1sp_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_o$data, file='richness_1sp_data.tsv', quote=FALSE, sep='\t')

p_an <-plot_anova_diversity_pval(gio_sp_2, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="Tissue_2",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_an_o+scale_colour_manual(values = set_8)
p1 <- p_an_o+scale_colour_manual(values = set_8)
ggsave(p1, file="Richness_2sp.pdf", width = 30, height = 20, units = "cm")
write.table(p_an_pvalues, file='richness_2sp_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_o$data, file='richness_2sp_data.tsv', quote=FALSE, sep='\t')

p_an <-plot_anova_diversity_pval(gio_sp3_a, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="Description_3",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_an_o+scale_colour_manual(values = set_3_rb)
p1 <- p_an_o+scale_colour_manual(values = set_3_rb)
ggsave(p1, file="Richness_3sp_a.pdf", width = 30, height = 20, units = "cm")
write.table(p_an_pvalues, file='richness_3spa_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_o$data, file='richness_3spa_data.tsv', quote=FALSE, sep='\t')

p_an <-plot_anova_diversity_pval(gio_sp3_c, method = c("shannon", "simpson", "evenness"),
                                 grouping_column ="Description_3",pValueCutoff=0.05)
p_an_o <- p_an[[1]]; p_an_pvalues <- p_an[[2]]
p_an_o <- p_an_o + ggtitle("Diversity index + evenness \n * = p<0.05, ** = p<0.01, *** = p<0.001") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_an_o+scale_colour_manual(values = set_3_rb)
p1 <- p_an_o+scale_colour_manual(values = set_3_rb)
ggsave(p1, file="Richness_3sp_c.pdf", width = 30, height = 20, units = "cm")
write.table(p_an_pvalues, file='richness_3spc_pval.tsv', quote=FALSE, sep='\t')
write.table(p_an_o$data, file='richness_3spc_data.tsv', quote=FALSE, sep='\t')

###-----------Metacoder----------

#Metacoder provides heat trees, which colour the branches of the tree that 
# has any significant change in abudance between two treatments. 

tmp <- tax_glom(giovanni_2, taxrank = "Genus")

temp <- gio_g_2
tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)

obj_bas <- subset_taxa(temp, Phylum=="p__Basidiomycota")
tax_table(obj_bas) <- tax_table(obj_bas)[,2:7]
head(tax_table(obj_bas))

obj_asco <- subset_taxa(temp, Phylum=="p__Ascomycota")
tax_table(obj_asco) <- tax_table(obj_asco)[,2:7]
head(tax_table(obj_asco))

#parse phyloseq object giovanni
obj_all <- parse_phyloseq(temp)
obj_asco <- parse_phyloseq(obj_asco)
obj_bas <- parse_phyloseq(obj_bas)

#Remove hashtag in front of whichever subset you wanna run the rest of this next bit of code on.
obj <- obj_all
#obj <- obj_asco
#obj <- obj_bas

#Use:
#tissuegroup <- obj$data$sample_data$trunk_canes_1
tissuegroup <- obj$data$sample_data$Tissue_2
#tissuegroup <- obj$data$sample_data$Description_3

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
#Plot p-values
hist(obj$data$diff_table$wilcox_p_value); obj$data$diff_table

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
#          node_color_range = set_3_rb,
          node_color_range = c("purple", "maroon", "orange"),
          node_size_axis_label = "species count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
#          title = "Red = Higher Abundance in Syptomatic tissue",
          title = "Orange = Higher Abundance in Canes\n Purple = Higher Abundance in Trunks",
          title_size = 0.05)

#set any difference that is not significant to zero 
test <- within(obj$data$diff_table, log2_median_ratio[log2_median_ratio == "Inf"] <- 0 )
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test
obj$data$diff_table

##Create heat tree matrix - if you use obj$data$diff_table, treatment1 corresponds to higher abundance in last specified colour. 

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          node_color_range = c("purple", "grey95", "orange"),
          node_size_axis_label = "species count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
          title = "Orange Nodes = Higher Abundance in Canes, p<0.05 \n Purple Nodes = Higher Abundance in Trunks, p<0.05 ",
          title_size = 0.05,
          edge_label_size_range = c(1,5) )


heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
 #                node_color_range = set_3_rb,
                  node_color_range = c("purple", "gray95", "orange"),
                 node_size_axis_label = "species count",
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_color_axis_label = "Log2 ratio median proportions",
                 edge_label_color_range = set_3_rb)
#pt 3

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_range = c("blue", "gray95", "red"),
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 edge_label_color_range = c("purple", "gray90", "orange"))
#

###------------------Heat tree pt. 1---------------

temp <- tax_glom(gio_sp1_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)

#parse phyloseq object giovanni
obj_all <- parse_phyloseq(temp)

obj <- obj_all
tissuegroup <- obj$data$sample_data$trunk_canes_1

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          node_color_range = c("purple", "grey95", "orange"),
          node_size_axis_label = "species count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
          title = "Orange Nodes = Higher Abundance in Canes \n Purple Nodes = Higher Abundance in Trunks",
          title_size = 0.05,
          edge_label_size_range = c(1,5),
          output_file = 'heattree_1.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_1.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pt. 2---------------

temp <- tax_glom(gio_sp2_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)


obj_bas <- subset_taxa(temp, Phylum=="p__Basidiomycota")
tax_table(obj_bas) <- tax_table(obj_bas)[,2:7]
head(tax_table(obj_bas))

obj_asco <- subset_taxa(temp, Phylum=="p__Ascomycota")
tax_table(obj_asco) <- tax_table(obj_asco)[,2:7]
head(tax_table(obj_asco))

#parse phyloseq object giovanni
obj_all <- parse_phyloseq(temp)
obj_asco <- parse_phyloseq(obj_asco)
obj_bas <- parse_phyloseq(obj_bas)

obj <- obj_all
tissuegroup <- obj$data$sample_data$Tissue_2

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

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
                 node_color_range = c("darkblue", "royalblue", "gray95", "orange", "darkorange2"),
                 #title = "Heat Tree Matrix, pt. 2, tissue groups",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_2.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_2.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pt. 2, basidiomycota---------------

obj <- obj_bas
tissuegroup <- obj$data$sample_data$Tissue_2

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

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
                 node_color_range = c("darkblue", "royalblue", "gray95", "orange", "darkorange2"),
                 #title = "Heat Tree Matrix, pt. 2, basidiomycota",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_2bas.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_2bas.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pt. 2, ascomycota---------------

obj <- obj_asco
tissuegroup <- obj$data$sample_data$Tissue_2

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")


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
                 node_color_range = c("darkblue", "royalblue", "gray95", "orange", "darkorange2"),
                 #title = "Heat Tree Matrix, pt. 2, ascomycota",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_2asco.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_2asco.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pt. 3 arms---------------

temp <- tax_glom(gio_sp3a_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)

obj_all <- parse_phyloseq(temp)
obj <- obj_all
tissuegroup <- obj$data$sample_data$Description_3

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

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
                 node_color_range = c("blue", "gray95", "red"),
                 title = "Heat Tree Matrix, pt. 3, Arms",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_3a.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_3a.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pt. 3 canes---------------

temp <- tax_glom(gio_sp3c_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)
#parse phyloseq object giovanni
obj_all <- parse_phyloseq(temp)
obj <- obj_all
tissuegroup <- obj$data$sample_data$Description_3

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")

heat_tree_matrix(obj,
                 data = "diff_table",
                 node_size = n_obs,
                 node_label = taxon_names,
                 node_color = log2_median_ratio,
                 node_color_trans = "linear",
                 node_color_interval = c(-3, 3),
                 edge_color_interval = c(-3, 3),
                 node_color_range = c("blue", "gray95", "red"),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 title = "Heat Tree Matrix, pt. 3, Canes",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_3c.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_3c.tsv', quote=FALSE, sep='\t')


###

###-----------------HEAT TREES, P-VALUES---------


###------------------Heat tree pval pt. 1---------------

temp <- tax_glom(gio_sp1_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)

#parse phyloseq object giovanni
obj_all <- parse_phyloseq(temp)

obj <- obj_all
tissuegroup <- obj$data$sample_data$trunk_canes_1

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

#Test p-values
obj$data$diff_table$wilcox_p_value <- p.adjust(obj$data$diff_table$wilcox_p_value, method = "fdr")
test <- within(obj$data$diff_table, log2_median_ratio[wilcox_p_value > 0.05] <- 0 )
obj$data$diff_table <- test

heat_tree(obj, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = log2_median_ratio, 
          node_color_interval = c(-2, 2),
          node_color_range = c("purple", "grey95", "orange"),
          node_size_axis_label = "species count",
          node_color_axis_label = "Log 2 ratio of median proportions",
          layout = "davidson-harel",
          initial_layout = "reingold-tilford",
          title = "Orange Nodes = Higher Abundance in Canes \n Purple Nodes = Higher Abundance in Trunks",
          title_size = 0.05,
          edge_label_size_range = c(1,5),
          output_file = 'heattree_1_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_1_pval.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pval pt. 2---------------

temp <- tax_glom(gio_sp2_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)


obj_bas <- subset_taxa(temp, Phylum=="p__Basidiomycota")
tax_table(obj_bas) <- tax_table(obj_bas)[,2:7]
head(tax_table(obj_bas))

obj_asco <- subset_taxa(temp, Phylum=="p__Ascomycota")
tax_table(obj_asco) <- tax_table(obj_asco)[,2:7]
head(tax_table(obj_asco))

#parse phyloseq object giovanni
obj_all <- parse_phyloseq(temp)
obj_asco <- parse_phyloseq(obj_asco)
obj_bas <- parse_phyloseq(obj_bas)

obj <- obj_all
tissuegroup <- obj$data$sample_data$Tissue_2

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

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
                 node_color_range = c("darkblue", "royalblue", "gray95", "orange", "darkorange2"),
                 #title = "Heat Tree Matrix, pt. 2, tissue groups",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_2_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_2_pval.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pval pt. 2, basidiomycota---------------

obj <- obj_bas
tissuegroup <- obj$data$sample_data$Tissue_2

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

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
                 node_color_range = c("darkblue", "royalblue", "gray95", "orange", "darkorange2"),
                 #title = "Heat Tree Matrix, pt. 2, basidiomycota",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_2bas_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_2bas_pval.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pval pt. 2, ascomycota---------------

obj <- obj_asco
tissuegroup <- obj$data$sample_data$Tissue_2

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

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
                 node_color_range = c("darkblue", "royalblue", "gray95", "orange", "darkorange2"),
                 #title = "Heat Tree Matrix, pt. 2, ascomycota",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_2asco_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_2asco_pval.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pval pt. 3 arms---------------

temp <- tax_glom(gio_sp3a_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)

obj_all <- parse_phyloseq(temp)
obj <- obj_all
tissuegroup <- obj$data$sample_data$Description_3

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

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
                 node_color_range = c("blue", "gray95", "red"),
                 title = "Heat Tree Matrix, pt. 3, Arms",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_3a_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_3a_pval.tsv', quote=FALSE, sep='\t')


###------------------Heat tree pval pt. 3 canes---------------

temp <- tax_glom(gio_sp3c_filter1, taxrank = "Genus")

tax_table(temp) <- gsub("s__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("g__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("f__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("o__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("c__unidentified", "", tax_table(temp)); 
tax_table(temp) <- gsub("p__unidentified", "", tax_table(temp)); 
tax_table(temp)
#parse phyloseq object giovanni
obj_all <- parse_phyloseq(temp)
obj <- obj_all
tissuegroup <- obj$data$sample_data$Description_3

# Convert counts to proportions
obj$data$otu_table <- calc_obs_props(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

# Calculate per-taxon proportions 
obj$data$tax_abund <- calc_taxon_abund(obj, data = "otu_table", cols = obj$data$sample_data$sample_id)

#Calculate per-type occurence
obj$data$tax_occ <- calc_n_samples(obj, "tax_abund", groups = tissuegroup)

#Calculate difference between treatments
obj$data$diff_table <- compare_groups(obj, data = "tax_abund", cols = obj$data$sample_data$sample_id, groups = tissuegroup)

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
                 node_color_range = c("blue", "gray95", "red"),
                 layout = "davidson-harel",
                 initial_layout = "reingold-tilford",
                 title = "Heat Tree Matrix, pt. 3, Canes",
                 title_size = 0.05,
                 node_size_axis_label = "Number of species",
                 node_color_axis_label = "Log2 ratio median proportions",
                 output_file = 'heattree_3c_pval.pdf')

test <- obj$data$diff_table %>%
  mutate(taxon_names = taxon_names(obj)[taxon_id], 
         taxon_ranks = taxon_ranks(obj)[taxon_id]) %>%
  arrange(desc(abs(obj$data$diff_table$log2_median_ratio)))
write.table(test, file='heattree_3c_pval.tsv', quote=FALSE, sep='\t')


###


###-----------------SUBSETTING------------------

X_Canes <- subset_samples(giovanni_1, trunk_canes_1=="canes")
X_Trunks <- subset_samples(giovanni_1, trunk_canes_1=="trunk")

test <- subset_taxa(giovanni_1, Genus=="g__Fomitiporia")

#subset transposed phyloseq object into sample sites
ID<- as.vector(unique(sample_data(giovanni_2)$Tissue_2))
ID2 <- as.vector(c("canes_2", "a","b","c","d","e","f","g"))
ID3 <- as.vector(c(canes_2,a,b,c,d,e,f,g))
for (i in 1:length(ID)){ 
  temp <- subset_samples(giovanni_2, Tissue_2==ID[i])
  name <- paste(ID2[i])
  assign(name, temp)
}

rm(ID, ID2, ID3, i)

###----------------------Merge samples----------------


test = prune_taxa(taxa_sums(gio_sp_1) > 0, gio_sp_1)

groups_gio = c("canes", "trunk")
sample_data(test)$group_gio <- get_variable(test, "trunk_canes_1") %in% groups_gio

merged_sp_1 = merge_samples(test, "trunk_canes_1")

merged_sp_1
sample_names(test)
sample_names(merged_sp_1)

sample_data(merged_sp_1)$trunk_canes_1 = sample_names(merged_sp_1)
sample_data(merged_sp_1)$group_gio = sample_names(merged_sp_1) %in% groups_gio
merged_sp_1

#With merged unassigned

test = prune_taxa(taxa_sums(merged_taxunassigned1) > 0, merged_taxunassigned1)

groups_gio = c("canes", "trunk")
sample_data(test)$group_gio <- get_variable(test, "trunk_canes_1") %in% groups_gio

merged_sp_1u = merge_samples(test, "trunk_canes_1")

merged_sp_1u
sample_names(test)
sample_names(merged_sp_1u)

sample_data(merged_sp_1u)$trunk_canes_1 = sample_names(merged_sp_1u)
sample_data(merged_sp_1u)$group_gio = sample_names(merged_sp_1u) %in% groups_gio
merged_sp_1u



#Pt. 2 

test <- gio_sp_2
test

dim(sample_data(test))
sample_data(test)
dim(otu_table(test))
dim(tax_table(test))

test = prune_taxa(taxa_sums(gio_sp_2) > 0, gio_sp_2)

groups_gio = c("Canes", "(a)_Graft_Union", "(b)_Trunk", "(c)_Upper_trunk", "(d)_Arm_1", "(e)_Spur_1", "(f)_Arm_2", "(g)_Spur_2")
sample_data(test)$group_gio <- get_variable(test, "Tissue_2") %in% groups_gio

merged_sp = merge_samples(test, "Tissue_2")

merged_sp
sample_names(test)
sample_names(merged_sp)

sample_data(merged_sp)$Tissue_2 = sample_names(merged_sp)
sample_data(merged_sp)$group_gio = sample_names(merged_sp) %in% groups_gio
merged_sp_2 <- merged_sp

#Part 3

test = prune_taxa(taxa_sums(gio_sp_3) > 0, gio_sp_3)

groups_gio = c("Asymptomatic_cane_in_symptomatic_plant", "Symptomatic_cane", "Asymptomatic_cane_in_asymptomatic_plant", 
"Asymptomatic_arm_in_asymptomatic_plant", "Asymptomatic_arm_in_symptomatic_plant", "Symptomatic_arm")
sample_data(test)$group_gio <- get_variable(test, "Description_3") %in% groups_gio

merg_sp_3 = merge_samples(test, "Description_3")
sample_names(test)
sample_names(merg_sp_3)

sample_data(merg_sp_3)$Description_3 = sample_names(merg_sp_3)
sample_data(merg_sp_3)$group_gio = sample_names(merg_sp_3) %in% groups_gio
merged_sp_3 

merged_sp3c <- subset_samples(merged_sp_3, trunk_canes_1==1)
merged_sp3c = prune_taxa(taxa_sums(merged_sp3c) > 0, merged_sp3c); merged_sp3c
merged_sp3a <- subset_samples(merged_sp_3, trunk_canes_1==2)
merged_sp3a = prune_taxa(taxa_sums(merged_sp3a) > 0, merged_sp3a); merged_sp3a

#Arms
test = prune_taxa(taxa_sums(gio_g3_a) > 0, gio_g3_a)

groups_gio = c("Asymptomatic_arm_in_asymptomatic_plant", "Asymptomatic_arm_in_symptomatic_plant", "Symptomatic_arm")
sample_data(test)$group_gio <- get_variable(test, "Description_3") %in% groups_gio
merged_g3a = merge_samples(test, "Description_3")

sample_names(test); sample_names(merged_g3a)

sample_data(merged_g3a)$Description_3 = sample_names(merged_g3a)
sample_data(merged_g3a)$group_gio = sample_names(merged_g3a) %in% groups_gio
merged_g3a

#Canes
test = prune_taxa(taxa_sums(gio_g3_c) > 0, gio_g3_c)

groups_gio = c("Asymptomatic_cane_in_asymptomatic_plant", "Asymptomatic_cane_in_symptomatic_plant", "Symptomatic_cane")
sample_data(test)$group_gio <- get_variable(test, "Description_3") %in% groups_gio
merged_g3c = merge_samples(test, "Description_3")

sample_names(test); sample_names(merged_g3c)

sample_data(merged_g3c)$Description_3 = sample_names(merged_g3c)
sample_data(merged_g3c)$group_gio = sample_names(merged_g3c) %in% groups_gio
merged_g3c

#/

###Merge unassigned taxa

badtaxa = c("Unassigned", "k__Fungi", "o__Capnodiales", "o__Diaporthales", "o__Dothideales", 
            "o__Filobasidiales", "o__Helotiales", "o__Hypocreales", "o__Pleosporales", "o__Sporidiobolales", 
            "p__Ascomycota", "p__Basidiomycota", "c__Agaricomycetes", "c__Dothideomycetes", "c__Eurotiomycetes", 
            "c__Sordariomycetes", "c__Tremellomycetes")

merged_taxunassigned1 <- merge_taxa(gio_sp_1, badtaxa)
taxa_names(merged_taxunassigned1) <-  gsub("k__Fungi", "Unassigned", taxa_names(merged_taxunassigned1)) 

merged_sp2_un <- merge_taxa(merged_sp_2, badtaxa)
taxa_names(merged_sp2_un) <-  gsub("k__Fungi", "Unassigned", taxa_names(merged_sp2_un))

badtaxa = c("k__Fungi", "o__Capnodiales", "o__Diaporthales", "o__Dothideales", 
            "o__Filobasidiales", "o__Helotiales", "o__Hypocreales", "o__Pleosporales", 
            "p__Ascomycota", "p__Basidiomycota", "c__Agaricomycetes", "c__Dothideomycetes", "c__Eurotiomycetes", 
            "c__Sordariomycetes", "c__Tremellomycetes")

merged_sp3_un <- merge_taxa(merged_sp_3, badtaxa)
taxa_names(merged_sp3_un) <-  gsub("k__Fungi", "Unassigned", taxa_names(merged_sp3_un))

merged_sp3c_un <- subset_samples(merged_sp3_un, trunk_canes_1==1)
merged_sp3c_un = prune_taxa(taxa_sums(merged_sp3c_un) > 0, merged_sp3c_un); merged_sp3c_un
merged_sp3a_un <- subset_samples(merged_sp3_un, trunk_canes_1==2)
merged_sp3a_un = prune_taxa(taxa_sums(merged_sp3a_un) > 0, merged_sp3a_un); merged_sp3a_un


###----------------------Core biome----------------

core.taxa.standard <- core_members(merged_sp_1, detection = 0, prevalence = 95/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_mycobiome_1_sp.tsv', quote=FALSE, sep='\t')

core.taxa.standard <- core_members(merged_sp_2, detection = 0, prevalence = 95/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_mycobiome_2_sp.tsv', quote=FALSE, sep='\t')

core.taxa.standard <- core_members(merged_sp_3, detection = 0, prevalence = 95/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_mycobiome_3_sp.tsv', quote=FALSE, sep='\t')

core.taxa.standard <- core_members(merged_sp3a, detection = 0, prevalence = 95/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_mycobiome_3_sp_arms.tsv', quote=FALSE, sep='\t')

core.taxa.standard <- core_members(merged_sp3c, detection = 0, prevalence = 95/100)
core.taxa.standard
write.table(core.taxa.standard, file='core_mycobiome_3_sp_canes.tsv', quote=FALSE, sep='\t')

###-----------Uniques, pt. 1--------------------

test2 <- as.data.frame(t(otu_table(merged_sp_1u)))
test2

temp <- as.vector(rownames(subset(test2, canes==0 )))
write.table(temp, file='unique1_trunks.tsv', quote=FALSE, sep='\t')
u1_trunks <- length(temp)

temp <- as.vector(rownames(subset(test2, trunk==0 )))
write.table(temp, file='unique1_canes.tsv', quote=FALSE, sep='\t')
u1_canes <- length(temp)

###-----------Uniques, pt. 2--------------------

test2 <- as.data.frame(t(otu_table(merged_sp_2)))

colnames(test2) = c("a_Graft_Union", "b_Trunk", "c_Upper_trunk", "d_Arm_1", "e_Spur_1", "f_Arm_2", "g_Spur_2", "Canes")  
colnames(test2)

#temp <- as.vector(rownames(subset(test2, Canes==0 & a_Graft_Union==0 & b_Trunk==0 & c_Upper_trunk==0 &  
# d_Arm_1==0 & e_Spur_1==0 & f_Arm_2==0 & g_Spur_2==0)))
#write.table(temp, file='unique2_.tsv', quote=FALSE, sep='\t')
#u2_ <- length(temp)

temp <- as.vector(rownames(subset(test2, Canes==0 & a_Graft_Union==0 & b_Trunk==0 & c_Upper_trunk==0 & 
                                    d_Arm_1==0 & e_Spur_1==0 & f_Arm_2==0 )))
write.table(temp, file='unique2_g_Spur_2.tsv', quote=FALSE, sep='\t')
u2_spur2 <- length(temp)

temp <- as.vector(rownames(subset(test2, Canes==0 & a_Graft_Union==0 & b_Trunk==0 & c_Upper_trunk==0 &  
                                    d_Arm_1==0 & e_Spur_1==0 & g_Spur_2==0)))
write.table(temp, file='unique2_f_Arm_2.tsv', quote=FALSE, sep='\t')
u2_arm2 <- length(temp)

temp <- as.vector(rownames(subset(test2, Canes==0 & a_Graft_Union==0 & b_Trunk==0 & c_Upper_trunk==0 & 
                                    d_Arm_1==0  & f_Arm_2==0 & g_Spur_2==0)))
write.table(temp, file='unique2_e_Spur_1.tsv', quote=FALSE, sep='\t')
u2_Spur1 <- length(temp)

temp <- as.vector(rownames(subset(test2, Canes==0 & a_Graft_Union==0 & b_Trunk==0 & c_Upper_trunk==0 &  
                                     e_Spur_1==0 & f_Arm_2==0 & g_Spur_2==0)))
write.table(temp, file='unique2_d_Arm_1.tsv', quote=FALSE, sep='\t')
u2_Arm1 <- length(temp)

temp <- as.vector(rownames(subset(test2, Canes==0 & a_Graft_Union==0 & b_Trunk==0  & 
                                    d_Arm_1==0 & e_Spur_1==0 & f_Arm_2==0 & g_Spur_2==0)))
write.table(temp, file='unique2_c_Upper_trunk.tsv', quote=FALSE, sep='\t')
u2_UpperTrunk <- length(temp)

temp <- as.vector(rownames(subset(test2, Canes==0 & a_Graft_Union==0 & c_Upper_trunk==0 & 
                                    d_Arm_1==0 & e_Spur_1==0 & f_Arm_2==0 & g_Spur_2==0)))
write.table(temp, file='unique2_b_Trunk.tsv', quote=FALSE, sep='\t')
u2_Trunk <- length(temp)

temp <- as.vector(rownames(subset(test2, Canes==0  & b_Trunk==0 & c_Upper_trunk==0 & 
                                    d_Arm_1==0 & e_Spur_1==0 & f_Arm_2==0 & g_Spur_2==0)))
write.table(temp, file='unique2_a_Graft_Union.tsv', quote=FALSE, sep='\t')
u2_GraftUnion <- length(temp)

temp <- as.vector(rownames(subset(test2, a_Graft_Union==0 & b_Trunk==0 & c_Upper_trunk==0 &  
                                    d_Arm_1==0 & e_Spur_1==0 & f_Arm_2==0 & g_Spur_2==0)))
write.table(temp, file='unique2_Canes.tsv', quote=FALSE, sep='\t')
u2_Canes <- length(temp)

#---------between-group biomes
test2 <- as.data.frame(t(otu_table(merged_sp_2)))
colnames(test2) = c("a_Graft_Union", "b_Trunk", "c_Upper_trunk", "d_Arm_1", "e_Spur_1", "f_Arm_2", "g_Spur_2", "Canes")  
colnames(test2)
#temp <- as.vector(rownames(subset(test2,Canes==0&a_Graft_Union==0&b_Trunk==0&c_Upper_trunk==0
#&d_Arm_1==0&e_Spur_1==0&f_Arm_2==0&g_Spur_2==0)))
#write.table(temp, file='unique2_.tsv', quote=FALSE, sep='\t')
#u2_ <- length(temp)

temp <- as.vector(rownames(subset(test2,Canes==0&c_Upper_trunk==0&d_Arm_1==0&e_Spur_1==0&f_Arm_2==0&g_Spur_2==0)))
#write.table(temp, file='unique2_.tsv', quote=FALSE, sep='\t')
a_b <- length(temp)

temp <- as.vector(rownames(subset(test2,Canes==0&a_Graft_Union==0&d_Arm_1==0&e_Spur_1==0&f_Arm_2==0&g_Spur_2==0)))
#write.table(temp, file='unique2_.tsv', quote=FALSE, sep='\t')
b_c <- length(temp)

temp <- as.vector(rownames(subset(test2,Canes==0&a_Graft_Union==0&b_Trunk==0&e_Spur_1==0&f_Arm_2==0&g_Spur_2==0)))
#write.table(temp, file='unique2_.tsv', quote=FALSE, sep='\t')
c_d <- length(temp)

temp <- as.vector(rownames(subset(test2,Canes==0&a_Graft_Union==0&b_Trunk==0&c_Upper_trunk==0&f_Arm_2==0&g_Spur_2==0)))
#write.table(temp, file='unique2_.tsv', quote=FALSE, sep='\t')
d_e <- length(temp)

temp <- as.vector(rownames(subset(test2,Canes==0&a_Graft_Union==0&b_Trunk==0&c_Upper_trunk==0&d_Arm_1==0&e_Spur_1==0)))
#write.table(temp, file='unique2_.tsv', quote=FALSE, sep='\t')
f_g <- length(temp)

###-----------Uniques, pt. 3--------------------

#Arms
test2 <- as.data.frame(t(otu_table(merged_g3a)))

temp <- as.vector(rownames(subset(test2, Asymptomatic_arm_in_asymptomatic_plant==0 & 
                                    Asymptomatic_arm_in_symptomatic_plant==0 & Symptomatic_arm ==0)))
write.table(temp, file='unique3_.tsv', quote=FALSE, sep='\t')
u3_ <- length(temp)

temp <- as.vector(rownames(subset(test2, Asymptomatic_arm_in_asymptomatic_plant==0 & Asymptomatic_arm_in_symptomatic_plant==0)))
write.table(temp, file='unique3_sym_arm.tsv', quote=FALSE, sep='\t')
u3_sa <- length(temp)
temp <- as.vector(rownames(subset(test2, Asymptomatic_arm_in_symptomatic_plant==0 & Symptomatic_arm ==0)))
write.table(temp, file='unique3_Asym_arm_asym_plant.tsv', quote=FALSE, sep='\t')
u3_aaap <- length(temp)
temp <- as.vector(rownames(subset(test2, Asymptomatic_arm_in_asymptomatic_plant==0& Symptomatic_arm ==0)))
write.table(temp, file='unique3_Asym_arm_sym_plant.tsv', quote=FALSE, sep='\t')
u3_aasp <- length(temp)

#Canes
test2 <- as.data.frame(t(otu_table(merged_g3c)))

temp <- as.vector(rownames(subset(test2, Asymptomatic_cane_in_asymptomatic_plant==0 & 
                                    Asymptomatic_cane_in_symptomatic_plant==0 & Symptomatic_cane ==0)))
write.table(temp, file='unique3_.tsv', quote=FALSE, sep='\t')

temp <- as.vector(rownames(subset(test2, Asymptomatic_cane_in_asymptomatic_plant==0 & Asymptomatic_cane_in_symptomatic_plant==0)))
write.table(temp, file='unique3_sym_cane.tsv', quote=FALSE, sep='\t')
u3_sc <- length(temp)
temp <- as.vector(rownames(subset(test2, Asymptomatic_cane_in_symptomatic_plant==0 & Symptomatic_cane ==0)))
write.table(temp, file='unique3_Asym_cane_asym_plant.tsv', quote=FALSE, sep='\t')
u3_acap <- length(temp)
temp <- as.vector(rownames(subset(test2, Asymptomatic_cane_in_asymptomatic_plant==0& Symptomatic_cane ==0)))
write.table(temp, file='unique3_Asym_cane_sym_plant.tsv', quote=FALSE, sep='\t')
u3_acsp <- length(temp)

###----------------------Test------------------------

plot_bar(merged, "Class", fill = "Genus", facet_grid = ~Tissue_2)

allTaxa = taxa_names(merged)
minuscore <- allTaxa[!(allTaxa %in% core.taxa.standard)]
Minuscore = prune_taxa(minuscore, merged)

physeq <- normalise_data(merged, norm.method = "relative")
p <- plot_taxa(merged, grouping_column = "Tissue_2", method = "hellinger", number.taxa = 21, 
               filename = NULL)
print(p)



###----------------------PERMANOVA---------------------------

pseq.rel <- microbiome::transform(gio_sp2_filter1, "compositional")
otu <- microbiome::abundances(pseq.rel)
meta <- microbiome::meta(pseq.rel)

p <- microbiome::plot_landscape(pseq.rel, method = "NMDS", distance = "jaccard", col = "Tissue_2", size = 3)
print(p)

data <- as(sample_data(gio_sp3a_filter1), "data.frame")
permanova <- vegan::adonis(distance(gio_sp2_filter01, method="jaccard") ~ Tissue_2,
       data = data)

permanova <- vegan::adonis(t(otu)~ Tissue_2 + Description_3,
                    data = meta, permutations=999, method = "jaccard")
write.table(as.data.frame(permanova$aov.tab), file='permanova_3c_jaccard.tsv', quote=FALSE, sep='\t')

print(as.data.frame(permanova$aov.tab))

dist <- vegan::vegdist(t(otu))
temp <- anova(vegan::betadisper(dist, meta$Tissue_2))
write.table(temp, file='anova_2.tsv', quote=FALSE, sep='\t')

#Under construction
coef <- coefficients(permanova)
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")


###----------------------Heatmap---------------

#temp <- tax_glom(merg_sp_3, taxrank = "Genus")
temp <- merged_taxunassigned1
tmp <- prune_taxa(names(sort(taxa_sums(temp),TRUE)[1:40]), temp)
p <- plot_heatmap(tmp, method="NMDS", distance="bray", sample.label="trunk_canes_1", 
                  sample.order = "trunk_canes_1", title="40 Most Abundant Taxa") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank())
print(p)

temp <- merged_sp2_un
tmp <- prune_taxa(names(sort(taxa_sums(temp),TRUE)[1:40]), temp)
p <- plot_heatmap(tmp, method="NMDS", distance="jaccard", sample.label="Tissue_2", 
                  sample.order = "Tissue_2", title="40 Most Abundant Taxa") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank())
print(p)

temp <- gio_sp3_a
tmp <- prune_taxa(names(sort(taxa_sums(temp),TRUE)[1:40]), temp)
p <- plot_heatmap(tmp, method="NMDS", distance="jaccard", sample.label="Description_3", 
                  sample.order = "Description_3", title="40 Most Abundant Taxa") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank())
print(p)

temp <- gio_sp3_c
tmp <- prune_taxa(names(sort(taxa_sums(temp),TRUE)[1:20]), temp)
p <- plot_heatmap(tmp, method="NMDS", distance="jaccard", sample.label="Description_3", 
                  sample.order = "Description_3", title="20 Most Abundant Taxa") +
  theme(axis.title.y=element_blank(), axis.title.x =element_blank())
print(p)

###----------------------Ordination--------------------

ord <- ordinate(gio_sp3_c, "NMDS", "jaccard")

p <- plot_ordination(gio_sp3_c, ord, type="sample", color="Description_3", title="NMDS+Jaccard ordination")+
  scale_color_manual(values=set_3_rb) +
  stat_ellipse(aes(group=Description_3),type = "norm") +
  geom_point(size=2)
p
p + facet_wrap(~Description_3, 3)

print(p)

dist1 <- as.matrix(phyloseq::distance(gio_sp_1, "bray", type="sample"))
write.table(dist1, file='bray_dist1.tsv', quote=FALSE, sep='\t')
dist1 <- as.matrix(phyloseq::distance(gio_sp_2, "bray", type="sample"))
write.table(dist1, file='bray_dist2.tsv', quote=FALSE, sep='\t')
dist1 <- as.matrix(phyloseq::distance(gio_sp3_a, "bray", type="sample"))
write.table(dist1, file='bray_dist3a.tsv', quote=FALSE, sep='\t')
dist1 <- as.matrix(phyloseq::distance(gio_sp3_c, "bray", type="sample"))
write.table(dist1, file='bray_dist3c.tsv', quote=FALSE, sep='\t')



###-----------------------LCBD-------------------------

temp <- normalise_data(merged_sp_1u, norm.method = "relative")
p <- plot_taxa(temp, grouping_column="trunk_canes_1", method="ab.jaccard",number.taxa=21,filename=NULL)
print(p + ggtitle("Relative abundance, 20 most common taxa (Trunks vs. Canes)"))

temp <- normalise_data(merged_sp2_un, norm.method = "relative")
p <- plot_taxa(temp, grouping_column="Tissue_2", method="ab.jaccard",number.taxa=21,filename=NULL)
print(p + ggtitle("Relative abundance, 20 most common taxa (Tissue Groups)"))

temp <- normalise_data(merged_sp3a_un, norm.method = "relative")
p <- plot_taxa(temp, grouping_column="Description_3", method="ab.jaccard",number.taxa=21,filename=NULL)
print(p + ggtitle("Relative abundance, 20 most common taxa (Arms)"))

temp <- normalise_data(merged_sp3c_un, norm.method = "relative")
p <- plot_taxa(temp, grouping_column="Description_3", method="ab.jaccard",number.taxa=21,filename=NULL)
print(p + ggtitle("Relative abundance, 20 most common taxa (Canes)"))


###-----------------------Co-occurence Networks----------------------------

temp <- gio_sp3_c

#transpose object
OTU <- t(otu_table(temp, taxa_are_rows = TRUE))
TAX = as.matrix(tax_table(temp))
met <- as.data.frame(sample_data(temp))
temp <- phyloseq(OTU, TAX)
MBS_physeq<-merge_phyloseq(temp, met)
rm(OTU,TAX,met,temp)

tmp <- taxa_level(MBS_physeq, "Genus")

co_occr<- co_occurence_network(tmp, grouping_column = "Description_3", rhos = 0.35, method="cor", 
                                   qval_threshold=0.05, select.condition = "Symptomatic_cane", scale.vertex.size=4, scale.edge.width=10, 
                                   plotNetwork=T, plotBetweennessEeigenvalue=T)

#Tissue_2 groups: Canes, (a)_Graft_Union, (b)_Trunk, (c)_Upper_trunk, (d)_Arm_1, (e)_Spur_1, 
# (f)_Arm_2, (g)_Spur_2
#Description_3 groups: Asymptomatic_cane_in_symptomatic_plant, Symptomatic_cane, Asymptomatic_cane_in_asymptomatic_plant, 
# Asymptomatic_arm_in_asymptomatic_plant, Asymptomatic_arm_in_symptomatic_plant, Symptomatic_arm

require(visNetwork)
g <- co_occr$net$graph

data <- toVisNetworkData(g)
visNetwork(nodes = data$nodes, edges = data$edges, main="Giovanni dataset, canes. ", width = 900)%>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)



###-----------------------Don't use! Subset uniques and core-----------------

#In order to determine which are present in canes but not trunks, we have to determine presence in each. 

unique_trunks = prune_taxa(taxa_sums(Canes)==0, Canes)
unique_canes = prune_taxa(taxa_sums(Trunks)==0, Trunks)

all_canes = prune_taxa(taxa_sums(Canes)!=0, Canes)
all_trunks = prune_taxa(taxa_sums(Trunks)!=0, Trunks)

rownames_canes <- rownames(otu_table(unique_canes))
rownames_trunks <- rownames(otu_table(unique_trunks))

unique_canes = prune_taxa(rownames_canes, giovanni)
unique_trunks = prune_taxa(rownames_trunks, giovanni)


# Define the taxa you don't want like this:
#badTaxa = c("bad1", "bad2", "bad3")
allTaxa = taxa_names(giovanni)
myTaxa <- allTaxa[!(allTaxa %in% rownames_canes)]
tmp = prune_taxa(myTaxa, giovanni)

allTaxa = taxa_names(tmp)
myTaxa <- allTaxa[!(allTaxa %in% rownames_trunks)]
core_biome = prune_taxa(myTaxa,tmp)
rm(tmp, allTaxa, myTaxa)

write.table(tax_table(core_biome), file='core_mycobiome_1.tsv', quote=FALSE, sep='\t')
write.table(tax_table(unique_trunks), file='unique_trunks_1.tsv', quote=FALSE, sep='\t')
write.table(tax_table(unique_canes), file='unique_canes_1.tsv', quote=FALSE, sep='\t')


###-----------------Pie chart and venn diagram-----------

## data input (number of reads mapped to each category)
core <- as.numeric(nrow(tax_table(core_biome)))
canes_uni <- as.numeric(nrow(tax_table(unique_canes)))
trunks_uni <- as.numeric(nrow(tax_table(unique_trunks)))
trunks_pie <- as.numeric(nrow(tax_table(all_trunks)))
canes_pie <- as.numeric(nrow(tax_table(all_canes)))
total <- as.numeric(nrow(tax_table(giovanni)))

core+canes_uni+trunks_uni
total

piedf <- data.frame(
  group = c("Trunks", "Canes", "Present in both"),
  value = c(trunks_uni,canes_uni,core)
)
head(piedf)



bp<- ggplot(piedf, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp

pie <- bp + coord_polar("y", start=1000)
pie

pie_values <- c(trunks_uni,canes_uni,core)

pie(
  pie_values,
  c(trunks_uni,canes_uni,core), 
  col = c("purple","orange","maroon"), 
  main = "Species distribution", 
  border=NA, 
  init.angle=45)
legend("topright", 
       c("Unique to Trunks","Unique to Canes","Present in Both"), 
       cex = 1, 
       fill = c("purple", "orange", "maroon"))

###--------------VennPieR--------

## data input (number of reads mapped to each category)
total
trunks_pie
canes_pie
trunks_uni
canes_uni
core
both=core
rest=total

# parameter for pie chart
iniR=0.2 # initial radius
colors=list(NO='white',
            #total='black', 
            #rest='#e5f5e0', 
            Present_in_Both='maroon',
            trunks_uni='violet',
            canes_uni='gold',
            Unique_to_Trunks='purple',
            Unique_to_Canes='orange'
                    )
colors_legend=list(NO='white',
            Present_in_Both='maroon',
            Unique_to_Trunks='violet',
            Unique_to_Canes='orange'
            
)

library('plotrix')

# from outer circle to inner circle
#0 circle: blank
pie(1, radius=iniR, init.angle=90, col=c('white'), border = NA, labels='', main ="Species Distribution")

#2 circle: for transparent colours, uniques only
floating.pie(0,0,c(both, canes_uni, trunks_uni),radius=3*iniR, 
             startpos=pi/4, col=as.character(colors[c('NO','canes_uni', 'trunks_uni')]),border=NA)

#1 circle: for opaque colors and both
floating.pie(0,0, c(both, canes_uni, trunks_uni), radius=2*iniR, 
             startpos=pi/4, col=as.character(colors[c('Present_in_Both','Unique_to_Canes','Unique_to_Trunks')]), border = NA) 

pieangles <- floating.pie(0,0, c(both, canes_uni, trunks_uni), 
                          radius=2*iniR, 
                          startpos=pi/4, 
                          col=as.character(colors[c('Present_in_Both','Unique_to_Canes','Unique_to_Trunks')]), 
                          border = NA) 

pie.labels(x=0,y=0,pieangles,c(both,canes_uni,trunks_uni),radius=.5,bg="white",border=TRUE,
           minangle=NA,boxed=FALSE,explode=0, cex=.80)

legend(.5, 5*iniR, gsub("_"," ",names(colors_legend)[-1]), col=as.character(colors_legend[-1]), pch=19,bty='n', ncol=1)

## or, in one column with reads count and %
#names=gsub("_"," ",names(colors)[-1])
#values = sapply(names(colors)[-1], get)
#percent=format(100*values/total, digits=2, trim=T)
#values = format(values, big.mark=",", scientific=FALSE, trim=T)
#cl=as.character(colors[-1])
#pchs=rep(19, length(cl)); pchs[1]=1;
#legend(0, 5*iniR, paste(names," (",values,", ", percent,"%)", sep=""), col=cl, pch=pchs,bty='n', ncol=1, cex=0.6)

###-------------------------Venn diagrams---------------------

library(VennDiagram)

grid.newpage()
draw.pairwise.venn(area1 = trunks_pie, area2 = canes_pie, cross.area = both, category = c("Trunks", "Canes"), lty = rep("blank", 2), 
                   fill = c("purple", "orange"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

###-------------------------Don't use - Likert plots--------------


p_an_p <- plot_anova_diversity_pval(test.log, method = c("shannon", "evenness"),grouping_column ="Tissue_2",pValueCutoff=0.05)

p <- plot_taxa(giovanni, grouping_column="Tissue_2", method="ab.sorensen",number.taxa=20,filename=NULL)
print(p + ggtitle("Relative abundance, 20 most common taxa (Trunks vs. Canes)"))

most_abund <- as.character(unique(p$data$Taxa))

allTaxa = taxa_names(giovanni)
myTaxa <- allTaxa[(allTaxa %in% most_abund)]
most_abund = prune_taxa(myTaxa,giovanni)
rm(allTaxa, myTaxa)

most_abund

MA_Canes <- subset_samples(most_abund, trunk_canes_1=="canes")
MA_Trunks <- subset_samples(most_abund, trunk_canes_1=="trunk")

ma_canes <- rowSums(otu_table(MA_Canes))
ma_canes_neg <- ma_canes*-1
ma_trunks <- rowSums(otu_table(MA_Trunks))
ma_overall <- rowSums(otu_table(most_abund))
ma_matrix <- cbind(ma_overall,ma_canes,ma_trunks)
ma_dt <- as.data.frame(ma_matrix[order(ma_matrix[,1], decreasing = TRUE),])
ma_matrix <- as.matrix(ma_dt)
#rownames_ma <- rownames(ma_dt)

#ma_c <- as.data.frame(ma_dt$ma_canes)
#ma_t <- as.data.frame(ma_dt$ma_trunks)

rownames_ma <- rownames(ma_matrix)

ma_c <- as.data.frame(ma_matrix[,2])
ma_t <- as.data.frame(ma_matrix[,3])

ggplot() + geom_bar(data=ma_t, aes(x = rownames_ma, y=ma_t$`ma_matrix[, 3]`, fill = "purple"), position="stack", stat="identity") +
  geom_bar(data=ma_c, aes(x = rownames_ma, y=-ma_c$`ma_matrix[, 2]`, fill = "orange"), position="stack", stat="identity") +
  geom_hline(yintercept = 0, color =c("white")) +
  scale_fill_identity("", labels = c("trunks", "canes"), breaks=c("purple", "orange"), guide="legend") + 
  coord_flip() +
  labs(title="Relative Abundance, 20 most common ", y="",x="") +
  theme(plot.title = element_text(size=14, hjust=0.5)) +
  theme(axis.text.y = element_text(hjust=0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
#  theme(legend.position = "bottom") 

giovanni_1 <- giovanni

###----------------Uniques---------------------

#Subsets of tissue groups

ID<- as.vector(unique(sample_data(giovanni_2)$Tissue_2))
for (i in 1:length(ID)){ 
  temp <- subset_samples(giovanni_2, Tissue_2==ID[i])
  name <- paste(ID[i])
  assign(name, temp)
}


#Unique taxa
core_biome_2
unique_graft_union
unique_trunk_2
unique_upper_trunk
unique_arm1
unique_spur1
unique_arm2
unique_spur2

# Define the taxa you don't want like this:
#badTaxa = c("bad1", "bad2", "bad3")
allTaxa = taxa_names(giovanni)
myTaxa <- allTaxa[!(allTaxa %in% rownames_canes)]
tmp = prune_taxa(myTaxa, giovanni)

allTaxa = taxa_names(tmp)
myTaxa <- allTaxa[!(allTaxa %in% rownames_trunks)]
core_biome = prune_taxa(myTaxa,tmp)
rm(tmp, allTaxa, myTaxa)



p_an_p <- plot_anova_diversity_pval(raw_data_2, method = c("shannon", "evenness"),grouping_column ="Tissue_2",pValueCutoff=0.05)




unique_trunks = prune_taxa(taxa_sums(Canes)==0, Canes)
unique_canes = prune_taxa(taxa_sums(Trunks)==0, Trunks)

all_canes = prune_taxa(taxa_sums(Canes)!=0, Canes)
all_trunks = prune_taxa(taxa_sums(Trunks)!=0, Trunks)

rownames_canes <- rownames(otu_table(unique_canes))
rownames_trunks <- rownames(otu_table(unique_trunks))

unique_canes = prune_taxa(rownames_canes, giovanni)
unique_trunks = prune_taxa(rownames_trunks, giovanni)

# Define the taxa you don't want like this:
#badTaxa = c("bad1", "bad2", "bad3")
allTaxa = taxa_names(giovanni)
myTaxa <- allTaxa[!(allTaxa %in% rownames_canes)]
tmp = prune_taxa(myTaxa, giovanni)

allTaxa = taxa_names(tmp)
myTaxa <- allTaxa[!(allTaxa %in% rownames_trunks)]
core_biome = prune_taxa(myTaxa,tmp)
rm(tmp, allTaxa, myTaxa)

write.table(tax_table(core_biome), file='core_mycobiome_1.tsv', quote=FALSE, sep='\t')
write.table(tax_table(unique_trunks), file='unique_trunks_1.tsv', quote=FALSE, sep='\t')
write.table(tax_table(unique_canes), file='unique_canes_1.tsv', quote=FALSE, sep='\t')

