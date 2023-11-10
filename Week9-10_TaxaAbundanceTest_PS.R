# Paisley Samuel
# November 2023
# ____________________________________________ #


### Testing Abundances of taxa between groups using phyloseq object
#### You start with a phyloseq object but end with a dataframe

# Packages used 
library(microbiome)
library(microbiomeutilities)
library(vegan)
library(phyloseq)
library(tidyverse)
library(pgirmess)
library(multcompView)

#Setting Working directory and Seed
setwd("insert path to folder containing the data files")
set.seed(1998) # choose a number to ensure you get the same results each time you run the code

## Creating abundance table by loading in feature and metadata tables
dat<-read.csv("feature_Y123_ADJUSTED.csv", header=TRUE, row.names = 1)
dat <- t(dat)  # rows should be samples so we have to transpose the dataframe (df); transposing converts dataframe to matrix
row.names(dat) # row names should now be the sample names after transposing
metadata <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1)
dat <- as.data.frame(dat)  # need to convert back into a df to complete the rest of the steps

#this chunk of code matches the samples in the feature table to the samples in
#the metadata (since there were some samples that were filtered out after data cleaning)
common.rownames <- intersect(rownames(dat), rownames(metadata)) #finds the common row names (sample names in this case) between the two dfs
dat <- dat[common.rownames,]              #subsets the data to only include the samples from the common.rownames list
metadata <- metadata[common.rownames,]
all.equal(rownames(dat),rownames(metadata)) #run this code to ensure that both dfs have the same samples (should read TRUE if all samples match)

#finding the singletons and doubletons (ASVs that appear only 1 or 2 times)
otu.abund<-which(colSums(dat)>2)

library(vegan)
dat.dom<-dat[,otu.abund] #subsets the data to only include the ASVs that were NOT singletons or doubletons
dat.pa<-decostand(dat.dom, method ="pa")   #converting dominant taxa to presence/absence abundance
dat.otus.01per<-which(colSums(dat.pa) > (0.01*nrow(dat.pa)))  #finding the taxa that appear at least 1% of the time
dat.01per<-dat.dom[,dat.otus.01per]  #subsetting the dominant taxa abundance table to include only the taxa that appear at least 1% of the time


## Making phyloseq objects (using dat.01per abundance)
library(phyloseq)
library(microbiome)
library(microbiomeutilities)
asvdat <- as.data.frame(t(dat.01per)) #taxa has to be rows now so the df was transposed and converted back to a df
taxdat <- read.csv("taxonomy_Y123_edited&cleaned.csv", header = TRUE, row.names = 1) #loading in taxonomy table
meta <- read.csv("Metadata-Diversity_BATCH.csv", header = TRUE, row.names = 1) #loading in metadata (or you can just use the metadata variable that is loaded)
asvmat <- data.matrix(asvdat)
taxmat <- as.matrix(taxdat) # use as.matrix NOT as.data.matrix as the latter will convert the data into numbers
ASV <- otu_table(asvmat, taxa_are_rows = TRUE)
TAX <- tax_table(taxmat)
META <- sample_data(meta)
pseq <- phyloseq(ASV,TAX,META) #creating phyloseq object that contains the metadata, taxonomy, and abundances


#Agglomerating the data at the genus level (or whatever level you'll be testing your abundances)
ps_genus <- tax_glom(pseq_rel, taxrank= "Genus")
#Subsetting taxa to the specific taxa that you'll be testing for differences between groups
## EXAMPLE: subset out Microcystis and Bdellovibrio
ps_genusMB <- subset_taxa(ps_genus, Genus %in% c("Microcystis_PCC-7914", "Bdellovibrio"))
MB_df <- psmelt(ps_genusMB) #converting into a df in order to create a single "Abundance" column
names(MB_df) #returns the names of the columns
#making Year a factor so its read as factor and not as numeric
MB_df$Year <- as_factor(MB_df$Year)

#Plotting just to see the differences visually (OPTIONAL)
MB_plot <- ggplot(MB_df, aes(x=Year, y=Abundance, fill=Genus))+
  geom_bar(position="fill", stat="identity")+
  guides(fill=guide_legend(title="Genus"))

##EXAMPLE 1
#Testing to see if there is sig. dif. between years in terms of Microcystis abundance
#Created a separate df to contain only the microcystis rows
M_sub <- MB_df[MB_df$Genus == "Microcystis_PCC-7914", ]

#testing each level of the predictor variable for normality
shapiro.test(subset(M_sub, Year=="1")$Abundance)
shapiro.test(subset(M_sub, Year=="2")$Abundance)
shapiro.test(subset(M_sub, Year=="3")$Abundance)
bartlett.test(M_sub$Abundance~M_sub$Year)

#NOT normally distributed and variances are not homogeneous, so we use a non-parametric test (Kruskal-Wallis)
kruskal.test(M_sub$Abundance~M_sub$Year) #p = 3.002e-07 so there is a difference between means 

#conducting a post-hoc test to see where the differences lie
kmc <- kruskalmc(M_sub$Abundance~M_sub$Year) # multiple-comparison test
kmc # comparisons TRUE= significantly different or FALSE= not significantly different
# To look for homogeneous groups, and give each group a code (letter):
test <- kmc$dif.com$stat.signif # select logical vector
names(test) <- row.names(kmc$dif.com)# add comparison names
# create a list with "homogeneous groups" coded by letter
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let # significant letters for the multiple comparison test
#   1   2   3 
# "a" "b" "c"   -> if all letters are different = each year is significantly different from each other

#calculating and comparing the means between each year (abundance)
year_abmean <- by(M_sub$Abundance, M_sub$Year, mean)
year_abmean
#Year 1 = 0.01460767
#Year 2 = 0.02511472
#Year 3 = 0.03536371

#Looking into the raw Abundance totals by year for Microcystis
M_sub %>%
  group_by(Year) %>%
  summarise(total = sum(Abundance))


##Example 2
#Testing to see if there is sig. dif. between years in terms of Bdellovibrio abundance
B_sub <- MB_df[MB_df$Genus == "Bdellovibrio", ]

#testing each level of the predictor variable for normality
shapiro.test(subset(B_sub, Year=="1")$Abundance)
shapiro.test(subset(B_sub, Year=="2")$Abundance)
shapiro.test(subset(B_sub, Year=="3")$Abundance)
bartlett.test(B_sub$Abundance~B_sub$Year)

#NOT normally distributed and variances are barely homogeneous, but we still use a non-parametric test (kruskal-wallis)
kruskal.test(B_sub$Abundance~B_sub$Year) #p = 0.00356 so there is a difference between means 

#conducting a post-hoc test
kmc <- kruskalmc(B_sub$Abundance~B_sub$Year)
kmc
test <- kmc$dif.com$stat.signif
names(test) <- row.names(kmc$dif.com)
let <- multcompLetters(test, compare="<", threshold=0.05,
                       Letters=c(letters, LETTERS, "."),reversed = FALSE)
let
#   1   2   3 
# "ab" "a" "b"   year 2 and 3 are significantly different from each other

#calculating and comparing the totals between each year (abundance)
year_bdmean <- by(B_sub$Abundance, B_sub$Year, mean)
year_bdmean
#Year 1 = 0.001026737
#Year 2 = 0.000898819
#Year 3 = 0.001127009

#Looking into the raw Abundance totals by year for Bdellovibrio
B_sub %>%
  group_by(Year) %>%
  summarise(total = sum(Abundance))
