#' ---
#' title: "BINF 6210 - Geography and Evolutionary Diversification of the Corvus genus"
#' author: "Thomas Papp-Simon"
#' date: "October 28th, 2022"
#' ---

##### Assignment 2 - BINF6210 Software Tools
##### Thomas Papp-Simon (1219370)
##### 28/10/2022

##### Geography and Evolutionary Diversification of the Corvus genus.

##### Introduction
# The members of the Corvus genus are one of the most interesting creatures in the world. Their remarkable intelligence sets them apart from other bird groups. With over 46 known extant species, it is no surprise that members of this genus are found all around the world (Garcia-Porta et al., 2022). Different species of Corvus have adapted to various geographical climates as well, which makes them all the more remarkable. There are many ecological factors to consider in the context of diversification. It is also important to realize the impact that these factors can potentially have on allopatric speciation within a genus. Therefore, I thought it would be interesting to explore this topic with a diverse genus such as Corvus. In this project, the objective is to analyze sequence and geographical data for the Corvus genus in order to investigate the relationship between geography and the diversification of this genus. Do the species within the Corvus genus love within the same geographic region, or is there clear evidence of allopatric divergence?

# Set working directory for myself
#setwd("C:/Users/thoma/Documents/BINF6210/6210-2")

# Installing and loading necessary packages

# install.packages("tidyverse")
# install.packages("rentrez")
# install.packages("seqinr")
# install.packages("rgbif")
# BiocManager::install(c("Biostrings", "muscle", "msa", "DECIPHER"))
# install.packages("stringi")
# install.packages("ape")
# install.packages("RSQLite")
# install.packages("phangorn")
# install.packages("phytools")

library(tidyverse)
library(rentrez)
library(seqinr)
library(rgbif)
library(Biostrings)
library(stringi)
library(ape)
library(RSQLite)
library(muscle)
library(DECIPHER)
library(phangorn)
library(phytools)

##### 1 - Data Acquisition, Exploration, Filtering, and Quality Control.

##### Sequence data acquisition from NCBI using the rentrez package.

# Using the entrez_db_searchable() function to give me a list of search fields for the Nucleotide database from NCBI.
entrez_db_searchable(db = "nuccore")

# Using the entrez_search() function to search the Nucleotide database for DNA sequence data on my taxonomic group of interest (Ara). I need to specify the genus name (Ara[ORGN]) and the gene of interest (COI[Gene]). I have also specified the retmax parameter to 300 to give me at least 300 hits for my search.
COI_search <- entrez_search(db = "nuccore", term = "Corvus[ORGN] AND COI[Gene] AND 400:700[SLEN]")
COI_search # only 20 hits, need to alter the retmax parameter to the max number of hits 

COI_maxHits <- COI_search$count
COI_search_max <- entrez_search(db = "nuccore", term = "Corvus[ORGN] AND COI[Gene] AND 400:700[SLEN]", retmax = COI_maxHits)
COI_search_max # Now we have the right amount of IDs for our search

# Using entrez_fetch() to get DNA sequence data in FASTA format and write it to a file.
COI_fetch <- entrez_fetch(db = "nuccore", id = COI_search_max$ids, rettype = "fasta")
write(COI_fetch, "COI_fetch_Corvus.fasta", sep = "\n")

# Convert to DNAStringSet object
COI_stringSet <- readDNAStringSet("COI_fetch_Corvus.fasta")

# Putting the sequences from the FASTA file into a dataframe so that we can work with it later
dfCOI <- data.frame(COI_Title = names(COI_stringSet), COI_Sequence = paste(COI_stringSet))



##### Cleaning/filtering the sequence data from NCBI.

view(dfCOI)

# Create new column called "Species_Name" in the dfCOI dataframe. Using the word() function to extract the 2nd and 3rd word from the existing COI_Title variable to get the species name for each record.
dfCOI$Species_Name <- word(dfCOI$COI_Title, 2L, 3L) 
length(unique(dfCOI$Species_Name)) # 15 unique species
unique(dfCOI$Species_Name) 

# Checking how much sequence data I have per species in my dataset.
dfSeqperSpecies <- dfCOI %>%
  group_by(Species_Name) %>%
  dplyr::count()
dfSeqperSpecies # Looks like I only have 1 sequence for a few species within my dataset, and some other species have a higher representation than others. Later on I will randomly sample 1 sequence from the data to represent each species, and then perform the alignment and clustering based off of that.



# Filtering nucleotide data by removing NA's, unspecified nucleotides from the ends of the sequences, and any missing nucleotides (dashes). Also filtering for sequences that have less than 1% of their nucleotides as unidentified (N's) (I want to preferably keep high quality sequences).
dfCOI_Corvus_Filtered <- dfCOI %>%
  filter(!is.na(COI_Sequence)) %>%
  mutate(Nucleotides = str_remove_all(COI_Sequence, "^N+|N+$|-")) %>%
  filter(str_count(Nucleotides, "N") <= (0.01* str_count(Nucleotides)))


# Making sure that my data filtering steps worked.
#class(dfCOI_Corvus_Filtered$Nucleotides)

sum(is.na(dfCOI_Corvus_Filtered$Nucleotides)) # 0 missing values
sum(str_count(dfCOI_Corvus_Filtered$Nucleotides, "-")) # No missing nucleotides

# Summary stats of the filtered data
summary(str_count(dfCOI_Corvus_Filtered$Nucleotides))

# Looking at the unique species in our dataset
unique(dfCOI_Corvus_Filtered$Species_Name) 


# Figure 1: 
# After filtering my data, I can use a histogram to show the distribution of sequence length. I expect the majority of the sequences to be around 600-650 base pairs since the literature states that the usual gene length for COI is around 650 base pairs (Yang et al., 2019).

# Put sequence length into a variable to be used in the generation of my histogram
seq_len <- str_count(dfCOI_Corvus_Filtered$Nucleotides)

# Use ggplot and geom_histogram to generate a histogram for the distribution of sequence lengths
ggplot(dfCOI_Corvus_Filtered, mapping = aes(x = seq_len)) + 
  geom_histogram(bins = 15, colour = "black", fill = "#009e73") + # Used a colour-blind friendly colour for the bars
  scale_y_continuous(expand = c(0,0), limits = c(0, 50)) + # This is to align the bars to the x-axis of the histogram (so there's no gap)
  ggtitle("COI sequence lengths for the Corvus genus") + # Title of histogram
  labs(x = "Sequence Length (base pairs)", y = "Frequency") + # x and y axis labels
  geom_vline(xintercept=mean(seq_len), lwd=1.5, linetype=1, color="magenta") + # mean line
  theme( 
    panel.background = element_blank(), # Remove background for a cleaner look
    axis.line.y = element_line(size = 0.75), # Make axes a bit thicker to make it look nicer
    axis.line.x = element_line(size = 0.75),
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", colour = "black", size = 20, hjust = 0.5) # Formatting title of plot
  ) 

##### Figure 1: Histogram showing the distribution of COI sequence lengths for the Corvus genus (n = 124). 



##### Geographical data acquisition from GBIF using the rgbif package.

# Let's look at the species in the sequence dataset from NCBI.
unique(dfCOI$Species_Name) # So I will get 100 occurrence records for each of these species from GBIF.

# Create a variable that contains these unique species names
species_list = unique(dfCOI$Species_Name)
species_list

# Using the occ_data function from the rgbif package to acquire occurrence data from GBIF for my taxonomic group. This function is faster than the occ_search function as it only gives me the occurrence data (which is what I'm interested in), whereas the occ_search function will also return data on taxon hierarchies and media (source: https://docs.ropensci.org/rgbif/reference/occ_data.html#occ-data-vs-occ-search).

# I'm going to acquire the same number of records (100) for occurrence data for each unique species in the species_list variable. This will ensure that I have the geographical data for every species in my phylogeny, and that there will be equal representation for each species. I also want to make sure that the data contains coordinate information by setting the 'hasCoordinate' argument to TRUE. 
GBIF_Corvus_occ <- occ_data(scientificName = species_list, hasCoordinate = TRUE, limit = 10)

# The results of the occ_data function gave me a list of 2 things for each species called 'meta' and 'data'. I only need the 'data' as it contains the information in the form of a dataframe that I need for my analysis. I need to combine all these dataframes into one using the bind_rows function from the dplyr package. I needed to specify that I'm only extracting the 'data' object for each species using the lapply function and the "[[" is used as an identifier for that dataframe. 
df_GBIF_Corvus <- bind_rows(lapply(GBIF_Corvus_occ, "[[", "data"))

# Checking the unique species for this data
unique(df_GBIF_Corvus$species) 

# I noticed that I did not get any data for one of the species: Corvus albicollis.
# It turns out that the 'Corvus albicollis' is spelt differently from NCBI  than in the GBIF database, which is why I didn't get any occurrence data for this particular species. So I need to get occurrence data for this species separately and then merge it to the rest of the data.
albicollis <- occ_data(scientificName = "Corvus albicollis", hasCoordinate = TRUE, limit = 10)
albicollis_data <- albicollis$data
df_GBIF_Corvus_merged <- bind_rows(df_GBIF_Corvus, albicollis_data)

# Checking if I have successfully merged the data frames
unique(df_GBIF_Corvus_merged$species) # Yes it worked, all the species are present!


# Write the data to a tsv file and then read it in so I don't have to re-download it every time. (I commented it out so it doesn't run every time)
#write_tsv(df_GBIF_Corvus_merged, "CorvusGBIFOccData.tsv") 
dfGBIF <- read_tsv(file = "CorvusGBIFOccData.tsv")



##### Cleaning/filtering occurrence data from GBIF.

# Looking at the dataframe to see what variables I would need for my analysis
view(dfGBIF)
names(dfGBIF)

# Create a subset of that dataframe with only the variables that we need for the analysis
dfGBIF_subset <- dfGBIF %>%
  select(genus, species, gbifID, country, countryCode, decimalLatitude, decimalLongitude)
view(dfGBIF_subset) 

# Compare the species names in the data from GBIF and NCBI.
unique(dfGBIF_subset$species) 
length(unique(dfGBIF_subset$species)) # 15 species in GBIF data

unique(dfCOI$Species_Name) 
length(unique(dfCOI$Species_Name))# 15 unique species in data from NCBI.


# Filter out the NA's from our variables.
dfGBIF_subset <- dfGBIF_subset %>%
  filter(!is.na(species)) %>%
  filter(!is.na(gbifID)) %>%
  filter(!is.na(country)) %>%
  filter(!is.na(decimalLatitude)) %>%
  filter(!is.na(decimalLongitude))

# What are the unique countries in the GBIF dataset?
unique(dfGBIF_subset$country)
# How many unique countries?
length(unique(dfGBIF_subset$country)) # 28 unique countries

# Counting how many times each species is recorded in each country
dfSpecies_by_Country <- dfGBIF_subset %>%
  group_by(country,species) %>%
  dplyr::count() # I had to specify to use the count function from the dplyr package because it wouldn't work otherwise for some reason.
view(dfSpecies_by_Country)

# There are some interesting things I noticed just by looking at the resulting dataframe. For instance, all records for Corvus coronoides, kubaryi, albicollis are from Australia, Northern Mariana Islands, and South Africa respectively. Perhaps this is indicative of geographic isolation for these particular species? 



##### 2 - Main Analysis

##### Sequence Alignment and Clustering

# Making sure my data is in the right format (dataframe) for functions from the Bioconductor package that are going to be used for the sequence alignment.
class(dfCOI_Corvus_Filtered) # already a data frame

# In order to be able to identify each unique sequence, I need to give each sequence a "name". In this case I will extract the accession number for each record and assign it to each sequence.
dfCOI_Corvus_Filtered$Accession_number <- word(dfCOI_Corvus_Filtered$COI_Title, 1L)

# Decided to remove this sequence because it was giving me problems in my alignment.
dfCOI_Corvus_Filtered2 <- dfCOI_Corvus_Filtered %>%
  filter(!Accession_number == "AY030179.1")

# Set the seed before I sample so we get the same result every time.
set.seed(900)

# Creating a subset of my sequence data to randomly select 1 sequence per species. These will be the sequences that are going to be used in my alignment and in clustering/the generation of my phylogenetic tree.
dfCOI_Species_Subset <- dfCOI_Corvus_Filtered2 %>%
  group_by(Species_Name) %>%
  sample_n(1)

# Making sure my data is in the right format for the subsequent analysis.
dfCOI_Species_Subset <- as.data.frame(dfCOI_Species_Subset)

# Convert the Nucleotides variable into a DNAStringSet object
dfCOI_Species_Subset$Nucleotides <- DNAStringSet(dfCOI_Species_Subset$Nucleotides)

# I created a column with the shortened version of each species name that I will use for the visualization of the phylogeny later on (formatted "C. speciesname").
dfCOI_Species_Subset$Short_SpeciesName <- word(dfCOI_Species_Subset$Species_Name, 2L)
dfCOI_Species_Subset$Short_SpeciesName <- paste("C.", dfCOI_Species_Subset$Short_SpeciesName, sep = " ")
# I realized that the species name of C. albicollis was still misspelled, which is why the labelling on my geophylogeny was not working. This ended up fixing the issue.
dfCOI_Species_Subset <- dfCOI_Species_Subset %>%
  mutate(Short_SpeciesName = str_replace(Short_SpeciesName,"C. albicolis", "C. albicollis"))

# Assigning the respective species name to each sequence using the names() function.
names(dfCOI_Species_Subset$Nucleotides) <- dfCOI_Species_Subset$Short_SpeciesName

# Using the muscle algorithm from the muscle package to perform the sequence alignment. I have set the gapopen argument to -1000, which means that it will penalize gaps within sequences in the context of the alignment score. As a result, we shouldn't see many gaps in the middle of my sequences after alignment. I also set the use.names argument to TRUE to preserve the names I assigned to each sequence beforehand (The species name).
dfCOI_Corvus_alignment <- DNAStringSet(muscle::muscle(dfCOI_Species_Subset$Nucleotides, gapopen = -1000), use.names = TRUE)

# Looking at the sequences after alignment in my browser.
BrowseSeqs(dfCOI_Corvus_alignment)

# Getting the mean of the number of gaps within the sequences after aligning them.
dfCOI_Corvus_alignment %>%
  lapply(str_count, "-") %>%
  unlist %>%
  mean # ~71 gaps per sequence on average.

# Write the results of the alignment to a fasta file in case I need to look at it in other software like MEGA.
#writeXStringSet(dfCOI_Corvus_alignment, file = "COI_Corvus_alignment.fasta")



# Convert the results of my alignment to a DNAbin object.
dnaBin_COI_Corvus <- as.DNAbin(dfCOI_Corvus_alignment)
class(dnaBin_COI_Corvus)



# Creating a distance matrix using the dist.dna function.
distanceMatrix <- dist.dna(dnaBin_COI_Corvus, model = "K80", as.matrix = TRUE, pairwise.deletion = TRUE)

# Convert to distance matrix specifically
distanceMatrix <- as.dist(distanceMatrix) 

# Generating clusters using the hclust function. I messed around with different methods and found that I got very similar results for every method. 
clusters.COI <- hclust(distanceMatrix,method="average")


# Transform the clusters (hclust object) into a phylo object that can be user to generate my geophylogeny later on.
phylo_tree <- as.phylo(clusters.COI)


##### Creating a geophylogeny for the Corvus genus
# Create a subset dataframe with the variables required for the geophylogeny. I need coordinate data and the species name to assign to each GPS coordinate.
GBIF_coord <- dfGBIF_subset %>%
  select(species, decimalLatitude, decimalLongitude)

# Formatting the species name so that it matches the labels on the tips of the phylogenetic tree I generated beforehand.
GBIF_coord$short_species <- word(GBIF_coord$species, 2L)
GBIF_coord$short_species <- paste("C.", GBIF_coord$short_species, sep = " ")

# Create another subset with the formatted/altered species names.
GBIF_coord2 <- GBIF_coord %>%
  select(short_species ,decimalLatitude, decimalLongitude)

# Make sure it's a data frame.
dfGBIF_coord <- as.data.frame(GBIF_coord2)
class(dfGBIF_coord)

# I need to convert the dataframe with the coordinate data into a matrix.
GBIF_matrix <- as.matrix(GBIF_coord2)

# Assign names to each row of the matrix (species name) and then removing the first row so that only the GPS coordinates remain in my matrix.
rownames(GBIF_matrix) <- dfGBIF_coord[,1]
GBIF_matrix <- GBIF_matrix[,-c(1)]

# Figure 2: 
# Generating the geophylogeny using the phylo.to.map function from phytools. I also wanted to make the lines for each species different by adding a rainbow-colour palette.
phylo_colors <- setNames(sample(rainbow(n=Ntip(phylo_tree))), phylo_tree$tip.label)
phylo.to.map(phylo_tree, coords = GBIF_matrix, colors = phylo_colors, ftype = "i", fsize = 0.8, psize = 0.7, lty = 1)

##### Figure 2: Geophylogeny of the Corvus genus. 10 data points obtained from GBIF were used for each species (n = 150).


##### Results & Discussion
# As mentioned in the introduction, there are at least 46 extant members of the Corvus genus. This genus also exhibits a widespread distribution across the world. As a result, I wanted to investigate the link between geography and diversification for this genus. For this analysis, I chose cytochrome c oxidase subunit I (COI) as the marker of choice for analyzing sequence data for the Corvus genus. This mitochondrial gene is commonly used in phylogenetic analyses, and has also been used in the study of phylogenetic relationships for members of this genus (Mansha et al., 2021). During my analysis, I was only able to acquire COI sequence data for only 15 of the 46 known species of this genus from NCBI. Therefore, the results of this analysis do not encompass the entirety of the Corvus genus. However, I was still able to achieve promising results by the end of my analysis. I used the Nucleotide database from NCBI to get the COI sequence data for this genus. During my search, I chose to filter the sequence length parameter so that it would only return sequences of 400 to 700 base pairs. COI is usually near 650 base pairs in length, so I thought that filtering for sequence length will  give me more accurate hits. After filtering the sequence data, I generated a histogram to showcase the distribution of COI sequence lengths for this genus (Figure 1). After this, I was able to successfuly align the sequences, and then use the results of my alignment to generate clusters for the members of this genus. Different models were tested during the generation of these clusters, and I found that they all gave similar results. 

# The other half of my analysis consisted of gathering occurrence data for the Corvus genus from GBIF. I was able to successfully gather this geographical data for all of the species in my sequence dataset. For this analysis, I wanted to make sure that I had an equal amount of geographical data for each species. I chose to use 10 data points for each species. At first I tried using more data points, but the results were unclear. With this data, I was able to generate a geophylogeny for the Corvus genus (Figure 2). This geophylogeny allowed me to visualize the widespread distribution and diversity of this genus. Just by looking at the geophylogeny, it is clear that there are distinct clusters that exist for certain species of this genus. For instance, Corvus albicollis is clearly clustered within the region of South Africa. However, the distribution of some other species seems to be more dispersed. This is the case for Corvus dauuricus, where its geographical distribution seems to range from India all the way to Russia and even Japan! I also noticed that there are no data points within South America. I thought that it was very interesting to see how species within similar clusters are so geographically dispersed from each other. Garcia-Porta et al. (2022) studied the different factors that lead to the global radiation of crows and ravens, and found that numerous ecological factors need to be taken into account when analyzing the divergence of this genus. For instance, the divergence of this genus could have been a result of their ability to survive new environments, even if the conditions are sub-optimal for their survival. In the context of my analysis, I would agree that there is evidence of allopatric divergence for this genus. Of course, one of the limitations of my analysis was the inability to acquire data for every species of the genus, which would have given me the opportunity to accurately analyze the genus in its entirety. If I were to research this topic in the future, I would make sure to have data for all species of the genus. I would also dive deeper to research the link between ecological factors and the divergence of this genus.


##### Acknowledgements
# I would like to thank Jesse Wolf and Nishita Sharif for their helpful comments and feedback for this project.


##### References

# as.phylo: Conversion between "phylo" and "hclust" trees: https://www.rdocumentation.org/packages/ape/versions/1.2-7/topics/as.phylo

# Combine two data frames by rows (rbind) when they have different sets of columns: https://stackoverflow.com/questions/3402371/combine-two-data-frames-by-rows-rbind-when-they-have-different-sets-of-columns

# Combining occurrence data from rgbif::occ_data/occ_search: https://discuss.ropensci.org/t/combining-occurrence-data-from-rgbif-occ-data-occ-search/1667/1

# Filter rows which contain a certain string: https://stackoverflow.com/questions/22850026/filter-rows-which-contain-a-certain-string

# Garcia-Porta, J., Sol, D., Pennell, M., Sayol, F., Kaliontzopoulou, A., & Botero, C. A. (2022). Niche expansion and adaptive divergence in the global radiation of crows and ravens. Nature Communications, 13(1), Article 1. https://doi.org/10.1038/s41467-022-29707-5

# ggplot2 histogram plot : Quick start guide - R software and data visualization: http://www.sthda.com/english/wiki/ggplot2-histogram-plot-quick-start-guide-r-software-and-data-visualization

# How to remove rows and columns of a matrix in R: https://www.educative.io/answers/how-to-remove-rows-and-columns-of-a-matrix-in-r

# Mansha, M., Khan, M. A., & Hussain, T. (2021). Molecular Identification and Phylogenetic relationship of Corvus splendens (Common House Crow) by Cytochrome c oxidase subunit I gene [Preprint]. In Review. https://doi.org/10.21203/rs.3.rs-237824/v1

# Projecting a phylogenetic tree onto a map with multiple geographic points per taxon: http://blog.phytools.org/2019/03/projecting-phylogenetic-tree-onto-map.html

# Some thoughts on working with rgbif occurrence data, including mapping: https://discuss.ropensci.org/t/some-thoughts-on-working-with-rgbif-occurrence-data-including-mapping/1105

# Yang, C.-H., Wu, K.-C., Chuang, L.-Y., & Chang, H.-W. (2019). Decision Theory-Based COI-SNP Tagging Approach for 126 Scombriformes Species Tagging. Frontiers in Genetics, 10. https://www.frontiersin.org/articles/10.3389/fgene.2019.00259
