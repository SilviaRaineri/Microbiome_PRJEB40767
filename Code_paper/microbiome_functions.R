##### Microbiome Functions ####
library(tidyverse)
library(ggpubr)
library(phyloseq)
library(vegan)
library(ggrepel)
options(scipen=999)

### Color schemes 

my_colors <- list(Treatment = c('#b3b3b3','#a8c0fc','#ff9b9b', '#fccd76', '#55d2cc', '#d43a5b','#ffa227','#359c7a', '#1b2152'),
                  Group = c('#a8c0fc','#ff9b9b', '#fccd76', '#55d2cc', '#d43a5b','#ffa227','#359c7a', '#1b2152'),
         Day = c('#02eb6b','#32468d'), # Pre / Day_42
         groups_4 =  c('#f95e59', '#a851db', '#fca311', '#3a75c4', '#83b5d1'), # high_f_gain/ high_f_loss/ low_f_gain/ low_f_loss
         weight = c('#32468d', '#ff686d','grey'), # gain/ loss/ stable
         adiv_classification =c('#159598', '#ffab8f', '#d74a86'))

names(my_colors$Treatment) <- c('Pre', 'Ctrl', 'RDCOB', 'Bupropion', 'Naltrexone','RDCOB_Bupropion','RDCOB_Naltrexone', 'Bupropion_Naltrexone', 'Sibutramine')
names(my_colors$Day) <- c('Pre', 'Day_42')
names(my_colors$weight) <- c('gain', 'loss', 'stable')
names(my_colors$Group) <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
names(my_colors$groups_4) <- c('high_food_weight_gain', 'high_food_weight_loss', 'low_food_weight_gain', 'low_food_weight_loss', 'n.a.')
names(my_colors$adiv_classification) <- c('high', 'low', 'medium')

viz <- c('#ffd700','#ffb14e','#fa8775','#ea5f94', '#cd34b5','#9d02d7','#0000ff')
### Species table  and map setup ---------------------------------------------------------

setup_species_table <- function(path) {
  species <- readxl::read_xlsx(path)
  species1 <- apply(species, 2, as.numeric)
  rownames(species1) <- species$ID
  species1 <- species1[,-1]
  colnames(species1)
  return(species1)
}


#species <- setup_species_table('./Data/species_c1rc2 copy.xlsx')

setup_map <- function(path) {
  map <- readxl::read_xlsx(path)
  rownames(map) <- map$Sample
  setdiff(rownames(map), colnames(species)) # This needs to be done because we have metadata about samples which have not been sequenced
  a <- which(rownames(map) == "50")
  b <- which(rownames(map) == "50H")
  c <- which(rownames(map) == "54")
  d <- which(rownames(map) == "48")
  map <- map[-c(a,b,c,d),]
  rownames(map) <- map$Sample
  return(map)
}

setup_tax <- function(path, rownames = species_table) {
  tax.table <- read_delim(path, delim = '\t') %>% as.matrix() # NB better to edit tax table on Excel or Spreadsheet, then import it into R
  rownames(tax.table) <- rownames(rownames)
  return(tax.table)
}


#map <- setup_map('./Data/cohort12_map.xlsx')
### Phyloseq object setup -----------------------------------------------------------------

make_phyloseq <- function(otu, tax, map) {
  OTU <- otu  %>% otu_table(otu, taxa_are_rows = T)
  TAX <- tax_table(tax) 
  rownames(TAX) <- rownames(OTU)
  DATA <- sample_data(map)
  physeq = phyloseq(OTU, TAX, DATA)
  return(physeq)
}


### PCA for Beta Diversity  ---------------------------------------------------------------
plotPCA <- function(data, color = c('Timepoint', 'Treatment', 'weight_loss')) {
  # Calculate distance matrix
  dist.matrix <- vegdist(t(data), method = 'bray')
  # Calculate PCA
  pc <- cmdscale(dist.matrix,k=2)
  myPCA<- princomp(pc); myPCA
  # Prepare dataframe for plot
  my.df <- data.frame(V1 =pc[,1], V2= pc[,2], ID = map$Sample, treatment = map$Treatment, timepoint = map$Day, weight_loss = map$weight_loss)
  # Plot
  if (color == 'Timepoint') {
    map$Day <- factor(map$Day, levels = c('Pre', 'Day_42'))
    myplot <- ggplot(my.df, aes(x = V1, y = V2, color = map$Day)) +
      geom_point(stat= 'identity', size = 2.5) +
      scale_color_manual(values = color.Day) +
      labs(x= 'PC1', y = 'PC2', shape = 'Timepoint',color = 'Timepoint')
  } else if (color == 'Treatment') {
    map$Treatment <- factor(map$Treatment, levels = c('Pre', 'Ctrl', 'RDCOB', 'Bupropion', 'Naltrexone','RDCOB_Bupropion','RDCOB_Naltrexone', 'Bupropion_Naltrexone', 'Sibutramine'))
    myplot <- ggplot(my.df, aes(x = V1, y = V2, color = map$Treatment)) +
      geom_point(stat= 'identity', size = 2.5) +
      scale_color_manual(values = color.Treatment) +
      labs(x= 'PC1', y = 'PC2',  color = 'Treatment')
} else if (color == 'weight_loss') {
  map$weight_loss <- factor(map$weight_loss, levels = c('gain', 'loss', 'stable'))
  myplot <- ggplot(my.df, aes(x = V1, y = V2, color = map$weight_loss)) +
    geom_point(stat= 'identity', size = 2.5) +
    scale_color_manual(values = color.weight.loss) +
    labs(x= 'PC1', y = 'PC2',  color = 'Weight Loss')} 
  return(list(myplot, myPCA))
    }


beta_diversity <- function(data) {
  # Calculate distance matrix
  t.data <- t(data)
  dist.matrix <- vegdist(t.data, method = 'bray')
  # Calculate PCA
  pc <- cmdscale(dist.matrix,k=2)
  myPCA<- princomp(pc)
  # Prepare dataframe for plot
  my.df <- data.frame(V1 =pc[,1], V2= pc[,2], ID = map$Sample[match(rownames(t.data), map$Sample)], treatment = map$Treatment[match(rownames(t.data), map$Sample)], timepoint = map$Day[match(rownames(t.data), map$Sample)], weight_loss = map$weight_loss[match(rownames(t.data), map$Sample)], group = map$Group[match(rownames(t.data), map$Sample)])
  return(my.df)
}

beta_diversity.variance <- function(data) {
  # Calculate distance matrix
  t.data <- t(data)
  dist.matrix <- vegdist(t.data, method = 'bray')
  # Calculate PCA
  pc <- cmdscale(dist.matrix,k=2)
  myPCA<- princomp(pc)
  return(summary(myPCA))
}



plot_Beta_div <- function(data, color = variable_color, color_scheme = my_colors) {
  col <- which(colnames(data) == color)
  a <-  ggplot(data, aes(x = V1, y = V2, color = data[, col])) +
    geom_point(stat='identity', size = 2.5, aes(shape = group)) +
    scale_color_manual(values = color_scheme) +
    labs(x= 'PC1', y = 'PC2',  title = 'Beta Diversity', color = color) +
    theme_bw()
return(a)
}


### Alpha Diversity -----------------------------------------------------------------------

alpha.diversity <- function(data, measure = c('Shannon', 'Simpson', 'InvSimpson', 'Chao1')) {
  otu_table(data) <- otu_table(round(as((otu_table(data)), "matrix")), taxa_are_rows(data))
  alpha.data <- estimate_richness(data, measures = measure)
  rownames(alpha.data) <- gsub("X", "", rownames(alpha.data))
  alpha.df <- data.frame(alpha.data, weight = map$weight_loss[match(rownames(alpha.data), map$Sample)], Group = map$Group[match(rownames(alpha.data), map$Sample)], Day = map$Day[match(rownames(alpha.data), map$Sample)],
                         Weight = map$Weight[match(rownames(alpha.data), map$Sample)], cohort = map$cohort[match(rownames(alpha.data), map$Sample)], number = map$number[match(rownames(alpha.data), map$Sample)])
  
  return(alpha.df)
}

alpha.boxplot <- function(data, variable_x, variable_y, color = variable_fill, wrap = variable_wrap, ref_group = variable_ref_group, color_scheme = my_colors) {
  if (missing(wrap)) {
    ex <- which(colnames(data) == variable_x)
    why <- which(colnames(data) == variable_y)
    fil <- which(colnames(data) == color)
    p <-  ggboxplot(data, x = colnames(data)[ex] , y = colnames(data)[why], color = colnames(data)[fil]) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=14)) +
      stat_compare_means( method = 't.test', label = 'p.format', ref.group = ref_group) +
      # geom_text_repel(aes(label = as.character(number))) +
      geom_line(aes(group = number), color = "#ddd1db") +
      geom_point(aes(color = Group)) +
      scale_color_manual(values = color_scheme, labels = c("Ctrl", "Sib"))+
      labs(title = 'Alpha Diversity ')
  } else {
  ex <- which(colnames(data) == variable_x)
  why <- which(colnames(data) == variable_y)
  w <- which(colnames(data) == wrap)
  fil <- which(colnames(data) == color)
  w2 <- noquote(colnames(data)[w])
  a <-  ggboxplot(data, x = colnames(data)[ex] , y = colnames(data)[why], color = colnames(data)[fil])  +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=14)) +
    #geom_text_repel(aes(label = as.character(number))) +
    geom_point(aes(color = Group)) +
    stat_compare_means( method = 't.test', label = 'p.format', ref.group = ref_group) +
    scale_color_manual(values = color_scheme,  labels = c("Ctrl", "Sib"))+
    labs(title = 'Alpha Diversity ')
  p <- facet(a, facet.by = w2)
  }
  return(p)
}


### Bacteroidetes/Firmicutes Ratio --------------------------------------------------------

bact_firmi_ratio <- function(physeq) {
  firmicutes <- subset_taxa(physeq, Phylum == 'Firmicutes')
  bacteroidetes <- subset_taxa(physeq, Phylum == 'Bacteroidetes')
  ratio <- colSums(otu_table(bacteroidetes))/colSums(otu_table(firmicutes))
  ratio.df <- data.frame(Samples = sample_data(physeq)$Sample, Treatment = sample_data(physeq)$Treatment, Ratio_Bacteroidetes_to_Firmicutes = ratio,
                         group = sample_data(physeq)$Group, weight_loss = sample_data(physeq)$weight_loss, day= sample_data(physeq)$Day)
  return(ratio.df)  
}


### Setup for Anova
setup.species.anova <- function(physeq, metadata) {
  # This function prepares species table and metadata table for use
  # with anova sets of functions. 
  species <- as.data.frame(otu_table(physeq,taxa_are_rows = T))
  # 2. Delete rows whose rowSum is 0
  species <- species[-which(rowSums(species) == 0),]
  # 3. log10 transformation 
  taxonomy_relab <- log10(species + 0.0001)
  return(taxonomy_relab)}

#### Gene by species composition plot ---------------------------
get.gene.by.species.top10.df <- function(gene.of.interest, main.table, species.of.interest) {# NB!!! NAs turned into zeros, make sure to double check if value 100 in "other" at the end is because of complete absence of any species in the sample or just absence of top10 species
  # Subset main table for gene of interest
  goi <- subset(main.table, gene_name==gene.of.interest) 
  # Take out k_unassigned row
  goi <- goi[-grep('k__unassigned', rownames(goi)),]
  # keep nomenclature as rownames
  nom <- goi$nomenclature
  # transform into numeric class
  goi <- apply(goi[,-c(19:22)], 2, as.numeric)
  # add rownames again
  rownames(goi) <- nom
  # transform into percentage counts
  perc <- apply(goi, 2, function(x) x/sum(x)*100)
  perc[is.na(perc)] <- 0 # NB!!! NAs turned into zeros, make sure to double check if value 100 in "other" at the end is because of complete absence of any species in the sample or just absence of top10 species
  # get rid of rows with sum = 0
  perc <- perc[-which(rowSums(perc)==0),]
  # order perc table by rowSums 
  top10 <- order(rowSums(perc), decreasing = TRUE)
  # subset main into top 10 + species of interest
  if(missing(species.of.interest)) {
    
    perc10 <- perc[top10[1:10],]
    # Calculate new colSums
    c.sum <- colSums(perc10)
    # add "other"column to get total colSum = 100
    perc10 <- rbind(perc10, other = 100-c.sum)
    return(perc10)
  } else {
    perc10 <- perc[c(top10[1:9], which(rownames(perc)==species.of.interest)),]
    # Calculate new colSums
    c.sum <- colSums(perc10)
    # add "other"column to get total colSum = 100
    perc10 <- rbind(perc10, other = 100-c.sum)
    return(perc10)}
}
