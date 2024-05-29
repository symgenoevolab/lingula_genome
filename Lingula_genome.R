########## R scripts used to produce plots for
########## Lewin, Liao and Luo 2024
########## Brachiopod genome and the evolution of BMP signalling


#### Initially created: 02/06/2023 (Thomas D. Lewin)
#### Last edited: 23/05/2024 (Thomas D. Lewin)

############################# Load general packages ############################

library(dplyr)
library(ggplot2)
require(RColorBrewer)

############### Create heatmap of BMP gene repertoires (Fig. 1c) ###############

# Load specific packages
require(ComplexHeatmap)

# Read in data
bmp_table <- read.table("BMP_table_for_R.txt", sep="\t", header = TRUE)

# Remove 'Gene' column and make it into the row name
rownames(bmp_table) <- bmp_table[, 1]
bmp_table$Gene <- NULL

# Manually set colours to use for heatmap
set_cols <-  c("gray90", "#A6CEE3", "#2078B4","#F8D238", "#FF7F00",
               "red","black", "black","black","black","black","#FB5ABE")
colors = structure(set_cols, names = as.character(0:11))

# Convert table into matrix
bmp_table <- as.matrix(bmp_table)

# plot heatmap using ComplexHeatmap package
Heatmap(bmp_table,
        rect_gp = gpar(col = "white", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        col = colors)

################ Calculate transcriptome age index (Fig. 5f,g) #################

# Install myTAI package if needed
devtools::install_github("drostlab/myTAI")

# Load specific packages
library(myTAI)

# First, read in the GenEra output
gene_ages <- read.table("gene_ages.tsv", sep="\t", header = FALSE)

# Adjust gene names, if necessary
gene_ages$V1 <- sapply(strsplit(gene_ages$V1, "_"), function(x) paste(x[1], x[2], sep = "_"))

# Load expression matrices
manipulation <- read.table("BMP_manipulation_replicates_expression_matrix.tsv", sep="\t", header = TRUE)

## For each row in embryo, get the corresponding Phylostratum

# First, create a new column to store the phylostrata
manipulation$Phylostratum <- NA

# Second, loop through the rows of the 'manipulation' dataset
for (i in 1:nrow(manipulation)) {
  # Find the index of the row in 'gene_ages' where the values match
  match_index <- match(manipulation$X[i], gene_ages$V1)
  
  # Check if a matching row was found, then add the phylostrata value to the new column
  if (!is.na(match_index)) {
    manipulation$Phylostratum[i] <- gene_ages$V3[match_index]
  }
}

# Change column names and re-order
colnames(manipulation)[colnames(manipulation) == "X"] <- "GeneID"
manipulation <- cbind(manipulation[length(manipulation)], manipulation[-length(manipulation)])

## Check NAs

# NAs represent genes where GenEra classifies them as contamination or HGT 
# and therefore does not assign them a Phylostratum.

# what percentage of rows with NAs do we have?
(sum(!complete.cases(manipulation))/nrow(manipulation))*100
# around 1%

# remove NA rows
manipulation <- na.omit(manipulation)

#Calculate TAI with myTAI
manipulationTAI <- TAI(manipulation)

# Write TAI to output csv file
write.csv(manipulationTAI, "TAI.csv")

# Plot TAI
PlotSignature( ExpressionSet = manipulation,
               measure       = "TAI", 
               TestStatistic = "FlatLineTest",
               xlab          = "Condition", 
               ylab          = "TAI" )

################## Produce BUSCO plot (Extended Data Fig. 1) ###################

# Read in data
busco <- data.frame(read.table("buscos.txt", sep="\t", header = TRUE))

# Set order of BUSCO variable
busco$BUSCO <- factor(busco$BUSCO, levels = c( "Missing", "Fragmented",
                                               "Duplicated","Single-copy"))

# Reorder by percent of single-copy orthologues
ordered_busco <- busco %>%
  arrange(desc(BUSCO == "Single-copy"), desc(Percent))

# Create a sorting order variable based on "Single-copy" percentage
order_variable <- ordered_busco$Species[ordered_busco$BUSCO == "Single-copy"]

# Manually set colours
custom_colors <- c("Single-copy" = "#A6CEE3", "Duplicated" = "#2078B4", "Fragmented" = "#FDBF6F", "Missing" = "#FF7F00")

# Create stacked bar plot
ggplot(data = ordered_busco, aes(x = reorder(Species, match(Species, order_variable)), y = Percent, fill = BUSCO)) +  geom_bar(stat='identity') +
  ylab("BUSCO genes (%; n = 954)") +
  xlab("Species") +
  scale_fill_manual(values = custom_colors) +
  theme_classic() +
  theme(text = element_text(size = 16))


######### Produce gene and repeat density plots (Extended Data Fig. 3) #########

# load specific packages
require(RIdeogram)
# vignette: https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html

## First, plot gene density 

# Load in gene gtf from Lingula with the Lingula karyotype file made using SyntenyFinder.
# SyntenyFinder Github: https://github.com/symgenoevolab/SyntenyFinder
# Window of 100 kb
lan_gene_density <- GFFex(input = "Lan.gtf", karyotype = "lan_karyotype.txt", feature = "gene", window = 100000)

# Load in karyotype file
lan_karyotype <- read.table("lan_karyotype.txt", sep="\t", header = TRUE, stringsAsFactors = F)

# Read in BMP gene coordinates to add to plot
bmp <- read.table("BMP_gene_coordinates.txt", header = TRUE, sep = "\t")

# Create plot of gene density and BMP gene locations
ideogram(karyotype = lan_karyotype, overlaid = lan_gene_density,  label = bmp, label_type = "marker",output = "gene_density.svg", width = 100)

## Next, plot repeat density

# Load in repeat gtf produced by RepeatMasker
# Window of 100 kb
lan_repeat_density <- GFFex(input = "Lan_repeats.gtf", karyotype = "lan_karyotype.txt", feature = "dispersed_repeat", window = 100000)

# Create plot of repeat density
ideogram(karyotype = lan_karyotype, overlaid = lan_repeat_density, output = "repeat_density.svg", width = 100)

### rRNA density

# Load in rRNA density gtf from barrnap
# Window of 100 kb
lan_rRNA_density <- GFFex(input = "Lan_rRNAs.gff", karyotype = "lan_karyotype.txt", feature = "rRNA", window = 100000)

#Create plot of rRNA density
ideogram(karyotype = lan_karyotype, overlaid = lan_rRNA_density, output = "rRNA_density.svg", width = 100, colorset1 = c("white", "#ba0404", "#801b1b"))

## Plot repeat (%) by chr length

# Load in data for genome stats
scaf <- data.frame(read.table("scaffold_stats.txt", sep="\t", header = TRUE))

# Set variables
y <- scaf$Masked
x <- scaf$Length

# Make plot 
plot(x, y, xlim = c(20, 75), pch = 16, cex = 1.5,
     xlab = "Chromosome Length (Mb)",
     ylab = "Repeat content (%)")

# Fit 5 linear models
lm1 <- lm(y~x, data=scaf)
lm2 <- lm(y~poly(x,2,raw=TRUE), data=scaf)
lm3 <- lm(y~poly(x,3,raw=TRUE), data=scaf)
lm4 <- lm(y~poly(x,4,raw=TRUE), data=scaf)
lm5 <- lm(y~poly(x,5,raw=TRUE), data=scaf)

# Calculate a range of 100 no.s starting at 20 and ending at 75
x_axis <- seq(20,75,length=100)

# plot lm on the 
lines(x_axis, predict(lm1, data.frame(x=x_axis)), col='red')
lines(x_axis, predict(lm2, data.frame(x=x_axis)), col='yellow')
lines(x_axis, predict(lm3, data.frame(x=x_axis)), col='purple')
lines(x_axis, predict(lm4, data.frame(x=x_axis)), col='blue')
lines(x_axis, predict(lm5, data.frame(x=x_axis)), col='orange')

# calculate r-squared values of each lm
summary(lm1)$adj.r.squared
summary(lm2)$adj.r.squared
summary(lm3)$adj.r.squared
summary(lm4)$adj.r.squared
summary(lm5)$adj.r.squared

# Get summary statistics
summary(lm2)

# ANOVA to test difference between R2 values of lm1 and lm2
anova(lm1, lm2)

########### Produce gene expression heatmaps (Extended Data Fig. 8) ############

# load packages
require(ComplexHeatmap)

# Read in log2(TPM+1) values
expression <- data.frame(read.table("Gene_expression.txt", sep="\t", header = TRUE))

# Remove 'Gene' column and make it into the row name
rownames(expression) <- expression[, 1]
expression$Gene <- NULL

# Convert to matrix
expression <-as.matrix(expression)

# Set colour matrix 
my_palette <- colorRampPalette(c("royal blue", "white", "red"))(n = 64)

# Make heatmap with ComplexHeatmap
Heatmap(expression, cluster_rows= TRUE, cluster_columns = FALSE ,col = my_palette)


################ Genome statistic plots (Supplementary Fig. 1) #################

# Read in data for genome stats
stats <- data.frame(read.table("genome_stats.txt", sep="\t", header = TRUE))

## GC content

# Arrange data by GC content 
ordered_GC<- stats[order(-stats$GC), ]

# Make dot plot
ggplot(ordered_GC, aes(x=reorder(Species, -GC), y=GC)) +
  geom_point(color="black", size = 3.5) + 
  ylab("GC content (%)") +
  xlab("Species") +
  theme_classic() +
  theme(text = element_text(size = 16))

## Genome size 

# Arrange data by Genome size 
ordered_Size <- stats[order(-stats$Size_Mb) , ]

# Make lollipop plot
ggplot(ordered_Size, aes(x=reorder(Species, -Size_Mb), y=Size_Mb)) +
  geom_point(color="black", size = 3.5) + 
  geom_segment( aes(x=Species, xend=Species, y=0, yend=Size_Mb), color="black")+
  ylab("Genome size (Mb)") +
  xlab("Species") +
  theme_classic() +
  theme(text = element_text(size = 16))

## Scaffold N50

# Arrange data by Scaffold N50
ordered_N50 <- stats[order(-stats$Scaf_N50_Mb) , ]

# Remove Sme from this one
ordered_N50 <- ordered_N50[ordered_N50$Species != "Sme", ]

# Make lollipop plot
ggplot(ordered_N50, aes(x=reorder(Species, -Scaf_N50_Mb), y=Scaf_N50_Mb)) +
  geom_point(color="black", size = 3.5) + 
  geom_segment( aes(x=Species, xend=Species, y=0, yend=Scaf_N50_Mb), color="black")+
  ylab("Scaffold N50 (Mb)") +
  xlab("Species") +
  theme_classic() +
  theme(text = element_text(size = 16))+
  ylim(c(0,1))

################################ End of script #################################
