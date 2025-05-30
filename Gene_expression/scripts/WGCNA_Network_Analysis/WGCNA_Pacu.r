
############# GWENA (Gene Whole co-Expression Network Analysis) #############

## load libraries
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("GWENA")
library("WGCNA")              #BiocManager::install("WGCNA", force = TRUE)
library("flashClust")         #install.packages("flashClust")
library("pheatmap")  
library("clusterProfiler")    #BiocManager::install("clusterProfiler")
library("simplifyEnrichment") #BiocManager::install("simplifyEnrichment")
library("genefilter")           #BiocManager::install("genefilter") 
library("DESeq2")               #BiocManager::install("DESeq2")
library("factoextra")           #install.packages("factoextra")
library("NbClust")              #install.packages("NbClust")
library("ComplexHeatmap")       #BiocManager::install("ComplexHeatmap", force = TRUE)
library("tidyverse")            
library("RColorBrewer")
library("ggplot2")              
library("goseq")                #BiocManager::install("goseq")
library("gridExtra")            #install.packages("gridExtra")
#library("VennDiagram")          #install.packages("VennDiagram")
library("patchwork")            #install.packages("patchwork")
library("dplyr")
library(GWENA)
library(magrittr) # Not mandatory, we use the pipe `%>%` to ease readability.


#treatment information
treatmentinfo <- read.csv("RNAseq_Pacu_data.csv", header = TRUE, sep = ";")

#gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix_noIso_Pacu.csv", row.names="gene_id"), colClasses = double)

#remove samples below 5 million reads
gcount <- gcount %>% dplyr::select(-B1, -F10)
treatmentinfo <- treatmentinfo[!(treatmentinfo$sample_id %in% c("B1", "F10")), ]

# Reorder treatmentinfo to match the columns of gcount_filt
treatmentinfo <- treatmentinfo[match(colnames(gcount), treatmentinfo$sample_id), ]

## Quality-filter gene counts
#Set filter values for PoverA: smallest sample size per treat is 2, so 2/21 (21 samples) is 0.1
#This means that 2 out of 21 (0.1) samples need to have counts over 10.
#So P=10 percent of the samples have counts over A=10. 
filt <- filterfun(pOverA(0.10,10))

#create filter for the counts data
gfilt <- genefilter(gcount, filt)

#identify genes to keep by count filter
gkeep <- gcount[gfilt,]

#identify gene lists
gn.keep <- rownames(gkeep)

#gene count data filtered in PoverA, P percent of the samples have counts over A
gcount_filt <- as.data.frame(gcount[which(rownames(gcount) %in% gn.keep),])

#How many rows do we have before and after filtering?
nrow(gcount) #Before                
# [1] 33730

nrow(gcount_filt) #After
# [1] 25423

### Read normalization
# Normalize our read counts using VST-normalization in DESeq2
# Construct the DESeq2 dataset

treatmentinfo$temp <- factor(treatmentinfo$temp)

#Create a DESeqDataSet design from gene count matrix and labels. Here we set the design to look at 
#any differences in gene expression across samples attributed to temp.

#Set DESeq2 design
gdds <- DESeqDataSetFromMatrix(countData = gcount_filt,
                               colData = treatmentinfo,
                               design = ~temp)

# Log-transform the count data using a variance stabilizing transformation (vst). 

SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst  to transform our data. Size factors should be less than for to use vst
print(sizeFactors(SF.gdds)) #View size factors

gvst <- vst(gdds, blind=FALSE) #apply a variance stabilizing transformation to minimize effects of small counts and normalize to library size

### Principal component plot of samples
gPCAdata <- plotPCA(gvst, intgroup = c("temp"), returnData=TRUE)
percentVar <- round(100*attr(gPCAdata, "percentVar")) #plot PCA of samples with all data
PCA <- ggplot(gPCAdata, aes(PC1, PC2, color=temp)) + 
  geom_point(size=3) +
  geom_text(aes(label = name), vjust = -1, size = 3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #scale_color_manual(labels = c("adult_meso","adult_shal","planu_meso", "planu_shal"), values = c("adult_meso"="blue", "adult_shal"="indianred3", "planu_meso"="deepskyblue", "planu_shal"="orange")) +
  coord_fixed() +
  theme_bw() + #Set background color
  theme(panel.border = element_blank(), # Set border
        #panel.grid.major = element_blank(), #Set major gridlines
        #panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background
  theme(legend.position = ("top")); PCA #set title attributes`

ggsave(file = "PCA_all_vst.png", PCA)


#### Compile GWENA Dataset

# GWENA support expression matrix data coming from either RNA-seq or microarray experiments. 
# Expression data have to be stored as text or spreadsheet files and formatted with genes as columns 
# and samples as rows 

#Transpose the filtered gene count matrix so that the gene IDs are rows and the sample IDs are columns.
datExpr <- as.data.frame(t(assay(gvst))) #transpose to output to a new data frame with the column names as row names. And make all data numeric

#Check for genes and samples with too many missing values with goodSamplesGenes. There shouldn't be any because we performed pre-filtering
gsg = goodSamplesGenes(datExpr, verbose = 3)
# [1] allOK is TRUE

# Cluster the samples to look for obvious outliers
#Look for outliers by examining the sample tree:
sampleTree = hclust(dist(datExpr), method = "average")

# Plot the sample tree
pdf(paste0('sampleTree','.pdf'))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


# Number of genes
ncol(datExpr)
#> [1] 25423
# Number of samples
nrow(datExpr)
#> [1] 21

# Checking expression data set is correctly defined
is_data_expr(datExpr)
# $bool
# [1] TRUE

## Network building
# Gene co-expression networks are an ensemble of genes (nodes) linked to each other (edges) 
# according to the strength of their relation.
#GWENA uses Spearman correlation by default. It is less sensitive to outliers which are 
#frequent in transcriptomics datasets 

# I will use a signed network because we have a relatively high softPower, according 
# to >12 (https://peterlangfelder.com/2018/11/25/__trashed/). 
# Moreover, in expression data where you are interested in when expression on one gene increases or decreases 
# with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).


threads_to_use <- 12
threads_to_use <- 14
# net <- build_net(datExpr, cor_func = "spearman", 
#                  n_threads = threads_to_use)

net <- build_net(datExpr, 
                 cor_func = "spearman", 
                 n_threads = threads_to_use, 
                 network_type = "signed") 

net <- build_net(datExpr, 
                 cor_func = "pearson", 
                 n_threads = threads_to_use, 
                 network_type = "signed") 


# Power selected :
net$metadata$power
#18

# Fit of the power law to data ($R^2$) :
fit_power_table <- net$metadata$fit_power_table
fit_power_table[fit_power_table$Power == net$metadata$power, "SFT.R.sq"]
#0.7342994

# Plot R^2 (SFT.R.sq) vs Power
library(ggplot2)
fit_power_table <- net$metadata$fit_power_table
ggplot(fit_power_table, aes(x = Power, y = SFT.R.sq)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  labs(title = "Scale-Free Topology Fit (R²) vs Power",
       x = "Soft-thresholding Power",
       y = "Scale-Free Topology Fit (R²)") +
  theme_minimal()

####  Modules detection
# At this point, the network is a complete graph: all nodes are connected to all other nodes with different strengths. 
# Because gene co-expression networks have a scale free property, groups of genes are strongly linked with one another. 
# In co-expression networks these groups are called modules and assumed to be representative of genes working together to a 
# common set of functions.

modules <- detect_modules(datExpr, 
                            net$network, 
                            detailled_result = TRUE,
                            merge_threshold = 0.25) #merge modules with >75% similarity 

#Important: Module 0 contains all genes that did not fit into any modules.

# Since this operation tends to create multiple smaller modules with highly similar expression profile 
# (based on the eigengene of each), they are usually merged into one.

# Number of modules before merging :
length(unique(modules$modules_premerge))
#> [1] 36
# Number of modules after merging: 
length(unique(modules$modules))
#> [1] 4

layout_mod_merge <- plot_modules_merge(
  modules_premerge = modules$modules_premerge, 
  modules_merged = modules$modules)

#Resulting modules contain more genes whose repartition can be seen by a simple barplot.
ggplot2::ggplot(data.frame(modules$modules %>% stack), 
                ggplot2::aes(x = ind)) + ggplot2::stat_count() +
  ggplot2::ylab("Number of genes") +
  ggplot2::xlab("Module")

#Each of the modules presents a distinct profile, which can be plotted in two figures to separate the positive (+ facet) and negative (- facet) correlations profile. As a summary of this profile, the eigengene (red line) is displayed to act as a signature.
plot_expression_profiles(datExpr, modules$modules)

















########## WGCNA

### Network construction and consensus module detection
# Choosing a soft-thresholding power: Analysis of network topology 
# The soft thresholding power is the number to which the co-expression similarity is raised to calculate adjacency. 

#Choose a set of soft-thresholding powers
powers <- c(seq(from = 1, to=19, by=2), c(21:30)) #Create a string of numbers from 1 through 10, and even numbers from 10 through 20
powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2))

#Call the network topology analysis function
sft <-pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
#Scale-free topology fit index as a function of the soft-thresholding power
pdf(paste0('network','.pdf'))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# # # this line corresponds to using an R^2 cut-off
abline(h=0.8,col="red")
# # # Mean connectivity as a function of the soft-thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#The lowest scale-free topology fit index R^2 recommended by Langfelder and Horvath is 0.8. 
#From the graph, it appears that our soft thresholding power is 5 because it is the lowest 
#power before the R^2=0.8 threshold that maximizes with model fit (5 is the number right above the red line).

### Network construction and module detection:
# Co-expression adjacency and topological overlap matrix similarity
# Co-expression similarity and adjacency, using the soft thresholding power 5 and translate the adjacency into topological overlap matrix to calculate 
# the corresponding dissimilarity. 
#I will use a signed network because in expression data where you are interested in when expression on one gene increases or 
#decreases with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).

options(stringsAsFactors = FALSE)
enableWGCNAThreads() #Allow multi-threading within WGCNA

#Run analysis
softPower=5 #Set softPower to 5
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, TOM, dissTOM, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/adjTOM.RData")
save(dissTOM, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/dissTOM.RData") 

# Clustering using TOM
#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/dissTOMClustering.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()


## Modules identification
# Module identification is cutting the branches off the tree in the dendrogram above. We want large modules, so we set the minimum module size 
# relatively high (minimum size = 30).

minModuleSize = 30 #default value used most often
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize)
table(dynamicMods) #list modules and respective sizes
save(dynamicMods, geneTree, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/dyMod_geneTree.RData")

dyMod_geneTree
## [1] "dynamicMods" "geneTree"

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
  #  blue     brown     green turquoise    yellow 
  #  7849       482       147     16569       376 

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/dissTOMColorClustering.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar or choose not to merge
# Plot module similarity based on eigengene value

#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = 5)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster again and plot the results
METree = flashClust(as.dist(MEDiss), method = "average")

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >85% eigengene similarity (most studies use 80-90% similarity)

MEDissThres= 0.20 #merge modules that are 80% similar

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/eigengeneClustering2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/mergedClusters.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
ncol(MEs) 
# [1] 4
#Instead of 103 modules, we now have 56, a much more reasonable number.

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to temp, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the temp for numeric values

treatmentinfo$temp <- factor(paste0("t", as.character(treatmentinfo$temp)))
#allTraits <- names(treatmentinfo$temp)
allTraits <-levels(treatmentinfo$temp)
allTraits$tt12 <- c(0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0)
allTraits$tt18 <- c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0)
allTraits$tt25 <- c(0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0)
allTraits$tt26.8<-c(1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0)
allTraits$tt30 <- c(0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1)
allTraits$tt35 <- c(0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0)

datTraits <- as.data.frame(allTraits)
datTraits <- datTraits[ , -(1:6)]
dim(datTraits)
## [1] 21 6

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)
#       tt12 tt18 tt25 tt26.8 tt30 tt35
# B10    0    0    0      1    0    0
# B6     1    0    0      0    0    0
# B7     0    0    0      0    0    1
# B8     0    0    1      0    0    0
# B9     0    0    0      0    1    0
# F6     1    0    0      0    0    0
# F7     0    0    0      0    0    1
# F8     0    0    1      0    0    0
# F9     0    0    0      0    1    0
# G1     0    1    0      0    0    0
# G10    0    0    0      1    0    0
# G6     1    0    0      0    0    0
# G7     0    0    0      0    0    1
# G8     0    0    1      0    0    0
# G9     0    0    0      0    1    0
# H1     0    1    0      0    0    0
# H10    0    0    0      1    0    0
# H6     1    0    0      0    0    0
# H7     0    0    0      0    0    1
# H8     0    0    1      0    0    0
# H9     0    0    0      0    1    0

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=5)$eigengenes
MEs = orderMEs(MEs0)
names(MEs)
#[1] "MEblue"      "MEturquoise" "MEbrown"     "MEyellow"

## Correlations of traits and genes with eigengenes
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/Temp clustering based on module-trait correlation.pdf")
plot(moduleTraitTree, main = "Group clustering based on module-trait correlation", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datExpr)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

# Plot as clustered Heatmap
#add bold sigignificant p-values, dendrogram with WGCNA MEtree cut-off, module clusters

#Create list of pvalues for eigengene correlation with specific temperature
heatmappval <- signif(moduleTraitPvalue, 1)

#Make list of heatmap row colors
htmap.colors <- names(MEs)
htmap.colors <- gsub("ME", "", htmap.colors)

library(dendsort)
row_dend = dendsort(hclust(dist(moduleTraitCor)))
col_dend = dendsort(hclust(dist(t(moduleTraitCor))))

pdf(file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pacu/Module-trait-relationship-heatmap3_Pacu.pdf", height = 11.5, width = 8)
ht=Heatmap(moduleTraitCor, name = "Module-Trait Eigengene Correlation", 
        col = blueWhiteRed(50), 
        row_names_side = "left", row_dend_side = "left",
        #width = unit(4, "in"), height = unit(8.5, "in"), 
        column_order = 1:6, column_dend_reorder = FALSE, cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), column_split = 3, column_dend_height = unit(0.2, "in"),
        cluster_rows = METree, row_split = 10, row_gap = unit(2.5, "mm"), border = TRUE,
        cell_fun = function(j, i, x, y, w, h, col) {
        if(heatmappval[i, j] <= 0.05) {
            grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "bold"))
        }
        else {
            grid.text(sprintf("%s", heatmappval[i, j]), x, y, gp = gpar(fontsize = 8, fontface = "plain"))
        }},
        column_names_gp =  gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10, alpha = 0.75, border = TRUE, fill = htmap.colors))
draw(ht)
dev.off()


## Make dataframe for Strader plots

# View module eigengene data
head(MEs)
#          MEblue MEturquoise    MEbrown   MEyellow
# B10 -0.11059365  0.08986028 -0.2647869 -0.1181387
# B6  -0.11890986  0.16004599 -0.2013546 -0.1321223
# B7   0.44978874 -0.44155162 -0.3418433 -0.2462229
# B8  -0.11619856  0.10449661 -0.2554843 -0.1280705
# B9  -0.08325844  0.08042010 -0.2286536 -0.1171563
# F6  -0.14366328  0.16916203 -0.2635769 -0.1386889

names(MEs)
#[1] "MEblue"      "MEturquoise" "MEbrown"     "MEyellow"  

Strader_MEs <- MEs
Strader_MEs$temp <- treatmentinfo$temp
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)
#         MEblue MEturquoise    MEbrown   MEyellow   temp sample_id
# B10 -0.11059365  0.08986028 -0.2647869 -0.1181387 tt26.8       B10
# B6  -0.11890986  0.16004599 -0.2013546 -0.1321223   tt12        B6
# B7   0.44978874 -0.44155162 -0.3418433 -0.2462229   tt35        B7
# B8  -0.11619856  0.10449661 -0.2554843 -0.1280705   tt25        B8
# B9  -0.08325844  0.08042010 -0.2286536 -0.1171563   tt30        B9
# F6  -0.14366328  0.16916203 -0.2635769 -0.1386889   tt12        F6

# Rename columns 1 to 4
colnames(Strader_MEs)[1:4] <- c("Module1", "Module2", "Module3", "Module4")

head(Strader_MEs)
#        Module1     Module2    Module3    Module4   temp sample_id
# B10 -0.11059365  0.08986028 -0.2647869 -0.1181387 tt26.8       B10
# B6  -0.11890986  0.16004599 -0.2013546 -0.1321223   tt12        B6
# B7   0.44978874 -0.44155162 -0.3418433 -0.2462229   tt35        B7
# B8  -0.11619856  0.10449661 -0.2554843 -0.1280705   tt25        B8
# B9  -0.08325844  0.08042010 -0.2286536 -0.1171563   tt30        B9
# F6  -0.14366328  0.16916203 -0.2635769 -0.1386889   tt12        F6

library(tidyr)
library(ggplot2)

# Convert to long format for ggplot
Strader_MEs_long <- Strader_MEs %>%
  pivot_longer(
    cols = starts_with("Module"),
    names_to = "Module",
    values_to = "Eigengene"
  )

# Boxplot with dots, faceted by Module
ggplot(Strader_MEs_long, aes(x = temp, y = Eigengene)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(aes(color = sample_id), width = 0.2, size = 2, show.legend = FALSE) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_bw() +
  labs(title = "Eigengene Expression by Module and Temperature",
       x = "Temperature Group",
       y = "Module Eigengene Value")


# Filter for Module1 and Module2 only
Strader_MEs_long_subset <- Strader_MEs_long %>%
  filter(Module %in% c("Module1", "Module2"))

# Boxplot with dots, faceted by Module (only Module1 and Module2)
ggplot(Strader_MEs_long_subset, aes(x = temp, y = Eigengene)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(aes(color = sample_id), width = 0.2, size = 2, show.legend = FALSE) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_bw() +
  labs(title = "Eigengene Expression by Module and Temperature (Module1 & Module2)",
       x = "Temperature Group",
       y = "Module Eigengene Value")

# Define ordered temperature levels
Strader_MEs_long_subset$temp <- factor(Strader_MEs_long_subset$temp, 
                                       levels = c("tt12", "tt18", "tt25", "tt26.8", "tt30", "tt35"))

# Create a palette of cyan shades (light to dark)
cyan_palette <- colorRampPalette(c("#76e4f7", "#00bfff", "#0057b7"))(length(levels(Strader_MEs_long_subset$temp)))

# Plot with custom colors for boxplots and dots, and no grid
ggplot(Strader_MEs_long_subset, aes(x = temp, y = Eigengene, group = temp)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp), width = 0.2, size = 2, show.legend = FALSE) +
  facet_wrap(~ Module, scales = "free_y") +
  scale_color_manual(values = cyan_palette) +
  scale_x_discrete(labels = c("tt12" = "12", "tt18" = "18", "tt25" = "25", "tt26.8" = "26.8", "tt30" = "30", "tt35" = "35")) +
  theme_bw() +
  labs(x = "Temperature",
       y = "Module Eigengene Value") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )

### Use R package segmented
#segmented is designed for regression models where you want to detect one 
# or more breakpoints (change-points) in the relationship between a numeric 
# predictor (here, temperature) and a response (eigengene value).

library(segmented)

# Make sure temp is numeric (remove "tt" if needed)
Strader_MEs_long$temp_num <- as.numeric(sub("^tt", "", Strader_MEs_long$temp))

# Filter for the module of interest (e.g., Module1)
df_mod1 <- Strader_MEs_long %>% filter(Module == "Module1")

# Fit linear model to all replicates
fit <- lm(Eigengene ~ temp_num, data = df_mod1)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit <- segmented(fit, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit)
plot(seg_fit)
# 	***Regression Model with Segmented Relationship(s)***
# Call: 
# segmented.lm(obj = fit, seg.Z = ~temp_num, npsi = 1)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 29.596  0.185
# Coefficients of the linear terms:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.1272243  0.0244934  -5.194 7.31e-05 ***
# temp_num     0.0005576  0.0011518   0.484    0.634    
# U1.temp_num  0.1024257  0.0038244  26.782       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.02579 on 17 degrees of freedom
# Multiple R-Squared: 0.9887,  Adjusted R-squared: 0.9867 
# Boot restarting based on 6 samples. Last fit:
# Convergence attained in 2 iterations (rel. change 1.3154e-14)

###### Output Section	- Meaning
# Estimated Break-Point -	Temperature where the trend changes (here, ~29.6°C)
# temp_num	- Slope before breakpoint (not significant here)
# U1.temp_num	- Change in slope after breakpoint (significant, large increase)
# R-squared	- Proportion of variance explained by the model (very high here)

# Before 29.6°C: Eigengene value changes very little with temperature.
# After 29.6°C: Eigengene value increases sharply with temperature.
# The breakpoint at ~29.6°C is where the trend in eigengene value changes.

##  Get fitted values and confidence intervals
# Get predicted values and confidence intervals
library(segmented)
library(boot)

# Set up the bootstrap function
# Function to fit segmented and predict at grid points
boot_seg <- function(data, indices, grid_temp) {
  d <- data[indices, ]
  fit <- lm(Eigengene ~ temp_num, data = d)
  seg_fit <- try(segmented(fit, seg.Z = ~temp_num, npsi = 1), silent = TRUE)
  if (inherits(seg_fit, "try-error")) {
    return(rep(NA, length(grid_temp)))
  }
  pred <- predict(seg_fit, newdata = data.frame(temp_num = grid_temp))
  return(pred)
}

#Run the bootstrap
set.seed(123)
grid_temp <- seq(min(df_mod1$temp_num), max(df_mod1$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod1,
  statistic = function(data, indices) boot_seg(data, indices, grid_temp),
  R = 1000 # Number of bootstrap replicates
)

#Calculate confidence intervals
# Each row is a bootstrap, each column is a grid point
boot_preds <- boot_out$t

# Calculate 2.5% and 97.5% quantiles for each grid point
ci_lower <- apply(boot_preds, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(boot_preds, 2, quantile, probs = 0.975, na.rm = TRUE)
fit_median <- apply(boot_preds, 2, median, na.rm = TRUE)

# Prepare for plotting
pred_df1_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod1$temp_num <- as.numeric(sub("^tt", "", df_mod1$temp))

df_mod1$temp_num_f <- factor(df_mod1$temp_num)
ggplot(df_mod1, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df1_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df1_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = cyan_palette) +
  scale_x_continuous(
    breaks = unique(df_mod1$temp_num),
    labels = gsub("^t", "", unique(df_mod1$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module1: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


# ##  Get fitted values and confidence intervals
# # Get predicted values and confidence intervals
# pred_df <- data.frame(
#   temp_num = seq(min(df_mod1$temp_num), max(df_mod1$temp_num), length.out = 200)
# )
# pred <- predict(seg_fit, newdata = pred_df, se.fit = TRUE)

# pred_df$fit <- pred$fit
# pred_df$se.fit <- pred$se.fit
# pred_df$lower <- pred$fit - 1.96 * pred$se.fit
# pred_df$upper <- pred$fit + 1.96 * pred$se.fit

# library(ggplot2)

# ggplot(df_mod1, aes(x = temp_num, y = Eigengene)) +
#   geom_point(alpha = 0.6) +
#   geom_line(data = pred_df, aes(x = temp_num, y = fit), color = "blue", linewidth = 1) +
#   geom_ribbon(
#     data = pred_df,
#     aes(x = temp_num, ymin = lower, ymax = upper),
#     alpha = 0.2, fill = "blue",
#     inherit.aes = FALSE
#   ) +
#   labs(x = "Temperature (°C)", y = "Module Eigengene Value",
#        title = "Segmented Regression with 95% Confidence Interval") +
#   theme_bw()

# ##Overlay on your Strader plot
# df_mod1$temp_num <- as.numeric(sub("^tt", "", df_mod1$temp))

# df_mod1$temp_num_f <- factor(df_mod1$temp_num)
# ggplot(df_mod1, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
#   geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
#   geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
#   geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
#   geom_line(data = pred_df, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
#   geom_ribbon(data = pred_df, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
#   scale_color_manual(values = cyan_palette) +
#   scale_x_continuous(
#     breaks = unique(df_mod1$temp_num),
#     labels = gsub("^tt", "", unique(df_mod1$temp))
#   ) +
#   theme_bw() +
#   labs(x = "Temperature (°C)", y = "Module Eigengene Value",
#        title = "Module1: Strader Plot with Segmented Regression") +
#   theme(
#     strip.background = element_rect(fill = "#e0ffff", color = NA),
#     panel.grid = element_blank()
#   )


# Filter for the module of interest (e.g., Module2)
df_mod2 <- Strader_MEs_long %>% filter(Module == "Module2")

# Fit linear model to all replicates
fit <- lm(Eigengene ~ temp_num, data = df_mod2)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit <- segmented(fit, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit)
plot(seg_fit)
# 	***Regression Model with Segmented Relationship(s)***
# Call: 
# segmented.lm(obj = fit, seg.Z = ~temp_num, npsi = 1)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 29.783  0.285
# Coefficients of the linear terms:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.182871   0.036387   5.026 0.000104 ***
# temp_num    -0.003278   0.001711  -1.915 0.072414 .  
# U1.temp_num -0.097599   0.005681 -17.179       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.03831 on 17 degrees of freedom
# Multiple R-Squared: 0.9751,  Adjusted R-squared: 0.9706 
# Boot restarting based on 6 samples. Last fit:
# Convergence attained in 2 iterations (rel. change 7.2697e-13)

##  Get fitted values and confidence intervals
# Get predicted values and confidence intervals
library(segmented)
library(boot)

# Set up the bootstrap function
# Function to fit segmented and predict at grid points
boot_seg <- function(data, indices, grid_temp) {
  d <- data[indices, ]
  fit <- lm(Eigengene ~ temp_num, data = d)
  seg_fit <- try(segmented(fit, seg.Z = ~temp_num, npsi = 1), silent = TRUE)
  if (inherits(seg_fit, "try-error")) {
    return(rep(NA, length(grid_temp)))
  }
  pred <- predict(seg_fit, newdata = data.frame(temp_num = grid_temp))
  return(pred)
}

#Run the bootstrap
set.seed(123)
grid_temp <- seq(min(df_mod2$temp_num), max(df_mod2$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod2,
  statistic = function(data, indices) boot_seg(data, indices, grid_temp),
  R = 1000 # Number of bootstrap replicates
)

#Calculate confidence intervals
# Each row is a bootstrap, each column is a grid point
boot_preds <- boot_out$t

# Calculate 2.5% and 97.5% quantiles for each grid point
ci_lower <- apply(boot_preds, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(boot_preds, 2, quantile, probs = 0.975, na.rm = TRUE)
fit_median <- apply(boot_preds, 2, median, na.rm = TRUE)

# Prepare for plotting
pred_df2_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod2$temp_num <- as.numeric(sub("^tt", "", df_mod2$temp))

df_mod2$temp_num_f <- factor(df_mod2$temp_num)
ggplot(df_mod2, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df2_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df2_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = cyan_palette) +
  scale_x_continuous(
    breaks = unique(df_mod2$temp_num),
    labels = gsub("^t", "", unique(df_mod2$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module2: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


# ##  Get fitted values and confidence intervals
# # Get predicted values and confidence intervals
# pred_df <- data.frame(
#   temp_num = seq(min(df_mod2$temp_num), max(df_mod2$temp_num), length.out = 200)
# )
# pred <- predict(seg_fit, newdata = pred_df, se.fit = TRUE)

# pred_df$fit <- pred$fit
# pred_df$se.fit <- pred$se.fit
# pred_df$lower <- pred$fit - 1.96 * pred$se.fit
# pred_df$upper <- pred$fit + 1.96 * pred$se.fit

# library(ggplot2)

# ggplot(df_mod2, aes(x = temp_num, y = Eigengene)) +
#   geom_point(alpha = 0.6) +
#   geom_line(data = pred_df, aes(x = temp_num, y = fit), color = "blue", linewidth = 1) +
#   geom_ribbon(
#     data = pred_df,
#     aes(x = temp_num, ymin = lower, ymax = upper),
#     alpha = 0.2, fill = "blue",
#     inherit.aes = FALSE
#   ) +
#   labs(x = "Temperature (°C)", y = "Module 2 Eigengene Value",
#        title = "Segmented Regression with 95% Confidence Interval") +
#   theme_bw()


# ##Overlay on your Strader plot
# df_mod2$temp_num <- as.numeric(sub("^tt", "", df_mod2$temp))

# df_mod2$temp_num_f <- factor(df_mod2$temp_num)
# ggplot(df_mod2, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
#   geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
#   geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
#   geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
#   geom_line(data = pred_df, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
#   geom_ribbon(data = pred_df, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
#   scale_color_manual(values = cyan_palette) +
#   scale_x_continuous(
#     breaks = unique(df_mod2$temp_num),
#     labels = gsub("^tt", "", unique(df_mod2$temp))
#   ) +
#   theme_bw() +
#   labs(x = "Temperature (°C)", y = "Module Eigengene Value",
#        title = "Module2: Strader Plot with Segmented Regression") +
#   theme(
#     strip.background = element_rect(fill = "#e0ffff", color = NA),
#     panel.grid = element_blank()
#   )


##Save only the essentials for Network Analysis

save(
  datExpr,     
  moduleColors,
  control_samples,
  treat30_samples,
  treat35_samples,
  edge_threshold,
  file = "network_analysis_essentials_Pacu.RData"
)







############## Find Hub Genes and Network Analysis

# Export the network of a specific module ("blue")
blue_genes <- colnames(datExpr)[moduleColors == "blue"] ##7849
# Subset your expression matrix to blue module genes
datExpr_blue <- datExpr[, blue_genes]

control_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "tt26.8"]
treat30_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "tt30"]
treat35_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "tt35"]

##Save only the essentials for Network Analysis

save(
  datExpr,     
  moduleColors,
  control_samples,
  treat30_samples,
  treat35_samples,
  file = "network_analysis_essentials_Pacu.RData"
)

# Subset expression data for each group
datExpr_control_blue <- datExpr[control_samples, blue_genes ]
datExpr_treat30_blue<- datExpr[treat30_samples, blue_genes ]
datExpr_treat35_blue<- datExpr[treat35_samples, blue_genes]

# Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# kWithin: sum of connection strengths with other module genes
#add genes names to connectivity vectors
kWithin_control <- softConnectivity(datExpr_control_blue)
names(kWithin_control) <- colnames(datExpr_control_blue)

kWithin_treat30 <- softConnectivity(datExpr_treat30_blue)
names(kWithin_treat30) <- colnames(datExpr_treat30_blue)

kWithin_treat35 <- softConnectivity(datExpr_treat35_blue)
names(kWithin_treat35) <- colnames(datExpr_treat35_blue)

# Number of top 10% hub genes for each group
#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module
n_hubs_control <- ceiling(0.10 * length(kWithin_control))
n_hubs_treat30 <- ceiling(0.10 * length(kWithin_control))
n_hubs_treat35 <- ceiling(0.10 * length(kWithin_control))

# Get top 10% hub genes for each group
top10pct_hubs_control <- names(sort(kWithin_control, decreasing = TRUE))[1:n_hubs_control]
top10pct_hubs_treat30 <- names(sort(kWithin_treat30, decreasing = TRUE))[1:n_hubs_treat30]
top10pct_hubs_treat35 <- names(sort(kWithin_treat35, decreasing = TRUE))[1:n_hubs_treat35]

cat("Top 10% hub genes (control):\n")
print(top10pct_hubs_control)
cat("Top 10% hub genes (treat30):\n")
print(top10pct_hubs_treat30)
cat("Top 10% hub genes (treat35):\n")
print(top10pct_hubs_treat35)

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_blue.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_blue.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_blue.csv", row.names = FALSE)

#### Identify the top Hub Gene per module and temp

# Number of top 10% hub genes for each group
n_hubs <- ceiling(0.10 * length(blue_genes))

# Get top 10% hub genes for each group
top10pct_hubs_control <- names(sort(kWithin_control, decreasing = TRUE))[1:n_hubs]
top10pct_hubs_treat30 <- names(sort(kWithin_treat30, decreasing = TRUE))[1:n_hubs]
top10pct_hubs_treat35 <- names(sort(kWithin_treat35, decreasing = TRUE))[1:n_hubs]

# Find the top hub gene (highest connectivity) within the top 10% for each group
top_hub_control <- top10pct_hubs_control[1]
top_hub_treat30 <- top10pct_hubs_treat30[1]
top_hub_treat35 <- top10pct_hubs_treat35[1]

cat("Top hub gene (control, top 10%):", top_hub_control, "\n")
#Top hub gene (control, top 10%): Pocillopora_acuta_HIv2___TS.g18650.t1 
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Top hub gene (treat30, top 10%): Pocillopora_acuta_HIv2___RNAseq.g11463.t1 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Top hub gene (treat35, top 10%): Pocillopora_acuta_HIv2___RNAseq.g19582.t1 


####### Compare module preservation across conditions 
####### Highlight Nodes That Change Hub Status
# Highlight Top 10% Hubs
# Hub changing status between control and temps at break-points per (30 and 35) for all modules

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)

#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module

# Select only the top 10 hub genes for each condition (less computationally heavy for plotting)
top10_genes_control <- top10pct_hubs_control[1:10]
top10_genes_treat30 <- top10pct_hubs_treat30[1:10]
top10_genes_treat35 <- top10pct_hubs_treat35[1:10]

# Subset expression data to these genes
datExpr_control_top10 <- datExpr_control_blue[, top10_genes_control, drop = FALSE]
datExpr_treat30_top10 <- datExpr_treat30_blue[, top10_genes_treat30, drop = FALSE]
datExpr_treat35_top10 <- datExpr_treat35_blue[, top10_genes_treat35, drop = FALSE]

# Build adjacency/correlation matrices
adj_control_top10 <- abs(cor(datExpr_control_top10, method = "pearson"))
adj_treat30_top10 <- abs(cor(datExpr_treat30_top10, method = "pearson"))
adj_treat35_top10 <- abs(cor(datExpr_treat35_top10, method = "pearson"))

# Create igraph objects
g_control_top10 <- graph_from_adjacency_matrix(adj_control_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
g_treat30_top10 <- graph_from_adjacency_matrix(adj_treat30_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
g_treat35_top10 <- graph_from_adjacency_matrix(adj_treat35_top10, mode = "undirected", weighted = TRUE, diag = FALSE)

# Highlight the top hub gene (by kWithin) in red, others in gray
V(g_control_top10)$color <- ifelse(names(V(g_control_top10)) == top_hub_control, "red", "gray")
V(g_treat30_top10)$color <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, "red", "gray")
V(g_treat35_top10)$color <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, "red", "gray")

V(g_control_top10)$label <- ifelse(names(V(g_control_top10)) == top_hub_control, top_hub_control, "")
V(g_treat30_top10)$label <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, top_hub_treat30, "")
V(g_treat35_top10)$label <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, top_hub_treat35, "")

# Plot networks
pdf("network_top10genes_blue.pdf", width = 12, height = 8)
par(mfrow = c(1, 3))
plot(g_control_top10, main = "Control: Top 10 Hubs", vertex.label = V(g_control_top10)$label, vertex.size = 8)
plot(g_treat30_top10, main = "Treat30: Top 10 Hubs", vertex.label = V(g_treat30_top10)$label, vertex.size = 8)
plot(g_treat35_top10, main = "Treat35: Top 10 Hubs", vertex.label = V(g_treat35_top10)$label, vertex.size = 8)
par(mfrow = c(1, 1))
dev.off()

#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 5
all_genes <- blue_genes

# Build membership matrix
membership <- data.frame(
  Gene = all_genes,
  Control = as.integer(all_genes %in% top10pct_hubs_control),
  Treat30 = as.integer(all_genes %in% top10pct_hubs_treat30),
  Treat35 = as.integer(all_genes %in% top10pct_hubs_treat35)
)

# Filter to keep only genes that are a hub in at least one condition
membership <- membership %>%
  filter(Control == 1 | Treat30 == 1 | Treat35 == 1)

# Summarize overlap patterns
membership_summary <- membership %>%
  group_by(Control, Treat30, Treat35) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  mutate(Cluster = "Blue module Genes")

# Create a label for each interaction pattern
membership_summary$pattern <- interaction(membership_summary$Control, membership_summary$Treat30, membership_summary$Treat35)
pattern_counts <- setNames(membership_summary$Freq, membership_summary$pattern)
pattern_labels <- paste0(names(pattern_counts), " (n=", pattern_counts, ")")
names(pattern_labels) <- names(pattern_counts)

# Plot as before
ggplot(membership_summary,
       aes(axis1 = Cluster,
           axis2 = factor(Control, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           axis3 = factor(Treat30, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           axis4 = factor(Treat35, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           y = Freq)) +
  geom_alluvium(aes(fill = pattern), width = 1/8, alpha = 0.8) +
  geom_stratum(width = 1/6, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("Cluster", "Control", "Treat30", "Treat35"),
                   labels = c("Cluster" = "Blue Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Blue: Overlap of Top 10% Hub Genes Across Conditions",
       y = "Number of Genes", x = "") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

## bar plot
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Assume membership_summary and pattern_labels are already created as in your alluvial plot code

# Calculate total number of unique hub genes across all conditions
unique_hubs <- unique(c(top10pct_hubs_control, top10pct_hubs_treat30, top10pct_hubs_treat35))
total_hubs <- length(unique_hubs)

# Add percent column to membership_summary
membership_summary <- membership_summary %>%
  mutate(Percent = 100 * Freq / total_hubs)

# For plotting, order patterns by frequency (optional)
membership_summary$pattern <- factor(
  membership_summary$pattern,
  levels = membership_summary$pattern[order(-membership_summary$Freq)]
)

# Assign colors (same as alluvial)
bar_colors <- setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels))

barBlue <- ggplot(membership_summary, aes(x = "Blue", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Blue: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )




####### repeat for turquoise module
# Export the network of a specific module ("turquoise")
turquoise_genes <- colnames(datExpr)[moduleColors == "turquoise"] ##16569
# Subset your expression matrix to turquoise module genes
datExpr_turquoise <- datExpr[, turquoise_genes]

# Subset expression data for each group
datExpr_control_turquoise <- datExpr[control_samples, turquoise_genes ]
datExpr_treat30_turquoise<- datExpr[treat30_samples, turquoise_genes ]
datExpr_treat35_turquoise<- datExpr[treat35_samples, turquoise_genes]

# Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# kWithin: sum of connection strengths with other module genes
#add genes names to connectivity vectors
kWithin_control <- softConnectivity(datExpr_control_turquoise)
names(kWithin_control) <- colnames(datExpr_control_turquoise)

kWithin_treat30 <- softConnectivity(datExpr_treat30_turquoise)
names(kWithin_treat30) <- colnames(datExpr_treat30_turquoise)

kWithin_treat35 <- softConnectivity(datExpr_treat35_turquoise)
names(kWithin_treat35) <- colnames(datExpr_treat35_turquoise)

# Number of top 10% hub genes for each group
#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module
n_hubs_control <- ceiling(0.10 * length(kWithin_control))
n_hubs_treat30 <- ceiling(0.10 * length(kWithin_control))
n_hubs_treat35 <- ceiling(0.10 * length(kWithin_control))

# Get top 10% hub genes for each group
top10pct_hubs_control <- names(sort(kWithin_control, decreasing = TRUE))[1:n_hubs_control]
top10pct_hubs_treat30 <- names(sort(kWithin_treat30, decreasing = TRUE))[1:n_hubs_treat30]
top10pct_hubs_treat35 <- names(sort(kWithin_treat35, decreasing = TRUE))[1:n_hubs_treat35]

cat("Top 10% hub genes (control):\n")
print(top10pct_hubs_control)
cat("Top 10% hub genes (treat30):\n")
print(top10pct_hubs_treat30)
cat("Top 10% hub genes (treat35):\n")
print(top10pct_hubs_treat35)

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_turquoise.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_turquoise.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_turquoise.csv", row.names = FALSE)

#### Identify the top Hub Gene per module and temp

# Number of top 10% hub genes for each group
n_hubs <- ceiling(0.10 * length(turquoise_genes))

# Get top 10% hub genes for each group
top10pct_hubs_control <- names(sort(kWithin_control, decreasing = TRUE))[1:n_hubs]
top10pct_hubs_treat30 <- names(sort(kWithin_treat30, decreasing = TRUE))[1:n_hubs]
top10pct_hubs_treat35 <- names(sort(kWithin_treat35, decreasing = TRUE))[1:n_hubs]

# Find the top hub gene (highest connectivity) within the top 10% for each group
top_hub_control <- top10pct_hubs_control[1]
top_hub_treat30 <- top10pct_hubs_treat30[1]
top_hub_treat35 <- top10pct_hubs_treat35[1]

cat("Top hub gene (control, top 10%):", top_hub_control, "\n")
#Top hub gene (control, top 10%): Pocillopora_acuta_HIv2___RNAseq.g12020.t1
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Top hub gene (treat30, top 10%): Pocillopora_acuta_HIv2___RNAseq.g8212.t1 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Top hub gene (treat35, top 10%): Pocillopora_acuta_HIv2___RNAseq.g19245.t1 


####### Compare module preservation across conditions 
####### Highlight Nodes That Change Hub Status
# Highlight Top 10% Hubs
# Hub changing status between control and temps at break-points per (30 and 35) for all modules

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)

#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module

# Select only the top 10 hub genes for each condition (less computationally heavy for plotting)
top10_genes_control <- top10pct_hubs_control[1:10]
top10_genes_treat30 <- top10pct_hubs_treat30[1:10]
top10_genes_treat35 <- top10pct_hubs_treat35[1:10]

# Subset expression data to these genes
datExpr_control_top10 <- datExpr_control_turquoise[, top10_genes_control, drop = FALSE]
datExpr_treat30_top10 <- datExpr_treat30_turquoise[, top10_genes_treat30, drop = FALSE]
datExpr_treat35_top10 <- datExpr_treat35_turquoise[, top10_genes_treat35, drop = FALSE]

# Build adjacency/correlation matrices
adj_control_top10 <- abs(cor(datExpr_control_top10, method = "pearson"))
adj_treat30_top10 <- abs(cor(datExpr_treat30_top10, method = "pearson"))
adj_treat35_top10 <- abs(cor(datExpr_treat35_top10, method = "pearson"))

# Create igraph objects
g_control_top10 <- graph_from_adjacency_matrix(adj_control_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
g_treat30_top10 <- graph_from_adjacency_matrix(adj_treat30_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
g_treat35_top10 <- graph_from_adjacency_matrix(adj_treat35_top10, mode = "undirected", weighted = TRUE, diag = FALSE)

# Highlight the top hub gene (by kWithin) in red, others in gray
V(g_control_top10)$color <- ifelse(names(V(g_control_top10)) == top_hub_control, "red", "gray")
V(g_treat30_top10)$color <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, "red", "gray")
V(g_treat35_top10)$color <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, "red", "gray")

V(g_control_top10)$label <- ifelse(names(V(g_control_top10)) == top_hub_control, top_hub_control, "")
V(g_treat30_top10)$label <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, top_hub_treat30, "")
V(g_treat35_top10)$label <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, top_hub_treat35, "")

# Plot networks
pdf("network_top10genes_turquoise.pdf", width = 12, height = 8)
par(mfrow = c(1, 3))
plot(g_control_top10, main = "Control: Top 10 Hubs", vertex.label = V(g_control_top10)$label, vertex.size = 8)
plot(g_treat30_top10, main = "Treat30: Top 10 Hubs", vertex.label = V(g_treat30_top10)$label, vertex.size = 8)
plot(g_treat35_top10, main = "Treat35: Top 10 Hubs", vertex.label = V(g_treat35_top10)$label, vertex.size = 8)
par(mfrow = c(1, 1))
dev.off()

#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 5
all_genes <- turquoise_genes

# Build membership matrix
membership <- data.frame(
  Gene = all_genes,
  Control = as.integer(all_genes %in% top10pct_hubs_control),
  Treat30 = as.integer(all_genes %in% top10pct_hubs_treat30),
  Treat35 = as.integer(all_genes %in% top10pct_hubs_treat35)
)

# Filter to keep only genes that are a hub in at least one condition
membership <- membership %>%
  filter(Control == 1 | Treat30 == 1 | Treat35 == 1)

# Summarize overlap patterns
membership_summary <- membership %>%
  group_by(Control, Treat30, Treat35) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  mutate(Cluster = "Cluster 5 Genes")

# Create a label for each interaction pattern
membership_summary$pattern <- interaction(membership_summary$Control, membership_summary$Treat30, membership_summary$Treat35)
pattern_counts <- setNames(membership_summary$Freq, membership_summary$pattern)
pattern_labels <- paste0(names(pattern_counts), " (n=", pattern_counts, ")")
names(pattern_labels) <- names(pattern_counts)

# Plot as before
ggplot(membership_summary,
       aes(axis1 = Cluster,
           axis2 = factor(Control, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           axis3 = factor(Treat30, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           axis4 = factor(Treat35, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           y = Freq)) +
  geom_alluvium(aes(fill = pattern), width = 1/8, alpha = 0.8) +
  geom_stratum(width = 1/6, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("Cluster", "Control", "Treat30", "Treat35"),
                   labels = c("Cluster" = "Turquoise Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Turquoise: Overlap of Top 10% Hub Genes Across Conditions",
       y = "Number of Genes", x = "") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

## bar plot
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Assume membership_summary and pattern_labels are already created as in your alluvial plot code

# Calculate total number of unique hub genes across all conditions
unique_hubs <- unique(c(top10pct_hubs_control, top10pct_hubs_treat30, top10pct_hubs_treat35))
total_hubs <- length(unique_hubs)

# Add percent column to membership_summary
membership_summary <- membership_summary %>%
  mutate(Percent = 100 * Freq / total_hubs)

# For plotting, order patterns by frequency (optional)
membership_summary$pattern <- factor(
  membership_summary$pattern,
  levels = membership_summary$pattern[order(-membership_summary$Freq)]
)

# Assign colors (same as alluvial)
bar_colors <- setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels))

barTurquoise <- ggplot(membership_summary, aes(x = "Turquoise", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Turquoise: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


##### combine all bar plots into one figure
library(patchwork)
barBlue + barTurquoise 








########### Quantify network stability or rewiring across conditions:
# Edge Overlap / Jaccard Index
# Definition:
# The Jaccard index measures the similarity between two sets (here, edge sets).
# It is defined as:
# Jaccard = (Number of shared edges) / (Total number of unique edges in both networks)
# To calculate the Jaccard Index and percent edge change for the entire set of genes in the blue and turquoise modules (not just the top 30), 
# use all genes in the module for your adjacency/correlation matrices and edge comparisons.

# Use all genes in the blue and turquoise modules
all_blue_genes <- colnames(datExpr)[moduleColors == "blue"]
all_turquoise_genes <- colnames(datExpr)[moduleColors == "turquoise"]

# Subset expression data for each condition
datExpr_control_all <- datExpr[control_samples, all_blue_genes]
datExpr_treat30_all <- datExpr[treat30_samples, all_blue_genes]
datExpr_treat35_all <- datExpr[treat35_samples, all_blue_genes]

datExpr_control_all <- datExpr[control_samples, all_turquoise_genes]
datExpr_treat30_all <- datExpr[treat30_samples, all_turquoise_genes]
datExpr_treat35_all <- datExpr[treat35_samples, all_turquoise_genes]

#  Build Adjacency/Correlation Matrices
adj_control_all <- abs(cor(datExpr_control_all, method = "pearson"))
adj_treat30_all <- abs(cor(datExpr_treat30_all, method = "pearson"))
adj_treat35_all <- abs(cor(datExpr_treat35_all, method = "pearson"))

# Threshold and Remove Diagonal
edge_threshold <- 0.5
adj1_all <- (adj_control_all > edge_threshold)
adj2_all <- (adj_treat30_all > edge_threshold)
diag(adj1_all) <- 0
diag(adj2_all) <- 0

#or
adj1_all <- (adj_control_all > edge_threshold)
adj2_all <- (adj_treat35_all > edge_threshold) #### run separately for 30 and 35
diag(adj1_all) <- 0
diag(adj2_all) <- 0

# Get Edge Indices
edges_control_all <- which(adj1_all & upper.tri(adj1_all), arr.ind = TRUE)
edges_treat30_all <- which(adj2_all & upper.tri(adj2_all), arr.ind = TRUE)
#or
edges_treat35_all <- which(adj2_all & upper.tri(adj2_all), arr.ind = TRUE)

# Calculate Jaccard Index and Percent Change
# Convert to character for set operations
edges_control_set <- paste(edges_control_all[,1], edges_control_all[,2], sep = "-")
edges_treat30_set <- paste(edges_treat30_all[,1], edges_treat30_all[,2], sep = "-")
#or
edges_treat35_set <- paste(edges_treat35_all[,1], edges_treat35_all[,2], sep = "-")

# Jaccard Index
shared_edges <- intersect(edges_control_set, edges_treat30_set)
all_edges <- union(edges_control_set, edges_treat30_set)
jaccard_index <- length(shared_edges) / length(all_edges)
cat("Jaccard index (all blue genes):", jaccard_index, "\n")
#Jaccard index (all blue genes): 0.4356688 ##30
# Interpretation:
# 1 = identical networks
# 0 = no shared edges
cat("Jaccard index (all turquoise genes):", jaccard_index, "\n")
#Jaccard index (all turquoise genes): 0.4298702 ##30

shared_edges <- intersect(edges_control_set, edges_treat35_set)
all_edges <- union(edges_control_set, edges_treat35_set)
jaccard_index <- length(shared_edges) / length(all_edges)
cat("Jaccard index (all blue genes):", jaccard_index, "\n")
#Jaccard index (all blue genes): 0.4331329 ##35
cat("Jaccard index (all turquoise genes):", jaccard_index, "\n")
#Jaccard index (all turquoise genes): 0.3083096 ##35




########## Global Network Properties
# Compare Network density

# Interpretation of Gene Network Density:
# High Density: Indicates a tightly connected network where most genes interact with each other.
# May suggest strong coordination and communication between genes.
# Could imply functional relationships or co-regulation of genes.
# Low Density: Indicates a more sparse network with fewer interactions.
# Might suggest more specialized functions or less interconnectedness between genes.
# Could indicate a network with more independent modules or pathways. 
# A denser network may be more resilient to perturbations because of alternative pathways and redundancies.
# Changes in network density can be associated with development, disease, or responses to environmental stimuli. 

#To statistically compare network density (or other global properties) across temperatures using permutations, you can use a permutation 
#(randomization) test. This approach is robust and publishable, especially when sample sizes are small.

# Workflow for permutation testing of network density:

library(igraph)
# 1. Calculate Observed Density Difference
# Suppose you want to compare control (tt26.8) vs treat30 (tt30):
# After thresholding and before creating the graph:
adj1_all <- (adj_control_all > edge_threshold)    #### select first from above the module of interest before running this (line 1288)
adj2_all <- (adj_treat30_all > edge_threshold)
diag(adj1_all) <- 0
diag(adj2_all) <- 0

adj1_all <- (adj_control_all > edge_threshold)    #### select first from above the module of interest before running this (line 1288)
adj2_all <- (adj_treat35_all > edge_threshold)
diag(adj1_all) <- 0
diag(adj2_all) <- 0

# Force symmetry by taking the maximum (logical OR) and explicitly set lower triangle
adj1_all <- adj1_all | t(adj1_all)
adj1_all[lower.tri(adj1_all)] <- t(adj1_all)[lower.tri(adj1_all)]
adj2_all <- adj2_all | t(adj2_all)
adj2_all[lower.tri(adj2_all)] <- t(adj2_all)[lower.tri(adj2_all)]

adj1_all[is.na(adj1_all)] <- FALSE
adj2_all[is.na(adj2_all)] <- FALSE

# Now create igraph objects
g_control <- graph_from_adjacency_matrix(adj1_all, mode = "undirected", diag = FALSE)
g_treat30 <- graph_from_adjacency_matrix(adj2_all, mode = "undirected", diag = FALSE)
g_treat35 <- graph_from_adjacency_matrix(adj2_all, mode = "undirected", diag = FALSE)

# Calculate observed densities
density_control <- edge_density(g_control)
density_treat30 <- edge_density(g_treat30)
density_treat35 <- edge_density(g_treat35)

obs_diff <- density_treat30 - density_control
obs_diff <- density_treat35 - density_control

# 2. Permutation Testing

library(parallel)
n_perm <- 200

all_samples <- c(control_samples, treat30_samples)
all_samples <- c(control_samples, treat35_samples)

n_control <- length(control_samples)
n_treat30 <- length(treat30_samples)
n_treat35 <- length(treat35_samples)

# Set up a cluster (use the number of physical cores you have)
cl <- makeCluster(4) # or detectCores() for all available

# Export needed variables and libraries to the cluster
clusterExport(cl, varlist = c("datExpr", "all_blue_genes", "n_control", "n_treat30", "all_samples", "edge_threshold"))
clusterExport(cl, varlist = c("datExpr", "all_blue_genes", "n_control", "n_treat35", "all_samples", "edge_threshold"))
clusterExport(cl, varlist = c("datExpr", "all_turquoise_genes", "n_control", "n_treat30", "all_samples", "edge_threshold"))
clusterExport(cl, varlist = c("datExpr", "all_turquoise_genes", "n_control", "n_treat35", "all_samples", "edge_threshold"))

clusterEvalQ(cl, library(igraph))

# Remove columns (genes) with zero variance in either group
zero_var_control <- apply(datExpr_control_all, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var_treat30 <- apply(datExpr_treat30_all, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var <- zero_var_control | zero_var_treat30
if (any(zero_var)) {
  datExpr_control_all <- datExpr_control_all[, !zero_var, drop = FALSE]
  datExpr_treat30_all <- datExpr_treat30_all[, !zero_var, drop = FALSE]
  all_turquoise_genes <- all_turquoise_genes[!zero_var]
}

perm_diffs <- parLapply(cl, 1:n_perm, function(i) {
  permuted <- sample(all_samples)
  perm_control <- permuted[1:n_control]
  perm_treat30 <- permuted[(n_control+1):(n_control+n_treat30)]
  
  datExpr_control_perm <- datExpr[perm_control, all_turquoise_genes]
  datExpr_treat30_perm <- datExpr[perm_treat30, all_turquoise_genes]
  
zero_var_control <- apply(datExpr_control_perm, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var_treat30 <- apply(datExpr_treat30_perm, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var <- zero_var_control | zero_var_treat30
if (any(zero_var)) {
  datExpr_control_perm <- datExpr_control_perm[, !zero_var, drop = FALSE]
  datExpr_treat30_perm <- datExpr_treat30_perm[, !zero_var, drop = FALSE]
}

  adj_control_perm <- abs(cor(datExpr_control_perm, method = "pearson"))
  adj_treat30_perm <- abs(cor(datExpr_treat30_perm, method = "pearson"))
  
  adj1_perm <- (adj_control_perm > edge_threshold)
  adj2_perm <- (adj_treat30_perm > edge_threshold)
  diag(adj1_perm) <- 0
  diag(adj2_perm) <- 0
  
  # Only force symmetry if needed
  if (!isSymmetric(adj1_perm)) adj1_perm <- adj1_perm | t(adj1_perm)
  if (!isSymmetric(adj2_perm)) adj2_perm <- adj2_perm | t(adj2_perm)
  adj1_perm[is.na(adj1_perm)] <- FALSE
  adj2_perm[is.na(adj2_perm)] <- FALSE
  
  g_control_perm <- graph_from_adjacency_matrix(adj1_perm, mode = "undirected", diag = FALSE)
  g_treat30_perm <- graph_from_adjacency_matrix(adj2_perm, mode = "undirected", diag = FALSE)
  
  density_control_perm <- edge_density(g_control_perm)
  density_treat30_perm <- edge_density(g_treat30_perm)
  
  density_treat30_perm - density_control_perm
})

stopCluster(cl)
perm_diffs <- unlist(perm_diffs)

# Calculate p-value (two-sided)
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
cat("Observed density difference_turquoise30:", obs_diff, "\n")  ###30
#Observed density difference: -0.09917594 
# The observed density difference is negative, meaning the network density is lower in the treatment 
# compared to control.
cat("Permutation p-value_turquoise30:", p_value, "\n")   ###30
# #Permutation p-value: 0.745 
# # 3. Interpretation
# # If p < 0.05: The difference in network density between conditions is statistically significant.
# # If p ≥ 0.05: The observed difference could be due to chance.


###35

# Remove columns (genes) with zero variance in either group
zero_var_control <- apply(datExpr_control_all, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var_treat35 <- apply(datExpr_treat35_all, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var <- zero_var_control | zero_var_treat35
if (any(zero_var)) {
  datExpr_control_all <- datExpr_control_all[, !zero_var, drop = FALSE]
  datExpr_treat35_all <- datExpr_treat35_all[, !zero_var, drop = FALSE]
  all_turquoise_genes <- all_turquoise_genes[!zero_var]
}

perm_diffs <- parLapply(cl, 1:n_perm, function(i) {
  permuted <- sample(all_samples)
  perm_control <- permuted[1:n_control]
  perm_treat35<- permuted[(n_control+1):(n_control+n_treat35)]
  
  datExpr_control_perm <- datExpr[perm_control, all_turquoise_genes]
  datExpr_treat35_perm <- datExpr[perm_treat35, all_turquoise_genes]

zero_var_control <- apply(datExpr_control_perm, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var_treat35<- apply(datExpr_treat35_perm, 2, function(x) sd(x, na.rm = TRUE) == 0)
zero_var <- zero_var_control | zero_var_treat35
if (any(zero_var)) {
  datExpr_control_perm <- datExpr_control_perm[, !zero_var, drop = FALSE]
  datExpr_treat35_perm <- datExpr_treat35_perm[, !zero_var, drop = FALSE]
}
  
  adj_control_perm <- abs(cor(datExpr_control_perm, method = "pearson"))
  adj_treat35_perm <- abs(cor(datExpr_treat35_perm, method = "pearson"))
  
  adj1_perm <- (adj_control_perm > edge_threshold)
  adj2_perm <- (adj_treat35_perm > edge_threshold)
  diag(adj1_perm) <- 0
  diag(adj2_perm) <- 0
  
  # Only force symmetry if needed
  if (!isSymmetric(adj1_perm)) adj1_perm <- adj1_perm | t(adj1_perm)
  if (!isSymmetric(adj2_perm)) adj2_perm <- adj2_perm | t(adj2_perm)
  adj1_perm[is.na(adj1_perm)] <- FALSE
  adj2_perm[is.na(adj2_perm)] <- FALSE
  
  g_control_perm <- graph_from_adjacency_matrix(adj1_perm, mode = "undirected", diag = FALSE)
  g_treat35_perm <- graph_from_adjacency_matrix(adj2_perm, mode = "undirected", diag = FALSE)
  
  density_control_perm <- edge_density(g_control_perm)
  density_treat35_perm <- edge_density(g_treat35_perm)
  
  density_treat35_perm - density_control_perm
})

stopCluster(cl)
perm_diffs <- unlist(perm_diffs)

# Calculate p-value (two-sided)
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
cat("Observed density difference_turquoise35:", obs_diff, "\n")  ###35
#Observed density difference:
# The observed density difference is negative, meaning the network density is lower in the treatment 
# compared to control.
cat("Permutation p-value_turquoise35:", p_value, "\n")   ###35
#Permutation p-value:
# 3. Interpretation
# If p < 0.05: The difference in network density between conditions is statistically significant.
# If p ≥ 0.05: The observed difference could be due to chance.