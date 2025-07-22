
############# WGCNA (Weighted Gene Co-expression Network Analysis) #############
############# Pcom (Porites compressa) RNAseq data #############
# Federica Scucchia June 2025

## load libraries
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

#treatment information
treatmentinfo <- read.csv("RNAseq_Pcom_data.csv", header = TRUE, sep = ";")

#gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix_noIso_Pcom2.csv", row.names="gene_id"), colClasses = double)

#remove samples below 5 million reads
gcount <- gcount %>% dplyr::select(-B7)
treatmentinfo <- treatmentinfo[!(treatmentinfo$sample_id %in% c("B7")), ]

# Reorder treatmentinfo to match the columns of gcount_filt
treatmentinfo <- treatmentinfo[match(colnames(gcount), treatmentinfo$sample_id), ]

## Quality-filter gene counts
#Set filter values for PoverA: smallest sample size per treat is 3, so 3/22 (22 samples) is 0.14
#This means that 3 out of 22 (0.14) samples need to have counts over 10.
#So P=14 percent of the samples have counts over A=10. 
filt <- filterfun(pOverA(0.14,10))

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
# [1] 44130

nrow(gcount_filt) #After
# [1] 28378


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

SF.gdds <- estimateSizeFactors( gdds ) #estimate size factors to determine if we can use vst to transform our data. Size factors should be less than four to use vst
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

ggsave(file = "PCA_all_vst_Pcom.png", PCA)


#### Compile WGCNA Dataset

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
#> [1] 28378
# Number of samples
nrow(datExpr)
#> [1] 22

# Checking expression data set is correctly defined
is_data_expr(datExpr)
# $bool
# [1] TRUE

#### Save all genes for Viseago (enrichment analysis)
# Extract gene names from column headers
gene_names <- colnames(datExpr)
# Save as background.txt (one gene per line, no header)
writeLines(gene_names, "background.txt")

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
#From the graph, it appears that our soft thresholding power is 7 because it is the lowest 
#power before the R^2=0.8 threshold that maximizes with model fit (7 is the number right above the red line).

### Network construction and module detection:
# Co-expression adjacency and topological overlap matrix similarity
# Co-expression similarity and adjacency, using the soft thresholding power 7and translate the adjacency into topological overlap matrix to calculate 
# the corresponding dissimilarity. 
#I will use a signed network because in expression data where you are interested in when expression on one gene increases or 
#decreases with expression level of another you would use a signed network (when you are interested in the direction of change, correlation and anti-correlation, you use a signed network).

options(stringsAsFactors = FALSE)
enableWGCNAThreads() #Allow multi-threading within WGCNA

#Run analysis
softPower=7 #Set softPower to 7
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, TOM, dissTOM, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/adjTOM.RData")
save(dissTOM, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/dissTOM.RData") 

# Clustering using TOM
#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/dissTOMClustering.pdf", width=20, height=20)
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
save(dynamicMods, geneTree, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/dyMod_geneTree.RData")

dyMod_geneTree
## [1] "dynamicMods" "geneTree"

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
# dynamicColors
# antiquewhite4         bisque4           black            blue           brown          brown4 
#              75             127             780            4910            3114             129 
#          coral1          coral2            cyan       darkgreen        darkgrey     darkmagenta 
#              75              75             342             265             250             170 
#  darkolivegreen      darkorange     darkorange2         darkred   darkseagreen4   darkslateblue 
#             171             213             134             274              78             120 

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/dissTOMColorClustering.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar or choose not to merge
# Plot module similarity based on eigengene value

#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors, softPower = 7)
MEs = MEList$eigengenes

#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

#Cluster again and plot the results
METree = flashClust(as.dist(MEDiss), method = "average")

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >80% eigengene similarity (most studies use 80-90% similarity)

MEDissThres= 0.20 #merge modules that are 80% similar

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/eigengeneClustering2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/mergedClusters.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
ncol(MEs) 
# [1] 38


# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to temp, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the temp for numeric values

treatmentinfo$temp <- factor(paste0("t", as.character(treatmentinfo$temp)))
#allTraits <- names(treatmentinfo$temp)
allTraits <-levels(treatmentinfo$temp)
allTraits$tt12 <- c(0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0)
allTraits$tt18 <- c(1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0)
allTraits$tt25 <- c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1)
allTraits$tt26.8<-c(0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0)
allTraits$tt30 <- c(0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0)
allTraits$tt35 <- c(0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0)

datTraits <- as.data.frame(allTraits)
datTraits <- datTraits[ , -(1:6)]
dim(datTraits)
## [1] 22 6

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)
#     tt12 tt18 tt25 tt26.8 tt30 tt35
# B1     0    1    0      0    0    0
# B10    0    0    0      1    0    0
# B6     1    0    0      0    0    0
# B8     0    0    1      0    0    0
# B9     0    0    0      0    1    0
# F10    0    0    0      1    0    0
# F7     0    0    0      0    0    1
# F9     0    0    0      0    1    0
# G10    0    0    0      1    0    0
# G6     1    0    0      0    0    0
# G7     0    0    0      0    0    1
# G9     0    0    0      0    1    0
# H1     0    1    0      0    0    0
# H10    0    0    0      1    0    0
# H6     1    0    0      0    0    0
# H7     0    0    0      0    0    1
# H8     0    0    1      0    0    0
# H9     0    0    0      0    1    0
# F6     1    0    0      0    0    0
# F8     0    0    1      0    0    0
# G1     0    1    0      0    0    0
# G8     0    0    1      0    0    0

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors,softPower=5)$eigengenes
MEs = orderMEs(MEs0)
names(MEs) #head
# [1] "MEdarkred"        "MEturquoise"      "MEblack"          "MEpink"           "MElightgreen"    
#  [6] "MEmediumpurple3"  "MEdarkturquoise"  "MEyellow4"        "MEskyblue2"       "MEfloralwhite"   
# [11] "MEcoral1"         "MEcyan"           "MEblue"           "MEgrey60"         "MEpaleturquoise" 

## Correlations of traits and genes with eigengenes
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/Temp clustering based on module-trait correlation.pdf")
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

pdf(file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Pcom/Module-trait-relationship-heatmap_Pcom.pdf", height = 9, width = 8)
ht=Heatmap(moduleTraitCor, name = "Module-Trait Eigengene Correlation", 
        col = blueWhiteRed(50), 
        row_names_side = "left", row_dend_side = "left",
        #width = unit(4, "in"), height = unit(8.5, "in"), 
        column_order = 1:6, column_dend_reorder = FALSE, cluster_columns = hclust(dist(t(moduleTraitCor)), method = "average"), column_split = 6, column_dend_height = unit(0.2, "in"),
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
#      MEdarkred MEturquoise    MEblack      MEpink MElightgreen MEmediumpurple3 MEdarkturquoise
# B1  -0.22664724 -0.02375949 -0.1633693 -0.11379582   -0.2374829     -0.25648779     -0.34484138
# B10  0.13565047  0.09591296 -0.1056327 -0.04978991   -0.1173814      0.15293867      0.08895811
# B6   0.06967086  0.16416343 -0.1582944 -0.04771290   -0.1613230     -0.13123134     -0.11420461
# B8   0.17911387  0.11708019 -0.0927854  0.01983302   -0.1011754      0.30128045      0.09824477

names(MEs)
#  [1] "MEdarkred"        "MEturquoise"      "MEblack"          "MEpink"           "MElightgreen"    
#  [6] "MEmediumpurple3"  "MEdarkturquoise"  "MEyellow4"        "MEskyblue2"       "MEfloralwhite"   

Strader_MEs <- MEs
Strader_MEs$temp <- treatmentinfo$temp
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)
#       MEdarkred MEturquoise    MEblack      MEpink MElightgreen MEmediumpurple3 MEdarkturquoise
# B1  -0.22664724 -0.02375949 -0.1633693 -0.11379582   -0.2374829     -0.25648779     -0.34484138
# B10  0.13565047  0.09591296 -0.1056327 -0.04978991   -0.1173814      0.15293867      0.08895811

# Calculate 10 over-arching expression patterns using mean eigengene for each module in a cluster
C1_Strader_MEs <- select(Strader_MEs, MEdarkred:MEturquoise)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)
C2_Strader_MEs <- select(Strader_MEs, MEblack:MEyellow4)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)
C3_Strader_MEs <- select(Strader_MEs, MEskyblue2:MEpaleturquoise)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)
C4_Strader_MEs <- select(Strader_MEs, MEhoneydew1)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)
C5_Strader_MEs <- select(Strader_MEs, MEantiquewhite4:MEyellow)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)
C6_Strader_MEs <- select(Strader_MEs, MEgreen:MEivory)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)
C7_Strader_MEs <- select(Strader_MEs, MEmediumorchid:MEindianred4)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)
C8_Strader_MEs <- select(Strader_MEs, MEmediumpurple2:MEnavajowhite2)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)
C9_Strader_MEs <- select(Strader_MEs, MEbrown4:MEdarkgreen)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)
C10_Strader_MEs <- select(Strader_MEs, MEsalmon4:MElightpink4)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

Strader_MEs$temp <- as.character(Strader_MEs$temp)
expressionProfile_data <- as.data.frame(cbind(temp = Strader_MEs$temp, cluster1= C1_Strader_MEs$Mean, cluster2 = C2_Strader_MEs$Mean, 
                          cluster3 = C3_Strader_MEs$Mean, cluster4 = C4_Strader_MEs$Mean,cluster5 = C5_Strader_MEs$Mean,cluster6 = C6_Strader_MEs$Mean,
                          cluster7 = C7_Strader_MEs$Mean,cluster8 = C8_Strader_MEs$Mean,cluster9 = C9_Strader_MEs$Mean,cluster10 = C10_Strader_MEs$Mean))

head(expressionProfile_data)
#    temp           cluster1            cluster2            cluster3            cluster4
# 1   t18 -0.125203368871824  -0.238283336890726  -0.139233175300969   0.189039412050515
# 2 t26.8  0.115781719170205  -0.004922996031784 -0.0953968799033358   0.484090743299738
# 3   t12  0.116917146569611  -0.125645469846288  -0.156813727864529   -0.15329821040734
# 4   t25     0.148097029453  0.0391409392243317  -0.101096903260113 -0.0308790113589991

# save data of selected clusters for GO enrichment analysis

### save mean eigengene values for cluster9
meanEigenClust9 <- expressionProfile_data$cluster9
write.csv(meanEigenClust9, file = "meanEigenClust9.csv")
write.csv(C9_Strader_MEs, file = "C9_Strader_MEs.csv")

### save mean eigengene values for cluster2
meanEigenClust2 <- expressionProfile_data$cluster2
write.csv(meanEigenClust2, file = "meanEigenClust2.csv")
write.csv(C2_Strader_MEs, file = "C2_Strader_MEs.csv")

### save mean eigengene values for cluster3
meanEigenClust3 <- expressionProfile_data$cluster3
write.csv(meanEigenClust3, file = "meanEigenClust3.csv")
write.csv(C3_Strader_MEs, file = "C3_Strader_MEs.csv")

### save mean eigengene values for cluster1
meanEigenClust1 <- expressionProfile_data$cluster1
write.csv(meanEigenClust1, file = "meanEigenClust1.csv")
write.csv(C1_Strader_MEs, file = "C1_Strader_MEs.csv")

### save mean eigengene values for cluster5
meanEigenClust4 <- expressionProfile_data$cluster4
write.csv(meanEigenClust4, file = "meanEigenClust4.csv")
write.csv(C5_Strader_MEs, file = "C4_Strader_MEs.csv")

### save mean eigengene values for cluster10
meanEigenClust10 <- expressionProfile_data$cluster10
write.csv(meanEigenClust10, file = "meanEigenClust10.csv")
write.csv(C10_Strader_MEs, file = "C10_Strader_MEs.csv")

### save mean eigengene values for cluster7
meanEigenClust7 <- expressionProfile_data$cluster7
write.csv(meanEigenClust7, file = "meanEigenClust7.csv")
write.csv(C7_Strader_MEs, file = "C7_Strader_MEs.csv")

### save mean eigengene values for cluster8
meanEigenClust8 <- expressionProfile_data$cluster8
write.csv(meanEigenClust8, file = "meanEigenClust8.csv")
write.csv(C8_Strader_MEs, file = "C8_Strader_MEs.csv")

expressionProfile_data$temp
#  [1] "t18"   "t26.8" "t12"   "t35"   "t25"   "t30"   "t26.8" "t18"   "t12"   "t35"   "t25"   "t30"  
# [13] "t18"   "t26.8" "t12"   "t35"   "t25"   "t30"   "t18"   "t26.8" "t12"   "t35"   "t25"   "t30"  

# Convert columns 2 to 11 with clusters to numeric
expressionProfile_data[ , 2:11] <- lapply(expressionProfile_data[ , 2:11], as.numeric)


### Plot mean module eigengene for each cluster
library(tidyr)
library(dplyr)
# Convert to long format for ggplot
expressionProfile_data_long <- expressionProfile_data %>%
  pivot_longer(
    cols = starts_with("cluster"),
    names_to = "cluster",
    values_to = "Eigengene"
  )

# Boxplot with dots, faceted by cluster
ggplot(expressionProfile_data_long, aes(x = temp, y = Eigengene)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(width = 0.2, size = 2, show.legend = FALSE) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_bw() +
  labs(title = "Eigengene Expression by cluster and Temperature",
       x = "Temperature Group",
       y = "cluster Eigengene Value")

# Filter to exclude clusters 5 and 6 (not significant correlations in any of the temperature groups)
expressionProfile_data_long_subset <- expressionProfile_data_long %>%
  filter(cluster %in% c("cluster1", "cluster2", "cluster3", "cluster4", "cluster10", "cluster7", "cluster8", "cluster9"))

# Boxplot with dots, faceted by Module (only Module1 and Module2)
ggplot(expressionProfile_data_long_subset, aes(x = temp, y = Eigengene)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(aes(color = temp), width = 0.2, size = 2, show.legend = FALSE) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_bw() +
  labs(title = "Eigengene Expression by cluster and Temperature",
       x = "Temperature Group",
       y = "Cluster mean Eigengene Value")

# Define ordered temperature levels
expressionProfile_data_long_subset$temp <- factor(expressionProfile_data_long_subset$temp, 
                                       levels = c("t12", "t18", "t25", "t26.8", "t30", "t35"))

# Create a palette of orange shades (light to dark)
orange_palette <- colorRampPalette(c("#ffe0b2", "#ff9800", "#e65100"))(length(levels(expressionProfile_data_long_subset$temp)))

# Plot with custom colors for boxplots and dots, and no grid
ggplot(expressionProfile_data_long_subset, aes(x = temp, y = Eigengene, group = temp)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp), width = 0.2, size = 2, show.legend = FALSE) +
  facet_wrap(~ cluster, scales = "free_y") +
  scale_color_manual(values = orange_palette) +
  scale_x_discrete(labels = c("t12" = "12", "t18" = "18", "t25" = "25", "t26.8" = "26.8", "t30" = "30", "t35" = "35")) +
  theme_bw() +
  labs(x = "Temperature",
       y = "Cluster mean eigengene value") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )



### Use R package segmented
#segmented is designed for regression models where you want to detect one 
# or more breakpoints (change-points) in the relationship between a numeric 
# predictor (here, temperature) and a response (eigengene value).
### selected modules: 1,2,3,4,8,9,10

library(segmented)

# Make sure temp is numeric (remove "tt" if needed)
expressionProfile_data_long_subset$temp_num <- as.numeric(sub("^t", "", expressionProfile_data_long_subset$temp))

# Filter for the module of interest (e.g., cluster1)
df_mod1 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster1")

# Fit linear model to all replicates
fit1 <- lm(Eigengene ~ temp_num, data = df_mod1)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit1 <- segmented(fit1, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit1)
plot(seg_fit1)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 29.799  0.689
# Coefficients of the linear terms:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.022809   0.097976  -0.233    0.819
# temp_num     0.004276   0.004561   0.938    0.361
# U1.temp_num -0.108656   0.016926  -6.419       NA
# Residual standard error: 0.1067 on 18 degrees of freedom
# Multiple R-Squared: 0.7679,  Adjusted R-squared: 0.7292 


#compare models
AIC(fit1, seg_fit1)
#          df        AIC
# fit1      3  -6.318496
# seg_fit1  5 -30.438172
#Lower AIC suggests the segmented model best fits the data


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
df_mod1$temp_num <- as.numeric(sub("^t", "", df_mod1$temp))
df_mod1$temp_num_f <- factor(df_mod1$temp_num)

ggplot(df_mod1, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df1_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df1_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = orange_palette) +
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

# Filter for the module of interest (e.g., cluster2)
df_mod2 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster2")

# Fit linear model to all replicates
fit2 <- lm(Eigengene ~ temp_num, data = df_mod2)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit2 <- segmented(fit2, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit2)
plot(seg_fit2)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 28.845   2.33
# Coefficients of the linear terms:
#              Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -0.293692   0.141350  -2.078   0.0523 .
# temp_num     0.014282   0.006581   2.170   0.0436 *
# U1.temp_num -0.050622   0.024420  -2.073       NA  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1539 on 18 degrees of freedom
# Multiple R-Squared: 0.2829,  Adjusted R-squared: 0.1633 


#compare models
AIC(fit2, seg_fit2)
#          df       AIC
# fit2      3 -12.16394
# seg_fit2  5 -14.31118
#Lower AIC suggests the segmented model best fits the data


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
df_mod2$temp_num <- as.numeric(sub("^t", "", df_mod2$temp))
df_mod2$temp_num_f <- factor(df_mod2$temp_num)

ggplot(df_mod2, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df2_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df2_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = orange_palette) +
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


# Filter for the module of interest (e.g., cluster3)
df_mod3 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster3")

# Fit linear model to all replicates
fit3 <- lm(Eigengene ~ temp_num, data = df_mod3)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit3 <- segmented(fit3, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit3)
plot(seg_fit3)
# Estimated Break-Point(s):
#                Est. St.Err
# psi1.temp_num   30  1.375
# Coefficients of the linear terms:
#              Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -0.187200   0.103360  -1.811   0.0868 .
# temp_num     0.006121   0.004812   1.272   0.2195  
# U1.temp_num  0.056564   0.017856   3.168       NA  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1126 on 18 degrees of freedom
# Multiple R-Squared: 0.6151,  Adjusted R-squared: 0.5509 


#compare models
AIC(fit3, seg_fit3)
#          df       AIC
# fit3      3 -21.57335
# seg_fit3  5 -28.08415
#Lower AIC suggests the segmented model best fits the data


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
# your current bootstrap code is recalculating the segmented fit for each bootstrap sample—that’s the whole point of the bootstrap: 
# to estimate the variability of the fit (and thus the confidence intervals) by refitting the model on resampled data.
# But: The red line you plot (fit in pred_df1_boot) is the median of the bootstrapped fits at each grid point, not the fit from your 
# original seg_fit1 object.

set.seed(123)
grid_temp <- seq(min(df_mod3$temp_num), max(df_mod3$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod3,
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
pred_df3_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod3$temp_num <- as.numeric(sub("^t", "", df_mod3$temp))
df_mod3$temp_num_f <- factor(df_mod3$temp_num)

ggplot(df_mod3, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df3_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df3_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = orange_palette) +
  scale_x_continuous(
    breaks = unique(df_mod3$temp_num),
    labels = gsub("^t", "", unique(df_mod3$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module3: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


# Filter for the module of interest (e.g., cluster4)
df_mod4 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster4")

# Fit linear model to all replicates
fit4 <- lm(Eigengene ~ temp_num, data = df_mod4)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit4 <- segmented(fit4, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit4)
plot(seg_fit4)
# Estimated Break-Point(s):
#                Est. St.Err
# psi1.temp_num 26.8  5.891
# Coefficients of the linear terms:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.296590   0.200592  -1.479    0.157
# temp_num     0.014739   0.009339   1.578    0.132
# U1.temp_num -0.036063   0.034654  -1.041       NA
# Residual standard error: 0.2185 on 18 degrees of freedom
# Multiple R-Squared: 0.1409,  Adjusted R-squared: -0.002335 

#compare models
AIC(fit4, seg_fit4)
#          df        AIC
# fit4      3 -0.2044966
# seg_fit4  5  1.0903499


# Filter for the module of interest (e.g., cluster7)
df_mod7 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster7")

# Fit linear model to all replicates
fit7 <- lm(Eigengene ~ temp_num, data = df_mod7)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit7 <- segmented(fit7, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit7)
plot(seg_fit7)
# Estimated Break-Point(s):
#                Est. St.Err
# psi1.temp_num   30 16.158
# Coefficients of the linear terms:
#              Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.150487   0.169169  -0.890    0.385
# temp_num     0.005976   0.007876   0.759    0.458
# U1.temp_num  0.007881   0.029225   0.270       NA
# Residual standard error: 0.1842 on 18 degrees of freedom
# Multiple R-Squared: 0.09334,  Adjusted R-squared: -0.05777 



# Filter for the module of interest (e.g., cluster8)
df_mod8 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster8")

# Fit linear model to all replicates
fit8 <- lm(Eigengene ~ temp_num, data = df_mod8)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit8 <- segmented(fit8, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit8)
plot(seg_fit8)
# Estimated Break-Point(s):
#                Est. St.Err
# psi1.temp_num   18  1.718

# Coefficients of the linear terms:
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -0.57703    0.24089  -2.395   0.0277 *
# temp_num     0.04626    0.01620   2.856   0.0105 *
# U1.temp_num -0.07408    0.01859  -3.985       NA  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.1273 on 18 degrees of freedom
# Multiple R-Squared: 0.5655,  Adjusted R-squared: 0.4931 

#compare models
AIC(fit8, seg_fit8)
#          df        AIC
# fit8      3 -12.77489
# seg_fit8  5 -22.69031
#Lower AIC suggests the segmented model best fits the data


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
# your current bootstrap code is recalculating the segmented fit for each bootstrap sample—that’s the whole point of the bootstrap: 
# to estimate the variability of the fit (and thus the confidence intervals) by refitting the model on resampled data.
# But: The red line you plot (fit in pred_df1_boot) is the median of the bootstrapped fits at each grid point, not the fit from your 
# original seg_fit1 object.

set.seed(123)
grid_temp <- seq(min(df_mod8$temp_num), max(df_mod8$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod8,
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
pred_df8_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod8$temp_num <- as.numeric(sub("^t", "", df_mod8$temp))
df_mod8$temp_num_f <- factor(df_mod8$temp_num)

ggplot(df_mod8, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df8_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df8_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = orange_palette) +
  scale_x_continuous(
    breaks = unique(df_mod8$temp_num),
    labels = gsub("^t", "", unique(df_mod8$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module8: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


# Filter for the module of interest (e.g., cluster9)
df_mod9 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster9")

# Fit linear model to all replicates
fit9 <- lm(Eigengene ~ temp_num, data = df_mod9)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit9 <- segmented(fit9, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit9)
plot(seg_fit9)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 28.612  1.623
# Coefficients of the linear terms:
#             Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  0.52977    0.14305   3.703  0.00163 **
# temp_num    -0.02531    0.00666  -3.800  0.00131 **
# U1.temp_num  0.07543    0.02471   3.052       NA   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1558 on 18 degrees of freedom
# Multiple R-Squared: 0.5136,  Adjusted R-squared: 0.4325 


#compare models
AIC(fit9, seg_fit9)
#          df        AIC
# fit9      3  -5.665488
# seg_fit9  5 -13.786105
#Lower AIC suggests the segmented model best fits the data


#Run the bootstrap
set.seed(123)
grid_temp <- seq(min(df_mod9$temp_num), max(df_mod9$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod9,
  statistic = function(data, indices) boot_seg(data, indices, grid_temp),
  R = 1000 # Number of bootstrap replicates
)
#colSums(!is.na(boot_preds))

#Calculate confidence intervals
# Each row is a bootstrap, each column is a grid point
boot_preds <- boot_out$t

# Calculate 2.5% and 97.5% quantiles for each grid point
ci_lower <- apply(boot_preds, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(boot_preds, 2, quantile, probs = 0.975, na.rm = TRUE)
fit_median <- apply(boot_preds, 2, median, na.rm = TRUE)

# Prepare for plotting
pred_df9_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod9$temp_num <- as.numeric(sub("^t", "", df_mod9$temp))
df_mod9$temp_num_f <- factor(df_mod9$temp_num)

ggplot(df_mod9, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df9_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df9_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = orange_palette) +
  scale_x_continuous(
    breaks = unique(df_mod9$temp_num),
    labels = gsub("^t", "", unique(df_mod9$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module9: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


# Filter for the module of interest (e.g., cluster10)
df_mod10 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster10")

# Fit linear model to all replicates
fit10 <- lm(Eigengene ~ temp_num, data = df_mod10)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit10 <- segmented(fit10, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit10)
plot(seg_fit10)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 26.801   6.56
# Coefficients of the linear terms:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.533234   0.091388   5.835 1.58e-05 ***
# temp_num    -0.022993   0.004255  -5.404 3.90e-05 ***
# U1.temp_num  0.014752   0.015788   0.934       NA    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.09953 on 18 degrees of freedom
# Multiple R-Squared: 0.7197,  Adjusted R-squared: 0.673 

#no break-point but The coefficient for temp_num is highly significant (t = -5.404, p = 3.9e-05), indicating a strong linear relationship between temperature and eigengene value across the full range.
#compare models
AIC(fit10, seg_fit10)
#           df       AIC
# fit10      3 -35.29462
# seg_fit10  5 -33.50060
#Lower AIC suggests the linera model best fits the data


#Run the bootstrap
set.seed(123)
grid_temp <- seq(min(df_mod9$temp_num), max(df_mod10$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod9,
  statistic = function(data, indices) boot_seg(data, indices, grid_temp),
  R = 1000 # Number of bootstrap replicates
)
#colSums(!is.na(boot_preds))

#Calculate confidence intervals
# Each row is a bootstrap, each column is a grid point
boot_preds <- boot_out$t

# Calculate 2.5% and 97.5% quantiles for each grid point
ci_lower <- apply(boot_preds, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(boot_preds, 2, quantile, probs = 0.975, na.rm = TRUE)
fit_median <- apply(boot_preds, 2, median, na.rm = TRUE)

# Prepare for plotting
pred_df10_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod10$temp_num <- as.numeric(sub("^t", "", df_mod10$temp))
df_mod10$temp_num_f <- factor(df_mod10$temp_num)

ggplot(df_mod10, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df10_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df10_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = orange_palette) +
  scale_x_continuous(
    breaks = unique(df_mod10$temp_num),
    labels = gsub("^t", "", unique(df_mod10$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module10: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


############## Find Hub Genes and Network Analysis

### export network of clusters with significant breakpoints and trends (clusters 1, 2, 3, 8, 9) - see analysis above

# Subset samples for each group
control_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "t26.8"]
treat30_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "t30"]
treat35_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "t35"]

#### Export the network of a specific cluster - cluster1
cluster1_modules <- c("darkred", "turquoise")  # replace with your actual module colors for cluster1
inCluster <- moduleColors %in% cluster1_modules #Select genes belonging to these modules
cluster1_genes <- probes[inCluster] ##8814

# Subset expression data for each group
datExpr_control_cluster1 <- datExpr[control_samples, cluster1_genes ]
datExpr_treat30_cluster1 <- datExpr[treat30_samples, cluster1_genes ]
datExpr_treat35_cluster1 <- datExpr[treat35_samples, cluster1_genes ]

# Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# kWithin: sum of connection strengths with other module genes
#add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster1)
# names(kWithin_control) <- colnames(datExpr_control_cluster1)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster1)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster1)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster1)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster1)

# Use sum of absolute correlations as a proxy for connectivity
# For control group
cor_mat_control <- cor(datExpr_control_cluster1, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_control) <- 0
hub_score_control <- rowSums(abs(cor_mat_control), na.rm = TRUE)
names(hub_score_control) <- colnames(datExpr_control_cluster1)
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
top10pct_hubs_control_1 <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]

# For treat30 group
cor_mat_treat30 <- cor(datExpr_treat30_cluster1, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat30) <- 0
hub_score_treat30 <- rowSums(abs(cor_mat_treat30), na.rm = TRUE)
names(hub_score_treat30) <- colnames(datExpr_treat30_cluster1)
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top10pct_hubs_treat30_1 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]

# For treat35 group 
cor_mat_treat35 <- cor(datExpr_treat35_cluster1, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat35) <- 0
hub_score_treat35 <- rowSums(abs(cor_mat_treat35), na.rm = TRUE)
names(hub_score_treat35) <- colnames(datExpr_treat35_cluster1)
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top10pct_hubs_treat35_1 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_cluster1.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster1.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster1.csv", row.names = FALSE)

#### Identify the top Hub Gene per cluster and temp

# For control group
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
# Find the top hub gene among the top 10%
top_hub_control <- top10pct_hubs_control[which.max(hub_score_control[top10pct_hubs_control])]

# For treat30 group
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top_hub_treat30 <- top10pct_hubs_treat30[which.max(hub_score_treat30[top10pct_hubs_treat30])]

# For treat35 group
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top_hub_treat35 <- top10pct_hubs_treat35[which.max(hub_score_treat35[top10pct_hubs_treat35])]

cat("Top hub gene (control, top 10%):", top_hub_control, "\n")
#Porites_compressa_HIv1___TS.g30032.t1 
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Porites_compressa_HIv1___RNAseq.33470_t 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Porites_compressa_HIv1___RNAseq.g5672.t1 


#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 1
all_genes <- cluster1_genes

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
  mutate(Cluster = "Cluster 1 Genes")

# Create a label for each interaction pattern
membership_summary$pattern <- interaction(membership_summary$Control, membership_summary$Treat30, membership_summary$Treat35)
pattern_counts <- setNames(membership_summary$Freq, membership_summary$pattern)
pattern_labels <- paste0(names(pattern_counts), " (n=", pattern_counts, ")")
names(pattern_labels) <- names(pattern_counts)

# Plot 
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
                   labels = c("Cluster" = "Cluster 1 Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Cluster 1: Overlap of Top 10% Hub Genes Across Conditions",
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


# Save membership summary data for UpSet plot analysis
cluster1_membership_summary <- membership_summary %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster1") %>%
  select(Species, Cluster, Control, Treat30, Treat35, Freq, Percent, pattern)

write.csv(cluster1_membership_summary, 
          file = "Pcom_Cluster1_membership_summary_for_upset.csv", 
          row.names = FALSE)

# Also save the individual gene membership data
cluster1_gene_membership <- membership %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster1") %>%
  select(Species, Cluster, Gene, Control, Treat30, Treat35)

write.csv(cluster1_gene_membership, 
          file = "Pcom_Cluster1_gene_membership_for_upset.csv", 
          row.names = FALSE)
     




# Assign colors (same as alluvial)
bar_colors <- setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels))

bar1 <- ggplot(membership_summary, aes(x = "Cluster 1", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Cluster 1: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Save hub genes data for UpSet plot analysis
# Create a comprehensive data frame for cluster 1
cluster1_hub_data <- data.frame(
  Gene = c(top10pct_hubs_control, top10pct_hubs_treat30, top10pct_hubs_treat35),
  Species = "Pcom",
  Cluster = "Cluster1",
  Temperature = c(rep("Control", length(top10pct_hubs_control)),
                  rep("Treat30", length(top10pct_hubs_treat30)),
                  rep("Treat35", length(top10pct_hubs_treat35))),
  Set_ID = c(paste0("Pcom_Cluster1_Control"),
             paste0("Pcom_Cluster1_Treat30"), 
             paste0("Pcom_Cluster1_Treat35")),
  stringsAsFactors = FALSE
)

# Remove duplicates (genes that appear in multiple conditions)
cluster1_hub_data <- cluster1_hub_data[!duplicated(cluster1_hub_data[c("Gene", "Temperature")]), ]

# Save to CSV for later UpSet analysis
write.csv(cluster1_hub_data, 
          file = "Pcom_Cluster1_hub_genes_for_upset.csv", 
          row.names = FALSE)







#### Export the network of a specific cluster - cluster2
cluster2_modules <- c("black", "pink", "lightgreen", "mediumpurple3", "darkturquoise", "yellow4")  # replace with your actual module colors for cluster2
inCluster <- moduleColors %in% cluster2_modules
cluster2_genes <- probes[inCluster] ##3561

# Subset expression data for each group
datExpr_control_cluster2 <- datExpr[control_samples, cluster2_genes ]
datExpr_treat30_cluster2 <- datExpr[treat30_samples, cluster2_genes ]
datExpr_treat35_cluster2 <- datExpr[treat35_samples, cluster2_genes ]

# Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# kWithin: sum of connection strengths with other module genes
#add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster1)
# names(kWithin_control) <- colnames(datExpr_control_cluster1)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster1)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster1)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster1)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster1)

# Use sum of absolute correlations as a proxy for connectivity
# For control group
cor_mat_control <- cor(datExpr_control_cluster2, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_control) <- 0
hub_score_control <- rowSums(abs(cor_mat_control), na.rm = TRUE)
names(hub_score_control) <- colnames(datExpr_control_cluster2)
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
top10pct_hubs_control_2 <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]

# For treat30 group
cor_mat_treat30 <- cor(datExpr_treat30_cluster2, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat30) <- 0
hub_score_treat30 <- rowSums(abs(cor_mat_treat30), na.rm = TRUE)
names(hub_score_treat30) <- colnames(datExpr_treat30_cluster2)
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top10pct_hubs_treat30_2 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]

# For treat35 group 
cor_mat_treat35 <- cor(datExpr_treat35_cluster2, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat35) <- 0
hub_score_treat35 <- rowSums(abs(cor_mat_treat35), na.rm = TRUE)
names(hub_score_treat35) <- colnames(datExpr_treat35_cluster2)
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top10pct_hubs_treat35_2 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_cluster2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster2.csv", row.names = FALSE)


#### Identify the top Hub Gene per cluster and temp

# For control group
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
# Find the top hub gene among the top 10%
top_hub_control <- top10pct_hubs_control[which.max(hub_score_control[top10pct_hubs_control])]

# For treat30 group
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top_hub_treat30 <- top10pct_hubs_treat30[which.max(hub_score_treat30[top10pct_hubs_treat30])]

# For treat35 group
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top_hub_treat35 <- top10pct_hubs_treat35[which.max(hub_score_treat35[top10pct_hubs_treat35])]

cat("Top hub gene (control, top 10%):", top_hub_control, "\n")
#Porites_compressa_HIv1___RNAseq.g9915.t1  
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Porites_compressa_HIv1___TS.g25224.t1 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Porites_compressa_HIv1___RNAseq.g40300.t1


#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 2
all_genes <- cluster2_genes

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
  mutate(Cluster = "Cluster 2 Genes")

# Create a label for each interaction pattern
membership_summary$pattern <- interaction(membership_summary$Control, membership_summary$Treat30, membership_summary$Treat35)
pattern_counts <- setNames(membership_summary$Freq, membership_summary$pattern)
pattern_labels <- paste0(names(pattern_counts), " (n=", pattern_counts, ")")
names(pattern_labels) <- names(pattern_counts)

# Plot 
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
                   labels = c("Cluster" = "Cluster 2 Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Cluster 21: Overlap of Top 10% Hub Genes Across Conditions",
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


# Save membership summary data for UpSet plot analysis
cluster2_membership_summary <- membership_summary %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster2") %>%
  select(Species, Cluster, Control, Treat30, Treat35, Freq, Percent, pattern)

write.csv(cluster2_membership_summary, 
          file = "Pcom_Cluster2_membership_summary_for_upset.csv", 
          row.names = FALSE)

# Also save the individual gene membership data
cluster2_gene_membership <- membership %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster2") %>%
  select(Species, Cluster, Gene, Control, Treat30, Treat35)

write.csv(cluster2_gene_membership, 
          file = "Pcom_Cluster2_gene_membership_for_upset.csv", 
          row.names = FALSE)



# Assign colors (same as alluvial)
bar_colors <- setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels))

bar2 <- ggplot(membership_summary, aes(x = "Cluster 2", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Cluster 2: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Save hub genes data for UpSet plot analysis
# Create a comprehensive data frame for cluster 2
cluster2_hub_data <- data.frame(
  Gene = c(top10pct_hubs_control, top10pct_hubs_treat30, top10pct_hubs_treat35),
  Species = "Pcom",
  Cluster = "Cluster2",
  Temperature = c(rep("Control", length(top10pct_hubs_control)),
                  rep("Treat30", length(top10pct_hubs_treat30)),
                  rep("Treat35", length(top10pct_hubs_treat35))),
  Set_ID = c(paste0("Pcom_Cluster2_Control"),
             paste0("Pcom_Cluster2_Treat30"), 
             paste0("Pcom_Cluster2_Treat35")),
  stringsAsFactors = FALSE
)

# Remove duplicates (genes that appear in multiple conditions)
cluster2_hub_data <- cluster2_hub_data[!duplicated(cluster2_hub_data[c("Gene", "Temperature")]), ]

# Save to CSV for later UpSet analysis
write.csv(cluster2_hub_data, 
          file = "Pcom_Cluster2_hub_genes_for_upset.csv", 
          row.names = FALSE)





#### Export the network of a specific cluster - cluster3
cluster3_modules <- c("skyblue2", "floralwhite", "coral1", "cyan", "blue", "grey60", "paleturquoise")  # replace with your actual module colors for cluster3
inCluster <- moduleColors %in% cluster3_modules
cluster3_genes <- probes[inCluster] ##6935

# Subset expression data for each group
datExpr_control_cluster3 <- datExpr[control_samples, cluster3_genes ]
datExpr_treat30_cluster3 <- datExpr[treat30_samples, cluster3_genes ]
datExpr_treat35_cluster3 <- datExpr[treat35_samples, cluster3_genes ]

# Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# kWithin: sum of connection strengths with other module genes
#add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster3)
# names(kWithin_control) <- colnames(datExpr_control_cluster3)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster3)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster3)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster3)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster3)

# Use sum of absolute correlations as a proxy for connectivity
# For control group
cor_mat_control <- cor(datExpr_control_cluster3, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_control) <- 0
hub_score_control <- rowSums(abs(cor_mat_control), na.rm = TRUE)
names(hub_score_control) <- colnames(datExpr_control_cluster3)
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
top10pct_hubs_control_3 <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]

# For treat30 group
cor_mat_treat30 <- cor(datExpr_treat30_cluster3, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat30) <- 0
hub_score_treat30 <- rowSums(abs(cor_mat_treat30), na.rm = TRUE)
names(hub_score_treat30) <- colnames(datExpr_treat30_cluster3)
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top10pct_hubs_treat30_3 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]

# For treat35 group 
cor_mat_treat35 <- cor(datExpr_treat35_cluster3, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat35) <- 0
hub_score_treat35 <- rowSums(abs(cor_mat_treat35), na.rm = TRUE)
names(hub_score_treat35) <- colnames(datExpr_treat35_cluster3)
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top10pct_hubs_treat35_3 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_cluster3.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster3.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster3.csv", row.names = FALSE)

#### Identify the top Hub Gene per cluster and temp

# For control group
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
# Find the top hub gene among the top 10%
top_hub_control <- top10pct_hubs_control[which.max(hub_score_control[top10pct_hubs_control])]

# For treat30 group
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top_hub_treat30 <- top10pct_hubs_treat30[which.max(hub_score_treat30[top10pct_hubs_treat30])]

# For treat35 group
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top_hub_treat35 <- top10pct_hubs_treat35[which.max(hub_score_treat35[top10pct_hubs_treat35])]

cat("Top hub gene (control, top 10%):", top_hub_control, "\n")
#Porites_compressa_HIv1___RNAseq.g19695.t1 
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Porites_compressa_HIv1___RNAseq.g3874.t1 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Porites_compressa_HIv1___TS.g14122.t1 


#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 3
all_genes <- cluster3_genes

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
  mutate(Cluster = "Cluster 3 Genes")

# Create a label for each interaction pattern
membership_summary$pattern <- interaction(membership_summary$Control, membership_summary$Treat30, membership_summary$Treat35)
pattern_counts <- setNames(membership_summary$Freq, membership_summary$pattern)
pattern_labels <- paste0(names(pattern_counts), " (n=", pattern_counts, ")")
names(pattern_labels) <- names(pattern_counts)

# Plot 
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
                   labels = c("Cluster" = "Cluster 3 Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Cluster 3: Overlap of Top 10% Hub Genes Across Conditions",
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


# Save membership summary data for UpSet plot analysis
cluster3_membership_summary <- membership_summary %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster3") %>%
  select(Species, Cluster, Control, Treat30, Treat35, Freq, Percent, pattern)

write.csv(cluster3_membership_summary, 
          file = "Pcom_Cluster3_membership_summary_for_upset.csv", 
          row.names = FALSE)

# Also save the individual gene membership data
cluster3_gene_membership <- membership %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster3") %>%
  select(Species, Cluster, Gene, Control, Treat30, Treat35)

write.csv(cluster3_gene_membership, 
          file = "Pcom_Cluster3_gene_membership_for_upset.csv", 
          row.names = FALSE)




# Assign colors (same as alluvial)
bar_colors <- setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels))

bar3 <- ggplot(membership_summary, aes(x = "Cluster 3", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Cluster 3: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Save hub genes data for UpSet plot analysis
# Create a comprehensive data frame for cluster 3
cluster3_hub_data <- data.frame(
  Gene = c(top10pct_hubs_control, top10pct_hubs_treat30, top10pct_hubs_treat35),
  Species = "Pcom",
  Cluster = "Cluster3",
  Temperature = c(rep("Control", length(top10pct_hubs_control)),
                  rep("Treat30", length(top10pct_hubs_treat30)),
                  rep("Treat35", length(top10pct_hubs_treat35))),
  Set_ID = c(paste0("Pcom_Cluster3_Control"),
             paste0("Pcom_Cluster3_Treat30"), 
             paste0("Pcom_Cluster3_Treat35")),
  stringsAsFactors = FALSE
)

# Remove duplicates (genes that appear in multiple conditions)
cluster3_hub_data <- cluster3_hub_data[!duplicated(cluster3_hub_data[c("Gene", "Temperature")]), ]

# Save to CSV for later UpSet analysis
write.csv(cluster3_hub_data, 
          file = "Pcom_Cluster3_hub_genes_for_upset.csv", 
          row.names = FALSE)





#### Export the network of a specific cluster - cluster9
cluster9_modules <- c("brown4", "darkgreen")  # replace with your actual module colors for cluster9
inCluster <- moduleColors %in% cluster9_modules
cluster9_genes <- probes[inCluster] ##536

# Subset expression data for each group
datExpr_control_cluster9 <- datExpr[control_samples, cluster9_genes ]
datExpr_treat30_cluster9 <- datExpr[treat30_samples, cluster9_genes ]
datExpr_treat35_cluster9 <- datExpr[treat35_samples, cluster9_genes ]

# Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# kWithin: sum of connection strengths with other module genes
#add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster9)
# names(kWithin_control) <- colnames(datExpr_control_cluster9)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster9)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster9)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster9)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster9)

# Use sum of absolute correlations as a proxy for connectivity
# For control group
cor_mat_control <- cor(datExpr_control_cluster9, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_control) <- 0
hub_score_control <- rowSums(abs(cor_mat_control), na.rm = TRUE)
names(hub_score_control) <- colnames(datExpr_control_cluster9)
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
top10pct_hubs_control_9 <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]

# For treat30 group
cor_mat_treat30 <- cor(datExpr_treat30_cluster9, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat30) <- 0
hub_score_treat30 <- rowSums(abs(cor_mat_treat30), na.rm = TRUE)
names(hub_score_treat30) <- colnames(datExpr_treat30_cluster9)
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top10pct_hubs_treat30_9 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]

# For treat35 group 
cor_mat_treat35 <- cor(datExpr_treat35_cluster9, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat35) <- 0
hub_score_treat35 <- rowSums(abs(cor_mat_treat35), na.rm = TRUE)
names(hub_score_treat35) <- colnames(datExpr_treat35_cluster9)
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top10pct_hubs_treat35_9 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_cluster9.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster9.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster9.csv", row.names = FALSE)

#### Identify the top Hub Gene per cluster and temp

# For control group
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
# Find the top hub gene among the top 10%
top_hub_control <- top10pct_hubs_control[which.max(hub_score_control[top10pct_hubs_control])]

# For treat30 group
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top_hub_treat30 <- top10pct_hubs_treat30[which.max(hub_score_treat30[top10pct_hubs_treat30])]

# For treat35 group
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top_hub_treat35 <- top10pct_hubs_treat35[which.max(hub_score_treat35[top10pct_hubs_treat35])]

cat("Top hub gene (control, top 10%):", top_hub_control, "\n")
#Porites_compressa_HIv1___RNAseq.g36475.t1 
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Porites_compressa_HIv1___RNAseq.38615_t 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Porites_compressa_HIv1___RNAseq.g13898.t1 


#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 9
all_genes <- cluster9_genes

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
  mutate(Cluster = "Cluster 9 Genes")

# Create a label for each interaction pattern
membership_summary$pattern <- interaction(membership_summary$Control, membership_summary$Treat30, membership_summary$Treat35)
pattern_counts <- setNames(membership_summary$Freq, membership_summary$pattern)
pattern_labels <- paste0(names(pattern_counts), " (n=", pattern_counts, ")")
names(pattern_labels) <- names(pattern_counts)

# Plot 
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
                   labels = c("Cluster" = "Cluster 9 Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Cluster 9: Overlap of Top 10% Hub Genes Across Conditions",
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


# Save membership summary data for UpSet plot analysis
cluster9_membership_summary <- membership_summary %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster9") %>%
  select(Species, Cluster, Control, Treat30, Treat35, Freq, Percent, pattern)

write.csv(cluster9_membership_summary, 
          file = "Pcom_Cluster9_membership_summary_for_upset.csv", 
          row.names = FALSE)

# Also save the individual gene membership data
cluster9_gene_membership <- membership %>%
  mutate(Species = "Pcom",
         Cluster = "Cluster9") %>%
  select(Species, Cluster, Gene, Control, Treat30, Treat35)

write.csv(cluster9_gene_membership, 
          file = "Pcom_Cluster9_gene_membership_for_upset.csv", 
          row.names = FALSE)



# Assign colors (same as alluvial)
bar_colors <- setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels))

bar9 <- ggplot(membership_summary, aes(x = "Cluster 9", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Cluster 9: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


# Save hub genes data for UpSet plot analysis
# Create a comprehensive data frame for cluster 9
cluster9_hub_data <- data.frame(
  Gene = c(top10pct_hubs_control, top10pct_hubs_treat30, top10pct_hubs_treat35),
  Species = "Pcom",
  Cluster = "Cluster9",
  Temperature = c(rep("Control", length(top10pct_hubs_control)),
                  rep("Treat30", length(top10pct_hubs_treat30)),
                  rep("Treat35", length(top10pct_hubs_treat35))),
  Set_ID = c(paste0("Pcom_Cluster9_Control"),
             paste0("Pcom_Cluster9_Treat30"), 
             paste0("Pcom_Cluster9_Treat35")),
  stringsAsFactors = FALSE
)

# Remove duplicates (genes that appear in multiple conditions)
cluster9_hub_data <- cluster9_hub_data[!duplicated(cluster9_hub_data[c("Gene", "Temperature")]), ]

# Save to CSV for later UpSet analysis
write.csv(cluster9_hub_data, 
          file = "Pcom_Cluster9_hub_genes_for_upset.csv", 
          row.names = FALSE)


##### combine all bar plots into one figure
library(patchwork)
bar1 + bar2 + bar3 + bar9


######## Make a global alluvial plot showing the overlap of all top 10% hub genes from all clusters across the three temperature groups

# Combine all top 10% hub genes from all clusters and all conditions
all_top_hubs <- unique(c(
  top10pct_hubs_control_1, top10pct_hubs_treat30_1, top10pct_hubs_treat35_1,
  top10pct_hubs_control_2, top10pct_hubs_treat30_2, top10pct_hubs_treat35_2,
  top10pct_hubs_control_3, top10pct_hubs_treat30_3, top10pct_hubs_treat35_3,
  top10pct_hubs_control_9, top10pct_hubs_treat30_9, top10pct_hubs_treat35_9
))

# For each gene, check if it is a top hub in each group (any cluster)
membership_global <- data.frame(
  Gene = all_top_hubs,
  Control = as.integer(all_top_hubs %in% c(
    top10pct_hubs_control_1, top10pct_hubs_control_2, top10pct_hubs_control_3, top10pct_hubs_control_9)),
  Treat30 = as.integer(all_top_hubs %in% c(
    top10pct_hubs_treat30_1, top10pct_hubs_treat30_2, top10pct_hubs_treat30_3, top10pct_hubs_treat30_9)),
  Treat35 = as.integer(all_top_hubs %in% c(
    top10pct_hubs_treat35_1, top10pct_hubs_treat35_2, top10pct_hubs_treat35_3, top10pct_hubs_treat35_9))
)


# Filter to keep only genes that are a hub in at least one condition
membership_global <- membership_global %>%
  filter(Control == 1 | Treat30 == 1 | Treat35 == 1)

library(dplyr)
library(ggalluvial)

membership_summary_global <- membership_global %>%
  group_by(Control, Treat30, Treat35) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  mutate(Global = "All Top 10% Hubs")

membership_summary_global$pattern <- interaction(
  membership_summary_global$Control,
  membership_summary_global$Treat30,
  membership_summary_global$Treat35
)
pattern_counts <- setNames(membership_summary_global$Freq, membership_summary_global$pattern)
pattern_labels <- paste0(names(pattern_counts), " (n=", pattern_counts, ")")
names(pattern_labels) <- names(pattern_counts)

ggplot(membership_summary_global,
       aes(axis1 = Global,
           axis2 = factor(Control, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           axis3 = factor(Treat30, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           axis4 = factor(Treat35, levels = c(0, 1), labels = c("Not Hub", "Hub")),
           y = Freq)) +
  geom_alluvium(aes(fill = pattern), width = 1/8, alpha = 0.8) +
  geom_stratum(width = 1/6, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("Global", "Control", "Treat30", "Treat35"),
                   labels = c("Global" = "All Top 10% Hubs",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_brewer(type = "qual", palette = "Set2", name = "Pattern", labels = pattern_labels) +
  theme_minimal() +
  labs(title = "Global Overlap of Top 10% Hub Genes Across Clusters 1,2,3,9 and Conditions",
       y = "Number of Genes", x = "") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())






########### Quantify network stability or rewiring across conditions:
# Edge Overlap / Jaccard Index
# Definition:
# The Jaccard index measures the similarity between two sets (here, edge sets).
# It is defined as:
# Jaccard = (Number of shared edges) / (Total number of unique edges in both networks)
# To calculate the Jaccard Index and percent edge change for the entire set of genes in the selected clusters (not just the top10% hubs), 
# use all genes in the module for your adjacency/correlation matrices and edge comparisons.

# Use all genes in the selected clusters

cluster1_modules <- c("darkred", "turquoise")  
inCluster <- moduleColors %in% cluster1_modules #Select genes belonging to these modules
cluster1_genes <- probes[inCluster] ##8814

cluster2_modules <- c("black", "pink", "lightgreen", "mediumpurple3", "darkturquoise", "yellow4")  
inCluster <- moduleColors %in% cluster2_modules
cluster2_genes <- probes[inCluster] ##3561

cluster3_modules <- c("skyblue2", "floralwhite", "coral1", "cyan", "blue", "grey60", "paleturquoise")  
inCluster <- moduleColors %in% cluster3_modules
cluster3_genes <- probes[inCluster] ##6935

cluster9_modules <- c("brown4", "darkgreen")  
inCluster <- moduleColors %in% cluster9_modules
cluster9_genes <- probes[inCluster] ##536

cluster10_modules <- c("salmon4", "plum", "orangered4", "red", "darkorange", "lightpink4")  
inCluster <- moduleColors %in% cluster10_modules
cluster10_genes <- probes[inCluster] ##1885

## Subset expression data
datExpr_control_all <- datExpr[control_samples, cluster1_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster1_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster1_genes]

#or
datExpr_control_all <- datExpr[control_samples, cluster2_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster2_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster2_genes]

#or
datExpr_control_all <- datExpr[control_samples, cluster3_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster3_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster3_genes]

#or
datExpr_control_all <- datExpr[control_samples, cluster9_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster9_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster9_genes]

#or
datExpr_control_all <- datExpr[control_samples, cluster10_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster10_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster10_genes]

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
cat("Jaccard index (all cluster genes):", jaccard_index, "\n")
#Jaccard index ........................ ##30
# Interpretation:
# 1 = identical networks
# 0 = no shared edges

shared_edges <- intersect(edges_control_set, edges_treat35_set)
all_edges <- union(edges_control_set, edges_treat35_set)
jaccard_index <- length(shared_edges) / length(all_edges)
cat("Jaccard index (all cluster genes):", jaccard_index, "\n")
#Jaccard index ........................ ##35





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
adj1_all <- (adj_control_all > edge_threshold)    #### select first from above the cluster of interest before running this (line 1680)
adj2_all <- (adj_treat30_all > edge_threshold)
diag(adj1_all) <- 0
diag(adj2_all) <- 0

# adj1_all <- (adj_control_all > edge_threshold)    #### select first from above the cluster of interest before running this (line 1680)
# adj2_all <- (adj_treat35_all > edge_threshold)
# diag(adj1_all) <- 0
# diag(adj2_all) <- 0

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

# 2. Permutation Test
# a. Combine all samples from both groups.
# b. Randomly assign samples to two groups (same sizes as original), recalculate density for each, and compute the difference.
# c. Repeat many times (e.g., 1000 permutations) to build a null distribution.
# d. Calculate the p-value as the proportion of permuted differences as or more extreme than observed.

library(parallel)
set.seed(123)
n_perm <- 200
all_samples <- c(control_samples, treat30_samples)
n_control <- length(control_samples)
n_treat30 <- length(treat30_samples)

perm_diffs <- mclapply(1:n_perm, function(i) {
  permuted <- sample(all_samples)
  perm_control <- permuted[1:n_control]
  perm_treat30 <- permuted[(n_control+1):(n_control+n_treat30)]
  
  # Subset expression data
  datExpr_control_perm <- datExpr[perm_control, all_blue_genes]
  datExpr_treat30_perm <- datExpr[perm_treat30, all_blue_genes]
  
  # Build adjacency matrices
  adj_control_perm <- abs(cor(datExpr_control_perm, method = "pearson"))
  adj_treat30_perm <- abs(cor(datExpr_treat30_perm, method = "pearson"))
  
  # Threshold and remove diagonal
  adj1_perm <- (adj_control_perm > edge_threshold)
  adj2_perm <- (adj_treat30_perm > edge_threshold)
  diag(adj1_perm) <- 0
  diag(adj2_perm) <- 0

  # Force symmetry and handle NAs
  adj1_perm <- adj1_perm | t(adj1_perm)
  adj2_perm <- adj2_perm | t(adj2_perm)
  adj1_perm[is.na(adj1_perm)] <- FALSE
  adj2_perm[is.na(adj2_perm)] <- FALSE

  # Build igraph objects
  g_control_perm <- graph_from_adjacency_matrix(adj1_perm, mode = "undirected", diag = FALSE)
  g_treat30_perm <- graph_from_adjacency_matrix(adj2_perm, mode = "undirected", diag = FALSE)
  
  # Calculate densities
  density_control_perm <- edge_density(g_control_perm)
  density_treat30_perm <- edge_density(g_treat30_perm)
  
  # Return difference
  density_treat30_perm - density_control_perm
}, mc.cores = 4) # Set to the number of CPU cores you want to use

perm_diffs <- unlist(perm_diffs)

# Calculate p-value (two-sided)
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
cat("Observed density difference:", obs_diff, "\n")
cat("Permutation p-value:", p_value, "\n")
# 3. Interpretation
# If p < 0.05: The difference in network density between conditions is statistically significant.
# If p ≥ 0.05: The observed difference could be due to chance.



########### Jaccard for all tempratures compared to control

#Define clusters and temperature groups
clusters <- list(
  cluster1 = cluster1_genes,
  cluster2 = cluster2_genes,
  cluster3 = cluster3_genes,
  cluster9 = cluster9_genes
)

temp_groups <- c("t12", "t18", "t25", "t30", "t35")
control_group <- "t26.8"

#Loop over clusters and temperatures
# Prepare results table
jaccard_results <- data.frame(
  Cluster = character(),
  Temp = character(),
  Jaccard = numeric(),
  Shared_Edges = integer(),
  All_Edges = integer(),
  stringsAsFactors = FALSE
)

edge_threshold <- 0.6 #0.3, 0.4, 0.5, 0.6
# The edge threshold (e.g., edge_threshold <- 0.5) determines which gene-gene correlations are considered "edges" in the network.
# A higher threshold means only very strong correlations are counted as edges, resulting in sparser networks.
# A lower threshold includes more edges, possibly capturing weaker but biologically relevant connections.
# The Jaccard index is sensitive to the number and identity of edges: changing the threshold can alter the overlap between networks at different temperatures.

for (clust_name in names(clusters)) {
  genes <- clusters[[clust_name]]
  
  # Subset expression data for control
  datExpr_control <- datExpr[which(treatmentinfo$temp == control_group), genes]
  adj_control <- abs(cor(datExpr_control, method = "pearson"))
  adj_control_bin <- (adj_control > edge_threshold)
  diag(adj_control_bin) <- 0
  adj_control_bin <- adj_control_bin | t(adj_control_bin)
  adj_control_bin[is.na(adj_control_bin)] <- FALSE
  edges_control <- which(adj_control_bin & upper.tri(adj_control_bin), arr.ind = TRUE)
  edges_control_set <- paste(edges_control[,1], edges_control[,2], sep = "-")
  
  for (temp in temp_groups) {
    datExpr_temp <- datExpr[which(treatmentinfo$temp == temp), genes]
    adj_temp <- abs(cor(datExpr_temp, method = "pearson"))
    adj_temp_bin <- (adj_temp > edge_threshold)
    diag(adj_temp_bin) <- 0
    adj_temp_bin <- adj_temp_bin | t(adj_temp_bin)
    adj_temp_bin[is.na(adj_temp_bin)] <- FALSE
    edges_temp <- which(adj_temp_bin & upper.tri(adj_temp_bin), arr.ind = TRUE)
    edges_temp_set <- paste(edges_temp[,1], edges_temp[,2], sep = "-")
    
    shared_edges <- intersect(edges_control_set, edges_temp_set)
    all_edges <- union(edges_control_set, edges_temp_set)
    jaccard_index <- if(length(all_edges) > 0) length(shared_edges) / length(all_edges) else NA
    
    # Add to results table
    jaccard_results <- rbind(jaccard_results, data.frame(
      Cluster = clust_name,
      Temp = temp,
      Jaccard = jaccard_index,
      Shared_Edges = length(shared_edges),
      All_Edges = length(all_edges)
    ))
  }
}

# Save results
write.csv(jaccard_results, "jaccard_index_all_clusters_temps_t06_Pcom.csv", row.names = FALSE)


































######## Extract the genes from clusters for GO enrichment analysis with ViSEAGO

# For cluster 1
# cluster1_genes: vector of gene names in cluster 1
# datExpr: samples x all genes (genes as columns)

# Subset datExpr to cluster 1 genes
cluster1_counts <- datExpr[, cluster1_genes, drop = FALSE]

# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster1_counts_t <- as.data.frame(t(cluster1_counts))
cluster1_counts_t$gene_id <- rownames(cluster1_counts_t)

# Move gene_id to first column
cluster1_counts_t <- cluster1_counts_t[, c(ncol(cluster1_counts_t), 1:(ncol(cluster1_counts_t)-1))]

# Save as CSV
write.csv(cluster1_counts_t, file = "cluster1_gene_counts.csv", row.names = FALSE)


# For cluster 2
# cluster1_genes: vector of gene names in cluster 2
# datExpr: samples x all genes (genes as columns)

# Subset datExpr to cluster 2 genes
cluster2_counts <- datExpr[, cluster2_genes, drop = FALSE]

# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster2_counts_t <- as.data.frame(t(cluster2_counts))
cluster2_counts_t$gene_id <- rownames(cluster2_counts_t)

# Move gene_id to first column
cluster2_counts_t <- cluster2_counts_t[, c(ncol(cluster2_counts_t), 1:(ncol(cluster2_counts_t)-1))]

# Save as CSV
write.csv(cluster2_counts_t, file = "cluster2_gene_counts.csv", row.names = FALSE)


# For cluster 3
# cluster1_genes: vector of gene names in cluster 3
# datExpr: samples x all genes (genes as columns)
# Subset datExpr to cluster 3 genes
cluster3_counts <- datExpr[, cluster3_genes, drop = FALSE]
# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster3_counts_t <- as.data.frame(t(cluster3_counts))
cluster3_counts_t$gene_id <- rownames(cluster3_counts_t)
# Move gene_id to first column
cluster3_counts_t <- cluster3_counts_t[, c(ncol(cluster3_counts_t), 1:(ncol(cluster3_counts_t)-1))]
# Save as CSV
write.csv(cluster3_counts_t, file = "cluster3_gene_counts.csv", row.names = FALSE)

# For cluster 9
# cluster1_genes: vector of gene names in cluster 9
# datExpr: samples x all genes (genes as columns)
# Subset datExpr to cluster 9 genes
cluster9_counts <- datExpr[, cluster9_genes, drop = FALSE]
# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster9_counts_t <- as.data.frame(t(cluster9_counts))
cluster9_counts_t$gene_id <- rownames(cluster9_counts_t)
# Move gene_id to first column
cluster9_counts_t <- cluster9_counts_t[, c(ncol(cluster9_counts_t), 1:(ncol(cluster9_counts_t)-1))]
# Save as CSV
write.csv(cluster9_counts_t, file = "cluster9_gene_counts.csv", row.names = FALSE)

# For cluster 10
# cluster1_genes: vector of gene names in cluster 10
# datExpr: samples x all genes (genes as columns)
# Subset datExpr to cluster 10 genes
cluster10_counts <- datExpr[, cluster10_genes, drop = FALSE]
# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster10_counts_t <- as.data.frame(t(cluster10_counts))
cluster10_counts_t$gene_id <- rownames(cluster10_counts_t)
# Move gene_id to first column
cluster10_counts_t <- cluster10_counts_t[, c(ncol(cluster10_counts_t), 1:(ncol(cluster10_counts_t)-1))]
# Save as CSV
write.csv(cluster10_counts_t, file = "cluster10_gene_counts.csv", row.names = FALSE)
