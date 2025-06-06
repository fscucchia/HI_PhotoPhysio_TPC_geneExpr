
############# WGCNA (Weighted Gene Co-expression Network Analysis) #############
############# MCAP (Montipora capitata) RNAseq data #############
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
treatmentinfo <- read.csv("RNAseq_Mcap_data.csv", header = TRUE, sep = ";")

#gene count matrix
gcount <- as.data.frame(read.csv("gene_count_matrix_noIso_Mcap2.csv", row.names="gene_id"), colClasses = double)

# Reorder treatmentinfo to match the columns of gcount_filt
treatmentinfo <- treatmentinfo[match(colnames(gcount), treatmentinfo$sample_id), ]

## Quality-filter gene counts
#Set filter values for PoverA: smallest sample size per treat is 4, so 4/24 (24 samples) is 0.17
#This means that 4 out of 24 (0.17) samples need to have counts over 10.
#So P=17 percent of the samples have counts over A=10. 
filt <- filterfun(pOverA(0.17,10))

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
# [1] 54384

nrow(gcount_filt) #After
# [1] 27540


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

ggsave(file = "PCA_all_vst_Mcap.png", PCA)


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
#> [1] 27540
# Number of samples
nrow(datExpr)
#> [1] 24

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
save(adjacency, TOM, dissTOM, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/adjTOM.RData")
save(dissTOM, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/dissTOM.RData") 

# Clustering using TOM
#Form distance matrix
geneTree= flashClust(as.dist(dissTOM), method="average")

#We will now plot a dendrogram of genes. Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/dissTOMClustering.pdf", width=20, height=20)
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
save(dynamicMods, geneTree, file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/dyMod_geneTree.RData")

dyMod_geneTree <- load(file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/dyMod_geneTree.RData")

dyMod_geneTree
## [1] "dynamicMods" "geneTree"

# Plot the module assignment under the gene dendrogram
dynamicColors = labels2colors(dynamicMods) # Convert numeric labels into colors
table(dynamicColors)
# dynamicColors
#         bisque4           black            blue           brown          brown4            cyan 
#             122             615            5774            2404             122             364 
#       darkgreen        darkgrey     darkmagenta  darkolivegreen      darkorange     darkorange2 
#             223             216             177             180             201             131 

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/dissTOMColorClustering.pdf")
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

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/eigengeneClustering1.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

#Merge modules with >80% eigengene similarity (most studies use 80-90% similarity)

MEDissThres= 0.20 #merge modules that are 80% similar

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/eigengeneClustering2.pdf", width = 20)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=MEDissThres, col="red")
dev.off()

merge= mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/mergedClusters.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

#Save new colors

moduleColors = mergedColors # Rename to moduleColors
colorOrder = c("grey", standardColors(50)); # Construct numerical labels corresponding to the colors
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
ncol(MEs) 
# [1] 34

# Plot new tree
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster again and plot the results
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/eigengeneClustering3.pdf")
METree = flashClust(as.dist(MEDiss), method = "average")
MEtreePlot = plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()

# Relating modules to temp, quantifying module–trait associations
#Prepare trait data. Data has to be numeric, so I replaced the temp for numeric values

treatmentinfo$temp <- factor(paste0("t", as.character(treatmentinfo$temp)))
#allTraits <- names(treatmentinfo$temp)
allTraits <-levels(treatmentinfo$temp)
allTraits$tt12 <- c(0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0)
allTraits$tt18 <- c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0)
allTraits$tt25 <- c(0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0)
allTraits$tt26.8<-c(0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0)
allTraits$tt30 <- c(0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1)
allTraits$tt35 <- c(0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0)

datTraits <- as.data.frame(allTraits)
datTraits <- datTraits[ , -(1:6)]
dim(datTraits)
## [1] 24 6

rownames(datTraits) <- treatmentinfo$sample_id
print(datTraits)
#     tt12 tt18 tt25 tt26.8 tt30 tt35
# B1     0    1    0      0    0    0
# B10    0    0    0      1    0    0
# B6     1    0    0      0    0    0
# B7     0    0    0      0    0    1
# B8     0    0    1      0    0    0
# B9     0    0    0      0    1    0
# F1     0    1    0      0    0    0
# F10    0    0    0      1    0    0
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
names(MEs) #head
# [1] "MEplum2"           "MEdarkmagenta"     "MEdarkorange"      "MEblack"          
#  [5] "MEivory"           "MEpurple"          "MEsalmon4"         "MEskyblue"        
#  [9] "MEwhite"           "MEyellowgreen"     "MEdarkolivegreen"  "MEthistle2"       
# [13] "MEgreenyellow"     "MEpink"            "MElightcyan"       "MEorange"         

## Correlations of traits and genes with eigengenes
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
Colors=sub("ME","",names(MEs))

moduleTraitTree = hclust(dist(t(moduleTraitCor)), method = "average");
pdf(file="/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/Temp clustering based on module-trait correlation.pdf")
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

pdf(file = "/work/pi_hputnam_uri_edu/fscucchia/20250424_ENCORE_HawaiiTPC_Federica/ENCORE_Hawaii_TPC_Rstudio/GWENA/Mcap/Module-trait-relationship-heatmap_Mcap.pdf", height = 9, width = 8)
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
#         MEplum2 MEdarkmagenta MEdarkorange     MEblack    MEivory   MEpurple   MEsalmon4   MEskyblue
# B1   0.06594039    -0.2993033   -0.2934318 -0.14503290 0.02064023 -0.2621171 -0.12212608 -0.07317448
# B10 -0.27386213    -0.3248712   -0.2425439 -0.08834224 0.09621475 -0.2144316 -0.08429681  0.01814720
# B6   0.02367071    -0.3033032   -0.2669567 -0.09628341 0.09555983 -0.2273501 -0.04835808 -0.02781063
# B7  -0.14404129    -0.2532995   -0.3583206 -0.13122981 0.01916209  0.1329105  0.30885306 -0.12650133

names(MEs)
#  [1] "MEplum2"           "MEdarkmagenta"     "MEdarkorange"      "MEblack"          
#  [5] "MEivory"           "MEpurple"          "MEsalmon4"         "MEskyblue"    

Strader_MEs <- MEs
Strader_MEs$temp <- treatmentinfo$temp
Strader_MEs$sample_id <- rownames(Strader_MEs)
head(Strader_MEs)
#        MEplum2 MEdarkmagenta MEdarkorange     MEblack    MEivory   MEpurple   MEsalmon4   MEskyblue
# B1   0.06594039    -0.2993033   -0.2934318 -0.14503290 0.02064023 -0.2621171 -0.12212608 -0.07317448
# B10 -0.27386213    -0.3248712   -0.2425439 -0.08834224 0.09621475 -0.2144316 -0.08429681  0.01814720
# B6   0.02367071    -0.3033032   -0.2669567 -0.09628341 0.09555983 -0.2273501 -0.04835808 -0.02781063

# Calculate 10 over-arching expression patterns using mean eigengene for each module in a cluster
C1_Strader_MEs <- select(Strader_MEs, MEplum2:MEdarkorange)
C1_Strader_MEs$Mean <- rowMeans(C1_Strader_MEs)
C2_Strader_MEs <- select(Strader_MEs, MEblack:MEsalmon4)
C2_Strader_MEs$Mean <- rowMeans(C2_Strader_MEs)
C3_Strader_MEs <- select(Strader_MEs, MEskyblue:MEthistle2)
C3_Strader_MEs$Mean <- rowMeans(C3_Strader_MEs)
C4_Strader_MEs <- select(Strader_MEs, MEgreenyellow:MEpink)
C4_Strader_MEs$Mean <- rowMeans(C4_Strader_MEs)
C5_Strader_MEs <- select(Strader_MEs, MElightcyan:MEsalmon)
C5_Strader_MEs$Mean <- rowMeans(C5_Strader_MEs)
C6_Strader_MEs <- select(Strader_MEs, MElightgreen:MEsienna3)
C6_Strader_MEs$Mean <- rowMeans(C6_Strader_MEs)
C7_Strader_MEs <- select(Strader_MEs, MElightcyan1:MElightyellow)
C7_Strader_MEs$Mean <- rowMeans(C7_Strader_MEs)
C8_Strader_MEs <- select(Strader_MEs, MEturquoise)
C8_Strader_MEs$Mean <- rowMeans(C8_Strader_MEs)
C9_Strader_MEs <- select(Strader_MEs, MEdarkgreen:MEgrey60)
C9_Strader_MEs$Mean <- rowMeans(C9_Strader_MEs)
C10_Strader_MEs <- select(Strader_MEs, MEbrown4:MEskyblue3)
C10_Strader_MEs$Mean <- rowMeans(C10_Strader_MEs)

Strader_MEs$temp <- as.character(Strader_MEs$temp)
expressionProfile_data <- as.data.frame(cbind(temp = Strader_MEs$temp, cluster1= C1_Strader_MEs$Mean, cluster2 = C2_Strader_MEs$Mean, 
                          cluster3 = C3_Strader_MEs$Mean, cluster4 = C4_Strader_MEs$Mean,cluster5 = C5_Strader_MEs$Mean,cluster6 = C6_Strader_MEs$Mean,
                          cluster7 = C7_Strader_MEs$Mean,cluster8 = C8_Strader_MEs$Mean,cluster9 = C9_Strader_MEs$Mean,cluster10 = C10_Strader_MEs$Mean))

head(expressionProfile_data)
#   temp           cluster1            cluster2           cluster3           cluster4
# 1   t18 -0.175598238178533  -0.127158968149933 -0.169637664501881  -0.17484121047813
# 2 t26.8 -0.280425757295429 -0.0727139820327856 -0.146264686461688  -0.19489540478974
# 3   t12 -0.182196382534518 -0.0691079364058859 -0.150357492317928 -0.194728300706314
# 4   t35 -0.251887119971089  0.0824239549992588 -0.279288812026189 -0.190274112963049

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

### save mean eigengene values for cluster5
meanEigenClust5 <- expressionProfile_data$cluster5
write.csv(meanEigenClust5, file = "meanEigenClust5.csv")
write.csv(C5_Strader_MEs, file = "C5_Strader_MEs.csv")

### save mean eigengene values for cluster6
meanEigenClust6 <- expressionProfile_data$cluster6
write.csv(meanEigenClust6, file = "meanEigenClust6.csv")
write.csv(C6_Strader_MEs, file = "C6_Strader_MEs.csv")

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

# Filter to exclude clusters 1,4,7 and 10 (no significant correlations in any of the temperature groups)
expressionProfile_data_long_subset <- expressionProfile_data_long %>%
  filter(cluster %in% c("cluster2", "cluster3", "cluster5", "cluster6", "cluster8", "cluster9"))

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

# Create a palette of green shades (light to dark)
green_palette <- colorRampPalette(c("#69f0ae", "#00e676", "#006400"))(length(levels(expressionProfile_data_long_subset$temp)))

# Plot with custom colors for boxplots and dots, and no grid
ggplot(expressionProfile_data_long_subset, aes(x = temp, y = Eigengene, group = temp)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp), width = 0.2, size = 2, show.legend = FALSE) +
  facet_wrap(~ cluster, scales = "free_y") +
  scale_color_manual(values = green_palette) +
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
### selected modules: 2,3,5,6,8,9

library(segmented)

# Make sure temp is numeric (remove "tt" if needed)
expressionProfile_data_long_subset$temp_num <- as.numeric(sub("^t", "", expressionProfile_data_long_subset$temp))

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
# psi1.temp_num 20.688  5.547
# Coefficients of the linear terms:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)  0.19971    0.30992   0.644    0.527
# temp_num    -0.01517    0.02026  -0.749    0.463
# U1.temp_num  0.03022    0.02321   1.302       NA
# Residual standard error: 0.1719 on 20 degrees of freedom
# Multiple R-Squared: 0.1168,  Adjusted R-squared: -0.01569 

#There is no significant breakpoint in this fit.
#The estimated breakpoint (psi1.temp_num) is at 20.7°C, but the standard error is large (5.5), 
# indicating high uncertainty. The slope before the breakpoint (temp_num) and the change in slope 
# after the breakpoint (U1.temp_num) both have high p-values (0.463 and NA), meaning they are not statistically 
# significant. R-squared is very low (0.12), and the adjusted R-squared is negative, indicating the model 
# does not explain the data well. There is no evidence for a significant breakpoint or trend in eigengene 
# values for this cluster across temperature. The data do not support a segmented (piecewise) relationship here.

# Filter for the module of interest (e.g., cluster3)
df_mod3 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster3")

# Fit linear model to all replicates
fit3 <- lm(Eigengene ~ temp_num, data = df_mod3)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit3 <- segmented(fit3, seg.Z = ~temp_num, npsi = 1)

summary(seg_fit3)
#Multiple R-Squared: 0.03849 - There is no significant breakpoint in this fit
plot(seg_fit3)



# Filter for the module of interest (e.g., cluster5)
df_mod5 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster5")

# Fit linear model to all replicates
fit5 <- lm(Eigengene ~ temp_num, data = df_mod5)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit5 <- segmented(fit5, seg.Z = ~temp_num, npsi = 2)

summary(seg_fit5)
#Multiple R-Squared: 0.03849 - There is no significant breakpoint in this fit
plot(seg_fit5)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 22.097  3.282
# psi2.temp_num 28.663  1.177
# Coefficients of the linear terms:
#             Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  0.31840    0.17215   1.850   0.0809 .
# temp_num    -0.01749    0.01125  -1.554   0.1376  
# U1.temp_num  0.05189    0.03916   1.325       NA  
# U2.temp_num -0.10649    0.03987  -2.671       NA  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.09549 on 18 degrees of freedom
# Multiple R-Squared: 0.7367,  Adjusted R-squared: 0.6636 

# Are the breakpoints significant?
# The model fit is good (high R²).
# The second change in slope (U2.temp_num) is likely meaningful, given the t value, but the 
# package does not provide a p-value for it. Standard errors for breakpoints are not huge, so the 
# locations are reasonably certain. Slopes before and after breakpoints: Only the last change 
# in slope (U2.temp_num) appears to be substantial (t value -2.67).


AIC(fit5, seg_fit5)
#          df       AIC
# fit5      3 -21.75410
# seg_fit5  7 -37.53489
#Lower AIC/BIC suggests the segmented model best fits the data



##  Get fitted values and confidence intervals
# Get predicted values and confidence intervals
library(segmented)
library(boot)

# Set up the bootstrap function
# Function to fit segmented and predict at grid points
boot_seg <- function(data, indices, grid_temp) {
  d <- data[indices, ]
  fit <- lm(Eigengene ~ temp_num, data = d)
  seg_fit <- try(segmented(fit, seg.Z = ~temp_num, npsi = 2), silent = TRUE)
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
grid_temp <- seq(min(df_mod5$temp_num), max(df_mod5$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod5,
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
pred_df5_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod5$temp_num <- as.numeric(sub("^t", "", df_mod5$temp))

df_mod5$temp_num_f <- factor(df_mod5$temp_num)
ggplot(df_mod5, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df5_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df5_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = green_palette) +
  scale_x_continuous(
    breaks = unique(df_mod5$temp_num),
    labels = gsub("^t", "", unique(df_mod5$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module5: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


# Filter for the module of interest (e.g., cluster6)
df_mod6 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster6")

# Fit linear model to all replicates
fit6 <- lm(Eigengene ~ temp_num, data = df_mod6)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit6 <- segmented(fit6, seg.Z = ~temp_num, npsi = 2)

summary(seg_fit6)
# Call: 
# segmented.lm(obj = fit6, seg.Z = ~temp_num, npsi = 2)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 18.502  3.603
# psi2.temp_num 29.804  2.254
# Coefficients of the linear terms:
#             Estimate Std. Error t value Pr(>|t|)   
# (Intercept)  0.43566    0.19847   2.195   0.0415 * 
# temp_num    -0.04136    0.01297  -3.188   0.0051 **
# U1.temp_num  0.09099    0.04515   2.015       NA   
# U2.temp_num -0.08096    0.04597  -1.761       NA   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# Residual standard error: 0.1101 on 18 degrees of freedom
# Multiple R-Squared: 0.7474,  Adjusted R-squared: 0.6772 
# Boot restarting based on 6 samples. Last fit:
# Convergence attained in 2 iterations (rel. change 7.9211e-11))

#compare models
AIC(fit6, seg_fit6)
#          df       AIC
# fit6      3 -15.09419
# seg_fit6  7 -30.70362
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
  seg_fit <- try(segmented(fit, seg.Z = ~temp_num, npsi = 2), silent = TRUE)
  if (inherits(seg_fit, "try-error")) {
    return(rep(NA, length(grid_temp)))
  }
  pred <- predict(seg_fit, newdata = data.frame(temp_num = grid_temp))
  return(pred)
}

#Run the bootstrap
set.seed(123)
grid_temp <- seq(min(df_mod6$temp_num), max(df_mod6$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod6,
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
pred_df6_boot <- data.frame(
  temp_num = grid_temp,
  fit = fit_median,
  lower = ci_lower,
  upper = ci_upper
)

##Overlay on your Strader plot
df_mod6$temp_num <- as.numeric(sub("^t", "", df_mod6$temp))

df_mod6$temp_num_f <- factor(df_mod6$temp_num)
ggplot(df_mod6, aes(x = temp_num, y = Eigengene, group = temp_num_f)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", size = 0.7) +
  geom_boxplot(aes(color = temp_num_f, group = temp_num_f), outlier.shape = NA, fill = NA, size = 0.5) +
  geom_jitter(aes(color = temp_num_f), width = 0.2, size = 2, show.legend = FALSE) +
  geom_line(data = pred_df6_boot, aes(x = temp_num, y = fit), color = "red", linewidth = 0.7, inherit.aes = FALSE) +
  geom_ribbon(data = pred_df6_boot, aes(x = temp_num, ymin = lower, ymax = upper), alpha = 0.1, fill = "red", inherit.aes = FALSE) +
  scale_color_manual(values = green_palette) +
  scale_x_continuous(
    breaks = unique(df_mod6$temp_num),
    labels = gsub("^t", "", unique(df_mod6$temp))
  ) +
  theme_bw() +
  labs(x = "Temperature (°C)", y = "Module Eigengene Value",
       title = "Module6: Strader Plot with Bootstrapped Segmented Regression CI") +
  theme(
    strip.background = element_rect(fill = "#e0ffff", color = NA),
    panel.grid = element_blank()
  )


# Filter for the module of interest (e.g., cluster8)
df_mod8 <- expressionProfile_data_long_subset %>% filter(cluster == "cluster8")

# Fit linear model to all replicates
fit8 <- lm(Eigengene ~ temp_num, data = df_mod8)

# Fit segmented regression (e.g., 1 breakpoint)
seg_fit8 <- segmented(fit8, seg.Z = ~temp_num, npsi = 2)

summary(seg_fit8)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 19.662 12.142
# psi2.temp_num 29.400  1.355
# Coefficients of the linear terms:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept) -0.23296    0.17064  -1.365    0.189
# temp_num     0.01020    0.01115   0.915    0.372
# U1.temp_num -0.01996    0.03882  -0.514       NA
# U2.temp_num  0.10673    0.03952   2.701       NA

# Residual standard error: 0.09465 on 18 degrees of freedom
# Multiple R-Squared: 0.8387,  Adjusted R-squared: 0.7939 

# Boot restarting based on 6 samples. Last fit:
# Convergence attained in 2 iterations (rel. change 8.2986e-12)

#compare models
AIC(fit8, seg_fit8)
#          df       AIC
# fit8      3 -12.11957
# seg_fit8  7 -37.95740
#Lower AIC/BIC suggests the segmented model best fits the data


##  Get fitted values and confidence intervals
# Get predicted values and confidence intervals
library(segmented)
library(boot)

# Set up the bootstrap function
# Function to fit segmented and predict at grid points
boot_seg <- function(data, indices, grid_temp) {
  d <- data[indices, ]
  fit <- lm(Eigengene ~ temp_num, data = d)
  seg_fit <- try(segmented(fit, seg.Z = ~temp_num, npsi = 2), silent = TRUE)
  if (inherits(seg_fit, "try-error")) {
    return(rep(NA, length(grid_temp)))
  }
  pred <- predict(seg_fit, newdata = data.frame(temp_num = grid_temp))
  return(pred)
}

#Run the bootstrap
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
  scale_color_manual(values = green_palette) +
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
seg_fit9 <- segmented(fit9, seg.Z = ~temp_num, npsi = 2)

summary(seg_fit9)
# segmented.lm(obj = fit9, seg.Z = ~temp_num, npsi = 2)
# Estimated Break-Point(s):
#                  Est. St.Err
# psi1.temp_num 19.830  1.747
# psi2.temp_num 27.483  1.048
# Coefficients of the linear terms:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.86594    0.17857  -4.849 0.000129 ***
# temp_num     0.06472    0.01167   5.544  2.9e-05 ***
# U1.temp_num -0.14189    0.04063  -3.493       NA    
# U2.temp_num  0.10089    0.04136   2.440       NA    
# Residual standard error: 0.09905 on 18 degrees of freedom
# Multiple R-Squared: 0.7391,  Adjusted R-squared: 0.6666 


#compare models
AIC(fit9, seg_fit9)
#          df       AIC
# fit9      3 -12.69290
# seg_fit9  7 -35.77522
#Lower AIC/BIC suggests the segmented model best fits the data


##  Get fitted values and confidence intervals
# Get predicted values and confidence intervals
library(segmented)
library(boot)

# Set up the bootstrap function
# Function to fit segmented and predict at grid points
boot_seg <- function(data, indices, grid_temp) {
  d <- data[indices, ]
  fit <- lm(Eigengene ~ temp_num, data = d)
  seg_fit <- try(segmented(fit, seg.Z = ~temp_num, npsi = 2), silent = TRUE)
  if (inherits(seg_fit, "try-error")) {
    return(rep(NA, length(grid_temp)))
  }
  pred <- predict(seg_fit, newdata = data.frame(temp_num = grid_temp))
  return(pred)
}

#Run the bootstrap
set.seed(123)
grid_temp <- seq(min(df_mod9$temp_num), max(df_mod9$temp_num), length.out = 200)
boot_out <- boot(
  data = df_mod9,
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
  scale_color_manual(values = green_palette) +
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






############## Find Hub Genes and Network Analysis

### export network of clusters with significant breakpoints and trends (cluster5, cluster6, cluster8, cluster9) - see analysis above

# Subset samples for each group
control_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "t26.8"]
treat30_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "t30"]
treat35_samples <- treatmentinfo$sample_id[treatmentinfo$temp == "t35"]

#### Export the network of a specific cluster - cluster5
cluster5_modules <- c("lightcyan", "orange", "blue", "salmon") #Define the module colors for your cluster of interest
inCluster5 <- moduleColors %in% cluster5_modules #Select genes belonging to these modules
cluster5_genes <- colnames(datExpr)[inCluster5]  ##10326

##Save only the essentials for Network Analysis
# save(
#   datExpr,     
#   moduleColors,
#   control_samples,
#   treat30_samples,
#   treat35_samples,
#   file = "network_analysis_essentials_Mcap.RData"
# )

# Subset expression data for each group
datExpr_control_cluster5 <- datExpr[control_samples, cluster5_genes ]
datExpr_treat30_cluster5 <- datExpr[treat30_samples, cluster5_genes ]
datExpr_treat35_cluster5 <- datExpr[treat35_samples, cluster5_genes ]

# Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# kWithin: sum of connection strengths with other module genes
#add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster5)
# names(kWithin_control) <- colnames(datExpr_control_cluster5)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster5)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster5)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster5)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster5)

# Number of top 10% hub genes for each group
#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module
# Use sum of absolute correlations as a proxy for connectivity (ranks genes by their overall co-expression with all other genes in the module)

# For control group
cor_mat_control <- cor(datExpr_control_cluster5, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_control) <- 0
hub_score_control <- rowSums(abs(cor_mat_control), na.rm = TRUE)
names(hub_score_control) <- colnames(datExpr_control_cluster5)
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
top10pct_hubs_control_5 <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]

# For treat30 group
cor_mat_treat30 <- cor(datExpr_treat30_cluster5, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat30) <- 0
hub_score_treat30 <- rowSums(abs(cor_mat_treat30), na.rm = TRUE)
names(hub_score_treat30) <- colnames(datExpr_treat30_cluster5)
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top10pct_hubs_treat30_5 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]

# For treat35 group 
cor_mat_treat35 <- cor(datExpr_treat35_cluster5, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat35) <- 0
hub_score_treat35 <- rowSums(abs(cor_mat_treat35), na.rm = TRUE)
names(hub_score_treat35) <- colnames(datExpr_treat35_cluster5)
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top10pct_hubs_treat35_5 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_cluster5_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster5_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster5_2.csv", row.names = FALSE)


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
#Montipora_capitata_HIv3___RNAseq.g26531.t1 
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Montipora_capitata_HIv3___RNAseq.g20492.t1 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Montipora_capitata_HIv3___RNAseq.g3459.t1 

# ####### Compare module preservation across conditions 
# ####### Highlight Nodes That Change Hub Status
# # Highlight Top 10% Hubs
# # Hub changing status between control and temps at break-points per (30 and 35) for all clusters

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(igraph)

# # Select only the top 10 hub genes for each condition (less computationally heavy for plotting)
# top10_genes_control <- top10pct_hubs_control[1:10]
# top10_genes_treat30 <- top10pct_hubs_treat30[1:10]
# top10_genes_treat35 <- top10pct_hubs_treat35[1:10]

# # Subset expression data to these genes
# datExpr_control_top10 <- datExpr_control_cluster5[, top10_genes_control, drop = FALSE]
# datExpr_treat30_top10 <- datExpr_treat30_cluster5[, top10_genes_treat30, drop = FALSE]
# datExpr_treat35_top10 <- datExpr_treat35_cluster5[, top10_genes_treat35, drop = FALSE]

# # Build adjacency/correlation matrices
# adj_control_top10 <- abs(cor(datExpr_control_top10, method = "pearson"))
# adj_treat30_top10 <- abs(cor(datExpr_treat30_top10, method = "pearson"))
# adj_treat35_top10 <- abs(cor(datExpr_treat35_top10, method = "pearson"))

# # Create igraph objects
# g_control_top10 <- graph_from_adjacency_matrix(adj_control_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat30_top10 <- graph_from_adjacency_matrix(adj_treat30_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat35_top10 <- graph_from_adjacency_matrix(adj_treat35_top10, mode = "undirected", weighted = TRUE, diag = FALSE)

# # Highlight the top hub gene (by kWithin) in red, others in gray
# V(g_control_top10)$color <- ifelse(names(V(g_control_top10)) == top_hub_control, "red", "gray")
# V(g_treat30_top10)$color <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, "red", "gray")
# V(g_treat35_top10)$color <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, "red", "gray")

# V(g_control_top10)$label <- ifelse(names(V(g_control_top10)) == top_hub_control, top_hub_control, "")
# V(g_treat30_top10)$label <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, top_hub_treat30, "")
# V(g_treat35_top10)$label <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, top_hub_treat35, "")

# # Plot networks
# pdf("network_top10genes_cluster5.pdf", width = 12, height = 8)
# par(mfrow = c(1, 3))
# plot(g_control_top10, main = "Control: Top 10 Hubs", vertex.label = V(g_control_top10)$label, vertex.size = 8)
# plot(g_treat30_top10, main = "Treat30: Top 10 Hubs", vertex.label = V(g_treat30_top10)$label, vertex.size = 8)
# plot(g_treat35_top10, main = "Treat35: Top 10 Hubs", vertex.label = V(g_treat35_top10)$label, vertex.size = 8)
# par(mfrow = c(1, 1))
# dev.off()

#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 5
all_genes <- cluster5_genes

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
                   labels = c("Cluster" = "Cluster 5 Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Cluster 5: Overlap of Top 10% Hub Genes Across Conditions",
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

bar5 <- ggplot(membership_summary, aes(x = "Cluster 5", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Cluster 5: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


#### Export the network of a specific cluster - cluster6
cluster6_modules <- c("lightgreen", "sienna3")  # replace with your actual module colors for cluster6
inCluster6 <- moduleColors %in% cluster6_modules #Select genes belonging to these modules
cluster6_genes <- colnames(datExpr)[inCluster6] ##640

# Subset expression data for each group
datExpr_control_cluster6 <- datExpr[control_samples, cluster6_genes]
datExpr_treat30_cluster6 <- datExpr[treat30_samples, cluster6_genes]
datExpr_treat35_cluster6 <- datExpr[treat35_samples, cluster6_genes]

# # Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# # kWithin: sum of connection strengths with other module genes
# #add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster6)
# names(kWithin_control) <- colnames(datExpr_control_cluster6)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster6)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster6)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster6)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster6)

# Number of top 10% hub genes for each group
#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module
# Use sum of absolute correlations as a proxy for connectivity (ranks genes by their overall co-expression with all other genes in the module)

# For control group
cor_mat_control <- cor(datExpr_control_cluster6, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_control) <- 0
hub_score_control <- rowSums(abs(cor_mat_control), na.rm = TRUE)
names(hub_score_control) <- colnames(datExpr_control_cluster6)
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
top10pct_hubs_control_6 <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]

# For treat30 group
cor_mat_treat30 <- cor(datExpr_treat30_cluster6, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat30) <- 0
hub_score_treat30 <- rowSums(abs(cor_mat_treat30), na.rm = TRUE)
names(hub_score_treat30) <- colnames(datExpr_treat30_cluster6)
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top10pct_hubs_treat30_6 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]

# For treat35 group 
cor_mat_treat35 <- cor(datExpr_treat35_cluster6, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat35) <- 0
hub_score_treat35 <- rowSums(abs(cor_mat_treat35), na.rm = TRUE)
names(hub_score_treat35) <- colnames(datExpr_treat35_cluster6)
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top10pct_hubs_treat35_6 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_cluster6_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster6_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster6_2.csv", row.names = FALSE)


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
#Montipora_capitata_HIv3___RNAseq.g454.t1
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Montipora_capitata_HIv3___TS.g9235.t1  
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
# Montipora_capitata_HIv3___TS.g46337.t1 

####### Compare module preservation across conditions 
####### Highlight Nodes That Change Hub Status
# Highlight Top 10% Hubs
# Hub changing status between control and temps at break-points per (30 and 35) for all clusters

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)

# Select only the top 10 hub genes for each condition (less computationally heavy for plotting)
# top10_genes_control <- top10pct_hubs_control[1:10]
# top10_genes_treat30 <- top10pct_hubs_treat30[1:10]
# top10_genes_treat35 <- top10pct_hubs_treat35[1:10]

# # Subset expression data to these genes
# datExpr_control_top10 <- datExpr_control_cluster6[, top10_genes_control, drop = FALSE]
# datExpr_treat30_top10 <- datExpr_treat30_cluster6[, top10_genes_treat30, drop = FALSE]
# datExpr_treat35_top10 <- datExpr_treat35_cluster6[, top10_genes_treat35, drop = FALSE]

# # Build adjacency/correlation matrices
# adj_control_top10 <- abs(cor(datExpr_control_top10, method = "pearson"))
# adj_treat30_top10 <- abs(cor(datExpr_treat30_top10, method = "pearson"))
# adj_treat35_top10 <- abs(cor(datExpr_treat35_top10, method = "pearson"))

# # Create igraph objects
# g_control_top10 <- graph_from_adjacency_matrix(adj_control_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat30_top10 <- graph_from_adjacency_matrix(adj_treat30_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat35_top10 <- graph_from_adjacency_matrix(adj_treat35_top10, mode = "undirected", weighted = TRUE, diag = FALSE)

# # Highlight the top hub gene (by kWithin) in red, others in gray
# V(g_control_top10)$color <- ifelse(names(V(g_control_top10)) == top_hub_control, "red", "gray")
# V(g_treat30_top10)$color <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, "red", "gray")
# V(g_treat35_top10)$color <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, "red", "gray")

# V(g_control_top10)$label <- ifelse(names(V(g_control_top10)) == top_hub_control, top_hub_control, "")
# V(g_treat30_top10)$label <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, top_hub_treat30, "")
# V(g_treat35_top10)$label <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, top_hub_treat35, "")

# # Plot networks
# pdf("network_top10genes_cluster6.pdf", width = 12, height = 8)
# par(mfrow = c(1, 3))
# plot(g_control_top10, main = "Control: Top 10 Hubs", vertex.label = V(g_control_top10)$label, vertex.size = 8)
# plot(g_treat30_top10, main = "Treat30: Top 10 Hubs", vertex.label = V(g_treat30_top10)$label, vertex.size = 8)
# plot(g_treat35_top10, main = "Treat35: Top 10 Hubs", vertex.label = V(g_treat35_top10)$label, vertex.size = 8)
# par(mfrow = c(1, 1))
# dev.off()

#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 6
all_genes <- cluster6_genes

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
  mutate(Cluster = "Cluster 6 Genes")

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
                   labels = c("Cluster" = "Cluster 6 Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Cluster 6: Overlap of Top 10% Hub Genes Across Conditions",
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

bar6 <- ggplot(membership_summary, aes(x = "Cluster 6", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Cluster 6: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


#### Export the network of a specific cluster - cluster8
cluster8_modules <- c("turquoise")  # replace with your actual module colors for cluster8
inCluster8 <- moduleColors %in% cluster8_modules #Select genes belonging to these modules
cluster8_genes <- colnames(datExpr)[inCluster8] ##6303

# Subset expression data for each group
datExpr_control_cluster8 <- datExpr[control_samples, cluster8_genes]
datExpr_treat30_cluster8 <- datExpr[treat30_samples, cluster8_genes]
datExpr_treat35_cluster8 <- datExpr[treat35_samples, cluster8_genes]

# # Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# # kWithin: sum of connection strengths with other module genes
# #add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster8)
# names(kWithin_control) <- colnames(datExpr_control_cluster8)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster8)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster8)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster8)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster8)

# Number of top 10% hub genes for each group
#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module
# Use sum of absolute correlations as a proxy for connectivity (ranks genes by their overall co-expression with all other genes in the module)

# For control group
cor_mat_control <- cor(datExpr_control_cluster8, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_control) <- 0
hub_score_control <- rowSums(abs(cor_mat_control), na.rm = TRUE)
names(hub_score_control) <- colnames(datExpr_control_cluster8)
top10pct_hubs_control <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]
top10pct_hubs_control_8 <- names(sort(hub_score_control, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_control))]

# For treat30 group
cor_mat_treat30 <- cor(datExpr_treat30_cluster8, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat30) <- 0
hub_score_treat30 <- rowSums(abs(cor_mat_treat30), na.rm = TRUE)
names(hub_score_treat30) <- colnames(datExpr_treat30_cluster8)
top10pct_hubs_treat30 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]
top10pct_hubs_treat30_8 <- names(sort(hub_score_treat30, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat30))]

# For treat35 group 
cor_mat_treat35 <- cor(datExpr_treat35_cluster8, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat_treat35) <- 0
hub_score_treat35 <- rowSums(abs(cor_mat_treat35), na.rm = TRUE)
names(hub_score_treat35) <- colnames(datExpr_treat35_cluster8)
top10pct_hubs_treat35 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]
top10pct_hubs_treat35_8 <- names(sort(hub_score_treat35, decreasing = TRUE))[1:ceiling(0.10 * length(hub_score_treat35))]

# Save the top 10% hub genes as CSV files
write.csv(data.frame(Gene = top10pct_hubs_control),
          file = "top10pct_hub_genes_control_cluster8_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster8_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster8_2.csv", row.names = FALSE)

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
#Montipora_capitata_HIv3___RNAseq.g6968.t1  
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Montipora_capitata_HIv3___RNAseq.g44936.t1 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Montipora_capitata_HIv3___TS.g27549.t3 

####### Compare module preservation across conditions 
####### Highlight Nodes That Change Hub Status
# Highlight Top 10% Hubs
# Hub changing status between control and temps at break-points per (30 and 35) for all clusters

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(igraph)

# # Select only the top 10 hub genes for each condition (less computationally heavy for plotting)
# top10_genes_control <- top10pct_hubs_control[1:10]
# top10_genes_treat30 <- top10pct_hubs_treat30[1:10]
# top10_genes_treat35 <- top10pct_hubs_treat35[1:10]

# # Subset expression data to these genes
# datExpr_control_top10 <- datExpr_control_cluster8[, top10_genes_control, drop = FALSE]
# datExpr_treat30_top10 <- datExpr_treat30_cluster8[, top10_genes_treat30, drop = FALSE]
# datExpr_treat35_top10 <- datExpr_treat35_cluster8[, top10_genes_treat35, drop = FALSE]

# # Build adjacency/correlation matrices
# adj_control_top10 <- abs(cor(datExpr_control_top10, method = "pearson"))
# adj_treat30_top10 <- abs(cor(datExpr_treat30_top10, method = "pearson"))
# adj_treat35_top10 <- abs(cor(datExpr_treat35_top10, method = "pearson"))

# # Create igraph objects
# g_control_top10 <- graph_from_adjacency_matrix(adj_control_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat30_top10 <- graph_from_adjacency_matrix(adj_treat30_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat35_top10 <- graph_from_adjacency_matrix(adj_treat35_top10, mode = "undirected", weighted = TRUE, diag = FALSE)

# # Highlight the top hub gene (by kWithin) in red, others in gray
# V(g_control_top10)$color <- ifelse(names(V(g_control_top10)) == top_hub_control, "red", "gray")
# V(g_treat30_top10)$color <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, "red", "gray")
# V(g_treat35_top10)$color <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, "red", "gray")

# V(g_control_top10)$label <- ifelse(names(V(g_control_top10)) == top_hub_control, top_hub_control, "")
# V(g_treat30_top10)$label <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, top_hub_treat30, "")
# V(g_treat35_top10)$label <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, top_hub_treat35, "")

# # Plot networks
# pdf("network_top10genes_cluster8.pdf", width = 12, height = 8)
# par(mfrow = c(1, 3))
# plot(g_control_top10, main = "Control: Top 10 Hubs", vertex.label = V(g_control_top10)$label, vertex.size = 8)
# plot(g_treat30_top10, main = "Treat30: Top 10 Hubs", vertex.label = V(g_treat30_top10)$label, vertex.size = 8)
# plot(g_treat35_top10, main = "Treat35: Top 10 Hubs", vertex.label = V(g_treat35_top10)$label, vertex.size = 8)
# par(mfrow = c(1, 1))
# dev.off()

#### Using the ggalluvial package to visualize the overlap of the top 10% hub genes in selected clusters for control, treat30, and 
#### treat35.

library(ggalluvial)
library(dplyr)

# All genes in cluster 8
all_genes <- cluster8_genes

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
  mutate(Cluster = "Cluster 8 Genes")

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
                   labels = c("Cluster" = "Cluster 8 Genes",
                              "Control" = "Control",
                              "Treat30" = "Treat30",
                              "Treat35" = "Treat35")) +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = setNames(RColorBrewer::brewer.pal(length(pattern_labels), "Set2"), names(pattern_labels)),
    labels = pattern_labels
  ) +
  theme_minimal() +
  labs(title = "Cluster 8: Overlap of Top 10% Hub Genes Across Conditions",
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

bar8 <- ggplot(membership_summary, aes(x = "Cluster 8", y = Percent, fill = pattern)) +
  geom_bar(stat = "identity", width = 0.1, color = "black") +
  scale_fill_manual(
    name = "Hub Pattern\n(Control.Treat30.Treat35)",
    values = bar_colors,
    labels = pattern_labels
  ) +
  labs(
    title = "Cluster 8: Distribution of Hub Gene Overlap Patterns",
    x = "",
    y = "Percent of Unique Hub Genes"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )




#### Export the network of a specific cluster - cluster9
cluster9_modules <- c("darkgreen", "yellow", "bisque4", "grey60")  # replace with your actual module colors for cluster9
inCluster9 <- moduleColors %in% cluster9_modules #Select genes belonging to these modules
cluster9_genes <- colnames(datExpr)[inCluster9] ##2057

# Subset expression data for each group
datExpr_control_cluster9 <- datExpr[control_samples, cluster9_genes]
datExpr_treat30_cluster9 <- datExpr[treat30_samples, cluster9_genes]
datExpr_treat35_cluster9 <- datExpr[treat35_samples, cluster9_genes]

# # Calculate intramodular connectivity for for all genes in a cluster/module using only the genes in that module and only the samples for each condition
# # kWithin: sum of connection strengths with other module genes
# #add genes names to connectivity vectors
# kWithin_control <- softConnectivity(datExpr_control_cluster9)
# names(kWithin_control) <- colnames(datExpr_control_cluster9)
# kWithin_treat30 <- softConnectivity(datExpr_treat30_cluster9)
# names(kWithin_treat30) <- colnames(datExpr_treat30_cluster9)
# kWithin_treat35 <- softConnectivity(datExpr_treat35_cluster9)
# names(kWithin_treat35) <- colnames(datExpr_treat35_cluster9)

# Number of top 10% hub genes for each group
#select the top 10% hub genes in each module based on intramodular connectivity (common threshold)
#Top 10% connectivity → ~100 hubs per 1000-gene module
# Use sum of absolute correlations as a proxy for connectivity (ranks genes by their overall co-expression with all other genes in the module)

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
          file = "top10pct_hub_genes_control_cluster9_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat30),
          file = "top10pct_hub_genes_treat30_cluster9_2.csv", row.names = FALSE)
write.csv(data.frame(Gene = top10pct_hubs_treat35),
          file = "top10pct_hub_genes_treat35_cluster9_2.csv", row.names = FALSE)

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
#Montipora_capitata_HIv3___RNAseq.g24713.t1  
cat("Top hub gene (treat30, top 10%):", top_hub_treat30, "\n")
#Montipora_capitata_HIv3___RNAseq.g15653.t1 
cat("Top hub gene (treat35, top 10%):", top_hub_treat35, "\n")
#Montipora_capitata_HIv3___RNAseq.g15049.t1  

####### Compare module preservation across conditions 
####### Highlight Nodes That Change Hub Status
# Highlight Top 10% Hubs
# Hub changing status between control and temps at break-points per (30 and 35) for all clusters

# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(igraph)

# # Select only the top 10 hub genes for each condition (less computationally heavy for plotting)
# top10_genes_control <- top10pct_hubs_control[1:10]
# top10_genes_treat30 <- top10pct_hubs_treat30[1:10]
# top10_genes_treat35 <- top10pct_hubs_treat35[1:10]

# # Subset expression data to these genes
# datExpr_control_top10 <- datExpr_control_cluster9[, top10_genes_control, drop = FALSE]
# datExpr_treat30_top10 <- datExpr_treat30_cluster9[, top10_genes_treat30, drop = FALSE]
# datExpr_treat35_top10 <- datExpr_treat35_cluster9[, top10_genes_treat35, drop = FALSE]

# # Build adjacency/correlation matrices
# adj_control_top10 <- abs(cor(datExpr_control_top10, method = "pearson"))
# adj_treat30_top10 <- abs(cor(datExpr_treat30_top10, method = "pearson"))
# adj_treat35_top10 <- abs(cor(datExpr_treat35_top10, method = "pearson"))

# # Create igraph objects
# g_control_top10 <- graph_from_adjacency_matrix(adj_control_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat30_top10 <- graph_from_adjacency_matrix(adj_treat30_top10, mode = "undirected", weighted = TRUE, diag = FALSE)
# g_treat35_top10 <- graph_from_adjacency_matrix(adj_treat35_top10, mode = "undirected", weighted = TRUE, diag = FALSE)

# # Highlight the top hub gene (by kWithin) in red, others in gray
# V(g_control_top10)$color <- ifelse(names(V(g_control_top10)) == top_hub_control, "red", "gray")
# V(g_treat30_top10)$color <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, "red", "gray")
# V(g_treat35_top10)$color <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, "red", "gray")

# V(g_control_top10)$label <- ifelse(names(V(g_control_top10)) == top_hub_control, top_hub_control, "")
# V(g_treat30_top10)$label <- ifelse(names(V(g_treat30_top10)) == top_hub_treat30, top_hub_treat30, "")
# V(g_treat35_top10)$label <- ifelse(names(V(g_treat35_top10)) == top_hub_treat35, top_hub_treat35, "")

# # Plot networks
# pdf("network_top10genes_cluster9.pdf", width = 12, height = 8)
# par(mfrow = c(1, 3))
# plot(g_control_top10, main = "Control: Top 10 Hubs", vertex.label = V(g_control_top10)$label, vertex.size = 8)
# plot(g_treat30_top10, main = "Treat30: Top 10 Hubs", vertex.label = V(g_treat30_top10)$label, vertex.size = 8)
# plot(g_treat35_top10, main = "Treat35: Top 10 Hubs", vertex.label = V(g_treat35_top10)$label, vertex.size = 8)
# par(mfrow = c(1, 1))
# dev.off()

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

##### combine all bar plots into one figure
library(patchwork)
bar5 + bar6 + bar8 + bar9


######## Make a global alluvial plot showing the overlap of all top 10% hub genes from clusters 5, 6, 8, and 9 across the three temperature groups

# Combine all top 10% hub genes from all clusters and all conditions
all_top_hubs <- unique(c(
  top10pct_hubs_control_5, top10pct_hubs_treat30_5, top10pct_hubs_treat35_5,
  top10pct_hubs_control_6, top10pct_hubs_treat30_6, top10pct_hubs_treat35_6,
  top10pct_hubs_control_8, top10pct_hubs_treat30_8, top10pct_hubs_treat35_8,
  top10pct_hubs_control_9, top10pct_hubs_treat30_9, top10pct_hubs_treat35_9
))

# For each gene, check if it is a top hub in each group (any cluster)
membership_global <- data.frame(
  Gene = all_top_hubs,
  Control = as.integer(all_top_hubs %in% c(top10pct_hubs_control_5, top10pct_hubs_control_6, top10pct_hubs_control_8, top10pct_hubs_control_9)),
  Treat30 = as.integer(all_top_hubs %in% c(top10pct_hubs_treat30_5, top10pct_hubs_treat30_6, top10pct_hubs_treat30_8, top10pct_hubs_treat30_9)),
  Treat35 = as.integer(all_top_hubs %in% c(top10pct_hubs_treat35_5, top10pct_hubs_treat35_6, top10pct_hubs_treat35_8, top10pct_hubs_treat35_9))
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

membership_summary_global$pattern <- interaction(membership_summary_global$Control, membership_summary_global$Treat30, membership_summary_global$Treat35)
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
  labs(title = "Global Overlap of Top 10% Hub Genes Across Clusters and Conditions",
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

cluster5_modules <- c("lightcyan", "orange", "blue", "salmon") #Define the module colors for your cluster of interest
inCluster5 <- moduleColors %in% cluster5_modules #Select genes belonging to these modules
cluster5_genes <- colnames(datExpr)[inCluster5]

cluster6_modules <- c("lightgreen", "sienna3")
inCluster6 <- moduleColors %in% cluster6_modules #Select genes belonging to these modules
cluster6_genes <- colnames(datExpr)[inCluster6]

cluster8_modules <- c("turquoise") 
inCluster8 <- moduleColors %in% cluster8_modules #Select genes belonging to these modules
cluster8_genes <- colnames(datExpr)[inCluster8]

cluster9_modules <- c("darkgreen", "yellow", "bisque4", "grey60") 
inCluster9 <- moduleColors %in% cluster9_modules #Select genes belonging to these modules
cluster9_genes <- colnames(datExpr)[inCluster9]

# Subset expression data
datExpr_control_all <- datExpr[control_samples, cluster5_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster5_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster5_genes]

#or
datExpr_control_all <- datExpr[control_samples, cluster6_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster6_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster6_genes]

#or
datExpr_control_all <- datExpr[control_samples, cluster8_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster8_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster8_genes]

#or
datExpr_control_all <- datExpr[control_samples, cluster9_genes]
datExpr_treat30_all <- datExpr[treat30_samples, cluster9_genes]
datExpr_treat35_all <- datExpr[treat35_samples, cluster9_genes]

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




######## Extract the genes from clusters for GO enrichment analysis with ViSEAGO

# For cluster 5
# cluster5_genes: vector of gene names in cluster 5
# datExpr: samples x all genes (genes as columns)

# Subset datExpr to cluster 5 genes
cluster5_counts <- datExpr[, cluster5_genes, drop = FALSE]

# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster5_counts_t <- as.data.frame(t(cluster5_counts))
cluster5_counts_t$gene_id <- rownames(cluster5_counts_t)

# Move gene_id to first column
cluster5_counts_t <- cluster5_counts_t[, c(ncol(cluster5_counts_t), 1:(ncol(cluster5_counts_t)-1))]

# Save as CSV
write.csv(cluster5_counts_t, file = "cluster5_gene_counts.csv", row.names = FALSE)

# For cluster 6
# cluster6_genes: vector of gene names in cluster 6
# Subset datExpr to cluster 6 genes
cluster6_counts <- datExpr[, cluster6_genes, drop = FALSE]
# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster6_counts_t <- as.data.frame(t(cluster6_counts))
cluster6_counts_t$gene_id <- rownames(cluster6_counts_t)
# Move gene_id to first column
cluster6_counts_t <- cluster6_counts_t[, c(ncol(cluster6_counts_t), 1:(ncol(cluster6_counts_t)-1))]
# Save as CSV
write.csv(cluster6_counts_t, file = "cluster6_gene_counts.csv", row.names = FALSE)

# For cluster 8
# cluster8_genes: vector of gene names in cluster 8
# Subset datExpr to cluster 8 genes
cluster8_counts <- datExpr[, cluster8_genes, drop = FALSE]
# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster8_counts_t <- as.data.frame(t(cluster8_counts))
cluster8_counts_t$gene_id <- rownames(cluster8_counts_t)
# Move gene_id to first column
cluster8_counts_t <- cluster8_counts_t[, c(ncol(cluster8_counts_t), 1:(ncol(cluster8_counts_t)-1))]
# Save as CSV
write.csv(cluster8_counts_t, file = "cluster8_gene_counts.csv", row.names = FALSE)

# For cluster 9
# cluster9_genes: vector of gene names in cluster 9
# Subset datExpr to cluster 9 genes
cluster9_counts <- datExpr[, cluster9_genes, drop = FALSE]
# Transpose so genes are rows, samples are columns (ViSEAGO expects this format)
cluster9_counts_t <- as.data.frame(t(cluster9_counts))
cluster9_counts_t$gene_id <- rownames(cluster9_counts_t)
# Move gene_id to first column
cluster9_counts_t <- cluster9_counts_t[, c(ncol(cluster9_counts_t), 1:(ncol(cluster9_counts_t)-1))]
# Save as CSV
write.csv(cluster9_counts_t, file = "cluster9_gene_counts.csv", row.names = FALSE)