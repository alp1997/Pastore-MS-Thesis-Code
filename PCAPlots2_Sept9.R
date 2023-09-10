#####GeneExpressionAnalysis just for making PCA and heatmap plots
rm(list = ls())
######################################################
# Point R to count data, set an output file prefix 
######################################################
# Load the libraries we'll need in the following code:
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")
library("ashr")
library("goseq")
library("biomaRt")
######
# create an object with the directory containing your counts:
directory <- "/Users/amandapastore/Library/Mobile Documents/com~apple~CloudDocs/Academics Folder/M.S. Research/Data and Info/counts/"

sampleFiles <- list.files(directory, pattern = ".*counts$")

### Remove the samples that don't have enough reads (there's 3 of them)
exclude<- c("Sample25_258_19_3_S13_R1_001.sam.sam.counts", "Sample48_pool6_12_6_S30_R1_001.sam.sam.counts", "Sample4R_Pool6_12_1_S36_R1_001.sam.sam.counts")
sampleFiles<- sampleFiles[!sampleFiles%in%exclude]
length(sampleFiles)
##

######################################################
# Read the count data into R along with treatment information
######################################################
# load in metadata table
meta <- read.csv("/Users/amandapastore/Library/Mobile Documents/com~apple~CloudDocs/Academics Folder/M.S. Research/Data and Info/CSV Sheets for R/rnainfo_full.csv") 


meta$trial =as.factor(meta$trial)
##should be 33 now
meta$sample[meta$sample%in%exclude]

meta2<- meta[!meta$sample%in%exclude, ]
nrow(meta2)
# ensure that sampleFiles and metadata table columns are in the same order
str(meta2)
str(sampleFiles)

all(str_remove(sampleFiles, ".counts") == meta2[,1])
###FALSE, so put them in the correct order
sampleFiles <- sort(sampleFiles)
meta2 <- meta2[order(meta2$sample),]
all(sampleFiles == meta2[,1])



length(sampleFiles)
length(meta2)
head(sampleFiles)
head(meta2)
# now create a data frame with sample names, file names and treatment information. 
sampleTable <- data.frame(
  sampleName = meta2$sample,
  fileName = sampleFiles,
  clonal.line = meta2$clonal.line,
  predator = meta2$predator,
  trial = meta2$trial,
  bin = meta2$bin
)
# look at the data frame to ensure it is what you expect:
sampleTable
## Creation of sampletable2 to make sample table 3 below
sampleTable2 <- merge(sampleTable, residmerge[,c("clonal.line", "plasticity")], all.x= TRUE, all.y= FALSE)
str(sampleTable2)
head(sampleTable2)

### DO NOT RUN THIS LINE TWICE!! IT WILL KEEP REARRANGING COLUMNS!!!
sampleTable3 <- sampleTable2[,c(2, 3, 1, 4:ncol(sampleTable2))]

str(sampleTable3)
################################################
################################################
########################################################################
sampleTable3$predator<- as.factor(sampleTable3$predator)
sampleTable3$clonal.line<- as.factor(sampleTable3$clonal.line)
str(sampleTable3)
# create the DESeq data objects
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable3, 
  directory = directory, 
  design = ~predator*clonal.line
)

######################################################
# Filter out genes with very low expression
######################################################
# sum counts for each gene across samples
sumcounts <- rowSums(counts(ddsHTSeq))
# get genes with summed counts greater than 20
keep <- sumcounts > 20

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq <- ddsHTSeq[keep,]

######################################################
# Run the statistical analysis
######################################################

dds <- DESeq(ddsHTSeq)

str(dds)
######################################################
# Which results do we want?
resobjects<- resultsNames(dds)
######################################################
# we are interested in:
# 1) genes responding to predator exposure in all lines (control vs exposed)

res2 <- results(dds, contrast=c("predator","P","N"))
res3<- results(dds, contrast=list(c(resobjects[2], resobjects[8])))
res4<- results(dds, contrast=list(c(resobjects[2], resobjects[9])))
res5<- results(dds, contrast=list(c(resobjects[2], resobjects[10])))
res6<- results(dds, contrast=list(c(resobjects[2], resobjects[11])))
res7<- results(dds, contrast=list(c(resobjects[2], resobjects[12])))
resall<- rbind(res2, res3, res4, res5, res6, res7)
######################################################
# Quickly summarize results
######################################################

# get a quick summary of the table
summary(res2)
summary(res3)
summary(res4)
summary(res5)
summary(res6)
summary(res7)
summary(resall)

# PCA plot
# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)
str(vsd)

temp<- as.data.frame(resall) %>% drop_na(padj) %>% filter(padj<0.1)
temp2<- as.data.frame(resall) %>% drop_na(padj) %>% filter(padj<0.001)
nrow(temp)
##Temp contains the info needed for plot
str(temp)
#write.csv(temp, file="table6_MSThesis.csv", quote=FALSE)
######
genes <- rownames(temp)
genes2 <- rownames(temp2)
##global search and replace to say to remove everything after the period so that gene names with decimal are kept
realtest<- gsub("\\..*", "", genes2)
length(unique(genes))
length(unique(genes2))
length(genes2)
# heatmap of DE genes
##pull list for 2 thru 6 as well to see the clonal line comparison, instead of contrast list it will jsut be 
#####PCA plots with genes now
vsd2 <- vsd[names(vsd)%in%genes]
vsd3 <- vsd[names(vsd)%in%realtest]
table5<-assay(vsd3)
str(table5)
tabletemp<- names(vsd)[names(vsd)%in%genes2]
!genes2%in%names(vsd)
genes2
test2%in%names(vsd)
test1%in%names(vsd)
###Try with clonal.line and predator AND with pool, clonal line, treatment
dat <-plotPCA(vsd2, returnData=TRUE, intgroup=c("clonal.line","predator"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=clonal.line, shape=predator))+ geom_point(size=3) +
xlab(paste("PC1: ", round(attr(dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(dat,"percentVar")[2],2)*100, "% variation explained", sep="")) +
  labs(title="PCA Plot for Clonal Lines") + labs(fill="Clonal Line")
plot(p)
##another better plot:
p <- ggplot(dat,aes(x=PC1,y=PC2,col=predator, shape=clonal.line))+ geom_point(size=3) +
  xlab(paste("PC1: ", round(attr(dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(dat,"percentVar")[2],2)*100, "% variation explained", sep="")) +
  labs(title="PCA Plot for each Treatment") + labs(fill="Clonal Line")
plot(p)

# create a metadata data frame to add to the heatmaps
df <- data.frame(colData(dds)[,c("plasticity", "clonal.line", "predator")])
rownames(df) <- colnames(dds)
colnames(df) <- c("plasticity", "clonal.line", "predator")

#Make the heatmap based on an analysis that includes predator and clonal line, change the rld 
# use regularized log-scaled counts
pheatmap(
  assay(vsd2), 
  cluster_rows=TRUE, 
  show_rownames=FALSE,
  cluster_cols=TRUE,
  annotation_col=df
)


##############

# plot counts for individual genes

plotCounts(dds, gene=lfcorder[9], intgroup=c("clonal.line","treatment"))
############
#Volcano plot

res_shrinkall <- lfcShrink(dds,type="ashr",contrast=list(c("predator_P_vs_N")))

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrinkall$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrinkall$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.5,
     col=(log_padj > 10)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes") +
     title(main="Volcano Plot")

