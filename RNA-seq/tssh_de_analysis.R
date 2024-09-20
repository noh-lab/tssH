# Analysis of RNA-seq data for "Symbiotic T6SS affects horizontal transmission of Paraburkholderia bonniea among Dictyostelium discoideum amoeba hosts" by Chen et al.
# Please contact Suegene Noh (suegene.noh@colby.edu) with any issues

manuscript_theme = ggplot2::theme_bw() + ggplot2::theme(text=element_text(size=11))

## Set up
setwd("/personal/snoh/tssH")

# load packages:
libraries <- c("DESeq2", "ggplot2", "RColorBrewer","GOstats", "GSEABase", "dplyr", "DEGreport", "gridExtra")
lapply(libraries, require, character.only=T)

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install()

## 1. Input count data 
# prepare to tell a subsequent function where your files are.
directory <- "."
sampleFiles <- grep("qs",list.files(directory),value=TRUE) 
sampleName <- sub("(*).star.count","\\1",sampleFiles) # count tables

# set up the metadata table, that includes the countfile names for the samples
# samples: 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          dicty = c(rep("qs4",8), rep("qs864",8)),
                          time =  c("090","090","360","360","negative","negative","positive","positive", "090","090","360","360","negative","negative","positive","positive"),
                          burk = c("delta", "wildtype", "delta", "wildtype", "delta", "wildtype", "delta", "wildtype", "delta", "wildtype", "delta", "wildtype", "delta", "wildtype", "delta", "wildtype"),
                          stringsAsFactors = T
)
# dicty = host strain, time = time point, burk = delta v wildtype v none
# can view sample table by highlighting the term and pressing command enter

# relevel burk to contrast wildtype (base) vs. delta, and time to contrast negative vs. other
levels(sampleTable$burk)
sampleTable$burk <- factor(sampleTable$burk, levels = c("wildtype", "delta") )
levels(sampleTable$burk)

levels(sampleTable$time)
sampleTable$time <- factor(sampleTable$time, levels = c("negative", "090", "360", "positive") )
levels(sampleTable$time)

# make a dds object
dds.dicty <- DESeqDataSetFromHTSeqCount(sampleTable= sampleTable,
                                        directory= directory,
                                        design= ~burk + dicty)

# keep rows in which at least 3 samples have 10 or more counts for a given gene, eliminate low gene (sample) counts.
keep <- rowSums(assay(dds.dicty) >= 10) >= 3
table(keep)
dds.dicty <- dds.dicty[keep,]


############ ############# WT vs. DELTA ANALYSIS ############# #############
## 2. Differential expression
# run a differential expression analysis based on the GLM model 
dds.dicty <- DESeq(dds.dicty)
resultsNames(dds.dicty)

# specify contrast
res.dicty <- results(dds.dicty, contrast=c("burk","delta","wildtype"), alpha=0.05)

# re-order results 
res.dicty <- res.dicty[order(res.dicty$padj),]

# summary of results
summary(res.dicty)

# look at results table directly
mcols(res.dicty)$description
head(res.dicty)

## 3. Exploratory visualization
# transform the data 
rld.dicty <- rlogTransformation(dds.dicty, blind=F) # takes longer
rld.dicty <- vst(dds.dicty, blind=F) # takes shorter

# PCA plot
plotPCA(rld.dicty, intgroup=c("burk"))  
plotPCA(rld.dicty, intgroup=c("time")) 
plotPCA(rld.dicty, intgroup=c("dicty"))

# check how well our samples correlate with each other, by making a color scheme 
corr.colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
dist.colors <- colorRampPalette(brewer.pal(9, "Oranges"))(255)
color.map <- function(time) { if (time=="negative") "red" else "greenyellow" } # change this based off results names
cond.colors <- unlist(lapply(sampleTable$time, color.map))

# correlation matrix and heatmap. Lightest color = highest correlation. Correlation plots the same, just color bar at the top change
corr.mat <- as.matrix(cor(assay(rld.dicty), method="spearman"))
heatmap(corr.mat, Colv=F, Rowv=NA, scale="none", ColSideColors=cond.colors, col=corr.colors, main = "Time Correlation")

# distance matrix and visualize results in a heatmap. Here, the smallest distance has the lightest colors.
dist.mat <- as.matrix(dist(t(assay(rld.dicty))))
heatmap(dist.mat, Colv=F, Rowv=NA, scale="none", ColSideColors=cond.colors, col=dist.colors, main = "Time Distance")

rm(dds.dicty)


############ ############# TIME SERIES ANALYSIS ############# #############
## change the GLM
dds.dicty2 <- DESeqDataSetFromHTSeqCount(sampleTable= sampleTable,
                                        directory= directory,
                                        design= ~time + dicty)

# keep rows in which at least 3 samples have 10 or more counts for a given gene, eliminate low gene (sample) counts.
keep <- rowSums(assay(dds.dicty2) >= 10) >= 3
table(keep)
dds.dicty2 <- dds.dicty2[keep,]


## 3. Differential expression
# now run
dds.dicty2 <- DESeq(dds.dicty2)
resultsNames(dds.dicty2)

# specify 3 pairwise contrasts
res.dicty.090 <- results(dds.dicty2, contrast=c("time","090","negative"), alpha=0.05)
res.dicty.360 <- results(dds.dicty2, contrast=c("time","360","negative"), alpha=0.05)
res.dicty.pos <- results(dds.dicty2, contrast=c("time","positive","negative"), alpha=0.05)

# re-order these results so that the most significant ones are listed first.
res.dicty.090 <- res.dicty.090[order(res.dicty.090$padj),]
res.dicty.360 <- res.dicty.360[order(res.dicty.360$padj),]
res.dicty.pos <- res.dicty.pos[order(res.dicty.pos$padj),]

# summary of results
summary(res.dicty.090)
summary(res.dicty.360)
summary(res.dicty.pos)

# Get up and downregulated gene lists for venn diagrams
DE090<-data.frame(res.dicty.090)
upreg_090<- row.names(DE090[which(DE090$log2FoldChange>0 & DE090$padj<0.05),])
downreg_090<- row.names(DE090[which(DE090$log2FoldChange<0 & DE090$padj<0.05),])

DE360<-data.frame(res.dicty.360)
upreg_360<- row.names(DE360[which(DE360$log2FoldChange>0 & DE360$padj<0.05),])
downreg_360<- row.names(DE360[which(DE360$log2FoldChange<0 & DE360$padj<0.05),])

DEpos<-data.frame(res.dicty.pos)
upreg_pos<- row.names(DEpos[which(DEpos$log2FoldChange>0 & DEpos$padj<0.05),])
downreg_pos<- row.names(DEpos[which(DEpos$log2FoldChange<0 & DEpos$padj<0.05),])

# create a summary table for venn diagram
col <- brewer.pal(4, "Greys")
venn_time <- as.data.frame(DE090[,1])
row.names(venn_time) <- row.names(DE090)
venn_time$upreg_090 <- row.names(venn_time) %in% upreg_090
venn_time$downreg_090 <- row.names(venn_time) %in% downreg_090
venn_time$upreg_360 <- row.names(venn_time) %in% upreg_360
venn_time$downreg_360 <- row.names(venn_time) %in% downreg_360
venn_time$upreg_pos <- row.names(venn_time) %in% upreg_pos
venn_time$downreg_pos <- row.names(venn_time) %in% downreg_pos
summary(venn_time)

# upreg
limma::vennDiagram(venn_time[,c(2,4,6)], circle.col=col[c(2,4,3)])

# downreg
limma::vennDiagram(venn_time[,c(3,5,7)], circle.col=col[c(2,4,3)])


## 4. Results visualization
rld.dicty2 <- rlogTransformation(dds.dicty2, blind=F) # takes longer
rld.dicty2 <- vst(dds.dicty2, blind=F) # takes shorter

# before shrinkage
plotMA(res.dicty.090, main = "090 MPI v negative before shrinkage")
plotMA(res.dicty.360, main = "360 MPI v negative before shrinkage")
plotMA(res.dicty.pos, main = "time pos v negative before shrinkage")

# after shrinkage:
## 090
lfc090.dicty <- lfcShrink(dds.dicty2, coef="time_090_vs_negative", type="apeglm") #CHANGED THE COEF BASED OFF RESULTSNAMES
plotMA(lfc090.dicty, main = "090 MPI v negative after shrinkage")

## 360
lfc360.dicty <- lfcShrink(dds.dicty2, coef="time_360_vs_negative", type="apeglm") #CHANGED THE COEF BASED OFF RESULTSNAMES
plotMA(lfc360.dicty, main = "360 MPI v negative after shrinkage")

## pos
lfcpos.dicty <- lfcShrink(dds.dicty2, coef="time_positive_vs_negative", type="apeglm") #CHANGED THE COEF BASED OFF RESULTSNAMES
plotMA(lfcpos.dicty, main = "time pos v negative after shrinkage")

# volcano plot
#par(mfrow=c(1,2))
tab090.dicty <- as.data.frame(lfc090.dicty)
tab090.dicty$sig <- tab090.dicty$padj
tab090.dicty$sig <- ifelse(tab090.dicty$padj<=0.05, "FDR 0.05","not sig.")
tab090.dicty$sig <- factor(ifelse(tab090.dicty$padj>0.05 & tab090.dicty$padj<=0.1, "FDR 0.1", tab090.dicty$sig))
ggplot(tab090.dicty, aes(x=log2FoldChange, y=-log10(padj), color=sig)) + geom_point(shape=16, size=1) + scale_color_manual(values=c("red","yellowgreen","grey"),name="") + ylab("-log10(Padj)") + ggtitle("090 MPI vs. negative") + theme(plot.title = element_text(hjust = 0.5)) + manuscript_theme

tab360.dicty <- as.data.frame(lfc360.dicty)
tab360.dicty$sig <- tab360.dicty$padj
tab360.dicty$sig <- ifelse(tab360.dicty$padj<=0.05, "FDR 0.05","not sig.")
tab360.dicty$sig <- factor(ifelse(tab360.dicty$padj>0.05 & tab360.dicty$padj<=0.1, "FDR 0.1", tab360.dicty$sig))
ggplot(tab360.dicty, aes(x=log2FoldChange, y=-log10(padj), color=sig)) + geom_point(shape=16, size=1) + scale_color_manual(values=c("red","yellowgreen","grey"),name="") + ylab("-log10(Padj)") + ggtitle("360 MPI vs. negative") + theme(plot.title = element_text(hjust = 0.5)) + manuscript_theme

tabpos.dicty <- as.data.frame(lfcpos.dicty)
tabpos.dicty$sig <- tabpos.dicty$padj
tabpos.dicty$sig <- ifelse(tabpos.dicty$padj<=0.05, "FDR 0.05","not sig.")
tabpos.dicty$sig <- factor(ifelse(tabpos.dicty$padj>0.05 & tabpos.dicty$padj<=0.1, "FDR 0.1", tabpos.dicty$sig))
ggplot(tabpos.dicty, aes(x=log2FoldChange, y=-log10(padj), color=sig)) + geom_point(shape=16, size=1) + scale_color_manual(values=c("red","yellowgreen","grey"),name="") + ylab("-log10(Padj)") + ggtitle("positive vs. negative") + theme(plot.title = element_text(hjust = 0.5)) + manuscript_theme

# heatmap of expression change
heat.colors <- colorRampPalette(brewer.pal(9, "RdBu"))(255)

color090.map <- function(time) { if (time=="090") "red" else "greenyellow" } # change this based off results names
cond090.colors <- unlist(lapply(sampleTable$time, color090.map))
select090 <- row.names(res.dicty.090)[1:50]
heatmap(assay(rld.dicty2)[select090,], ColSideColors=cond090.colors, cexRow=0.4, col=heat.colors)

color360.map <- function(time) { if (time=="360") "red" else "greenyellow" } # change this based off results names
cond360.colors <- unlist(lapply(sampleTable$time, color360.map))
select360 <- row.names(res.dicty.360)[1:50]
heatmap(assay(rld.dicty2)[select360,], ColSideColors=cond360.colors, cexRow=0.4, col=heat.colors)

colorpos.map <- function(time) { if (time=="positive") "red" else "greenyellow" } # change this based off results names
condpos.colors <- unlist(lapply(sampleTable$time, colorpos.map))
selectpos <- row.names(res.dicty.pos)[1:50]
heatmap(assay(rld.dicty2)[selectpos,], ColSideColors=condpos.colors, cexRow=0.4, col=heat.colors)

# save all results tables to file
write.csv(tab090.dicty, file="tssH.de_090.txt")
write.csv(tab360.dicty, file="tssH.de_360.txt")
write.csv(tabpos.dicty, file="tssH.de_pos.txt")

# export figures
postscript(file=paste("tssh_pca",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=3)
plotPCA(rld.dicty2, intgroup=c("time")) + manuscript_theme
dev.off()

postscript(file=paste("tssh_volc1",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=4, height=3)
ggplot(tab090.dicty, aes(x=log2FoldChange, y=-log10(padj), color=sig)) + geom_point(shape=16, size=1) + scale_color_manual(values=c("red","yellowgreen","grey"),name="") + ylab("-log10(Padj)") + ggtitle("090 MPI vs. negative") + theme(plot.title = element_text(hjust = 0.5)) + manuscript_theme
dev.off()
postscript(file=paste("tssh_volc2",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=4, height=3)
ggplot(tab360.dicty, aes(x=log2FoldChange, y=-log10(padj), color=sig)) + geom_point(shape=16, size=1) + scale_color_manual(values=c("red","yellowgreen","grey"),name="") + ylab("-log10(Padj)") + ggtitle("360 MPI vs. negative") + theme(plot.title = element_text(hjust = 0.5)) + manuscript_theme
dev.off()
postscript(file=paste("tssh_volc3",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=4, height=3)
ggplot(tabpos.dicty, aes(x=log2FoldChange, y=-log10(padj), color=sig)) + geom_point(shape=16, size=1) + scale_color_manual(values=c("red","yellowgreen","grey"),name="") + ylab("-log10(Padj)") + ggtitle("positive vs. negative") + theme(plot.title = element_text(hjust = 0.5)) + manuscript_theme
dev.off()

postscript(file=paste("tssh_vennup",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=6)
limma::vennDiagram(venn_time[,c(2,4,6)], circle.col=col[c(2,4,3)])
dev.off()
postscript(file=paste("tssh_venndown",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=6, height=6)
limma::vennDiagram(venn_time[,c(3,5,7)], circle.col=col[c(2,4,3)])
dev.off()



############ ############# GO ENRICHMENT ANALYSIS using GOStats ############ #############
## 5. enrichment analysis of time series
# read in results tables
df.res.090 <- read.csv("tssH.de_090.txt", row.names = "X")
up.090 <- row.names(subset(df.res.090, padj<=0.05 & log2FoldChange>0))
down.090 <- row.names(subset(df.res.090, padj<=0.05 & log2FoldChange<0))
universe.090 <- row.names(subset(df.res.090, is.na(padj)==FALSE))

df.res.360 <- read.csv("tssH.de_360.txt", row.names = "X")
up.360 <- row.names(subset(df.res.360, padj<=0.05 & log2FoldChange>0))
down.360 <- row.names(subset(df.res.360, padj<=0.05 & log2FoldChange<0))
universe.360 <- row.names(subset(df.res.360, is.na(padj)==FALSE))

df.res.pos <- read.csv("tssH.de_pos.txt", row.names = "X")
up.pos <- row.names(subset(df.res.pos, padj<=0.05 & log2FoldChange>0))
down.pos <- row.names(subset(df.res.pos, padj<=0.05 & log2FoldChange<0))
universe.pos <- row.names(subset(df.res.pos, is.na(padj)==FALSE))


# create GeneSetCollection object
gocurateData <- read.table("gene_association.dictybase.filter.gostats.txt",h=F)
goCurate <- GOFrame(gocurateData, organism="Dictyostelium discoideum")
goCurFrame <- GOAllFrame(goCurate)
gsc <- GeneSetCollection(goCurFrame, setType=GOCollection())

# build a GSEAGOHyperGParams object
bp_up090.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc, geneIds=up.090, universeGeneIds=universe.090, ontology="BP", pvalueCutoff=0.01, conditional=T, testDirection="over")
bp_down090.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc, geneIds=down.090, universeGeneIds=universe.090, ontology="BP", pvalueCutoff=0.01, conditional=T, testDirection="over")

bp_up360.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc, geneIds=up.360, universeGeneIds=universe.360, ontology="BP", pvalueCutoff=0.01, conditional=T, testDirection="over")
bp_down360.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc, geneIds=down.360, universeGeneIds=universe.360, ontology="BP", pvalueCutoff=0.01, conditional=T, testDirection="over")

bp_uppos.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc, geneIds=up.pos, universeGeneIds=universe.pos, ontology="BP", pvalueCutoff=0.01, conditional=T, testDirection="over")
bp_downpos.params <- GSEAGOHyperGParams(name="dicty GSEA", geneSetCollection=gsc, geneIds=down.pos, universeGeneIds=universe.pos, ontology="BP", pvalueCutoff=0.01, conditional=T, testDirection="over")

# run the test itself, and obtain the results based on the p-value cutoff 
bp.up.090 <- hyperGTest(bp_up090.params)
bp.down.090 <- hyperGTest(bp_down090.params)

bp.up.360 <- hyperGTest(bp_up360.params)
bp.down.360 <- hyperGTest(bp_down360.params)

bp.up.pos <- hyperGTest(bp_uppos.params)
bp.down.pos <- hyperGTest(bp_downpos.params)

# results table 
summary(bp.up.090, categorySize=10) # 19 terms
summary(bp.down.090, categorySize=10) # 22 terms

summary(bp.up.360, categorySize=10) #102 terms
summary(bp.down.360, categorySize=10) #56 terms

summary(bp.up.pos, categorySize=10) #42
summary(bp.down.pos, categorySize=10) #49


# format GOSTAT dataframes for tables with gene ids
temp <- subset(as.data.frame(summary(bp.up.090)), Size>=10 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.up.090 <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
temp <- subset(as.data.frame(summary(bp.down.090)), Size>=10 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.down.090 <- cbind(temp[,1:2], round(temp[,3:6], digits=3))

temp <- subset(as.data.frame(summary(bp.up.360)), Size>=10 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.up.360 <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
temp <- subset(as.data.frame(summary(bp.down.360)), Size>=10 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.down.360 <- cbind(temp[,1:2], round(temp[,3:6], digits=3))

temp <- subset(as.data.frame(summary(bp.up.pos)), Size>=10 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.up.pos <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
temp <- subset(as.data.frame(summary(bp.down.pos)), Size>=10 & Size<400)[,c(7,1,4,5,6,2)]
df.bp.down.pos <- cbind(temp[,1:2], round(temp[,3:6], digits=3))
#rm(temp,bp.over.c, mf.over.c, bp.down.c, mf.down.c, bp.up.c, mf.up.c)

# these data can be accessed as below
methods(class="GOHyperGResult")
geneIdsByCategory(bp.over.c)["GO:0006260"]
#$`GO:0006260`
#[1] "DDB_G0272760" "DDB_G0292958"

# get genes in each category: x is the data frame, y is the hypergtest object
get.members2 <- function(x,y) {
  for(i in as.numeric(rownames(x))) {
    print(geneIdsByCategory(y)[i])
  }
}

sink('tssH.up090_genes.txt')
get.members2(df.bp.up.090, bp.up.090)
sink()

sink('tssH.down090_genes.txt')
get.members2(df.bp.down.090, bp.down.090)
sink()

sink('tssH.up360_genes.txt')
get.members2(df.bp.up.360, bp.up.360)
sink()

sink('tssH.down360_genes.txt')
get.members2(df.bp.down.360, bp.down.360)
sink()

sink('tssH.uppos_genes.txt')
get.members2(df.bp.up.pos, bp.up.pos)
sink()

sink('tssH.downpos_genes.txt')
get.members2(df.bp.down.pos, bp.down.pos)
sink()


# manually save results files for go-figure; make sure terms are in the next set of rows
go.up.090 <- summary(bp.up.090, categorySize=10) # 19 terms
go.down.090 <- summary(bp.down.090, categorySize=10) # 22 terms
go.up.360 <- summary(bp.up.360, categorySize=10) #102 terms
go.down.360 <- summary(bp.down.360, categorySize=10) #56 terms
go.up.pos <- summary(bp.up.pos, categorySize=10) #42
go.down.pos <- summary(bp.down.pos, categorySize=10) #49

write.table(go.up.090, file="tssH.GO_up090.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(go.down.090, file="tssH.GO_down090.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(go.up.360, file="tssH.GO_up360.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(go.down.360, file="tssH.GO_down360.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(go.up.pos, file="tssH.GO_uppos.txt", quote=F, sep="\t", row.names=T, col.names=T)
write.table(go.down.pos, file="tssH.GO_downpos.txt", quote=F, sep="\t", row.names=T, col.names=T)

# run go-figure on desktop
#./gofigure -i tssH.GO_up090.txt -o tssH.up090 -j gostats -n bpo -m 50 -q svg -e 40 -f small -n bpo -p viridis
#./gofigure -i tssH.GO_down090.txt -o tssH.down090 -j gostats -n bpo -m 50 -q svg -e 40 -f small -n bpo -p viridis
#./gofigure -i tssH.GO_up360.txt -o tssH.up360 -j gostats -n bpo -m 50 -q svg -e 40 -f small -n bpo -p viridis
#./gofigure -i tssH.GO_down360.txt -o tssH.down360 -j gostats -n bpo -m 50 -q svg -e 40 -f small -n bpo -p viridis
#./gofigure -i tssH.GO_uppos.txt -o tssH.uppos -j gostats -n bpo -m 50 -q svg -e 40 -f small -n bpo -p viridis
#./gofigure -i tssH.GO_downpos.txt -o tssH.downpos -j gostats -n bpo -m 50 -q svg -e 40 -f small -n bpo -p viridis


# venn diagram of go terms
length(unique(sort(c(go.up.090[,1], go.down.090[,1], go.up.360[,1], go.down.360[,1], go.up.pos[,1], go.down.pos[,1]))))

venn_go <- data.frame(matrix(NA, nrow=254, ncol=1))
row.names(venn_go) <- unique(sort(c(go.up.090[,1], go.down.090[,1], go.up.360[,1], go.down.360[,1], go.up.pos[,1], go.down.pos[,1])))
venn_go$go.up.090 <- row.names(venn_go) %in% go.up.090[,1]
venn_go$go.down.090 <- row.names(venn_go) %in% go.down.090[,1]
venn_go$go.up.360 <- row.names(venn_go) %in% go.up.360[,1]
venn_go$go.down.360 <- row.names(venn_go) %in% go.down.360[,1]
venn_go$go.up.pos <- row.names(venn_go) %in% go.up.pos[,1]
venn_go$go.down.pos <- row.names(venn_go) %in% go.down.pos[,1]
summary(venn_go)

# upreg
limma::vennDiagram(venn_go[,c(2,4,6)], circle.col=col[c(2,4,3)])

# downreg
limma::vennDiagram(venn_go[,c(3,5,7)], circle.col=col[c(2,4,3)])

# small overlap between 360 and pos
subset(go.up.360, GOBPID %in% row.names(subset(venn_go, go.up.360=="TRUE" & go.up.pos=="TRUE")))
subset(go.down.360, GOBPID %in% row.names(subset(venn_go, go.down.360=="TRUE" & go.down.pos=="TRUE")))
subset(go.up.360, GOBPID %in% row.names(subset(venn_go, go.up.360=="TRUE" & go.up.090=="TRUE")))
subset(go.down.360, GOBPID %in% row.names(subset(venn_go, go.down.360=="TRUE" & go.down.090=="TRUE")))


############ ############# PHAGOCYTOSIS GENE CLUSTERING ############ #############
## 6. cluster identification for genes of interest
# read in phagocytosis genes
endo <- read.table("endo_phago_path.tsv", sep="\t", h=T, as.is = T)

# read in results files again if necessary
DE090 <- read.csv(file="tssH.de_090.txt", row.names = 1)
DE360 <- read.csv(file="tssH.de_360.txt", row.names = 1)
DEpos <- read.csv(file="tssH.de_pos.txt", row.names = 1)

upreg_090<- row.names(DE090[which(DE090$log2FoldChange>0 & DE090$padj<0.05),])
downreg_090<- row.names(DE090[which(DE090$log2FoldChange<0 & DE090$padj<0.05),])

upreg_360<- row.names(DE360[which(DE360$log2FoldChange>0 & DE360$padj<0.05),])
downreg_360<- row.names(DE360[which(DE360$log2FoldChange<0 & DE360$padj<0.05),])

upreg_pos<- row.names(DEpos[which(DEpos$log2FoldChange>0 & DEpos$padj<0.05),])
downreg_pos<- row.names(DEpos[which(DEpos$log2FoldChange<0 & DEpos$padj<0.05),])

# first filter endo list for DE genes
upendo_090 <- subset(endo, ddb %in% upreg_090) # none
downendo_090 <- subset(endo, ddb %in% downreg_090) #1
downendo_090$sig <- c("downreg_090")
upendo_360 <- subset(endo, ddb %in% upreg_360) #29
upendo_360$sig <- c("upreg_360")
downendo_360 <- subset(endo, ddb %in% downreg_360) #4
downendo_360$sig <- c("downreg_360")
upendo_pos <- subset(endo, ddb %in% upreg_pos) #8
upendo_pos$sig <- c("upreg_pos")
downendo_pos <- subset(endo, ddb %in% downreg_pos) #10
downendo_pos$sig <- c("downreg_pos")
endo_sig <- rbind(downendo_090, upendo_360, downendo_360, upendo_pos, downendo_pos)

temp <- endo_sig[!duplicated(endo_sig$ddb),]
endo_sig <- temp

# for a single gene
d <- plotCounts(dds.dicty2, gene=downendo_090$ddb[1], intgroup="time", returnData=TRUE)
ggplot(d, aes(x=time, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + geom_line(aes(as.numeric(time), count), stat="smooth", method = "loess", span=1, se=F, alpha=0.3) + scale_y_log10(guide = guide_axis(angle = 45)) + scale_x_discrete(guide = guide_axis(angle = 45)) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 

# for a group of genes, one list item per gene
l <- vector("list", nrow(endo_sig))
i <- 1
while(i <= nrow(endo_sig)) {
l[[i]] <- plotCounts(dds.dicty2, gene=endo_sig$ddb[i], intgroup="time", returnData=TRUE)
i <- i+1
}

p <- vector("list", nrow(endo_sig))
i <- 1
while(i <= nrow(endo_sig)) {
p[[i]] <- ggplot(l[[i]], aes(x=time, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + geom_line(aes(as.numeric(time), count), stat="smooth", method = "loess", span=1, se=F, alpha=0.3, formula = 'y ~ x') + scale_y_log10(guide = guide_axis(angle = 45)) + scale_x_discrete(guide = guide_axis(angle = 45)) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) 
i <- i+1
}
grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]])
grid.arrange(p[[17]],p[[18]],p[[19]],p[[20]],p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]],p[[27]],p[[28]],p[[29]],p[[30]],p[[31]],p[[32]])
grid.arrange(p[[33]],p[[34]],p[[35]],p[[36]],p[[37]],p[[38]],p[[39]],p[[40]],p[[41]],p[[42]],p[[43]],p[[44]],p[[45]],p[[46]],p[[47]],p[[48]],p[[49]],p[[50]],p[[51]],p[[52]])
grid.arrange()


# get expression data, from objects above
exp_matrix <- assay(rlog(dds.dicty2))[endo_sig[,2],]
exp_meta <- sampleTable
rownames(exp_meta) <- sampleTable$sampleName


## run DEGreport::degPatterns 
# https://lpantano.github.io/DEGreport/reference/degPatterns.html
patt <- degPatterns(exp_matrix, exp_meta, time = "time", col = "dicty", minc = 1, plot = T)
patt$plot
exp_clust <- patt$df
table(exp_clust$cluster)

exp_clust1 <- subset(exp_clust, cluster=="1")
exp_clust2 <- subset(exp_clust, cluster=="2") 
exp_clust3 <- subset(exp_clust, cluster=="3") 
exp_clust4 <- subset(exp_clust, cluster=="4")
exp_clust5 <- subset(exp_clust, cluster=="5")
exp_clust6 <- subset(exp_clust, cluster=="6")

norm_clust <- patt[["normalized"]] # we want these summarized and normalized values instead of counts -> values, time, genes
temp <- norm_clust[,c(1,3,4,7)]
norm_clust <- temp

# take these clusters and make dataframes as above
d <- subset(norm_clust, genes==exp_clust1$genes[1])
rownames(d) <- d$sampleName
d$gene <- endo_sig[endo_sig$ddb==exp_clust1$genes[1],1]
temp <- d[,c(2,4,5)]
d <- temp 

df_clust1 <- data.frame(value=numeric(), time=character(), gene=factor())
j <- 1
while(j <= nrow(exp_clust1)) {
  d <- subset(norm_clust, genes==exp_clust1$genes[j])
  d$gene <- factor(endo_sig[endo_sig$ddb==exp_clust1$genes[j],1])
  temp1 <- d[,c(2,4,5)]
  temp2 <- rbind(df_clust1, temp1)
  df_clust1 <- temp2
  j <- j+1
}

df_clust2 <- data.frame(count=numeric(), time=character(), gene=factor())
j <- 1
while(j <= nrow(exp_clust2)) {
  d <- subset(norm_clust, genes==exp_clust2$genes[j])
  d$gene <- factor(endo_sig[endo_sig$ddb==exp_clust2$genes[j],1])
  temp1 <- d[,c(2,4,5)]
  temp2 <- rbind(df_clust2, temp1)
  df_clust2 <- temp2
  j <- j+1
}

df_clust3 <- data.frame(count=numeric(), time=character(), gene=factor())
j <- 1
while(j <= nrow(exp_clust3)) {
  d <- subset(norm_clust, genes==exp_clust3$genes[j])
  d$gene <- factor(endo_sig[endo_sig$ddb==exp_clust3$genes[j],1])
  temp1 <- d[,c(2,4,5)]
  temp2 <- rbind(df_clust3, temp1)
  df_clust3 <- temp2
  j <- j+1
}

df_clust4 <- data.frame(count=numeric(), time=character(), gene=factor())
j <- 1
while(j <= nrow(exp_clust4)) {
  d <- subset(norm_clust, genes==exp_clust4$genes[j])
  d$gene <- factor(endo_sig[endo_sig$ddb==exp_clust4$genes[j],1])
  temp1 <- d[,c(2,4,5)]
  temp2 <- rbind(df_clust4, temp1)
  df_clust4 <- temp2
  j <- j+1
}

df_clust5 <- data.frame(count=numeric(), time=character(), gene=factor())
j <- 1
while(j <= nrow(exp_clust5)) {
  d <- subset(norm_clust, genes==exp_clust5$genes[j])
  d$gene <- factor(endo_sig[endo_sig$ddb==exp_clust5$genes[j],1])
  temp1 <- d[,c(2,4,5)]
  temp2 <- rbind(df_clust5, temp1)
  df_clust5 <- temp2
  j <- j+1
}

df_clust6 <- data.frame(count=numeric(), time=character(), gene=factor())
j <- 1
while(j <= nrow(exp_clust6)) {
  d <- subset(norm_clust, genes==exp_clust6$genes[j])
  d$gene <- factor(endo_sig[endo_sig$ddb==exp_clust6$genes[j],1])
  temp1 <- d[,c(2,4,5)]
  temp2 <- rbind(df_clust6, temp1)
  df_clust6 <- temp2
  j <- j+1
}

# combine clusters 1 and 6
temp <- rbind(df_clust1, df_clust6)
df_clust1 <- temp

# move mcln and vacB from cluster 4 into 2
temp1 <- subset(df_clust4, gene=="mcln" | gene=="vacB")
temp2 <- rbind(df_clust2, temp1)
df_clust2 <- temp2
temp1 <- subset(df_clust4, gene!="mcln" & gene!="vacB")
df_clust4 <- temp1

# write to files
write.table(df_clust2, file="tssH.phago_clustA.txt", quote = F, sep="\t", row.names=T, col.names=T)
write.table(df_clust3, file="tssH.phago_clustB.txt", quote = F, sep="\t", row.names=T, col.names=T)
write.table(df_clust1, file="tssH.phago_clustC.txt", quote = F, sep="\t", row.names=T, col.names=T)
write.table(df_clust4, file="tssH.phago_clustD.txt", quote = F, sep="\t", row.names=T, col.names=T)
write.table(df_clust5, file="tssH.phago_clustE.txt", quote = F, sep="\t", row.names=T, col.names=T)


p_clust <- vector("list", max(exp_clust$cluster))
p_clust[[1]] <- ggplot(df_clust1, aes(x=time, y=value, group=gene, shape=gene)) + geom_point(position=position_jitter(w=0.1,h=0)) + geom_line(aes(as.numeric(time), value), stat="smooth", method = "loess", span=1, se=F, alpha=0.3, formula = 'y ~ x') + scale_x_discrete(guide = guide_axis(angle = 45), labels=c("negative" = "000", "090" = "090",  "360" = "360", "positive"="POS")) + scale_shape_manual(values=1:nlevels(df_clust1$gene)) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())  + theme(legend.title=element_blank(), legend.text = element_text(size = 7), legend.key.size = unit(0.5, 'lines')) + guides(shape=guide_legend(nrow=16)) 
p_clust[[2]] <- ggplot(df_clust2, aes(x=time, y=value, group=gene, shape=gene)) + geom_point(position=position_jitter(w=0.1,h=0)) + geom_line(aes(as.numeric(time), value), stat="smooth", method = "loess", span=1, se=F, alpha=0.3, formula = 'y ~ x') + scale_x_discrete(guide = guide_axis(angle = 45), labels=c("negative" = "000", "090" = "090",  "360" = "360", "positive"="POS")) + scale_shape_manual(values=1:nlevels(df_clust2$gene)) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())  + theme(legend.title=element_blank(), legend.text = element_text(size = 7), legend.key.size = unit(0.5, 'lines')) + guides(shape=guide_legend(nrow=16))   
p_clust[[3]] <- ggplot(df_clust3, aes(x=time, y=value, group=gene, shape=gene)) + geom_point(position=position_jitter(w=0.1,h=0)) + geom_line(aes(as.numeric(time), value), stat="smooth", method = "loess", span=1, se=F, alpha=0.3, formula = 'y ~ x') + scale_x_discrete(guide = guide_axis(angle = 45), labels=c("negative" = "000", "090" = "090",  "360" = "360", "positive"="POS")) + scale_shape_manual(values=1:nlevels(df_clust3$gene)) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())  + theme(legend.title=element_blank(), legend.text = element_text(size = 7), legend.key.size = unit(0.5, 'lines')) + guides(shape=guide_legend(nrow=16))  
p_clust[[4]] <- ggplot(df_clust4, aes(x=time, y=value, group=gene, shape=gene)) + geom_point(position=position_jitter(w=0.1,h=0)) + geom_line(aes(as.numeric(time), value), stat="smooth", method = "loess", span=1, se=F, alpha=0.3, formula = 'y ~ x') + scale_x_discrete(guide = guide_axis(angle = 45), labels=c("negative" = "000", "090" = "090",  "360" = "360", "positive"="POS")) + scale_shape_manual(values=1:nlevels(df_clust4$gene)) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())  + theme(legend.title=element_blank(), legend.text = element_text(size = 7), legend.key.size = unit(0.5, 'lines')) + guides(shape=guide_legend(nrow=16))  
p_clust[[5]] <- ggplot(df_clust5, aes(x=time, y=value, group=gene, shape=gene)) + geom_point(position=position_jitter(w=0.1,h=0)) + geom_line(aes(as.numeric(time), value), stat="smooth", method = "loess", span=1, se=F, alpha=0.3, formula = 'y ~ x') + scale_x_discrete(guide = guide_axis(angle = 45), labels=c("negative" = "000", "090" = "090",  "360" = "360", "positive"="POS")) + scale_shape_manual(values=1:nlevels(df_clust5$gene)) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + theme(legend.title=element_blank(), legend.text = element_text(size = 7), legend.key.size = unit(0.5, 'lines')) + guides(shape=guide_legend(nrow=16))  


grid.arrange(p_clust[[2]], p_clust[[3]], p_clust[[1]], p_clust[[4]], p_clust[[5]])

cairo_ps(filename=paste("/personal/snoh/tssH/tssh_cluster",format(Sys.time(),"%Y%m%d"),"eps",sep="."),onefile=F, width=7, height=6, fallback_resolution = 300, family="sans")
grid.arrange(p_clust[[2]], p_clust[[3]], p_clust[[1]], p_clust[[4]], p_clust[[5]])
dev.off()
