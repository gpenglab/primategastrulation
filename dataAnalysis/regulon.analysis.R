##regulon analysis of M.fas E18-E20 86samples by R--scenic
#Guizhong Cui
#R version 3.6.1 (2019-07-05) 

rm(list=ls())

#import pacakge
library(SingleCellExperiment)
library("AUCell") # ‘1.9.1’
library("RcisTarget") # ‘1.7.1’
library("GENIE3") # ‘1.8.0’
library("SCENIC") # ‘1.1.3’
library(BiocParallel)
library("foreach")
library(dplyr)
library("foreach")
library(Rtsne)


dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-10species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-10species.mc9nr.feather")


#load("data/sceMouseBrain.RData")

fn <-"Mfas.E16.E18.E20.sample.88.gene.17028.Stringtie.TPM"
exp<- read.table("Mfas.E16.E18.E20.sample.88.gene.17028.Stringtie.TPM.txt",
                 header = T,check.names = F)
# head(exp[,1:5])
# exp <-exp[-c(1:2),]
exp$gene <-rownames(exp)

name.c <- read.csv("mfas2human.genename.csv",header = T)

exp.h <- merge(name.c, exp,by="gene")
exp.h <-distinct(exp.h,rename, .keep_all = T)
rownames(exp.h)<-exp.h$rename
exprMat <- exp.h[,-c(1:2)]
# head(exprMat[,1:4])
write.table(exprMat,paste("Mfas2human.E18.E20.87sample.TPM",
                          dim(exprMat)[1],"gene","txt",sep = "."),quote=F,sep="\t")

#Cell info
cellInfo <- read.csv("sample.info.csv")
# head(cellInfo)
rownames(cellInfo)=cellInfo$sample
cellInfo=cellInfo[,-1]
cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)
colnames(cellInfo)[4]
colnames(cellInfo)[4] <-"CellType"
head(cellInfo)

dir.create("int")
## Warning in dir.create("int"): 'int' already exists
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c(
  "Epi"="forestgreen",
  "Mes"="lightslateblue",
  "En"="orange",
  "PS"="blueviolet",
  "Ex"="yellow4",
  "Na"="azure3",
  "Al"="grey"
))

colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
pdf("sector.legend.pdf",width = 7, height = 15)
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))
dev.off()

org="hgnc" # or hgnc, or dmel.mgi for mouse, hgnc for human, or dmel for fly
dbDir="/sibcb1/njinglab3/CGZ/refgenome/scenic/human.ref"# RcisTarget databases location
myDatasetTitle="Mfas_E18-20" # choose a name for your analysis
Options<- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=20)
Options@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
Options@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(Options, file="int/Options.Rds") 

###Co-expression network
#Gene filter/selection
exprMat<-as.matrix(exprMat)
summary(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions=Options,minCountsPerGene=0,minSamples=2) #

interestingGenes <- c("CER1", "TBXT", "OTX2","MESP1","SOX17")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

# This matrix is now ready for the co-expression analysis.
exprMat_filtered <- exprMat
dim(exprMat_filtered) #17019    87

###Correlation
runCorrelation <- function(exprMat_filtered,Options )
{
  corrMat <- cor(t(exprMat_filtered), method="spearman")
  saveRDS(corrMat, file=getIntName(Options, "corrMat"))
}
runCorrelation(exprMat_filtered, Options)

# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
exprMat_filtered=as.matrix(exprMat_filtered)
runGenie3(exprMat_filtered, Options)

Options  <- readRDS("int/Options.Rds")
Options @settings$verbose <- TRUE
Options @settings$nCores <- 20
Options @settings$seed <- 123

runSCENIC_1_coexNetwork2modules(Options)
runSCENIC_2_createRegulons(Options)
Options@settings$nCores <- 1

# scenicOptions@settings$defaultTsne$perpl <-30
runSCENIC_3_scoreCells(Options, exprMat_filtered)
str(Options )

#################################################### get AUC matrix
#The output folder contains several files that provide an overview of the results from each step.
loadInt(Options) 
aucell_regulonAUC <- loadInt(Options , "aucell_regulonAUC")
aucell.out <- getAUC(aucell_regulonAUC)#提前AUC matrix
nr <- nrow(aucell.out)
nc <- ncol(aucell.out)

write.table(aucell.out,paste("aucell_regulonAUC.Mfas2human.E18.E20",
                             nc,nr,"regulon.txt",sep = "."),quote=F,sep="\t")





