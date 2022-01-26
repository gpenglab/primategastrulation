### species/spatial domain- specificity score
### Guizhong Cui

#imput pacakges
library(readr)
library(data.table)
library(AUCell)
library(GSEABase)
library(tidyr)
library(dplyr)
library(philentropy)

##regulon list---------------------------------
dir = "../scenic.88s.mfas/regulon.genelist/"  
file_list = list.files(path = dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)  #获得csv文件列表     
 
for(i in 1:length(file_list)) {  
  
  df = fread(file = file_list[i])         
  df <- read.csv(file = paste0(file_list[i]),header = T,row.names = 1)
  ls.t[i] <- list(reg=df[,1])
  names(ls.t)[i] <- colnames(df)
}

saveRDS(ls.t,"Mfas.regulon.target.list.rds")

####get normalized expression matrix------------------
file1 <- "../Mfas_mus/Mfas.mus.high.orthology_matrix_Mfas86s-musDomain.111.gene.12538.FPKM.txt"

mm<- read.table(file1, header = T, check.names = F,sep = "\t")
mm$gene <- rownames(mm)

h2 <- read.table("../cross.species.3/human.gastruloid2.22.sample.28846.genes.UMI-corrected.transcripts.matrix.txt",header =T,check.names = F )
# head(h2[,1:5])
h2$gene <- rownames(h2)

s3 <- merge(mm, h2, by="gene")
rownames(s3) <- s3$gene
s3 <- s3[,-1]
# head(s3[,1:4])
write.table(s3,file =paste("Mfas86s.musDomain_hGast2.section.",dim(s3)[1],"gene",dim(s3)[2],"sample.FPKM.normalized.txt",sep = "."),
            quote = F,sep = "\t",row.names = T,col.names = T)

write.csv(colnames(s3),"s3.sample.name.csv")
s3.log <- log(s3+1)
head(s3.log[,1:4])
write.table(s3.log,file =paste("Mfas86s.musDomain_hGast2.section.",dim(s3)[1],"gene",dim(s3)[2],"sample.FPKM.log.normalized.txt",sep = "."),
            quote = F,sep = "\t",row.names = T,col.names = T)

#get scores of mfas regulon by aucell------------------------------------------------
fn <- "Mfas86s.musDomain_hGast2.mfasRegulon.aucell."
s3.log <- as.matrix(s3.log)
geneSets <- ls.t

cells_rankings <- AUCell_buildRankings(s3.log)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*1)
save(cells_AUC, file=paste0(fn,"cells_AUC.RData"))
aucell.out <- getAUC(cells_AUC)#提前AUC matrix
# head(aucell.out[,1:3])

write.table(aucell.out,file =paste("AUCell.homade.",fn,dim(aucell.out)[1],"gene",dim(aucell.out)[2],".txt",sep = "."),
            quote = F,sep = "\t",row.names = T,col.names = T)


##remove extended regulons, and caculate  Jensen–Shannon divergence (JSD)----------------- 
reg <- read.table("../AUCell.homade.noExtend..Mfas86s.musDomain_hGast2.mfasRegulon.aucell..490.gene.133..txt",header = T)
head(reg[,1:4])

rp <- read.csv("../ref.pattern.csv",header = T,row.names = 1)
head(rp[,1:5])

#creat jsd value matrix, row is regulon, col is referenc pattern
jsd.all <- data.frame()

summary(colnames(reg)==colnames(rp)) 
# colnames(reg) <- gsub("-",".",colnames(reg))

for (i in c(1:nrow(reg))) { #nrow(reg)
  # i=1
  rg <- reg[i,]
  # head(rg)
  
  for ( j in c(1:nrow(rp))) { #1:nrow(rp)
    # j=1
    rf <- rp[j,]
    da <- rbind(rf,rg)
    # jsd <- c()
    
    jsd.all[i,j] <- JSD(as.matrix(da),unit = "log2",est.prob = "empirical")
    rownames(jsd.all)[i] <- rownames(rg)
    colnames(jsd.all)[j] <- rownames(rf)
    
    
  }
  
}

write.csv(jsd.all, file = "jsd.value.regulonMfas.3species.pattern.csv" )

#creat p value matrix, row is regulon, col is referenc pattern····------------------------
p.all <- data.frame(matrix(0,nrow=nrow(reg),ncol=nrow(rp)))

summary(colnames(reg)==colnames(rp)) 

for (i in c(1:nrow(reg))) {
  # g=2
  rg <- reg[i,]
  # head(rg)
  
  for ( j in c(1:nrow(rp))) { #1:nrow(rp)
    # j=2
    rf <- rp[j,]
    da <- rbind(rf,rg)
    jsd <- c()
    
    jsd[1] <- JSD(as.matrix(da),unit = "log2",est.prob = "empirical")
    names(jsd)[1] <- rownames(da)[2]
    # head(jsd)
    
    for (s in 1:999) { #1000
      # s=1
      da[2,] <-  sample(rg)
      rownames(da)[2] <- paste0("random",s)
      jsd[s+1] <- JSD(as.matrix(da),unit = "log2",est.prob = "empirical")
      names(jsd)[s+1] <- rownames(da)[2]
    }
    od <- sort(jsd,decreasing=T)
    # head(od[1:5])
    pv <- which(names(od)==rownames(rg))/1000
    
    p.all[i,j] <- pv
    
    rownames(p.all)[i] <- rownames(rg)
    colnames(p.all)[j] <- rownames(rf)
  }
  
}

write.csv(p.all,file = paste0("p.value.1000.",dim(p.all)[1],".regulon.",dim(p.all)[2],".patterns.regulonMfas.3species.refPattern.csv"))


#########regulon specificity score (RSS)--species specificity score (SSS)--------------------------
jsd.all <- read.csv("jsd.value.regulonMfas.3species.pattern.csv",row.names = 1,header = T)
jsd.line <- read.csv("jsd.")
jsd.all <- as.matrix(jsd.all)
sss <- 1- sqrt(jsd.all)

write.csv(sss, "species.specificity.score_SSS.basedJSD.regulonMfas.3species.pattern.csv")

##########significat--------------------------------------------------------------------
##### pvalue of every pattern, commen
p.a <- read.csv("p.value.1000.regulonMfas.mfas.mus.lineage.csv",header = T,row.names = 1)
p.all <- read.csv("p.value.1000.490.regulon.19.patterns.regulonMfas.3species.refPattern.csv",header = T,check.names = F,row.names = 1)
pv <- cbind(p.a,p.all)
pval <- as.matrix(pv[,-ncol(pv)])
p.adj <- 1.001-pval
sigs <- (-log10(p.adj)) #significant score

write.csv(sigs, "significant.score_logPval.regulonMfas.3species.pattern.csv")





