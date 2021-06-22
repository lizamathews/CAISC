rm(list=ls())
load("exp.RData")
geneLoc<-read.table("mart_export_mouse.txt", header = TRUE)
genes<-table(geneLoc$MGIsymbol)
genes<-genes[genes==1]
genesWithExp<-intersect(names(genes), rownames(exp))  

####
geneLoc<-geneLoc[geneLoc$MGIsymbol %in% genesWithExp,];
geneLoc<-geneLoc[order(geneLoc[,"Chromosome"], geneLoc[,"Genestart"]),]
rownames(geneLoc)<-geneLoc$MGIsymbol
zz<-file("annoationGene.txt", "w")
for(gene in 1:dim(geneLoc)[1]){
  cat(geneLoc[gene, 1], file=zz);cat("\tchr", file = zz);
  cat(geneLoc[gene, 2], file=zz);cat("\t", file = zz);
  cat(geneLoc[gene, 3], file=zz);cat("\t", file = zz);
  cat(geneLoc[gene, 4], file=zz);cat("\n", file = zz);
}
close(zz)

####
exp<-exp[geneLoc$MGIsymbol,]
cells<-colnames(exp)
##Get cellInfo
load("cellInfo.RData")
rownames(cellInfo)<-cellInfo$cellname
cellInfo<-cellInfo[colnames(exp),]
cell_group4<-cellInfo[cellInfo[, "memb_pred"]==4,]

zz<-file("annoationCell.txt", "w")
for(cell in rownames(cellInfo)){
  if(cellInfo[cell, "memb_pred"]==4){
    cat(cell, file=zz);cat("\tTumor\n", file = zz);
  }else{
    cat(cell, file=zz);cat("\tTumor\n", file = zz);
  }
}
###Add group 4 as reference
for(cell in rownames(cell_group4)){
  cat(paste0("ref_", cell), file=zz);cat("\tHealthy\n", file = zz);
}
close(zz)
####
###group 1: half chr2 dup, half chr4 mono
###group 2: half chr2 dup
###group 3: half chr4 mono 
#quantile(geneLoc$Genestart[geneLoc$Chromosome==1])
#0%       25%       50%       75%      100% 
#3114224  52055421  91457511 144534992 181991129 
#> quantile(geneLoc$Genestart[geneLoc$Chromosome==4])
#0%       25%       50%       75%      100% 
#3167320  60438038 116714729 135415493 156331103 
cell_group1<-cellInfo[cellInfo$memb_pred==1,]; cell_group1<-cell_group1[1:15,]
cell_group2<-cellInfo[cellInfo$memb_pred==2,]; cell_group2<-cell_group2[1:15,]
cell_group3<-cellInfo[cellInfo$memb_pred==3,]; cell_group3<-cell_group3[1:15,]
cell_group4<-cellInfo[cellInfo$memb_pred==4,];

zz<-file("expressionmatrix.txt", "w")
zz1<-file("log_aneuploidy.txt", "w")
cat("gene", file=zz)
for(cell in colnames(exp)){
  cat("\t", file = zz);cat(cell, file=zz);
}
for(cell in rownames(cell_group4)){
  cat(paste0("\tref_", cell), file=zz);
}
cat("\n", file=zz)
iii<-0

for(gene in rownames(exp)){
  cat(geneLoc[gene, "Chromosome"]);cat("\t");cat(iii);cat("\n");iii<-iii+1;
  cat(gene, file=zz)
  for(cell in colnames(exp)){
    if(cell %in% rownames(cell_group1)){
      if(iii==1){ cat(cell, file=zz1);cat("\tgroup1\n", file=zz1)}
      if(geneLoc[gene, "Chromosome"]=="1" & geneLoc[gene, "Genestart"] < 164534992){##duplication
        cat("\t", file = zz);cat(exp[gene, cell]*4, file=zz); ###duplication
        cat(cell, file=zz1); cat("\t", file=zz1);cat(geneLoc[gene, "Chromosome"], file=zz1); cat("\t", file=zz1); cat(geneLoc[gene, "Genestart"], file=zz1); cat("\tI am here 68\n", file=zz1)
      }else{
        if(geneLoc[gene, "Chromosome"] %in% c("4","6") & geneLoc[gene, "Genestart"] < 156714729){##mono
          cat("\t", file = zz);cat(round(exp[gene, cell]/4), file=zz); ###mono
          #cat(geneLoc)
          cat(cell, file=zz1); cat("\t", file=zz1);cat(geneLoc[gene, "Chromosome"], file=zz1); cat("\t", file=zz1); cat(geneLoc[gene, "Genestart"], file=zz1); cat("I am here 73\n", file=zz1)
        }else{
          cat("\t", file = zz);cat(exp[gene, cell], file=zz); ###diploid
        }
      }
    }else{#if(cell %in% rownames(cell_group1)){
      if (cell %in% rownames(cell_group2)){
        if(iii==1){ cat(cell, file=zz1);cat("\tgroup2\n", file=zz1)}
        if(geneLoc[gene, "Chromosome"]=="1" & geneLoc[gene, "Genestart"] < 164534992){##duplication
          cat("\t", file = zz);cat(exp[gene, cell]*4, file=zz); ###duplication
          cat(cell, file=zz1); cat("\t", file=zz1);cat(geneLoc[gene, "Genestart"], file=zz1); cat("\t", file=zz1); cat(geneLoc[gene, "Genestart"], file=zz1); cat("I am here 82\n", file=zz1)
        }else{
          cat("\t", file = zz);cat(exp[gene, cell], file=zz); ###diploid
        }
      }else{##if (cell %in% rownames(cell_group2)){
        if(cell %in% rownames(cell_group3)){
          if(iii==1){ cat(cell, file=zz1);cat("\tgroup3\n", file=zz1)}
          if(geneLoc[gene, "Chromosome"] %in% c("4","6") & geneLoc[gene, "Genestart"] < 156714729){##mono
            cat("\t", file = zz);cat(round(exp[gene, cell]/4), file=zz); ###mono
            cat(cell, file=zz1); cat("\t", file=zz1);cat(geneLoc[gene, "Chromosome"], file=zz1); cat("\t", file=zz1); cat(geneLoc[gene, "Genestart"], file=zz1); cat("I am here 90\n", file=zz1)
          }else{
            cat("\t", file = zz);cat(exp[gene, cell], file=zz); ###diploid
          }
        }else{#if(cell %in% rownames(cell_group3)){
          if(iii==1){ cat(cell, file=zz1);cat("\tgroupn\n", file=zz1)}
          if(geneLoc[gene, "Chromosome"] %in% c("5", "9")){##dup
            if(cell %in% rownames(cell_group4)){
              cat("\t", file = zz);cat(round(exp[gene, cell]*4), file=zz); ###mono
            }else{
              cat("\t", file = zz);cat(exp[gene, cell], file=zz); ###diploid
            }
          }else{
            cat("\t", file = zz);cat(exp[gene, cell], file=zz); ###diploid
          }
        }
      }
    }
  }
  for(cell in rownames(cell_group4)){
    cat("\t", file = zz);cat(exp[gene, cell], file=zz); ###diploid
  }
  cat("\n", file=zz)
}
close(zz)
close(zz1)
