library(GSVA)
library(GEOquery)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(scales)
library(reshape2)
library(mogene10sttranscriptcluster.db)
library(mogene20sttranscriptcluster.db)
library(mouse4302.db)
library(cowplot)
library(emmeans)
library(nlme)

#########################################################################################################
## Imports Hacohen States

DEGList = readRDS("2019_07_06 T Cell State DEG List Unique States.rds")
DEG_Signatures = list()

for(i in 1:3){
  DEG_Signatures = c(DEG_Signatures,  Up = list(rownames(DEGList[[i]])[DEGList[[i]]$avg_logFC>0]))
  DEG_Signatures = c(DEG_Signatures,  Down = list(rownames(DEGList[[i]])[DEGList[[i]]$avg_logFC<0]))
}

names(DEG_Signatures) = paste(rep(c("Cluster1","Cluster3","Cluster5"),each = 2),names(DEG_Signatures),sep = "_")

Human_MouseConversion <- read_csv("Human_MouseConversion.xlsx.csv")
DEG_SignaturesConverted = list()
for(i in 1:6) DEG_SignaturesConverted = c(DEG_SignaturesConverted,list(toupper(merge(data.frame(DEG_Signatures[[i]]),Human_MouseConversion,by.x = 1, by.y = 2)[,2])))
names(DEG_SignaturesConverted) = names(DEG_Signatures)

Hacohen_Signatures = list(Cluster1 = list(Up = DEG_SignaturesConverted[[1]],Down = DEG_SignaturesConverted[[2]]),
                          Cluster3 = list(Up = DEG_SignaturesConverted[[3]],Down = DEG_SignaturesConverted[[4]]),
                          Cluster5 = list(Up = DEG_SignaturesConverted[[5]],Down = DEG_SignaturesConverted[[6]]))


rm(DEG_Signatures,DEGList,Human_MouseConversion,DEG_SignaturesConverted)

#########################################################################################################
## Functions for Analyzing Data

SignatureScore <- function(InputGeneMatrix,TestSignatures){
  GeneMatrixFiltered = merge(data.frame(TestSignatures),data.frame(InputGeneMatrix),by.x = 1,by.y = 0)
  
  GeneMatrixFiltered_Zscore = data.matrix(GeneMatrixFiltered[,2:ncol(GeneMatrixFiltered)])
  for(i in 1:nrow(GeneMatrixFiltered_Zscore)) GeneMatrixFiltered_Zscore[i,] = scale(GeneMatrixFiltered_Zscore[i,])
  
  MeanZ_Score = numeric(ncol(GeneMatrixFiltered_Zscore))
  
  ## Averages Z score of all genes in a signature (Using Geom mean will return an error for negative Z scores!!!)
  for(i in 1:ncol(GeneMatrixFiltered_Zscore)) MeanZ_Score[i] = mean(GeneMatrixFiltered_Zscore[,i])
  names(MeanZ_Score) = colnames(GeneMatrixFiltered_Zscore)
  
  MeanZ_Score
}

SignatureScore_UpMinusDown <- function(InputGeneMatrix,UpDownList){
  UpSignature = SignatureScore(InputGeneMatrix,UpDownList[["Up"]])
  DnSignature = SignatureScore(InputGeneMatrix,UpDownList[["Down"]])
  
  NetSignature = UpSignature - DnSignature
  NetSignature
}

SignatureScore_RowInd <- function(InputGeneMatrix,RowInds){
  GeneMatrixFiltered_Zscore = InputGeneMatrix[RowInds,]
  
  for(i in 1:nrow(GeneMatrixFiltered_Zscore)) GeneMatrixFiltered_Zscore[i,] = scale(GeneMatrixFiltered_Zscore[i,])
  
  MeanZ_Score = numeric(ncol(GeneMatrixFiltered_Zscore))
  
  ## Averages Z score of all genes in a signature (Using Geom mean will return an error for negative Z scores!!!)
  for(i in 1:ncol(GeneMatrixFiltered_Zscore)) MeanZ_Score[i] = mean(GeneMatrixFiltered_Zscore[,i])
  names(MeanZ_Score) = colnames(GeneMatrixFiltered_Zscore)
  
  MeanZ_Score
}

SignatureScore_Random <- function(InputGeneMatrix,UpDownList,nperms = 50){
  NetSignatures = matrix(nrow = nperms,ncol = ncol(InputGeneMatrix))
  colnames(NetSignatures) = colnames(InputGeneMatrix)
  for(i in 1:nperms){
    UpGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Up"]]))
    DnGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Down"]]))
    
    UpSignature = SignatureScore_RowInd(InputGeneMatrix,UpGenesRandom)
    DnSignature = SignatureScore_RowInd(InputGeneMatrix,DnGenesRandom)
    
    NetRandom = UpSignature - DnSignature
    NetSignatures[i,] = NetRandom  
  }
  colMeans(NetSignatures)
}

SignatureScorePlot <- function(InputGeneMatrix,SigList = Hacohen_Signatures){
  SignatureScores = matrix(nrow = length(SigList),ncol = ncol(InputGeneMatrix))
  colnames(SignatureScores) = colnames(InputGeneMatrix);rownames(SignatureScores) = names(SigList)
  RandomScores = SignatureScores
  for(i in 1:length(SigList)){
    SignatureScores[i,] = SignatureScore_UpMinusDown(InputGeneMatrix,SigList[[i]])
    RandomScores[i,] = SignatureScore_Random(InputGeneMatrix,SigList[[i]],nperms = 100)
  }
  
  SignatureScores = melt(SignatureScores,value.name = "ZScore", varnames=c('Cluster', 'Sample'))
  RandomScores = melt(RandomScores,value.name = "ZScore", varnames=c('Cluster', 'Sample'))
  SignatureScores$ExpBase = "Sig"
  RandomScores$ExpBase = "Base"
  
  Signature_Plot = rbind(SignatureScores,RandomScores)
  Signature_Plot$SampleName = sapply(strsplit(as.character(Signature_Plot$Sample),fixed = TRUE,split = "_"),function(x) paste(x[-length(x)],collapse = "_"))
  Signature_Plot$SampleName = factor(Signature_Plot$SampleName,levels = unique(Signature_Plot$SampleName))
  Signature_Plot$ExpBase = factor(Signature_Plot$ExpBase)
  
  ## LME Test for Significance - compares each point to the point before
  
  LME_Sig = Signature_Plot[Signature_Plot$ExpBase == "Sig",]
  LME_Sig$Sample = factor(sapply(strsplit(as.character(LME_Sig$Sample),fixed = TRUE,split = "_"),function(x) x[length(x)]))
  
  LME_Sig_PVals = matrix(nrow = length(SigList),ncol = length(levels(LME_Sig$SampleName))-1)
  
  for(Sigs in 1:length(SigList)){
    LME_Sigs = LME_Sig[LME_Sig$Cluster == names(SigList)[Sigs],]
    LME_Sig_Results = lme(ZScore~SampleName,random=~1|Sample,data = LME_Sigs)
    LME_Sig_Results = summary(emmeans(LME_Sig_Results,"pairwise"~SampleName,adjust = "none"))$contrasts
    
    
    PVals = LME_Sig_Results[nrow(LME_Sig_Results),6]
    SampleName = as.character(LME_Sig_Results[nrow(LME_Sig_Results),1])
    Index = nrow(LME_Sig_Results)
    
    for(i in 2:(length(unique(LME_Sig$SampleName))-1)){
      Index = Index - i
      PVals = c(PVals,LME_Sig_Results[Index,6])
      SampleName = c(SampleName,as.character(LME_Sig_Results[Index,1]))
    }
    
    PVals = rev(PVals)
    SampleName = rev(SampleName)
    
    LME_Sig_PVals[Sigs,] = PVals
  }
  
  colnames(LME_Sig_PVals) = SampleName
  rownames(LME_Sig_PVals) = names(SigList)
  
  print(LME_Sig_PVals)
  
  ## Plots Individual + Combined
  pal = scales::hue_pal()(5)[c(1,3,5)]
  
  Plots = list()
  CombinedPlot = NULL
  ## Individual Plots
  for(i in 1:length(names(SigList))){
    Sample = factor(levels(Signature_Plot$SampleName),levels(Signature_Plot$SampleName))
    ExpMean = numeric(length(names(SigList)))
    ExpSD = ExpMean; BaseMean = ExpMean; BaseSD = ExpMean
    for(j in 1:length(Sample)){
      Exp = Signature_Plot$ZScore[Signature_Plot$ExpBase=="Sig"&
                                    Signature_Plot$SampleName==Sample[j]&
                                    Signature_Plot$Cluster==names(SigList)[i]]
      Base = Signature_Plot$ZScore[Signature_Plot$ExpBase=="Base"&
                                     Signature_Plot$SampleName==Sample[j]&
                                     Signature_Plot$Cluster==names(SigList)[i]]
      
      ExpMean[j] = mean(Exp); ExpSD[j] = sd(Exp)
      BaseMean[j] = mean(Base); BaseSD[j] = sd(Base)
    }
    
    IndPlot = data.frame(Sample,ExpMean,ExpUpper = ExpMean + ExpSD,ExpLower = ExpMean - ExpSD,
                         BaseMean,BaseUpper = BaseMean + BaseSD,BaseLower = BaseMean - BaseSD)
    CombinedPlot = rbind(CombinedPlot,IndPlot)
    
    Plots[[i]] = ggplot(data = IndPlot,aes(x = Sample,y = BaseMean,group = 1)) +
      geom_ribbon(aes(ymin = BaseLower,ymax = BaseUpper),fill = "grey") +
      geom_line(color = "black") + 
      geom_line(aes(x = Sample,y = ExpMean),color = pal[i]) +
      geom_errorbar(aes(ymin = ExpLower,ymax = ExpUpper),width=0.4, size=0.5, color=pal[i]) +
      geom_point(aes(x = Sample,y = ExpMean),fill = pal[i],shape = 21,size = 4) +
      theme_cowplot() +
      theme(axis.text = element_blank(),axis.title = element_blank()) + 
      ylim(c(-1.9,1.9))
  }  
  
  CombinedPlot$Cluster = rep(names(SigList),each = length(levels(CombinedPlot$Sample)))
  
  
  Plots[[4]] = ggplot(data = CombinedPlot,aes(Sample,ExpMean,group = Cluster)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    geom_errorbar(aes(color = Cluster,ymin = ExpLower,ymax = ExpUpper),width=0.4,size = 0.5) +
    geom_line(aes(color = Cluster)) + 
    geom_point(aes(fill = Cluster),shape = 21,size = 4) + 
    theme_cowplot() +
    theme(legend.position = "none") +
    theme(axis.title = element_blank(),axis.text = element_blank()) + 
    ylim(c(-1.9,1.9)) 
  
  for(i in 1:2) Plots[[i]] = Plots[[i]] + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank())
  
  
  CombinedPlot = plot_grid(plotlist = Plots,align = "hv",nrow = 4)
  CombinedPlot
}

SignatureScorePlot_OneSampleRep <- function(InputGeneMatrix,SigList = Hacohen_Signatures){
  SignatureScores = matrix(nrow = length(SigList),ncol = ncol(InputGeneMatrix))
  colnames(SignatureScores) = colnames(InputGeneMatrix);rownames(SignatureScores) = names(SigList)
  RandomScores = SignatureScores
  for(i in 1:length(SigList)){
    SignatureScores[i,] = SignatureScore_UpMinusDown(InputGeneMatrix,SigList[[i]])
    RandomScores[i,] = SignatureScore_Random(InputGeneMatrix,SigList[[i]],nperms = 100)
  }
  
  SignatureScores = melt(SignatureScores,value.name = "ZScore", varnames=c('Cluster', 'Sample'))
  RandomScores = melt(RandomScores,value.name = "ZScore", varnames=c('Cluster', 'Sample'))
  SignatureScores$ExpBase = "Sig"
  RandomScores$ExpBase = "Base"
  
  
  Signature_Plot = rbind(SignatureScores,RandomScores)
  
  pal = scales::hue_pal()(5)[c(1,3,5)]
  
  Plots = list()
  ## Individual Plots
  for(i in 1:length(names(SigList))){
    IndPlot = Signature_Plot[Signature_Plot$Cluster==names(SigList)[i],]
    
    Plots[[i]] = ggplot(data = IndPlot[IndPlot$ExpBase=="Base",], aes(x = Sample, y = ZScore,group = 1)) + 
      geom_line(color = "black") + 
      geom_line(data = IndPlot[IndPlot$ExpBase=="Sig"&IndPlot$Sample!="CM"&IndPlot$Sample!="EM",],color = pal[i]) + 
      geom_point(data = IndPlot[IndPlot$ExpBase=="Sig",],fill = pal[i],shape = 21,size = 4) + 
      theme_cowplot() + 
      theme(axis.text = element_blank(),axis.title = element_blank()) + 
      ylim(-1.9,1.9)
  }
  
  Plots[[4]] = ggplot(data = Signature_Plot[Signature_Plot$ExpBase=="Sig",],aes(Sample,ZScore,group = Cluster)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    geom_point(aes(fill = Cluster),shape = 21,size = 1) + 
    geom_line(data = Signature_Plot[Signature_Plot$ExpBase == "Sig"&Signature_Plot$Sample!="CM"&Signature_Plot$Sample!="EM",],aes(color = Cluster)) + 
    geom_point(aes(fill = Cluster),shape = 21,size = 4) + 
    theme_cowplot() +
    theme(legend.position = "none") +
    theme(axis.title = element_blank(),axis.text = element_blank()) + 
    ylim(c(-1.9,1.9)) 
  
  for(i in 1:2) Plots[[i]] = Plots[[i]] + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank())
  
  CombinedPlot = plot_grid(plotlist = Plots,align = "hv",nrow = 4)
  CombinedPlot
}

SlopePlot_OneSampleRep <- function(InputGeneMatrix,SigList =Hacohen_Signatures){
  SignatureScores = matrix(nrow = length(SigList),ncol = ncol(InputGeneMatrix))
  colnames(SignatureScores) = colnames(InputGeneMatrix);rownames(SignatureScores) = names(SigList)
  RandomScores = SignatureScores
  for(i in 1:length(SigList)){
    SignatureScores[i,] = SignatureScore_UpMinusDown(InputGeneMatrix,SigList[[i]])
    RandomScores[i,] = SignatureScore_Random(InputGeneMatrix,SigList[[i]],nperms = 100)
  }
  
  SignatureScores = melt(SignatureScores,value.name = "ZScore", varnames=c('Cluster', 'Sample'))
  RandomScores = melt(RandomScores,value.name = "ZScore", varnames=c('Cluster', 'Sample'))
  SignatureScores$ExpBase = "Sig"
  RandomScores$ExpBase = "Base"
  
  
  Signature_Plot = rbind(SignatureScores,RandomScores)  
  
  levels(Signature_Plot$Sample) = c("Delete","Delete","0","5","10","15","20","25","30","45","60","90")
  Signature_Plot = Signature_Plot[Signature_Plot$Sample!="Delete",]
  Signature_Plot$Sample = as.numeric(as.character(Signature_Plot$Sample))
  Signature_Plot$ExpBase = factor(Signature_Plot$ExpBase,levels = c("Base","Sig"))
  
  pal = scales::hue_pal()(5)[c(1,3,5)]
  
  Plots = list()
  for(i in 1:length(names(SigList))){
    
    ## For LmFit % change we can try to inflate the ZScore such that there are no start points close to zero!
    ## LmFit = lm(ZScore + 3 ~ Sample*ExpBase,data = Signature_Plot[Signature_Plot$Cluster==names(SigList)[i],])
    
    LmFit = lm(ZScore ~ Sample*ExpBase,data = Signature_Plot[Signature_Plot$Cluster==names(SigList)[i],])
    
    em <- emtrends(LmFit,pairwise~ExpBase,var = 'Sample',adjust = "none")
    
    SlopeComparison = data.frame(em[[1]])
    
    # SlopeComparison[1,c(2,5,6)] = SlopeComparison[1,c(2,5,6)]*90/abs(Signature_Plot$ZScore[Signature_Plot$Cluster==names(SigList)[i]&Signature_Plot$Sample==0&Signature_Plot$ExpBase==SlopeComparison$ExpBase[1]])*100
    # SlopeComparison[2,c(2,5,6)] = SlopeComparison[2,c(2,5,6)]*90/abs(Signature_Plot$ZScore[Signature_Plot$Cluster==names(SigList)[i]&Signature_Plot$Sample==0&Signature_Plot$ExpBase==SlopeComparison$ExpBase[2]])*100
    
    SlopeComparison[1,c(2,5,6)] = SlopeComparison[1,c(2,5,6)]*90
    SlopeComparison[2,c(2,5,6)] = SlopeComparison[2,c(2,5,6)]*90
    
    Plots[[i]] = ggplot(data = SlopeComparison,aes(x = ExpBase,y = Sample.trend)) +
      scale_colour_manual(values = c("grey55",pal[i])) + 
      scale_fill_manual(values = c("grey55",pal[i])) +
      # geom_hline(yintercept = 0) +
      geom_errorbar(aes(ymin = lower.CL,ymax = upper.CL,color = ExpBase),width=0.4, size=0.5) + 
      geom_point(aes(fill = ExpBase),size = 4,shape = 22) +
      theme_cowplot() + 
      theme(legend.position = "none",axis.title = element_blank(),axis.text = element_blank()) +
      ylim(c(-2.5,2.5))
    print(names(SigList)[i])
    print(em[[2]])
  }
  
  Plots[[4]] = frame()
  for(i in 1:2) Plots[[i]] = Plots[[i]] + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank())
  
  CombinedPlot = plot_grid(plotlist = Plots,nrow = 4)
  CombinedPlot
}

#########################################################################################################
## Wherry Armstrong

Wherry = getGEO("GSE41867")
Wherry_Exprs = data.frame(exprs(Wherry[[1]]))

Annotated = data.frame(ACCNUM = sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse = ","),SYMBOL = sapply(contents(mogene10sttranscriptclusterSYMBOL),paste,collapse = ","))
Wherry_Exprs <- merge(Annotated, Wherry_Exprs,by=0,all=T)

Wherry_Exprs_Matrix = data.matrix(Wherry_Exprs[,4:ncol(Wherry_Exprs)])
rownames(Wherry_Exprs_Matrix) = toupper(Wherry_Exprs$SYMBOL)
colnames(Wherry_Exprs_Matrix) = c("N_1","N_2","N_3","N_4","Arm_D6_1","Arm_D6_2","Arm_D6_3","Arm_D6_4","Arm_D8_1","Arm_D8_2","Arm_D8_3","Arm_D8_4",
                                  "Arm_D15_1","Arm_D15_2","Arm_D15_3","Arm_D15_4","Arm_D30_1","Arm_D30_2","Arm_D30_3","Arm_D30_4","13_D6_1","13_D6_2",
                                  "13_D6_3","13_D6_4","13_D8_1","13_D8_2","13_D8_3","13_D8_4","13_D15_1","13_D15_2","13_D15_3","13_D15_4","13_D30_1",
                                  "13_D30_2","13_D30_3","13_D30_4")

Arm = Wherry_Exprs_Matrix[,1:20]
cl13 = Wherry_Exprs_Matrix[,c(1:4,21:36)]

boxplot(Arm)

rm(Annotated,Wherry,Wherry_Exprs,Wherry_Exprs_Matrix)

Arm = SignatureScorePlot(Arm)

#########################################################################################################
## Wherry cl. 13

cl13 = SignatureScorePlot(cl13)

#########################################################################################################
## Sarkar

Sarkar = getGEO("GSE10239")
Sarkar_Exprs = data.frame(exprs(Sarkar[[1]]))

## Annotation
x <- mouse4302SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
Symbols = data.frame(Symbol = unlist(xx))
Sarkar_Exprs = merge(Symbols,Sarkar_Exprs,by = 0)

Sarkar_Exprs_Matrix = data.matrix(Sarkar_Exprs[,3:ncol(Sarkar_Exprs)])
rownames(Sarkar_Exprs_Matrix) = toupper(Sarkar_Exprs$Symbol)
colnames(Sarkar_Exprs_Matrix) = c("N_1","N_2","N_3","Memory_1","Memory_2",
                                  "Memory_3","KLRG1int_1","KLRG1int_2","KLRG1int_3",
                                  "KLRG1hi_1","KLRG1hi_2","KLRG1hi_3")

Sarkar_Exprs_Matrix = Sarkar_Exprs_Matrix[,c(1:3,10:12,4:6)]
Sarkar_Exprs_Matrix = log2(Sarkar_Exprs_Matrix+1)
boxplot(Sarkar_Exprs_Matrix)

Sarkar_Arm = SignatureScorePlot(Sarkar_Exprs_Matrix)

rm(Sarkar_Exprs_Matrix,Sarkar_Exprs,Sarkar,x,xx,Symbols)

#########################################################################################################
## Scheitinger

GSE89307 <- read_delim("GSE89307_FPKM.txt","\t", escape_double = FALSE, trim_ws = TRUE)

Sheitinger_Matrix = data.matrix(GSE89307[,3:ncol(GSE89307)])
Sheitinger_Matrix = Sheitinger_Matrix[,-c(1:6,27:35)]
rownames(Sheitinger_Matrix) = GSE89307$Gene_Symbol
# Sheitinger_Matrix = Sheitinger_Matrix[rowSums(Sheitinger_Matrix)>0,]
Sheitinger_Matrix = log2((Sheitinger_Matrix+1))
Sheitinger_Matrix = Sheitinger_Matrix[rowSums(Sheitinger_Matrix)>0,]
rownames(Sheitinger_Matrix) = toupper(rownames(Sheitinger_Matrix))
Sheitinger_Matrix = Sheitinger_Matrix[,rev(1:ncol(Sheitinger_Matrix))]

Sheitinger_SV40TAG = SignatureScorePlot(Sheitinger_Matrix)

rm(GSE89307,Sheitinger_Matrix,mapped_probes)

#########################################################################################################
## Kupper

Kupper = getGEO("GSE79805")
Kupper_Exprs = data.frame(exprs(Kupper[[1]]))

Annotated = data.frame(ACCNUM = sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse = ","),SYMBOL = sapply(contents(mogene20sttranscriptclusterSYMBOL),paste,collapse = ","))
Kupper_Exprs <- merge(Annotated, Kupper_Exprs,by=0,all=T)

Kupper_Exprs_Matrix = data.matrix(Kupper_Exprs[,4:ncol(Kupper_Exprs)])
rownames(Kupper_Exprs_Matrix) = toupper(Kupper_Exprs$SYMBOL)
colnames(Kupper_Exprs_Matrix) = c("Naive", "CM", "EM", "Day5", "Day10", "Day15", "Day20", "Day25", "Day30", "Day45", "Day60", "Day90", "OT-I Fasbp4/5 dKO skin infiltrating cells 10 days post infection", "OT-I Fasbp4/5 dKO skin infiltrating cells 30 days post infection")

## Removes irrelevant timepoints
Kupper_Exprs_Matrix = Kupper_Exprs_Matrix[,-c(ncol(Kupper_Exprs_Matrix),ncol(Kupper_Exprs_Matrix)-1)]
Kupper_Exprs_Matrix = Kupper_Exprs_Matrix[,c(2,3,1,4:12)]

boxplot(Kupper_Exprs_Matrix)

set.seed(31415)
Kupper_VacV = SignatureScorePlot_OneSampleRep(Kupper_Exprs_Matrix)
set.seed(31415)
Kupper_VacV_Slopes = SlopePlot_OneSampleRep(Kupper_Exprs_Matrix)
rm(Annotated,Kupper_Exprs,Kupper_Exprs_Matrix,Kupper)

#########################################################################################################
## Combined Plot

CombinedPlot = plot_grid(Sarkar_Arm,Arm,cl13,Sheitinger_SV40TAG,Kupper_VacV,frame(),Kupper_VacV_Slopes,ncol = 7,
                         rel_widths = c(3,3,3,3,3,0.5,1))

print(CombinedPlot)

Scale = 2
pdf("BackScoring_RandomCondensed_v6.pdf",width = 6.55*Scale,height = 4.25*Scale)
print(CombinedPlot)
dev.off()
