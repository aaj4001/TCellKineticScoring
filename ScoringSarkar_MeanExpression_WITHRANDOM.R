library(GSVA)
library(GEOquery)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(mouse4302.db)
library(ggpubr)

DEGList = readRDS("2019_07_06 T Cell State DEG List Unique States.rds")
DEG_Signatures = list()

for(i in 1:3){
  DEG_Signatures = c(DEG_Signatures,  Up = list(rownames(DEGList[[i]])[DEGList[[i]]$avg_logFC>0]))
  DEG_Signatures = c(DEG_Signatures,  Down = list(rownames(DEGList[[i]])[DEGList[[i]]$avg_logFC<0]))
}

names(DEG_Signatures) = paste(rep(c("Cluster1","Cluster3","Cluster5"),each = 2),names(DEG_Signatures),sep = "_")

Human_MouseConversion <- read_csv("Human_MouseConversion.xlsx.csv")
DEG_SignaturesConverted = list()
for(i in 1:6) DEG_SignaturesConverted = c(DEG_SignaturesConverted,list(merge(data.frame(DEG_Signatures[[i]]),Human_MouseConversion,by.x = 1, by.y = 2)[,2]))
names(DEG_SignaturesConverted) = names(DEG_Signatures)

## Reformats Signature List to proper architecture

DEG_SignaturesConverted = list(Cluster1 = list(Up = DEG_SignaturesConverted[["Cluster1_Up"]], Down = DEG_SignaturesConverted[["Cluster1_Down"]]),
                               Cluster3 = list(Up = DEG_SignaturesConverted[["Cluster3_Up"]], Down = DEG_SignaturesConverted[["Cluster3_Down"]]),
                               Cluster5 = list(Up = DEG_SignaturesConverted[["Cluster5_Up"]], Down = DEG_SignaturesConverted[["Cluster5_Down"]]))

## Downloads Kupper Dataset and Annotates Genes
Sarkar = getGEO("GSE10239")
Sarkar_Exprs = data.frame(exprs(Sarkar[[1]]))

## Annotation
x <- mouse4302SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
Symbols = data.frame(Symbol = unlist(xx))
Sarkar_Exprs = merge(Symbols,Sarkar_Exprs,by = 0)

Sarkar_Exprs_Matrix = data.matrix(Sarkar_Exprs[,3:ncol(Sarkar_Exprs)])
rownames(Sarkar_Exprs_Matrix) = Sarkar_Exprs$Symbol
colnames(Sarkar_Exprs_Matrix) = c("N_1","N_2","N_3","Memory_1","Memory_2",
                                  "Memory_3","KLRG1int_1","KLRG1int_2","KLRG1int_3",
                                  "KLRG1hi_1","KLRG1hi_2","KLRG1hi_3")

Sarkar_Exprs_Matrix = Sarkar_Exprs_Matrix[,c(1:3,10:12,4:6)]
Sarkar_Exprs_Matrix = log2(Sarkar_Exprs_Matrix+1)
boxplot(Sarkar_Exprs_Matrix)

## Method for Scoring Signatures
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
  NetSignatures
}

## Does Random Scoring in Sarker
nperm = 100
RandomScores = NULL
for(i in 1:length(DEG_SignaturesConverted)) RandomScores = rbind(RandomScores,SignatureScore_Random(Sarkar_Exprs_Matrix,DEG_SignaturesConverted[[i]],nperms = nperm))
RandomScores = data.frame(RandomScores)
RandomScores$State = rep(names(DEG_SignaturesConverted),each = nperm)
RandomScores = reshape2::melt(RandomScores)
RandomScores$SigType = "Baseline"
ggdensity(RandomScores,"value",color ="variable",facet.by = "State")

RandomSigs = ggline(RandomScores,"variable","value",facet.by = "State",add = "mean_se")
ggpar(RandomSigs + theme(strip.background = element_blank(),strip.text = element_blank(),
                         axis.title = element_blank()),ylim = c(-1.9,1.9),legend = "none")

## Does Signature Scoring
GeneSetScores = matrix(nrow = length(DEG_SignaturesConverted),ncol = ncol(Sarkar_Exprs_Matrix))
for(i in 1:length(DEG_SignaturesConverted)) GeneSetScores[i,] = SignatureScore_UpMinusDown(Sarkar_Exprs_Matrix,DEG_SignaturesConverted[[i]])
GeneSetScores = data.frame(GeneSetScores)
colnames(GeneSetScores) = colnames(Sarkar_Exprs_Matrix)
GeneSetScores$State = names(DEG_SignaturesConverted)

GeneSetScores_Melt = reshape2::melt(GeneSetScores)
GeneSetScores_Melt$SigType = "Exp"

## Combines Scores and plots with t tests
CombineScores = rbind(RandomScores,GeneSetScores_Melt)
CombineScores$variable = factor(sapply(strsplit(as.character(CombineScores$variable),fixed = TRUE,split = "_"),function(x) x[1]),levels = c("N","KLRG1hi","Memory"))
CombineScores$SigType = factor(CombineScores$SigType)
CombineScores$State = factor(CombineScores$State)

pal = scales::hue_pal()(5)[c(1,3,5)]

# CombineScoresPlot = list()
# for(i in 1:3){
#   CombineScoresPlot[[i]] = ggline(CombineScores[CombineScores$State==levels(CombineScores$State)[i],],"variable","value",color = "SigType",add = "mean_se",palette = c("grey",pal[i])) +
#     stat_compare_means(aes(group = SigType),method = "t.test",label = "p.signif",label.y = 1.8)
# 
#   CombineScoresPlot[[i]] = ggpar(CombineScoresPlot[[i]] + theme(strip.background = element_blank(),strip.text = element_blank(),
#                                                               axis.title = element_blank()),ylim = c(-1.9,1.9),legend = "none")
# }
# 
# ggarrange(CombineScoresPlot[[1]],CombineScoresPlot[[2]],CombineScoresPlot[[3]],ncol = 3)

compare_means(value ~ SigType,data = CombineScores,ref.group = "Baseline",group.by = c("State","variable"),method = "t.test") -> A

Plots = list()
for(j in 1:3){
  StateValues = CombineScores[CombineScores$State==levels(CombineScores$State)[j]&CombineScores$SigType=="Exp",]
  Baselines = CombineScores[CombineScores$State==levels(CombineScores$State)[j]&CombineScores$SigType=="Baseline",]
  Color = pal[j]
  
  value = numeric(length(levels(StateValues$variable)))
  SD = value
  Upper = value
  Lower = value
  B_Value = value
  B_UpperCI = value
  B_LowerCI = value
  for(i in 1:length(levels(StateValues$variable))) {
    value[i] = mean(StateValues$value[as.numeric(StateValues$variable)==i])
    SD[i] = sd(StateValues$value[as.numeric(StateValues$variable)==i])
    B_Value[i] = mean(Baselines$value[as.numeric(Baselines$variable)==i])
    #B_UpperCI[i] = t.test(Baselines$value[as.numeric(Baselines$variable)==i])$conf.int[2]
    #B_LowerCI[i] = t.test(Baselines$value[as.numeric(Baselines$variable)==i])$conf.int[1]
    
    B_UpperCI[i] =  B_Value[i] + sd(Baselines$value[as.numeric(Baselines$variable)==i])
    B_LowerCI[i] = B_Value[i] - sd(Baselines$value[as.numeric(Baselines$variable)==i])
  }
  Upper = value + SD
  Lower = value - SD
  
  StatePlot = data.frame(variable = factor(levels(StateValues$variable),levels = levels(StateValues$variable)),
                         value,
                         Lower,
                         Upper)
  
  BaselinePlot = data.frame(variable = factor(levels(StateValues$variable),levels = levels(StateValues$variable)),
                            value = B_Value,
                            Upper = B_UpperCI,
                            Lower = B_LowerCI)
  
Plots[[j]] =  ggplot(BaselinePlot,aes(variable,value,group = 1))+
    geom_ribbon(aes(ymin = Lower, ymax = Upper),fill = "grey")+
    geom_line() + 
    geom_errorbar(data=StatePlot, mapping=aes(variable,ymin = Lower,ymax = Upper),width=0.4, size=0.5, color=Color) + 
    geom_point(data=StatePlot, mapping=aes(variable,value), size=4, shape=21,fill = Color)+
    theme_classic()+
    ylim(-1.9,1.9)+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

ggarrange(Plots[[1]],Plots[[2]],Plots[[3]],ncol = 3)

################################################################################################################################
## CombinedPlot

StatePlots = list()
for(j in 1:3){
  StateValues = CombineScores[CombineScores$State==levels(CombineScores$State)[j]&CombineScores$SigType=="Exp",]
  Baselines = CombineScores[CombineScores$State==levels(CombineScores$State)[j]&CombineScores$SigType=="Baseline",]
  Color = pal[j]
  
  value = numeric(length(levels(StateValues$variable)))
  SD = value
  Upper = value
  Lower = value
  B_Value = value
  B_UpperCI = value
  B_LowerCI = value
  for(i in 1:length(levels(StateValues$variable))) {
    value[i] = mean(StateValues$value[as.numeric(StateValues$variable)==i])
    SD[i] = sd(StateValues$value[as.numeric(StateValues$variable)==i])
    B_Value[i] = mean(Baselines$value[as.numeric(Baselines$variable)==i])
    B_UpperCI[i] = t.test(Baselines$value[as.numeric(Baselines$variable)==i])$conf.int[2]
    B_LowerCI[i] = t.test(Baselines$value[as.numeric(Baselines$variable)==i])$conf.int[1]
  }
  Upper = value + SD
  Lower = value - SD
  
  StatePlots[[levels(CombineScores$State)[j]]] = data.frame(variable = factor(levels(StateValues$variable),levels = levels(StateValues$variable)),
                         value,
                         Lower,
                         Upper)
}

CombinedPlot = ggplot(StatePlots[[1]],aes(variable,value,group = 1))+
  geom_point(size = 4, shape = 21, fill = pal[1]) + 
  geom_errorbar(data=StatePlots[[1]], mapping=aes(variable,ymin = Lower,ymax = Upper),width=0.4, size=0.5, color=pal[1]) + 
  geom_point(data = StatePlots[[2]],aes(variable,value,group = 1),size = 4, shape = 21, fill = pal[2]) + 
  geom_errorbar(data=StatePlots[[2]], mapping=aes(variable,ymin = Lower,ymax = Upper),width=0.4, size=0.5, color=pal[2]) + 
  geom_point(data = StatePlots[[3]],aes(variable,value,group = 1),size = 4, shape = 21, fill = pal[3]) + 
  geom_errorbar(data=StatePlots[[3]], mapping=aes(variable,ymin = Lower,ymax = Upper),width=0.4, size=0.5, color=pal[3]) + 
  theme_classic()+
  ylim(-1.9,1.9)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

saveRDS(ggarrange(Plots[[1]],Plots[[2]],Plots[[3]],CombinedPlot,nrow = 4),"Sarkar.rds")
