library(GSVA)
library(GEOquery)
library(mogene10sttranscriptcluster.db)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(scales)

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
Wherry = getGEO("GSE41867")
Wherry_Exprs = data.frame(exprs(Wherry[[1]]))

Annotated = data.frame(ACCNUM = sapply(contents(mogene10sttranscriptclusterACCNUM), paste, collapse = ","),SYMBOL = sapply(contents(mogene10sttranscriptclusterSYMBOL),paste,collapse = ","))
Wherry_Exprs <- merge(Annotated, Wherry_Exprs,by=0,all=T)

Wherry_Exprs_Matrix = data.matrix(Wherry_Exprs[,4:ncol(Wherry_Exprs)])
rownames(Wherry_Exprs_Matrix) = Wherry_Exprs$SYMBOL
colnames(Wherry_Exprs_Matrix) = c("N_1","N_2","N_3","N_4","Arm_D6_1","Arm_D6_2","Arm_D6_3","Arm_D6_4","Arm_D8_1","Arm_D8_2","Arm_D8_3","Arm_D8_4",
                                  "Arm_D15_1","Arm_D15_2","Arm_D15_3","Arm_D15_4","Arm_D30_1","Arm_D30_2","Arm_D30_3","Arm_D30_4","13_D6_1","13_D6_2",
                                  "13_D6_3","13_D6_4","13_D8_1","13_D8_2","13_D8_3","13_D8_4","13_D15_1","13_D15_2","13_D15_3","13_D15_4","13_D30_1",
                                  "13_D30_2","13_D30_3","13_D30_4")

Arm = Wherry_Exprs_Matrix[,1:20]
cl13 = Wherry_Exprs_Matrix[,c(1:4,21:36)]

boxplot(Arm)
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

## Does Random Scoring in Wherry Armstrong
nperm = 100
Arm_RandomScores = NULL
cl13_RandomScores = Arm_RandomScores
for(i in 1:length(DEG_SignaturesConverted)) {
  Arm_RandomScores = rbind(Arm_RandomScores,SignatureScore_Random(Arm,DEG_SignaturesConverted[[i]],nperms = nperm))
  cl13_RandomScores = rbind(cl13_RandomScores,SignatureScore_Random(cl13,DEG_SignaturesConverted[[i]],nperms = nperm))
}
Arm_RandomScores = data.frame(Arm_RandomScores)
Arm_RandomScores$State = rep(names(DEG_SignaturesConverted),each = nperm)
Arm_RandomScores = reshape2::melt(Arm_RandomScores)
Arm_RandomScores$SigType = "Baseline"

ggdensity(Arm_RandomScores,"value",color ="variable",facet.by = "State")
ggline(Arm_RandomScores,"variable","value",facet.by = "State",add = "mean_se")

## Does Signature Scoring in Wherry Armstrong
ArmScores = matrix(nrow = length(DEG_SignaturesConverted),ncol = ncol(Arm))

for(i in 1:length(DEG_SignaturesConverted)) ArmScores[i,] = SignatureScore_UpMinusDown(Arm,DEG_SignaturesConverted[[i]])
ArmScores = data.frame(ArmScores)
colnames(ArmScores) = colnames(Arm)
ArmScores$State = names(DEG_SignaturesConverted)
ArmScores = reshape2::melt(ArmScores)
ArmScores$SigType = "Exp"

# Combines Scores, plots, does stats
CombinedScores = rbind(Arm_RandomScores,ArmScores)
CombinedScores$variable = factor(sapply(strsplit(as.character(CombinedScores$variable),fixed = TRUE,split = "_"), function(x) x[length(x)-1]),
                                 levels = c("N","D6","D8","D15","D30"))
CombinedScores$SigType = factor(CombinedScores$SigType)
CombinedScores$State = factor(CombinedScores$State)

pal = scales::hue_pal()(5)[c(1,3,5)]

# CombinedScoresPlot = list()
# for(i in 1:3){
#   CombinedScoresPlot[[i]] = ggline(CombinedScores[CombinedScores$State==levels(CombinedScores$State)[i],],"variable","value",color = "SigType",add = "mean_se",palette = c("grey",pal[i])) +
#     stat_compare_means(aes(group = SigType),method = "t.test",label = "p.signif",label.y = 1.8)
# 
#   CombinedScoresPlot[[i]] = ggpar(CombinedScoresPlot[[i]] + theme(strip.background = element_blank(),strip.text = element_blank(),
#                                                                   axis.title = element_blank()),ylim = c(-1.9,1.9),legend = "none")
# }
# 
# ggarrange(CombinedScoresPlot[[1]],CombinedScoresPlot[[2]],CombinedScoresPlot[[3]],ncol = 3)

compare_means(value ~ SigType,data = CombinedScores,ref.group = "Baseline",group.by = c("State","variable"),method = "t.test")

CombineScores = CombinedScores
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
    geom_line(data=StatePlot, mapping=aes(variable,value),colour = Color)+
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
  geom_line(data=StatePlots[[1]], mapping=aes(variable,value),colour = pal[1])+
  geom_errorbar(data=StatePlots[[1]], mapping=aes(variable,ymin = Lower,ymax = Upper),width=0.4, size=0.5, color=pal[1]) + 
  geom_point(data = StatePlots[[2]],aes(variable,value,group = 1),size = 4, shape = 21, fill = pal[2]) + 
  geom_line(data=StatePlots[[2]], mapping=aes(variable,value),colour = pal[2])+
  geom_errorbar(data=StatePlots[[2]], mapping=aes(variable,ymin = Lower,ymax = Upper),width=0.4, size=0.5, color=pal[2]) + 
  geom_point(data = StatePlots[[3]],aes(variable,value,group = 1),size = 4, shape = 21, fill = pal[3]) + 
  geom_line(data=StatePlots[[3]], mapping=aes(variable,value),colour = pal[3])+
  geom_errorbar(data=StatePlots[[3]], mapping=aes(variable,ymin = Lower,ymax = Upper),width=0.4, size=0.5, color=pal[3]) + 
  theme_classic()+
  ylim(-1.9,1.9)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())

saveRDS(ggarrange(Plots[[1]],Plots[[2]],Plots[[3]],CombinedPlot,nrow = 4),"Wherry_Arm.rds")

