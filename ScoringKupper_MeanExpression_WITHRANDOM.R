library(GSVA)
library(GEOquery)
library(mogene20sttranscriptcluster.db)
library(readr)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(rstatix)
library(emmeans)

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
Kupper = getGEO("GSE79805")
Kupper_Exprs = data.frame(exprs(Kupper[[1]]))

Annotated = data.frame(ACCNUM = sapply(contents(mogene20sttranscriptclusterACCNUM), paste, collapse = ","),SYMBOL = sapply(contents(mogene20sttranscriptclusterSYMBOL),paste,collapse = ","))
Kupper_Exprs <- merge(Annotated, Kupper_Exprs,by=0,all=T)

Kupper_Exprs_Matrix = data.matrix(Kupper_Exprs[,4:ncol(Kupper_Exprs)])
rownames(Kupper_Exprs_Matrix) = Kupper_Exprs$SYMBOL
colnames(Kupper_Exprs_Matrix) = c("Naive", "CM", "EM", "Day 5", "Day 10", "Day 15", "Day 20", "Day 25", "Day 30", "Day 45", "Day 60", "Day 90", "OT-I Fasbp4/5 dKO skin infiltrating cells 10 days post infection", "OT-I Fasbp4/5 dKO skin infiltrating cells 30 days post infection")

## Removes irrelevant timepoints
Kupper_Exprs_Matrix = Kupper_Exprs_Matrix[,-c(ncol(Kupper_Exprs_Matrix),ncol(Kupper_Exprs_Matrix)-1)]


boxplot(Kupper_Exprs_Matrix)

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

## Does Random Scoring in Kupper
nperm = 100
RandomScores = NULL
for(i in 1:length(DEG_SignaturesConverted)) RandomScores = rbind(RandomScores,SignatureScore_Random(Kupper_Exprs_Matrix,DEG_SignaturesConverted[[i]],nperms = nperm))
RandomScores = data.frame(RandomScores)
RandomScores$State = rep(names(DEG_SignaturesConverted),each = nperm)
RandomScores = reshape2::melt(RandomScores)
RandomScores$SigType = "Baseline"
ggdensity(RandomScores,"value",color ="variable",facet.by = "State")

RandomSigs = ggline(RandomScores,"variable","value",facet.by = "State",add = "mean_se")
ggpar(RandomSigs + theme(strip.background = element_blank(),strip.text = element_blank(),
                                          axis.title = element_blank()),ylim = c(-1.9,1.9),legend = "none")

## Does Kupper Scoring
GeneSetScores = matrix(nrow = length(DEG_SignaturesConverted),ncol = ncol(Kupper_Exprs_Matrix))
for(i in 1:length(DEG_SignaturesConverted)) GeneSetScores[i,] = SignatureScore_UpMinusDown(Kupper_Exprs_Matrix,DEG_SignaturesConverted[[i]])
GeneSetScores = data.frame(GeneSetScores)
colnames(GeneSetScores) = colnames(Kupper_Exprs_Matrix)
GeneSetScores$State = names(DEG_SignaturesConverted)
GeneSetScores = reshape2::melt(GeneSetScores)
GeneSetScores$SigType = "Exp"

# Combines Scores, plots, does stats
levels(RandomScores$variable) = levels(GeneSetScores$variable)
CombinedScores = rbind(RandomScores,GeneSetScores)
CombinedScores$SigType = factor(CombinedScores$SigType)
CombinedScores$State = factor(CombinedScores$State)
CombinedScores$variable = factor(CombinedScores$variable,levels(CombinedScores$variable)[c(2,3,1,4:12)])

pal = scales::hue_pal()(5)[c(1,3,5)]

# CombinedScoresPlot = list()
# for(i in 1:3){
#   stat.test = NULL
#   for(j in 1:length(levels(CombinedScores$variable))){
#     value = CombinedScores$value[CombinedScores$State==levels(CombinedScores$State)[i]&CombinedScores$SigType=="Exp"&CombinedScores$variable==levels(CombinedScores$variable)[j]]
#     stat.test <- rbind(stat.test,CombinedScores[CombinedScores$State==levels(CombinedScores$State)[i]&CombinedScores$SigType=="Baseline"&CombinedScores$variable==levels(CombinedScores$variable)[j],] %>%
#                          group_by(variable) %>%
#                          t_test(value ~ 1, mu = value) %>%
#                          add_significance() %>%
#                          adjust_pvalue() %>%
#                          mutate(y.position = 1.8))
#   }
# 
#   CombinedScoresPlot[[i]] = ggline(CombinedScores[CombinedScores$State==levels(CombinedScores$State)[i],],"variable","value",color = "SigType",add = "mean_se",palette = c("grey",pal[i])) +
#     stat_pvalue_manual(stat.test, label = "p.signif", xmin = "variable", xmax = NULL)
#   CombinedScoresPlot[[i]] = ggpar(CombinedScoresPlot[[i]] + theme(strip.background = element_blank(),strip.text = element_blank(),
#                                                                   axis.title = element_blank()),ylim = c(-1.9,1.9),legend = "none")
# }
# 
# ggarrange(CombinedScoresPlot[[1]],CombinedScoresPlot[[2]],CombinedScoresPlot[[3]],ncol = 3)
# 
# # # Code for rstatix: https://github.com/kassambara/ggpubr/issues/79 for explanation on one-sample t test
# # stat.test = NULL
# # for(j in 1:length(levels(CombinedScores$variable))){
# #   value = CombinedScores$value[CombinedScores$State==levels(CombinedScores$State)[1]&CombinedScores$SigType=="Exp"&CombinedScores$variable==levels(CombinedScores$variable)[j]]
# #   stat.test <- rbind(stat.test,CombinedScores[CombinedScores$State==levels(CombinedScores$State)[1]&CombinedScores$SigType=="Baseline"&CombinedScores$variable==levels(CombinedScores$variable)[j],] %>%
# #     group_by(variable) %>%
# #     t_test(value ~ 1, mu = value) %>%
# #     add_significance() %>%
# #     adjust_pvalue() %>%
# #     mutate(y.position = 35))
# # }
# 
# ## Tests Stats Package with t.test function
# 
# Cluster_Searched = "Cluster1"
# Time_Searched = "Naive"
# 
# t.test(CombinedScores$value[CombinedScores$State==Cluster_Searched&CombinedScores$variable==Time_Searched&CombinedScores$SigType=="Baseline"],
#        mu = CombinedScores$value[CombinedScores$State==Cluster_Searched&CombinedScores$variable==Time_Searched&CombinedScores$SigType=="Exp"])

#############################################################################################################################################
## Plotting Function

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


#############################################################################################################################################
## Slope Analysis
# 
# ## Analyzes Slopes of Kupper Timecourse Enrichment
# 
# SlopeAnalysis = GeneSetScores
# levels(SlopeAnalysis$variable) = c("0","nope","nope","5","10","15","20","25","30","45","60","90")
# SlopeAnalysis = SlopeAnalysis[SlopeAnalysis$variable!="nope",]
# SlopeAnalysis$Time = as.numeric(as.character(SlopeAnalysis$variable))
# 
# LmFit = lm(value~Time*State,data = SlopeAnalysis)
# 
# em <- emtrends(LmFit,pairwise~State,var = 'Time',adhust = 'none')
# data.frame(em[[1]])
# ## Plots Points + Linear Regression
# 
# ggplot(data = SlopeAnalysis,aes(Time,value,color = State)) +
#   scale_color_manual(values = pal)+
#   geom_smooth(method = "lm")+
#   geom_point()+
#   theme_classic()

#############################################################################################################################################
## Slope Analysis v2
SlopeAnalysis_v2 = CombinedScores
levels(SlopeAnalysis_v2$variable) = c("0","nope","nope","5","10","15","20","25","30","45","60","90")
SlopeAnalysis_v2 = SlopeAnalysis_v2[SlopeAnalysis_v2$variable!="nope",]
SlopeAnalysis_v2$Time = as.numeric(as.character(SlopeAnalysis_v2$variable))

Clust = 3

ggplot(data = SlopeAnalysis_v2[as.numeric(SlopeAnalysis_v2$State)==Clust,],aes(Time,value,color = SigType)) +
  scale_color_manual(values = c("black",pal[Clust]))+
  geom_smooth(method = "lm") + 
  geom_point(data = SlopeAnalysis_v2[as.numeric(SlopeAnalysis_v2$State)==Clust&SlopeAnalysis_v2$SigType=="Exp",],fill = pal[Clust],size = 4,shape = 21,color = "black")+
  theme_classic()


LmFit = lm(value~Time*SigType,data = SlopeAnalysis_v2[as.numeric(SlopeAnalysis_v2$State)==Clust,])
em <- emtrends(LmFit,pairwise~SigType,var = 'Time',adjust = 'none')
em

ggplot(data.frame(em[[1]]),aes(SigType,Time.trend,group = 1,fill = SigType))+
  scale_fill_manual(values = c("black",pal[Clust]))+
  geom_errorbar(data = data.frame(em[[1]]),mapping = aes(SigType,ymin = lower.CL,ymax = upper.CL),width=0.4, size=0.5, color=c("black",pal[Clust])) +
  geom_point(shape = 21,size = 4) + 
  theme_classic()+
  theme(legend.position = "none")

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

saveRDS(ggarrange(Plots[[1]],Plots[[2]],Plots[[3]],CombinedPlot,nrow = 4),"Kupper.rds")
  