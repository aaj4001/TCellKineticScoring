library(RColorBrewer)
library(ggpubr)
library(scales)
library(reshape2)
library(cowplot)
library(emmeans)
library(nlme)
library(msigdbi)
library(dplyr)
library(rstatix)
library(magrittr)

## Imports Query Sets and Kinetic Sets  #########################################################################################################

DEGList = read.gmt("Hacohen_StateDEGs_Human.gmt")$genesets

DEGList_UpDown = list()
for(i in seq(1,6,by = 2)) DEGList_UpDown[[names(DEGList)[i]]] = list(Up = DEGList[[i]],Down = DEGList[[i+1]])
names(DEGList_UpDown) = sapply(strsplit(names(DEGList_UpDown),split = "_",fixed = TRUE),function(x) x[1])

palette = scales::hue_pal()(5)[c(1,3,5)]

load("KineticSets.rda")


## Backscoring Up Minus Down #########################################################################################################

SignatureScore <- function(InputGeneMatrix,TestSignatures, RowInds = NULL){
  ## If RowInds is passed, will ignore TestSignatures and select genes based on row index provided (use for random scoring)
  if(!is.null(RowInds)) {GeneMatrixFiltered = InputGeneMatrix[RowInds,]
  } else GeneMatrixFiltered = InputGeneMatrix[rownames(InputGeneMatrix) %in% TestSignatures,]
  
  GeneMatrixFiltered_Zscore = t(apply(GeneMatrixFiltered,1,scale)); colnames(GeneMatrixFiltered_Zscore) = colnames(GeneMatrixFiltered)

  Mean_ZScore = apply(GeneMatrixFiltered_Zscore,2,mean); names(Mean_ZScore) = colnames(GeneMatrixFiltered)
  Mean_ZScore
}

SignatureScore_UpMinusDown <- function(InputGeneMatrix,UpDownList){
  UpSignature = SignatureScore(InputGeneMatrix,UpDownList[["Up"]])
  if(!is.null(UpDownList[["Down"]])){
    DnSignature = SignatureScore(InputGeneMatrix,UpDownList[["Down"]])
    NetSignature = UpSignature - DnSignature
  } else NetSignature = UpSignature

  NetSignature
}

SignatureScore_Random <- function(InputGeneMatrix,UpDownList){
  UpGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Up"]]))
  UpSignature = SignatureScore(InputGeneMatrix,RowInds = UpGenesRandom)
  
  if(!is.null(UpDownList[["Down"]])){
    DnGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Down"]]))
    DnSignature = SignatureScore(InputGeneMatrix,RowInds = DnGenesRandom)
    NetRandom = UpSignature - DnSignature
  } else NetRandom = UpSignature

  NetRandom
}

## Function for analyzing expression data of signature scores, and conducting statistics if needed
## Will use this for bootstrapping of random scores as well

ExpressionAnalysisByTimepoint <- function(ExprsData, PlotType, OnlyMeans = F){
  if(PlotType == "Replicate_Line_LME"){
    Timepoints = unique(sapply(strsplit(names(ExprsData),split = "_"),function(y) paste(y[-length(y)],collapse = "_")))
    
    ## LME Model Calculation
    LME_Model = data.frame(Sample = names(ExprsData), Exp = ExprsData) %>%
      mutate(Rep = sapply(strsplit(Sample,split = "_"),function(y) y[length(y)]),
             Sample = sapply(strsplit(Sample,split = "_"),function(y) paste(y[-length(y)],collapse = "_")),
             Sample = factor(Sample, levels = Timepoints))
    
    LME_Sig_Results = lme(Exp~Sample,random=~1|Rep,data = LME_Model)
    LME_Sig_Results = summary(emmeans(LME_Sig_Results,"pairwise"~Sample,adjust = "none"))
    
    LME_Sig_Means = LME_Sig_Results[["emmeans"]] %>%
      mutate(lower.CL = emmean - (SE*1.96),upper.CL = emmean + (SE*1.96)) %>%
      select(c("Sample","emmean","SE","lower.CL","upper.CL")) %>%
      rename(Mean = emmean)
    
    LME_Sig_PVals = LME_Sig_Results$contrasts %>%
      mutate(group2 = factor(sapply(strsplit(contrast, split = " - "), function(x) x[2]), levels = Timepoints),
             group1 = factor(sapply(strsplit(contrast, split = " - "), function(x) x[1]), levels = Timepoints),
             PVal_Asterisk = as.character(symnum(p.value,cutpoints = c(-Inf,0.0001,0.001,0.01,0.05,Inf), symbols = c("****","***","**","*","")))) %>%
      subset(as.numeric(group2) - as.numeric(group1) == 1)
  } else if(PlotType=="Line_SlopeAnalysis"){
    Timepoints = factor(names(ExprsData),levels = names(ExprsData))
    LME_Sig_Means = data.frame(Sample = Timepoints,Mean = ExprsData,
                               Timepoints_v2 = c(-20,-10,0,5,10,15,20,25,30,45,60,90))
    
    LME_Sig_PVals = LME_Sig_Means %>% subset(Timepoints_v2 >=0 )
    LME_Sig_PVals = lm(Mean ~ Timepoints_v2,data = LME_Sig_PVals)
    
    LME_Sig_PVals = list(DeltaExpr = summary(LME_Sig_PVals)$coefficients[2,1]*90,
                         PVals = summary(LME_Sig_PVals)$coefficients[2,4])
    LME_Sig_PVals$PVals = as.character(symnum(LME_Sig_PVals$PVals,cutpoints = c(-Inf,0.0001,0.001,0.01,0.05,Inf), symbols = c("****","***","**","*","ns")))
    LME_Sig_PVals = paste0("DExprs_d0-d90: ",signif(LME_Sig_PVals$DeltaExpr,3),", Sig: ",LME_Sig_PVals$PVals)
  } else if(PlotType == "Replicate_Line_TTest"){
    Timepoints = unique(sapply(strsplit(names(ExprsData),split = "_"),function(y) paste(y[-length(y)],collapse = "_")))
    
    LME_Model = data.frame(Sample = names(ExprsData), Exp = ExprsData) %>%
      mutate(Rep = sapply(strsplit(Sample,split = "_"),function(y) y[length(y)]),
             Sample = sapply(strsplit(Sample,split = "_"),function(y) paste(y[-length(y)],collapse = "_")),
             Sample = factor(Sample, levels = Timepoints))
    
    LME_Sig_Means = LME_Model %>% dplyr::group_by(Sample) %>%
      summarize(Mean = mean(Exp), SE = plotrix::std.error(Exp),
                lower.CL = Mean - (SE*1.96),upper.CL = Mean + (SE*1.96))
    
    LME_Sig_PVals = LME_Model %>% pairwise_t_test(Exp ~ Sample, p.adjust.method = "none",pool.sd = T) %>%
      mutate(group1 = factor(group1, levels = Timepoints),
             group2 = factor(group2, levels = Timepoints),
             PVal_Asterisk = as.character(symnum(p,cutpoints = c(-Inf,0.0001,0.001,0.01,0.05,Inf), symbols = c("****","***","**","*","")))) %>%
      subset(as.numeric(group2) - as.numeric(group1) == 1)
  } else{LME_Sig_Means = NULL;LME_Sig_PVals = NULL}
  
  if(OnlyMeans){
    Output = LME_Sig_Means$Mean; names(Output) = LME_Sig_Means$Sample
  } 
  else Output = list(Means = LME_Sig_Means, PVals = LME_Sig_PVals)
  
  Output
}

# KineticSets = KineticLists;SigList = DEGList_UpDown;pal = rep("black",length(SigList));nperm = 10
# SameScale = T;PValsAnnot = T;PValsSize = 4.2

SignatureScorePlot <- function(KineticSets,SigList,pal = rep("black",length(SigList)),nperm = 100,
                               SameScale = T, PValsAnnot = T,PValsSize = 4.2){
  
    SignatureScoreResults = lapply(SigList, function(Sig) {
    KineticResult = lapply(KineticSets, function(Kinetic){
      RandScore = NULL
      for(i in 1:nperm) RandScore = rbind(RandScore, SignatureScore_Random(Kinetic$Exprs,Sig))
      
      return(list(SigScore = SignatureScore_UpMinusDown(Kinetic$Exprs,Sig),
                  RandScore = RandScore,
                  PlotType = Kinetic$PlotType))
    })
    return(KineticResult)})
    
    # Color Palette Assignment
    for(i in 1:length(SigList)) SignatureScoreResults[[i]] %<>% lapply(function(Kinetic){
      Kinetic$PlotColor = ifelse(!is.null(pal),pal[i],"black")
      return(Kinetic)
    })
    
    SignatureScoreResults %<>% unlist(recursive = F)
    
    ## Generates Signature Score Plots (plots signature score and random scores)
    ## Different plotting style for slope line and other plot types
    
    SignatureScorePlot = SignatureScoreResults %>% lapply(function(Result){
      SignatureScores = ExpressionAnalysisByTimepoint(Result$SigScore,Result$PlotType)
      
      RandomScores = apply(Result$RandScore,1,ExpressionAnalysisByTimepoint,PlotType = Result$PlotType,OnlyMeans = T)
      RandomScores = data.frame(Sample = factor(rownames(RandomScores), levels = rownames(RandomScores)), 
                                Mean = apply(RandomScores, 1, mean),
                                SE = apply(RandomScores,1,plotrix::std.error)) %>%
        mutate(lower.CL = Mean - SE, upper.CL = Mean + SE)
      
      SigScorePlot = ggplot(data = RandomScores,aes(x = Sample,y = Mean, group = 1)) + 
        geom_ribbon(data = RandomScores,aes(ymin = lower.CL,ymax = upper.CL,),fill = "gray") +
        geom_line(group = 0)
      
      if(Result$PlotType == "Line_SlopeAnalysis"){
        SigScorePlot = SigScorePlot + 
          geom_smooth(data = SignatureScores$Means %>% subset(!Sample %in% c("CM","EM")),method = "lm", formula = y~x, color = "transparent") +
          geom_line(data = SignatureScores$Means %>% subset(!Sample %in% c("CM","EM")),color = Result$PlotColor) + 
          geom_point(data = SignatureScores$Means,fill = Result$PlotColor,shape = 21,size = 4) +
          geom_smooth(data = SignatureScores$Means %>% subset(!Sample %in% c("CM","EM")),method = "lm",formula = y~x,color = colorspace::darken(Result$PlotColor,0.2),fill = "transparent") +
          theme_cowplot() +
          theme(axis.title = element_blank())
      } else{
        SigScorePlot = SigScorePlot + 
          geom_line(data = SignatureScores$Means,color = Result$PlotColor) + 
          geom_errorbar(data = SignatureScores$Means, aes(ymin = lower.CL,ymax = upper.CL),color = Result$PlotColor,width=0.4, size=0.5) + 
          geom_point(data = SignatureScores$Means,fill = Result$PlotColor,shape = 21,size = 4) +
          theme_cowplot() +
          theme(axis.title = element_blank())
      }

      list(Plot = SigScorePlot, AxisLimits = layer_scales(SigScorePlot)$y$get_limits(),
           PVals = SignatureScores$PVals, PlotType = Result$PlotType, PlotColor = Result$PlotColor)
    })
  
    ## Adjusts/Syncs Scales if User Specified (SameScale Variable)
    if(SameScale){
      AxisMin = min(sapply(SignatureScorePlot,function(x) x$AxisLimits[1]))
      AxisMax = max(sapply(SignatureScorePlot,function(x) x$AxisLimits[2]))
      
      SignatureScorePlot %<>% lapply(function(x){
        if(!PValsAnnot) x$Plot = x$Plot + 
            ylim(c(AxisMin,AxisMax))
        
        x$AxisLimits = c(AxisMin,AxisMax)
        x
      })
      
      NoAxisLabels = 1:length(SignatureScorePlot)
      for(i in NoAxisLabels[-(((length(SigList)-1)*length(KineticSets))+1)]) SignatureScorePlot[[i]]$Plot = SignatureScorePlot[[i]]$Plot + theme(axis.text.y = element_blank())
    }
    
    ## Adds P Values Calculated if User Specified (PValsAnnot Variable)
    if(PValsAnnot){
      SignatureScorePlot %<>% lapply(function(x){
        AxisLimit_PVal = (x$AxisLimits[2] - x$AxisLimits[1])*1.15 + x$AxisLimits[1]
        
        if(x$PlotType=="Line_SlopeAnalysis"){
          x$Plot = x$Plot + 
            annotate("text",label = x$PVals,x = 1, y = mean(c(x$AxisLimits[2],AxisLimit_PVal)),hjust = "left",color = colorspace::darken(x$PlotColor,0.2),vjust = "bottom") + 
            ylim(c(x$AxisLimits[1],AxisLimit_PVal))
        } else{
            x$Plot = x$Plot + 
              geom_text(data = x$PVals, mapping = aes(x = group2, y = x$AxisLimits[2], label = PVal_Asterisk),size = PValsSize,color = colorspace::darken(x$PlotColor,0.2),vjust = "bottom") +
              ylim(c(x$AxisLimits[1],AxisLimit_PVal))
        }
        
        return(x)})
    }
    
    ## Combines Plots
    CombinedPlot = lapply(SignatureScorePlot, function(x) x$Plot)
    
    ## Makes Pretty = removes x axis from top plots (just on bottom x axis), adjusts x axis labels
    PlotIndices = 1:length(CombinedPlot)
    for(i in 1:(7*(length(SigList)-1))) CombinedPlot[[i]] = CombinedPlot[[i]] + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
    for(i in (7*(length(SigList)-1)+1):length(CombinedPlot)) CombinedPlot[[i]] = CombinedPlot[[i]] + theme(axis.text.x = element_text(size = 6))
    
    
    plot_grid(plotlist= CombinedPlot, nrow = 3,align = "hv",rel_widths = c(0.7,1,1,1.2,1.4,0.5,0.7))
}

#########################################################################################################
## Combined Plot - plots just Cd8 clusters

SigScorePlot = SignatureScorePlot(KineticLists,DEGList_UpDown,pal = palette,nperm = 100)
SigScorePlot

Scale = 1.9
pdf("3D_Hacohen_BackScoring_v3 2 YF.pdf",width = 6.55*Scale,height = 3*Scale)
plot_grid(plotlist = CombinedPlot,align = "hv",ncol = 7,rel_widths = c(0.7,1,1,1.2,1.4,0.5,0.7))
dev.off()
