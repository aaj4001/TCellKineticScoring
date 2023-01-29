library(RColorBrewer)
library(ggpubr)
library(scales)
library(reshape2)
library(cowplot)
library(emmeans)
library(nlme)
library(msigdbi)

## Imports Query Sets and Kinetic Sets  #########################################################################################################

DEGList = read.gmt("Hacohen_StateDEGs_Human.gmt")$genesets

DEGList_UpDown = list()
for(i in seq(1,6,by = 2)) DEGList_UpDown[[names(DEGList)[i]]] = list(Up = DEGList[[i]],Down = DEGList[[i+1]])
names(DEGList_UpDown) = sapply(strsplit(names(DEGList_UpDown),split = "_",fixed = TRUE),function(x) x[1])

palette = scales::hue_pal()(5)[c(1,3,5)]

load("KineticSets.rda")


## Backscoring Up Minus Down #########################################################################################################

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
  if(!is.null(UpDownList[["Down"]])){
    DnSignature = SignatureScore(InputGeneMatrix,UpDownList[["Down"]])
    NetSignature = UpSignature - DnSignature
  } else{
    NetSignature = UpSignature
  }
  NetSignature
}

SignatureScore_RowInd <- function(InputGeneMatrix,RowInds){
  GeneMatrixFiltered = InputGeneMatrix[RowInds,]
  
  for(i in 1:nrow(GeneMatrixFiltered)) GeneMatrixFiltered[i,] = scale(GeneMatrixFiltered[i,])
  
  MeanZ_Score = colMeans(GeneMatrixFiltered)
  names(MeanZ_Score) = colnames(InputGeneMatrix)
  
  MeanZ_Score
}

SignatureScore_Random <- function(InputGeneMatrix,UpDownList){
  UpGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Up"]]))
  DnGenesRandom = sample(x = 1:nrow(InputGeneMatrix),length(UpDownList[["Down"]]))
  
  UpSignature = SignatureScore_RowInd(InputGeneMatrix,UpGenesRandom)
  if(!is.null(UpDownList[["Down"]])){
    DnSignature = SignatureScore_RowInd(InputGeneMatrix,DnGenesRandom)
    NetRandom = UpSignature - DnSignature
  } else{
    NetRandom = UpSignature
  }
  NetRandom
}

# KineticSets = YF_RNASeq;SigList = DEGList_UpDown;pal = palette;nperm = 100
SignatureScorePlot_TRMSlopeEval <- function(KineticSets,SigList,pal = rep("black",length(SigList)),nperm = 100){
  ScoringResults = list()
  for(i in 1:length(names(SigList))) {
    for(j in names(KineticSets)) {
      SigScore = SignatureScore_UpMinusDown(KineticSets[[j]]$Exprs,SigList[[i]])
      RandScore = matrix(nrow = nperm,ncol = ncol(KineticSets[[j]]$Exprs))
      for(Perms in 1:nperm) RandScore[Perms,] = SignatureScore_Random(KineticSets[[j]]$Exprs,SigList[[i]])
      colnames(RandScore) = colnames(KineticSets[[j]]$Exprs)
      
      ScoringResults[[paste(names(SigList)[i],j,sep = "_")]] = list(SigScore = SigScore,
                                                                    RandScore = RandScore,
                                                                    palette = pal[i],
                                                                    PlotType = KineticSets[[j]]$PlotType)
    }
  }
  ## Plotting Function
  # x = ScoringResults[[1]]
  ScorePlots = lapply(ScoringResults,function(x){
    ## If Experiment has replicates
    if(x$PlotType=="Replicate_Line"){
      Timepoints = unique(sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_")))
      
      ## Signature Score Calculation
      LME_Model = data.frame(Sample = names(x$SigScore),Exp = x$SigScore,ExpBase = "Exp")
      LME_Model$Rep = sapply(strsplit(names(x$SigScore),split = "_"),function(y) y[length(y)])
      LME_Model$Sample = sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_"))
      LME_Model$Sample = factor(LME_Model$Sample,levels = Timepoints)
      
      LME_Sig_Results = lme(Exp~Sample,random=~1|Rep,data = LME_Model)
      LME_Sig_Results = summary(emmeans(LME_Sig_Results,"pairwise"~Sample,adjust = "none"))
      PlottingMatrix_Sig = LME_Sig_Results[["emmeans"]]
      PlottingMatrix_Sig$lower.CL = PlottingMatrix_Sig$emmean - (PlottingMatrix_Sig$SE*1.96)
      PlottingMatrix_Sig$upper.CL = PlottingMatrix_Sig$emmean + (PlottingMatrix_Sig$SE*1.96)
      
      ######### P Value Calculation
      PValues = LME_Sig_Results$contrasts
      
      PVals = PValues[nrow(PValues),6]
      SampleName = as.character(PValues[nrow(PValues),1])
      Index = nrow(PValues)
      
      for(i in 2:(length(unique(LME_Model$Sample))-1)){
        Index = Index - i
        PVals = c(PVals,PValues[Index,6])
        SampleName = c(SampleName,as.character(PValues[Index,1]))
      } 
      
      PVals = c(PVals,1)
      SampleName = c(SampleName,"N")
      names(PVals) = SampleName
      PVals = rev(PVals)
      
      PVals = sapply(PVals,function(x){
        if(x < .0001){
          Sig = "****"
        } else if(x <.001){
          Sig = "***"
        } else if(x <.01){
          Sig = "**"
        } else if(x < .05){
          Sig = "*"
        }else{
          Sig = "ns"
        }
        Sig})
      
      ## Random Score Calculation
      RandScore = list()
      for(i in 1:nperm) RandScore[[i]] = x$RandScore[i,]
      LME_Rand_Results = sapply(RandScore,function(x) {
        LME_Model = data.frame(Sample = names(x),Exp = x,ExpBase = "Rand")
        LME_Model$Rep = sapply(strsplit(names(x),split = "_"),function(y) y[length(y)])
        LME_Model$Sample = sapply(strsplit(names(x),split = "_"),function(y) paste(y[-length(y)],collapse = "_"))
        LME_Model$Sample = factor(LME_Model$Sample,levels = Timepoints)
        
        LME_Rand_Results = lme(Exp~Sample,random=~1|Rep,data = LME_Model)
        LME_Rand_Results = summary(emmeans(LME_Rand_Results,"pairwise"~Sample,adjust = "none"))[["emmeans"]]$emmean
        names(LME_Rand_Results) = Timepoints
        LME_Rand_Results
      })
      RandMean = rowMeans(LME_Rand_Results)
      RandSE = numeric(length(RandMean)); for(i in 1:length(RandMean)) RandSE[i] = plotrix::std.error(LME_Rand_Results[i,])
      PlottingMatrix_Rand = data.frame(Sample = factor(Timepoints,Timepoints),emmean = RandMean,SE = RandSE, 
                                       lower.CL = RandMean-(RandSE*1.96), upper.CL = RandMean+(RandSE*1.96))
      
      ######### Plotting
      Plot = ggplot(data = PlottingMatrix_Rand,aes(x = Sample,y = emmean,ymin = lower.CL,ymax = upper.CL,group = 1)) +
        geom_ribbon(fill = "grey") +
        geom_line(color = "black") + 
        geom_line(data = PlottingMatrix_Sig,color = x$palette) +
        geom_errorbar(data = PlottingMatrix_Sig,color = x$palette,width=0.4, size=0.5) +
        geom_point(data = PlottingMatrix_Sig,fill = x$palette,shape = 21,size = 4) +
        theme_cowplot() +
        theme(axis.title = element_blank())
      
      
      Results = list(PVals = PVals,Plot = Plot)
      
      ## Else if there are no replicate samples (EG Kupper)  
    }else if(x$PlotType=="Slope_Line"){
      ######### P Values
      
      ## Signature Slope
      Signature_TRMCourse = x$SigScore[-(1:2)]
      Signature_TRMCourse = data.frame(ZScore = Signature_TRMCourse,
                                       Sample = names(Signature_TRMCourse))
      Signature_TRMCourse$Sample = as.numeric(c("0","5","10","15","20","25","30","45","60","90"))
      
      
      LmFit = lm(ZScore ~ Sample,data = Signature_TRMCourse)
      
      SlopeResults = list(DeltaExpr = summary(LmFit)$coefficients[2,1]*90,
                          PVals = summary(LmFit)$coefficients[2,4])
      
      ######### Plotting
      Timepoints = factor(names(x$SigScore),levels = names(x$SigScore))
      PlottingMatrix_Sig = data.frame(Timepoints,Means = x$SigScore,
                                      Timepoints_v2 = c(-20,-10,0,5,10,15,20,25,30,45,60,90))
      RandMean = colMeans(x$RandScore)
      RandSE = numeric(length(RandMean)); for(i in 1:length(RandSE)) RandSE[i] = plotrix::std.error(x$RandScore[,i])
      
      PlottingMatrix_Rand = data.frame(Timepoints,Means = RandMean, SE = RandSE,
                                       UpperCL = RandMean+(RandSE*1.96),LowerCL = RandMean-(RandSE*1.96),
                                       Timepoints_v2 = c(-20,-10,0,5,10,15,20,25,30,45,60,90))
      
      # Plot = ggplot(data = PlottingMatrix_Rand,aes(x = Timepoints_v2,y = Means,group = 1)) + 
      #   scale_x_continuous(breaks = c(-20,-10,0,5,10,15,20,25,30,45,60,90),
      #                      labels = names(x$SigScore)) +
      #   geom_ribbon(fill = "grey",mapping = aes(ymin = LowerCL, ymax = UpperCL)) +
      #   geom_line(color = "black",mapping = aes(ymin = LowerCL, ymax = UpperCL)) + 
      #   geom_smooth(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],method = "lm",color = x$palette) +
      #   geom_line(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],color = x$palette) +
      #   geom_point(data = PlottingMatrix_Sig,fill = x$palette,shape = 21,size = 4) +
      #   theme_cowplot() +
      #   theme(axis.title = element_blank())
      
      Plot = ggplot(data = PlottingMatrix_Rand,aes(x = Timepoints,y = Means,group = 1)) + 
        # scale_x_continuous(breaks = c(-20,-10,0,5,10,15,20,25,30,45,60,90),
        #                    labels = names(x$SigScore)) +
        geom_ribbon(fill = "grey",mapping = aes(ymin = LowerCL, ymax = UpperCL)) +
        geom_line(color = "black",mapping = aes(ymin = LowerCL, ymax = UpperCL)) + 
        geom_smooth(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],method = "lm",color = "transparent") +
        geom_line(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],color = x$palette) +
        geom_point(data = PlottingMatrix_Sig,fill = x$palette,shape = 21,size = 4) +
        geom_smooth(data = PlottingMatrix_Sig[!Timepoints%in%c("CM","EM"),],method = "lm",color = colorspace::darken(x$palette,0.2),fill = "transparent") +
        theme_cowplot() +
        theme(axis.title = element_blank())
      
      
      Results = list(Plot = Plot,PVals = SlopeResults)
    }else if(x$PlotType=="Replicate_Line_TTest"){
      Timepoints = unique(sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_")))
      
      ##########################
      ## Stats
      
      Timepoint_Stat = factor(sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_")),Timepoints)
      PVals = pairwise.t.test(x$SigScore,Timepoint_Stat,p.adj = "none")$p.value
      
      PVal = as.numeric(PVals)
      Names = character()
      for(i in 1:length(PVal)){
        ArrInd = which(PVals==PVals[i],arr.ind = TRUE)
        Names[i] = paste0(rownames(PVals)[ArrInd[1]]," - ",colnames(PVals)[ArrInd[2]])
      }
      names(PVal) = Names
      
      
      ##########################
      ## Random Calculation
      RandScore = list()
      for(i in 1:nperm) RandScore[[i]] = x$RandScore[i,]
      RandScore = sapply(RandScore,function(rand){
        RandMean = sapply(Timepoints,function(A){
          Means = mean(rand[grepl(A,names(rand))])
        })
      })
      
      Means = rowMeans(RandScore)
      SE = numeric(length(Means))
      for(i in 1:length(Means)) SE[i] = plotrix::std.error(RandScore[i,])
      
      RandScore = data.frame(Timepoints = factor(names(Means),Timepoints),Means = Means,SE = SE,
                             UpperCL = Means+(SE*1.96),LowerCL = Means-(SE*1.96))
      
      ##########################
      ## Signature Calculation
      
      PlotSig = data.frame(t(sapply(Timepoints,function(A){
        Means = mean(x$SigScore[grepl(A,names(x$SigScore))])
        SE = plotrix::std.error(x$SigScore[grepl(A,names(x$SigScore))])
        
        UpperCL = Means+(SE*1.96)
        LowerCL = Means-(SE*1.96)
        
        rbind(Means,SE,UpperCL,LowerCL)
      })))
      
      colnames(PlotSig) = cbind("Means","SE","UpperCL","LowerCL")
      PlotSig$Timepoints = factor(rownames(PlotSig),Timepoints)
      
      Plot = ggplot(data = PlotSig,aes(x = Timepoints,y = Means,ymin = LowerCL,ymax = UpperCL,group = 1)) + 
        geom_ribbon(data = RandScore,fill = "grey") +
        geom_line(data = RandScore,color = "black") + 
        geom_line(color = x$palette) +
        geom_point(fill = x$palette,shape = 21,size = 4) +
        geom_errorbar(color = x$palette,width=0.4, size=0.5) + 
        theme_cowplot() +
        theme(axis.title = element_blank())
      
      
      
      Results = list(Plot = Plot,PVals = PVal)
    }else{
      Timepoints = unique(sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_")))
      
      ## Signature Score Calculation
      LME_Model = data.frame(Sample = names(x$SigScore),Exp = x$SigScore,ExpBase = "Exp")
      LME_Model$Rep = sapply(strsplit(names(x$SigScore),split = "_"),function(y) y[length(y)])
      LME_Model$Sample = sapply(strsplit(names(x$SigScore),split = "_"),function(y) paste(y[-length(y)],collapse = "_"))
      LME_Model$Sample = factor(LME_Model$Sample,levels = Timepoints)
      
      LME_Sig_Results = lme(Exp~Sample,random=~1|Rep,data = LME_Model)
      LME_Sig_Results = summary(emmeans(LME_Sig_Results,"pairwise"~Sample,adjust = "none"))
      PlottingMatrix_Sig = LME_Sig_Results[["emmeans"]]
      PlottingMatrix_Sig$lower.CL = PlottingMatrix_Sig$emmean - (PlottingMatrix_Sig$SE*1.96)
      PlottingMatrix_Sig$upper.CL = PlottingMatrix_Sig$emmean + (PlottingMatrix_Sig$SE*1.96)      
      
      ######### P Value Calculation
      PValues = LME_Sig_Results$contrasts
      
      PVals = PValues[nrow(PValues),6];names(PVals) = levels(PValues$contrast)[1]
      
      PVals = sapply(PVals,function(x){
        if(x < .0001){
          Sig = "****"
        } else if(x <.001){
          Sig = "***"
        } else if(x <.01){
          Sig = "**"
        } else if(x < .05){
          Sig = "*"
        }else{
          Sig = "ns"
        }
        Sig})      
      
      ## Random Score Calculation
      RandScore = list()
      x$RandScore = x$RandScore[!is.nan(x$RandScore[,1]),]
      for(i in 1:nrow(x$RandScore)) RandScore[[i]] = x$RandScore[i,]
      LME_Rand_Results = sapply(RandScore,function(x) {
        LME_Model = data.frame(Sample = names(x),Exp = x,ExpBase = "Rand")
        LME_Model$Rep = sapply(strsplit(names(x),split = "_"),function(y) y[length(y)])
        LME_Model$Sample = sapply(strsplit(names(x),split = "_"),function(y) paste(y[-length(y)],collapse = "_"))
        LME_Model$Sample = factor(LME_Model$Sample,levels = Timepoints)
        
        LME_Rand_Results = lme(Exp~Sample,random=~1|Rep,data = LME_Model)
        LME_Rand_Results = summary(emmeans(LME_Rand_Results,"pairwise"~Sample,adjust = "none"))[["emmeans"]]$emmean
        names(LME_Rand_Results) = Timepoints
        LME_Rand_Results
      })
      RandMean = rowMeans(LME_Rand_Results)
      RandSE = numeric(length(RandMean)); for(i in 1:length(RandMean)) RandSE[i] = plotrix::std.error(LME_Rand_Results[i,])
      PlottingMatrix_Rand = data.frame(Sample = factor(Timepoints,Timepoints),emmean = RandMean,SE = RandSE, 
                                       lower.CL = RandMean-(RandSE*1.96), upper.CL = RandMean+(RandSE*1.96))
      
      ######### Plotting
      Plot = ggplot(data = PlottingMatrix_Rand,aes(x = Sample,y = emmean,ymin = lower.CL,ymax = upper.CL,group = 1)) +
        geom_ribbon(fill = "grey") +
        geom_line(color = "black") + 
        geom_line(data = PlottingMatrix_Sig,color = x$palette) +
        geom_errorbar(data = PlottingMatrix_Sig,color = x$palette,width=0.4, size=0.5) +
        geom_point(data = PlottingMatrix_Sig,fill = x$palette,shape = 21,size = 4) +
        theme_cowplot() +
        theme(axis.title = element_blank())
      
      Results = list(PVals = PVals,Plot = Plot)
      
    }})
  
  nKinetics = 1:length(KineticSets)
  PValues = lapply(ScorePlots,function(x) x$PVals)
  PVals = lapply(nKinetics,function(x) t(sapply(PValues[seq(x,length(ScoringResults),by = length(KineticSets))],function(x) x)))
  names(PVals) = names(KineticSets)
  
  print(PVals)
  
  
  
  ScorePlot = lapply(ScorePlots,function(x) x$Plot)
  
  ScorePlot
}

#########################################################################################################
## Combined Plot - plots just Cd8 clusters

SigScorePlot = SignatureScorePlot_TRMSlopeEval(KineticLists,DEGList_UpDown,pal = palette)
CombinedPlot = SigScorePlot

for(i in 1:(7*length(DEGList_UpDown))) CombinedPlot[[i]] = CombinedPlot[[i]] + ylim(c(-2,2.4)) + theme(axis.text = element_text(size = 6))
for(i in 1:(7*(length(DEGList_UpDown)-1))) CombinedPlot[[i]] = CombinedPlot[[i]] + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
for(i in seq(1,(7*length(DEGList_UpDown)))[-seq(1,(7*length(DEGList_UpDown)),7)]) CombinedPlot[[i]] = CombinedPlot[[i]] + theme(axis.text.y = element_blank())

Scale = 1.9
pdf("3D_Hacohen_BackScoring_v3 2 YF.pdf",width = 6.55*Scale,height = 3*Scale)
plot_grid(plotlist = CombinedPlot,align = "hv",ncol = 7,rel_widths = c(0.7,1,1,1.2,1.4,0.5,0.7))
dev.off()
