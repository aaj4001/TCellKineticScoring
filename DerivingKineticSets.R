## ---------------------------
##
## Purpose of script: Derives Kinetic Sets
##
## Author: Abhi Jaiswal
##
## Date Created: 2023-01-27
##
## Email: ajaiswal1995@gmail.com
##
## Packages:

library(GEOquery)
library(readr)
library(mogene10sttranscriptcluster.db)
library(mogene20sttranscriptcluster.db)
library(mouse4302.db)
library(biomaRt)
library(dplyr)
library(magrittr)


####################################################################################################################################
## Wherry

Wherry = getGEO("GSE41867")
Wherry_Exprs = data.frame(exprs(Wherry[[1]]))

## Standardizes probe calling for microarray datasets!!!
## Any probes that don't map to genes are given NA designation "---"
## For non-NA probes, any probes found mapping to the same gene will be compared
## Only the probe with the highest  median absolute deviation will be retained
## https://support.bioconductor.org/p/70133/

x <- mogene10sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])

Wherry_Exprs <- merge(xx, Wherry_Exprs,by.x = 1, by.y =0,all.y = T)

Wherry_Exprs = Wherry_Exprs %>% mutate(symbol = case_when(is.na(symbol) ~ "---",
                                                          !is.na(symbol) ~ symbol),
                                       MAD = apply(dplyr::select(Wherry_Exprs,starts_with("GSM")),1,mad))

Wherry_Exprs_NAs = Wherry_Exprs %>% subset(symbol == "---")
Wherry_Exprs_nonNAs = Wherry_Exprs %>% subset(symbol != "---") %>% group_by(symbol) %>%
  slice_max(MAD, with_ties = F)

Wherry_Exprs = rbind(Wherry_Exprs_nonNAs, Wherry_Exprs_NAs) %>% dplyr::select(symbol,starts_with("GSM"))


Wherry_Exprs_Matrix = data.matrix(Wherry_Exprs[,2:ncol(Wherry_Exprs)])
rownames(Wherry_Exprs_Matrix) = toupper(Wherry_Exprs$symbol)
colnames(Wherry_Exprs_Matrix) = c("N_1","N_2","N_3","N_4","Arm_D6_1","Arm_D6_2","Arm_D6_3","Arm_D6_4","Arm_D8_1","Arm_D8_2","Arm_D8_3","Arm_D8_4",
                                  "Arm_D15_1","Arm_D15_2","Arm_D15_3","Arm_D15_4","Arm_D30_1","Arm_D30_2","Arm_D30_3","Arm_D30_4","cl13_D6_1","cl13_D6_2",
                                  "cl13_D6_3","cl13_D6_4","cl13_D8_1","cl13_D8_2","cl13_D8_3","cl13_D8_4","cl13_D15_1","cl13_D15_2","cl13_D15_3","cl13_D15_4","cl13_D30_1",
                                  "cl13_D30_2","cl13_D30_3","cl13_D30_4")

Arm = Wherry_Exprs_Matrix[,1:20]
cl13 = Wherry_Exprs_Matrix[,c(1:4,21:36)]

boxplot(Arm)

rm(x,xx,mapped_probes,Wherry,Wherry_Exprs,Wherry_Exprs_Matrix, Wherry_Exprs_NAs, Wherry_Exprs_nonNAs)

####################################################################################################################################
## Sarkar

Sarkar = getGEO("GSE10239", getGPL = T)

ProbeIDs = Sarkar$GSE10239_series_matrix.txt.gz@featureData@data %>% dplyr::select(ID, `Gene Symbol`) %>%
  mutate(symbol = sapply(strsplit(`Gene Symbol`," /// "),function(x) x[1]))

Sarkar_Exprs = data.frame(exprs(Sarkar[[1]]))
Sarkar_Exprs = merge(ProbeIDs,Sarkar_Exprs,by.x = 1,by.y = 0, all.y = T)

## Standardizes probe calling for microarray datasets!!!
## Any probes that don't map to genes are given NA designation "---"
## For non-NA probes, any probes found mapping to the same gene will be compared
## Only the probe with the highest  median absolute deviation will be retained
## https://support.bioconductor.org/p/70133/

Sarkar_Exprs %<>% mutate(symbol = case_when(is.na(symbol) ~ "---",
                                            !is.na(symbol) ~ symbol),
                         MAD = apply(dplyr::select(Sarkar_Exprs,starts_with("GSM")),1,mad))

Sarkar_Exprs_NA = Sarkar_Exprs %>% subset(symbol == "---")
Sarkar_Exprs_nonNA = Sarkar_Exprs %>% subset(symbol != "---") %>% group_by(symbol) %>%
  slice_max(MAD, with_ties = F)

Sarkar_Exprs = rbind(Sarkar_Exprs_nonNA,Sarkar_Exprs_NA)
Sarkar_Exprs = Sarkar_Exprs %>% dplyr::select(symbol, starts_with("GSM"))

Sarkar_Arm = data.matrix(Sarkar_Exprs[,2:ncol(Sarkar_Exprs)])
rownames(Sarkar_Arm) = toupper(Sarkar_Exprs$symbol)
colnames(Sarkar_Arm) = c("N_1","N_2","N_3","Memory_1","Memory_2",
                         "Memory_3","KLRG1int_1","KLRG1int_2","KLRG1int_3",
                         "KLRG1hi_1","KLRG1hi_2","KLRG1hi_3")
Sarkar_Arm = Sarkar_Arm[,c(1:3,10:12,4:6)]

## Need to take log2 of microarray data!!
boxplot(Sarkar_Arm)
Sarkar_Arm = log2(Sarkar_Arm+1)
boxplot(Sarkar_Arm)

rm(Sarkar_Exprs,Sarkar,ProbeIDs, Sarkar_Exprs_nonNA,Sarkar_Exprs_NA)

####################################################################################################################################
## SV40-TAG

GSE89307 <- read_delim("DerivingKineticSets_RawData/GSE89307_FPKM.txt","\t", escape_double = FALSE, trim_ws = TRUE)

SV40 = data.matrix(GSE89307[,3:ncol(GSE89307)])
SV40 = SV40[,-c(1:6,27:35)]
rownames(SV40) = GSE89307$Gene_Symbol
colnames(SV40) = sapply(strsplit(colnames(SV40),split = "_"),function(x) paste(x[-1],collapse = "_"))
SV40 = log2((SV40+1))
SV40 = SV40[rowSums(SV40)>0,]
rownames(SV40) = toupper(rownames(SV40))
SV40 = SV40[,rev(1:ncol(SV40))]

rm(GSE89307)

####################################################################################################################################
## Kupper TRM

Kupper = getGEO("GSE79805")
Kupper_Exprs = data.frame(exprs(Kupper[[1]]))

x <- mogene20sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])

Kupper_Exprs <- merge(xx, Kupper_Exprs,by.x = 1, by.y =0,all.y = T)

## Standardizes probe calling for microarray datasets!!!
## Any probes that don't map to genes are given NA designation "---"
## For non-NA probes, any probes found mapping to the same gene will be compared
## Only the probe with the highest  median absolute deviation will be retained
## https://support.bioconductor.org/p/70133/

Kupper_Exprs %<>% mutate(symbol = case_when(is.na(symbol) ~ "---",
                                            !is.na(symbol) ~ symbol),
                         MAD = apply(dplyr::select(Kupper_Exprs,starts_with("GSM")),1,mad))

Kupper_Exprs_NA = Kupper_Exprs %>% subset(symbol == "---")
Kupper_Exprs_nonNA = Kupper_Exprs %>% subset(symbol != "---") %>% group_by(symbol) %>%
  slice_max(MAD, with_ties = F)

Kupper_Exprs = rbind(Kupper_Exprs_nonNA,Kupper_Exprs_NA)
Kupper_Exprs = Kupper_Exprs %>% dplyr::select(symbol, starts_with("GSM"))

Kupper_VacV = data.matrix(Kupper_Exprs[,2:ncol(Kupper_Exprs)])
rownames(Kupper_VacV) = toupper(Kupper_Exprs$symbol)
colnames(Kupper_VacV) = c("N", "CM", "EM", "d5", "d10", "d15", "d20", "d25", "d30", "d45", "d60", "d90", "OT-I Fasbp4/5 dKO skin infiltrating cells 10 days post infection", "OT-I Fasbp4/5 dKO skin infiltrating cells 30 days post infection")

## Removes irrelevant timepoints
Kupper_VacV = Kupper_VacV[,-c(ncol(Kupper_VacV),ncol(Kupper_VacV)-1)]
Kupper_VacV = Kupper_VacV[,c(2,3,1,4:12)]

boxplot(Kupper_VacV)
rm(Kupper,Kupper_Exprs,x,xx,mapped_probes, Kupper_Exprs_NA, Kupper_Exprs_nonNA)

####################################################################################################################################
## Human Yellow Fever Microarray

GSE26347 = getGEO("GSE26347",AnnotGPL = TRUE)

GSE26347_Eff = exprs(GSE26347[["GSE26347-GPL570_series_matrix.txt.gz"]])[,19:25]

ProbeIDs = GSE26347[["GSE26347-GPL570_series_matrix.txt.gz"]]@featureData@data %>%
  mutate(symbol = sapply(strsplit(`Gene symbol`,"///"),function(x) x[1])) %>%
  select(ID, symbol)

GSE26347_Eff = merge(ProbeIDs,GSE26347_Eff,by.x = 1,by.y = 0)

## Standardizes probe calling for microarray datasets!!!
## Any probes that don't map to genes are given NA designation "---"
## For non-NA probes, any probes found mapping to the same gene will be compared
## Only the probe with the highest  median absolute deviation will be retained
## https://support.bioconductor.org/p/70133/

GSE26347_Eff %<>% mutate(symbol = case_when(is.na(symbol) ~ "---",
                                            !is.na(symbol) ~ symbol),
                         MAD = apply(dplyr::select(GSE26347_Eff,starts_with("GSM")),1,mad))

GSE26347_Eff_NA = GSE26347_Eff %>% subset(symbol == "---")
GSE26347_Eff_nonNA = GSE26347_Eff %>% subset(symbol != "---") %>% group_by(symbol) %>%
  slice_max(MAD, with_ties = F)

GSE26347_Eff = rbind(GSE26347_Eff_nonNA,GSE26347_Eff_NA)
GSE26347_Eff = GSE26347_Eff %>% dplyr::select(symbol, starts_with("GSM"))

GSE26347_Effector = data.matrix(GSE26347_Eff[,-1])
rownames(GSE26347_Effector) = GSE26347_Eff$symbol
colnames(GSE26347_Effector) = c("Naive_BAW","Naive_MRR","Naive_EEB","Effector_BAW","Effector_MRR","Effector_EEB","Effector_BLS")

## Log2 transformation needed - see boxplot before and after
boxplot(GSE26347_Effector)
GSE26347_Effector = log2(GSE26347_Effector)
boxplot(GSE26347_Effector)


rm(GSE26347,GSE26347_Eff, ProbeIDs, GSE26347_Eff_NA, GSE26347_Eff_nonNA)

####################################################################################################################################
## Human Yellow Fever RNASeq

## Yellow Fever Effector FPKM
GSE100745 <- read_delim("./DerivingKineticSets_RawData/GSE100745_cufflinks_fpkm_gene.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

GeneSymbol = getBM(attributes = c("ensembl_gene_id_version","external_gene_name"),
                   filters = "ensembl_gene_id_version",
                   values = GSE100745$tracking_id,
                   mart = mart)

GSE100745 = merge(GeneSymbol,GSE100745,by = 1)

YF_RNASeq = data.matrix(GSE100745[,3:ncol(GSE100745)])
rownames(YF_RNASeq) = GSE100745$external_gene_name
colnames(YF_RNASeq) = c("Naive_AKJ", "Naive_YF653", "Naive_YF661", "Naive_RTA336", "Naive_MCB", "Naive_EVC422", "D14_YF679", 
                        "D14_YF661", "D14_YF631", "Mem_MCB", "Mem_EVC422", "Mem_E07", "Mem_RTA336", "Mem_AKJ")

YF_RNASeq = log2((YF_RNASeq+1))
YF_RNASeq = YF_RNASeq[rowSums(YF_RNASeq)>0,]
rownames(YF_RNASeq) = toupper(rownames(YF_RNASeq))

rm(GeneSymbol, GSE100745, mart)

####################################################################################################################################
## Puts together all Kinetic Sets into one rda file, includes PlotType instructions for how to handle stats and plotting

KineticLists = list(Sarkar_Arm = list(Exprs = Sarkar_Arm, PlotType = "Replicate_Line_LME"), 
                    Arm = list(Exprs = Arm, PlotType = "Replicate_Line_LME"), 
                    cl13 = list(Exprs = cl13, PlotType = "Replicate_Line_LME"), 
                    SV40 = list(Exprs = SV40, PlotType = "Replicate_Line_LME"),
                    Kupper_VacV = list(Exprs = Kupper_VacV, PlotType = "Line_SlopeAnalysis"),
                    YF_Effector = list(Exprs = GSE26347_Effector, PlotType = "Replicate_Line_LME"),
                    YF_RNASeq = list(Exprs = YF_RNASeq, PlotType = "Replicate_Line_TTest"))

save(KineticLists,file = "KineticSets.rda",compression_level = 9)
