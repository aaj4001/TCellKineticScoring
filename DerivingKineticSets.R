library(GEOquery)
library(readr)
library(mogene10sttranscriptcluster.db)
library(mogene20sttranscriptcluster.db)
library(mouse4302.db)

## Imports Datasets #########################################################################################################
## Wherry
Wherry = getGEO("GSE41867")
Wherry_Exprs = data.frame(exprs(Wherry[[1]]))

x <- mogene10sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])

Wherry_Exprs <- merge(xx, Wherry_Exprs,by.x = 1, by.y =0,all.y = T)

Wherry_Exprs_Matrix = data.matrix(Wherry_Exprs[,3:ncol(Wherry_Exprs)])
rownames(Wherry_Exprs_Matrix) = toupper(Wherry_Exprs$symbol)
colnames(Wherry_Exprs_Matrix) = c("N_1","N_2","N_3","N_4","Arm_D6_1","Arm_D6_2","Arm_D6_3","Arm_D6_4","Arm_D8_1","Arm_D8_2","Arm_D8_3","Arm_D8_4",
                                  "Arm_D15_1","Arm_D15_2","Arm_D15_3","Arm_D15_4","Arm_D30_1","Arm_D30_2","Arm_D30_3","Arm_D30_4","cl13_D6_1","cl13_D6_2",
                                  "cl13_D6_3","cl13_D6_4","cl13_D8_1","cl13_D8_2","cl13_D8_3","cl13_D8_4","cl13_D15_1","cl13_D15_2","cl13_D15_3","cl13_D15_4","cl13_D30_1",
                                  "cl13_D30_2","cl13_D30_3","cl13_D30_4")

Arm = Wherry_Exprs_Matrix[,1:20]
cl13 = Wherry_Exprs_Matrix[,c(1:4,21:36)]

boxplot(Arm)

rm(x,xx,mapped_probes,Wherry,Wherry_Exprs,Wherry_Exprs_Matrix)

## Sarkar
Sarkar = getGEO("GSE10239")
Sarkar_Exprs = data.frame(exprs(Sarkar[[1]]))

## Annotation
x <- mouse4302SYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])

Sarkar_Exprs = merge(xx,Sarkar_Exprs,by.x = 1,by.y = 0, all.y = T)

Sarkar_Arm = data.matrix(Sarkar_Exprs[,3:ncol(Sarkar_Exprs)])
rownames(Sarkar_Arm) = toupper(Sarkar_Exprs$symbol)
colnames(Sarkar_Arm) = c("N_1","N_2","N_3","Memory_1","Memory_2",
                         "Memory_3","KLRG1int_1","KLRG1int_2","KLRG1int_3",
                         "KLRG1hi_1","KLRG1hi_2","KLRG1hi_3")

Sarkar_Arm = Sarkar_Arm[,c(1:3,10:12,4:6)]
Sarkar_Arm = log2(Sarkar_Arm+1)
boxplot(Sarkar_Arm)

rm(Sarkar_Exprs,Sarkar,x,xx,mapped_probes)

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

## Kupper TRM
Kupper = getGEO("GSE79805")
Kupper_Exprs = data.frame(exprs(Kupper[[1]]))

x <- mogene20sttranscriptclusterSYMBOL
mapped_probes <- mappedkeys(x)
xx <- as.data.frame(x[mapped_probes])

Kupper_Exprs <- merge(xx, Kupper_Exprs,by.x = 1, by.y = 0,all.y=T)

Kupper_VacV = data.matrix(Kupper_Exprs[,3:ncol(Kupper_Exprs)])
rownames(Kupper_VacV) = toupper(Kupper_Exprs$symbol)
colnames(Kupper_VacV) = c("N", "CM", "EM", "d5", "d10", "d15", "d20", "d25", "d30", "d45", "d60", "d90", "OT-I Fasbp4/5 dKO skin infiltrating cells 10 days post infection", "OT-I Fasbp4/5 dKO skin infiltrating cells 30 days post infection")

## Removes irrelevant timepoints
Kupper_VacV = Kupper_VacV[,-c(ncol(Kupper_VacV),ncol(Kupper_VacV)-1)]
Kupper_VacV = Kupper_VacV[,c(2,3,1,4:12)]

boxplot(Kupper_VacV)
rm(Kupper,Kupper_Exprs,x,xx,mapped_probes)

KineticLists = list(Sarkar_Arm = Sarkar_Arm, 
                    Arm = Arm, 
                    cl13 = cl13, 
                    SV40 = SV40,
                    Kupper_VacV = Kupper_VacV)

save(KineticLists,file = "KineticSets.rda",compression_level = 9)
