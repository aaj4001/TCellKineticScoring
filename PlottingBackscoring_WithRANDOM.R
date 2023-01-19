library(ggplot2)
library(ggpubr)

ggarrange(readRDS("./Plotting_WITHRANDOM/Sarkar.rds"),
          readRDS("./Plotting_WITHRANDOM/Wherry_Arm.rds"),
          readRDS("./Plotting_WITHRANDOM/Wherry_Cl13.rds"),
          readRDS("./Plotting_WITHRANDOM/Schei.rds"),
          readRDS("./Plotting_WITHRANDOM/Kupper.rds"),
          ncol = 5)
