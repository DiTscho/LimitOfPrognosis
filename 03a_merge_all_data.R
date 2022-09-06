##### Load libraries ###############################################################################
library("tidyr")
library("dplyr")
library("purrr")
library("broom")
library("ipred")
library("pec")
library("GEOquery")
#### Load data #######################################################################################
load("data.Rdata")
##### Merge clinical and expression ######################################################################################
datasets = c("ExprImageRaw", "ExprImageFil", "METABRIC", "GSE11121", "GSE19615", "GSE22226A", "GSE22226B", "GSE26971", 
             "GSE96058", "GSE1456A","GSE1456B", "GSE42568", "GSE45255", "GSE4922A", "GSE4922B", "GSE7390", 
             "GSE7849", "GSE9195", "GSE9893", "GSE2603", "NKI", "TCGA")

df = tibble(dataset = datasets,
            clinical_data = list(
              data.ExprImage.raw,
              data.ExprImage.fil,
              data.MB.clin,
              data.X1.clin,
              data.X2.clin,
              data.X3A.clin,
              data.X3B.clin,
              data.X4.clin,
              data.X5.clin,
              data.X6A.clin,
              data.X6B.clin,
              data.X7.clin,
              data.X8.clin,
              data.X9A.clin,
              data.X9B.clin,
              data.X10.clin,
              data.X11.clin,
              data.X12.clin,
              data.X13.clin,
              data.X14.clin,
              data.X15.clin,
              data.X16.clin),
            expression_data = list(
              data.ExprImage.raw,
              data.ExprImage.fil,
              data.MB.expr,
              data.X1.expr,
              data.X2.expr,
              data.X3A.expr,
              data.X3B.expr,
              data.X4.expr,
              data.X5.expr,
              data.X6A.expr,
              data.X6B.expr,
              data.X7.expr,
              data.X8.expr,
              data.X9A.expr,
              data.X9B.expr,
              data.X10.expr,
              data.X11.expr,
              data.X12.expr,
              data.X13.expr,
              data.X14.expr,
              data.X15.expr,
              data.X16.expr)
            )

#### Save data #######################
saveRDS(df, "../database/platform.rds")
#df1 = readRDS("../database/platform.rds")





























