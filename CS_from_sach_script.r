#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(glmnet)
library(dplyr)


load("~/Documents/scripts/scripts_for_ingrid/SensuStricto_glmnet_classifier_object.R")
#load(args[2])
#df = read.table(args[1], header=TRUE)
df <- read.table("~/Documents/scripts/scripts_for_ingrid/final_dataframe.csv", header=TRUE)
df <- unique(df)
rownames(df) = make.names(df$Names, unique=TRUE)

df <- select(df, -c("tag", "tga", "taa","Names"))

df_filt = data.frame(
  row.names = rownames(df),
  df[,c(1:61)], 
  Length = df$GENELENGTH, 
  CAI = df$CAICLASS, 
  GC= df$GC_CONT, 
  GC3 = df$GC_CONT3,  
  Cost = df$Cost,
  perTM = (df$Exp_AA_in_transm)/(df$GENELENGTH/3),
  disord = (df$Disordered)/(df$GENELENGTH/3),
  Low_comp = (df$Complexity)/(df$GENELENGTH/3),
  Gravy = df$Gravy,
  Aromo = df$Aromo,
  Helix = df$H/(df$GENELENGTH/3),
  Sheet = df$E/(df$GENELENGTH/3),
  Aggregation = df$Aggregation/(df$GENELENGTH/3))	



df_filt$GC_GC3 = abs((df_filt$GC3-df_filt$GC))*abs((df_filt$GC3-0.5))

df_filt = as.data.frame(scale(df_filt[,!(names(df_filt) %in% c("label"))],
                                center=means_train,
                                scale=sd_train))

df_filt = df_filt[, !(names(df_filt) %in% c("GC", "GC3","Length"))]

results_pred_class = predict.cv.glmnet(glmcv__complete,
                                        as.matrix(df_filt),
                                        type="class",
                                        s="lambda.1se")

results_pred_score = predict.cv.glmnet(glmcv__complete, 
                                       as.matrix(df_filt),
                                       type="response",
                                       s="lambda.1se")

pdf("CS_histogram.pdf")
hist(results_pred_score, breaks=100)
dev.off()
