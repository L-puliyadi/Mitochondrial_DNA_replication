if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
if(!"igraph" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("igraph")
}
if(!"rstudioapi" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("rstudioapi")
}
library(RCy3)
library(igraph)
library(rstudioapi)

cytoscapePing()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.4"
l <- "50"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

    for(g in genes) {
      n <- n+1
    string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit= "0" newNetName=\"STRING',n,'\"', sep='')
    commandsRun(string.cmd)
    string2.cmd <- paste('string expand additionalNodes=\"',l,'\" network=\"CURRENT\"')
    commandsRun(string2.cmd)
    assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
    
    deleteTableColumn(column = 'stringdb::structures', table='node')
    deleteTableColumn(column = 'stringdb::sequence', table='node')
    deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
    deleteTableColumn(column = 'stringdb::STRING style', table='node')
    
    # apply visual style ()
    setVisualStyle('Heatpropogation')
    }

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17, data_18, data_19, data_20, data_21, data_22)

write.table(data_total, "leave_one_50_0.4.txt", sep="\t")