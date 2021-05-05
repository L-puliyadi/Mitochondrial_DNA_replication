library("BiocManager")

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

#1without POLG in all
name <- c('POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

    for(g in genes) {
      n <- n+1

    string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
    commandsRun(string.cmd)
    assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
   
    deleteTableColumn(column = 'stringdb::structures', table='node')
    deleteTableColumn(column = 'stringdb::sequence', table='node')
    deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
    deleteTableColumn(column = 'stringdb::STRING style', table='node')
    
    # apply visual style ()
    setVisualStyle('Heatpropogation')
    
    clearSelection()
    selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
    
    Input <- getSelectedNodes()
    setNodeShapeBypass(Input,'diamond')
    }

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_POLG_150_0.5.txt", sep="\t")

#2 WITHOUT POLG2
name <- c('POLG', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_POLG2_150_0.5.txt", sep="\t")

#3 WITHOUT MGME11
name <- c('POLG', 'POLG2', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_MGME1_150_0.5.txt", sep="\t")

#4 WITHOUT SSBP1
name <- c('POLG', 'POLG2', 'MGME1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_SSBP1_150_0.5.txt", sep="\t")

#5 WITHOUT C10ORF2
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_C10ORF2_150_0.5.txt", sep="\t")

#6 WITHOUT TOP1MT
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_TOP1MT_150_0.5.txt", sep="\t")

#7 WITHOUT RNASEH1
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_RNASEH1_150_0.5.txt", sep="\t")

#8 WITHOUT PRIMPOL
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_PRIMPOL_150_0.5.txt", sep="\t")

#9 WITHOUT EXOG
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_EXOG_150_0.5.txt", sep="\t")

#10 WITHOUT ATAD3A
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_ATAD3A_150_0.5.txt", sep="\t")

#11 WITHOUT RRM2B
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_RRM2B_150_0.5.txt", sep="\t")

#12 WITHOUT TK2
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_TK2_150_0.5.txt", sep="\t")

#13 WITHOUT DGUOK
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, POLRMT, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_DGUOK_150_0.5.txt", sep="\t")

#14 WITHOUT POLRMT
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'TFAM', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, TFAM, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17, 'MPV17')

write.table(data_total, "leavetwo_POLRMT_150_0.5.txt", sep="\t")

#15 WITHOUT TFAM
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TEFM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TEFM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_TFAM_150_0.5.txt", sep="\t")

#16 WITHOUT TEFM
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFB2M, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_TEFM_150_0.5.txt", sep="\t")

#17 WITHOUT TFB2M
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'MPV17')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, MPV17",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_TFB2M_150_0.5.txt", sep="\t")

#18without MPV17
name <- c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M')
genes <- c("POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, DGUOK, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, POLRMT, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, TFAM, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TEFM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TFB2M",
           "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM")

cutoffs <- "0.5"
limits <- "150"
n <- 0

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

for(g in genes) {
  n <- n+1
  
  string.cmd = paste('string protein query query=\"',g,'\" cutoff=\"',cutoffs,'\" species="Homo sapiens" limit=\"',limits,'\" newNetName=\"STRING',n,'\"', sep='')
  commandsRun(string.cmd)
  assign(paste("data_",n, sep = ""), getTableColumns('node', 'display name'))
  
  deleteTableColumn(column = 'stringdb::structures', table='node')
  deleteTableColumn(column = 'stringdb::sequence', table='node')
  deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
  deleteTableColumn(column = 'stringdb::STRING style', table='node')
  
  # apply visual style ()
  setVisualStyle('Heatpropogation')
  
  clearSelection()
  selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TFB2M', 'MPV17'), by.col='query term')
  
  Input <- getSelectedNodes()
  setNodeShapeBypass(Input,'diamond')
}

data_total <- cbind(data_1, data_2, data_3, data_4, data_5, data_6, data_7, data_8, data_9, data_10, data_11, data_12, data_13, data_14, data_15, data_16, data_17)

write.table(data_total, "leavetwo_MPV17_150_0.5.txt", sep="\t")