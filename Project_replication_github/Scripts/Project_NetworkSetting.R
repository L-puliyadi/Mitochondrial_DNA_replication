library("BiocManager")

if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
if(!"igraph" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("igraph")
}
if(!"clusterProfiler" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
library(RCy3)
library(igraph)
library(clusterProfiler)

cytoscapePing()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

genes <- "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17" 

cutoffs <- c("0.4", "0.5", "0.6", "0.7")
limits <- c("50","100","150","200")

#place individually in case it takes long time
#cutoffs <- c("0.5")
#limits <- c("150")

# load visual style from file
vizstyle.file <- file.path(getwd(), "Dynamics.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

#before changing command
#string.cmd = paste('string protein query query=\"',genes,'\" cutoff=\"',value,'\" species="Homo sapiens" limit=\"',l,'\" newNetName=\"STRING',l,'_',value,'\"', sep='')
#commandsRun(string.cmd)
#assign(paste("data_",l,'_',value, sep = ""), getTableColumns('node', c('display name', 'query term','compartment::mitochondrion')))

for(l in limits) {
  for (value in cutoffs) {
    string.cmd = paste('string protein query query=\"',genes,'\" cutoff=\"',value,'\" species="Homo sapiens" limit="0" newNetName=\"STRING',l,'_',value,'\"', sep='')
    commandsRun(string.cmd)
    string2.cmd <- paste('string expand additionalNodes=\"',l,'\" network=\"CURRENT\"')
    commandsRun(string2.cmd)
    assign(paste("data_",l,'_',value, sep = ""), getTableColumns('node', c('display name', 'query term','compartment::mitochondrion')))
    
    deleteTableColumn(column = 'stringdb::structures', table='node')
    deleteTableColumn(column = 'stringdb::sequence', table='node')
    deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
    deleteTableColumn(column = 'stringdb::STRING style', table='node')
    
    # apply visual style ()
    setVisualStyle('Dynamics')
    
    #Select input nodes
    clearSelection()
    selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M'), by.col='query term')
    Input <- getSelectedNodes()
    setNodeShapeBypass(Input,'diamond')
    
    #heat propagation
    #https://mkutmon.gitlab.io/tutorial-network-algorithms/NetworkPropagation.html
    diffusionBasic()
    clearSelection()
    setNodeColorMapping("diffusion_output_heat", c(0,0.5), c('#FFFFFF', '#FF0000'), style.name = "Dynamics")
  }
}
