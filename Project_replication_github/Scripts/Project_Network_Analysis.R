# Cytoscape version 3.8.2
# CyTargetLinker version 4.1.0
# stringApp version 1.6.0
# yFiles Organic Layout (app version 1.1).  

# R version 3.6.0
# Rstudio version 1.4.1106
# Rcy3 package version 2.4.0
# igraph package version 1.2.6
# rstudioapi package version 0.13
# Primer3plus version 2.4.2

#########################################################ENCODE#######################################################
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

#Put script and additional files in the same folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Create network
genes <- "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17" 
value <- "0.5"
l <- "150"

# load visual style from file
vizstyle.file <- file.path(getwd(), "Heatpropogation.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)


string.cmd = paste('string protein query query=\"',genes,'\" cutoff=\"',value,'\" species="Homo sapiens" limit="0" newNetName=\"STRING',l,'_',value,'\"', sep='')
commandsRun(string.cmd)
string2.cmd <- paste('string expand additionalNodes=\"',l,'\" network=\"CURRENT\"')
commandsRun(string2.cmd)

deleteTableColumn(column = 'stringdb::structures', table='node')
deleteTableColumn(column = 'stringdb::sequence', table='node')
deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
deleteTableColumn(column = 'stringdb::STRING style', table='node')

# apply visual style ()
setVisualStyle('Dynamics')
#renameNetwork('STRING_MitoCarta')

#Select input nodes
clearSelection()
selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17'), by.col='query term')

Input <- getSelectedNodes()
setNodeShapeBypass(Input,'diamond')

# heat propagation
# https://mkutmon.gitlab.io/tutorial-network-algorithms/NetworkPropagation.html
diffusionBasic()
clearSelection()
setNodeColorMapping("diffusion_output_heat", c(0,0.5), c('#FFFFFF', '#FF0000'), style.name = "Dynamics")
assign(paste("data_",l,'_',value, sep = ""), getTableColumns('node', c('display name', 'diffusion_output_heat', 'diffusion_output_rank')))

#Select one added protein
connections <- c("ACOT4", "AK2", "AK3", "AK4", "AK5", "ATAD3B", "AUH", "BTBD19", "C21orf91", "CCDC58", "CDA", "CMPK1", "CMPK2", "COLEC11", "COX4I1", "CRLF3", "CROT", "DARS2", "DCTD", "DCTPP1", "DDX28", "DHX30", "DNA2", "DTYMK", "DUT", "ENSP00000433415", "ENTPD1", "ENTPD3", "ENTPD8", "EPRS", "ERAL1", "EXO5", "FAM178A", "FASTKD5", "FEN1", "GABPA", "GABPB1", "GAL3ST3", "GARS", "GATC", "GKAP1", "GPR146", "GRSF1", "GSTK1", "GUK1", "IMMT", "INIP", "INTS3", "KIAA0020", "KLHL35", "LIG1", "LLPH", "LONP1", "LRPPRC", "LRRC37B", "MCM5", "MCM9", "MFN2", "MPV17L", "MPV17L2", "MRPL12", "MRPL21", "MT-ATP6", "MT-CO1", "MT-CO2", "MT-CYB", "MTERF1", "MTERF2", "MTERF3", "MTERF4", "MT-ND4", "MT-ND6", "MTPAP", "NDUFA5", "NDUFB8", "NDUFS3", "NIFK", "NIP7", "NKRF", "NME4", "NPTX2", "NRF1", "NSUN4", "NT5C", "NT5C1A", "NT5C1B", "NT5C2", "NT5C3A", "NT5C3B", "NT5M", "NTPCR", "NUDT19", "OXA1L", "PCNA", "PDE12", "PEAK1", "PEX10", "PEX12", "PEX13", "PEX14", "PEX2", "PIF1", "PNP", "POLDIP2", "POLE", "POLQ", "PPARGC1A", "PPRC1", "PRIM1", "PRIM2", "PTCD1", "PTCD2", "PTCD3", "RARS2", "RBM28", "REV3L", "RFWD3", "RHCG", "RNASEH2A", "RNASEH2C", "RPA1", "RPA2", "RPL30", "RRM2", "SAMHD1", "SLC25A4", "SSBP2", "SUCLA2", "SUCLG1", "SUPV3L1", "SURF1", "TACO1", "TBL3", "TFB1M", "TICRR", "TIMELESS", "TK1", "TMEM242", "TMEM70", "TOMM20", "TP53AIP1", "TRAPPC12", "TRMT10C", "TRMU", "TSSC1", "TYMP", "UCK1", "UCN", "UTP6", "ZUFSP")

for (c in connections){
selectNodes(c, by.col='display name')
  selectFirstNeighbors(direction = "any")
  assign(paste("degree_",c, sep = ""), getSelectedNodeCount())
  n <- 0
  count <- getSelectedNodes()
  
  
  if ("9606.ENSP00000268124" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000442563" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000366939" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000419665" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000309595" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000328835" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000313350" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000313816" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000287675" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000368030" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000251810" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000299697" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000264093" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000369383" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000465759" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000420588" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000462963" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000355471" %in% count) {
    n <- n+1
  }
  
  assign(paste("count_",c, sep = ""), n)
  
  clearSelection()
} 


Degree_0.5_150 <- rbind(degree_ACOT4, degree_AK2, degree_AK3, degree_AK4, degree_AK5, degree_ATAD3B, degree_AUH, degree_BTBD19, degree_C21orf91, degree_CCDC58, degree_CDA, degree_CMPK1, degree_CMPK2, degree_COLEC11, degree_COX4I1, degree_CRLF3, degree_CROT, degree_DARS2, degree_DCTD, degree_DCTPP1, degree_DDX28, degree_DHX30, degree_DNA2, degree_DTYMK, degree_DUT, degree_ENSP00000433415, degree_ENTPD1, degree_ENTPD3, degree_ENTPD8, degree_EPRS, degree_ERAL1, degree_EXO5, degree_FAM178A, degree_FASTKD5, degree_FEN1, degree_GABPA, degree_GABPB1, degree_GAL3ST3, degree_GARS, degree_GATC, degree_GKAP1, degree_GPR146, degree_GRSF1, degree_GSTK1, degree_GUK1, degree_IMMT, degree_INIP, degree_INTS3, degree_KIAA0020, degree_KLHL35, degree_LIG1, degree_LLPH, degree_LONP1, degree_LRPPRC, degree_LRRC37B, degree_MCM5, degree_MCM9, degree_MFN2, degree_MPV17L, degree_MPV17L2, degree_MRPL12, degree_MRPL21, degree_MTERF1, degree_MTERF2, degree_MTERF3, degree_MTERF4, degree_MTPAP, degree_NDUFA5, degree_NDUFB8, degree_NDUFS3, degree_NIFK, degree_NIP7, degree_NKRF, degree_NME4, degree_NPTX2, degree_NRF1, degree_NSUN4, degree_NT5C, degree_NT5C1A, degree_NT5C1B, degree_NT5C2, degree_NT5C3A, degree_NT5C3B, degree_NT5M, degree_NTPCR, degree_NUDT19, degree_OXA1L, degree_PCNA, degree_PDE12, degree_PEAK1, degree_PEX10, degree_PEX12, degree_PEX13, degree_PEX14, degree_PEX2, degree_PIF1, degree_PNP, degree_POLDIP2, degree_POLE, degree_POLQ, degree_PPARGC1A, degree_PPRC1, degree_PRIM1, degree_PRIM2, degree_PTCD1, degree_PTCD2, degree_PTCD3, degree_RARS2, degree_RBM28, degree_REV3L, degree_RFWD3, degree_RHCG, degree_RNASEH2A, degree_RNASEH2C, degree_RPA1, degree_RPA2, degree_RPL30, degree_RRM2, degree_SAMHD1, degree_SLC25A4, degree_SSBP2, degree_SUCLA2, degree_SUCLG1, degree_SUPV3L1, degree_SURF1, degree_TACO1, degree_TBL3, degree_TFB1M, degree_TICRR, degree_TIMELESS, degree_TK1, degree_TMEM242, degree_TMEM70, degree_TOMM20, degree_TP53AIP1, degree_TRAPPC12, degree_TRMT10C, degree_TRMU, degree_TSSC1, degree_TYMP, degree_UCK1, degree_UCN, degree_UTP6, degree_ZUFSP)
Count_0.5_150 <- rbind(count_ACOT4, count_AK2, count_AK3, count_AK4, count_AK5, count_ATAD3B, count_AUH, count_BTBD19, count_C21orf91, count_CCDC58, count_CDA, count_CMPK1, count_CMPK2, count_COLEC11, count_COX4I1, count_CRLF3, count_CROT, count_DARS2, count_DCTD, count_DCTPP1, count_DDX28, count_DHX30, count_DNA2, count_DTYMK, count_DUT, count_ENSP00000433415, count_ENTPD1, count_ENTPD3, count_ENTPD8, count_EPRS, count_ERAL1, count_EXO5, count_FAM178A, count_FASTKD5, count_FEN1, count_GABPA, count_GABPB1, count_GAL3ST3, count_GARS, count_GATC, count_GKAP1, count_GPR146, count_GRSF1, count_GSTK1, count_GUK1, count_IMMT, count_INIP, count_INTS3, count_KIAA0020, count_KLHL35, count_LIG1, count_LLPH, count_LONP1, count_LRPPRC, count_LRRC37B, count_MCM5, count_MCM9, count_MFN2, count_MPV17L, count_MPV17L2, count_MRPL12, count_MRPL21, count_MTERF1, count_MTERF2, count_MTERF3, count_MTERF4, count_MTPAP, count_NDUFA5, count_NDUFB8, count_NDUFS3, count_NIFK, count_NIP7, count_NKRF, count_NME4, count_NPTX2, count_NRF1, count_NSUN4, count_NT5C, count_NT5C1A, count_NT5C1B, count_NT5C2, count_NT5C3A, count_NT5C3B, count_NT5M, count_NTPCR, count_NUDT19, count_OXA1L, count_PCNA, count_PDE12, count_PEAK1, count_PEX10, count_PEX12, count_PEX13, count_PEX14, count_PEX2, count_PIF1, count_PNP, count_POLDIP2, count_POLE, count_POLQ, count_PPARGC1A, count_PPRC1, count_PRIM1, count_PRIM2, count_PTCD1, count_PTCD2, count_PTCD3, count_RARS2, count_RBM28, count_REV3L, count_RFWD3, count_RHCG, count_RNASEH2A, count_RNASEH2C, count_RPA1, count_RPA2, count_RPL30, count_RRM2, count_SAMHD1, count_SLC25A4, count_SSBP2, count_SUCLA2, count_SUCLG1, count_SUPV3L1, count_SURF1, count_TACO1, count_TBL3, count_TFB1M, count_TICRR, count_TIMELESS, count_TK1, count_TMEM242, count_TMEM70, count_TOMM20, count_TP53AIP1, count_TRAPPC12, count_TRMT10C, count_TRMU, count_TSSC1, count_TYMP, count_UCK1, count_UCN, count_UTP6, count_ZUFSP)
# need to perform the MT-ATP6, MT-CO2, MT-ND4 individually, also query proteins

write.table(Degree_0.5_100, "Degree_0.5_100.txt", sep = "\t")
write.table(Count_0.5_100, "Count_0.5_100.txt", sep = "\t")

#Transcription factor extention 
Proximal <- file.path(getwd(), "encode-proximal-2012.xgmml")
Distal <- file.path(getwd(), "encode-distal-2012.xgmml")
CTLextend.cmd = paste('cytargetlinker extend idAttribute="display name" linkSetFiles="', Proximal, ',', Distal, '" network=current direction=SOURCES', sep="")
commandsRun(CTLextend.cmd)
layoutNetwork()
renameNetwork('TF_extended_encode')

vizstyle.file <- file.path(getwd(), "CyTargetLinker.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)
setVisualStyle('CyTargetLinker')

#transcription factors
connectionsp <- c("BCL3", "RXRA", "SMARCB1", "CEBPB", "YY1", "TCF12", "SMC3", "E2F6", "EP300", "CTCF", "SIX5", "HNF4G", "SPI1", "RAD21", "ZBTB7A", "JUN", "HEY1", "ZBTB33", "POU2F2", "BATF", "IRF4", "MYC", "SRF", "JUNB", "SP1", "HNF4A", "HDAC2", "TFAP2C", "STAT3", "ZZZ3", "NFKB1", "FOSL2", "JUND", "TBP", "FOXA1", "EBF1", "BCL11A", "PAX5", "MAFK", "MAX", "TAF1", "USF1", "MEF2A", "NANOG", "BHLHE40", "TCF4", "EGR1", "ESRRA", "SMARCC1", "ZNF143", "TFAP2A", "FOXA2", "SMARCC2", "ATF3", "SREBF1", "CHD2", "E2F1", "ELF1", "E2F4", "NFYA", "ELK4", "PBX3", "NFE2", "STAT1", "FOS", "MAFF", "SP2", "SMARCA4", "TAL1", "ESR1", "GATA1", "STAT2", "CTCFL", "GTF2B", "SUZ12", "HSF1", "USF2", "NFYB", "SETDB1", "SIN3A", "REST", "NR3C1", "ZNF263", "BCLAF1", "ZNF274", "POLR3A", "BDP1", "GATA2", "NR2C2")

for (c in connectionsp){
  selectNodes(c, by.col='CTL.name')
  selectFirstNeighbors(direction = "outgoing")
  assign(paste("TF_",c, sep = ""), getSelectedNodes())
  clearSelection()
    
}

#TFP <- rbind(TF_BCL3, TF_RXRA, TF_SMARCB1, TF_CEBPB, TF_YY1, TF_TCF12, TF_SMC3, TF_E2F6, TF_EP300, TF_CTCF, TF_SIX5, TF_HNF4G, TF_SPI1, TF_RAD21, TF_ZBTB7A, TF_JUN, TF_HEY1, TF_ZBTB33, TF_POU2F2, TF_BATF, TF_IRF4, TF_MYC, TF_SRF, TF_JUNB, TF_SP1, TF_HNF4A, TF_HDAC2, TF_TFAP2C, TF_STAT3, TF_ZZZ3, TF_NFKB1, TF_FOSL2, TF_JUND, TF_TBP, TF_FOXA1, TF_EBF1, TF_BCL11A, TF_PAX5, TF_MAFK, TF_MAX, TF_TAF1, TF_USF1, TF_MEF2A, TF_NANOG, TF_BHLHE40, TF_TCF4, TF_EGR1, TF_ESRRA, TF_SMARCC1, TF_ZNF143, TF_TFAP2A, TF_FOXA2, TF_SMARCC2, TF_ATF3, TF_SREBF1, TF_CHD2, TF_E2F1, TF_ELF1, TF_E2F4, TF_NFYA, TF_ELK4, TF_PBX3, TF_NFE2, TF_STAT1, TF_FOS, TF_MAFF, TF_SP2, TF_SMARCA4, TF_TAL1, TF_ESR1, TF_GATA1, TF_STAT2, TF_CTCFL, TF_GTF2B, TF_SUZ12, TF_HSF1, TF_USF2, TF_NFYB, TF_SETDB1, TF_SIN3A, TF_REST, TF_NR3C1, TF_ZNF263, TF_BCLAF1, TF_ZNF274, TF_POLR3A, TF_BDP1, TF_GATA2, TF_NR2C2)

TF_BCL3
TF_RXRA
TF_SMARCB1
TF_CEBPB
TF_YY1
TF_TCF12
TF_SMC3
TF_E2F6
TF_EP300
TF_CTCF
TF_SIX5
TF_HNF4G
TF_SPI1
TF_RAD21
TF_ZBTB7A
TF_JUN
TF_HEY1
TF_ZBTB33
TF_POU2F2
TF_BATF
TF_IRF4
TF_MYC
TF_SRF
TF_JUNB
TF_SP1
TF_HNF4A
TF_HDAC2
TF_TFAP2C
TF_STAT3
TF_ZZZ3
TF_NFKB1
TF_FOSL2
TF_JUND
TF_TBP
TF_FOXA1
TF_EBF1
TF_BCL11A
TF_PAX5
TF_MAFK
TF_MAX
TF_TAF1
TF_USF1
TF_MEF2A
TF_NANOG
TF_BHLHE40
TF_TCF4
TF_EGR1
TF_ESRRA
TF_SMARCC1
TF_ZNF143
TF_TFAP2A
TF_FOXA2
TF_SMARCC2
TF_ATF3
TF_SREBF1
TF_CHD2
TF_E2F1
TF_ELF1
TF_E2F4
TF_NFYA
TF_ELK4
TF_PBX3
TF_NFE2
TF_STAT1
TF_FOS
TF_MAFF
TF_SP2
TF_SMARCA4
TF_TAL1
TF_ESR1
TF_GATA1
TF_STAT2
TF_CTCFL
TF_GTF2B
TF_SUZ12
TF_HSF1
TF_USF2
TF_NFYB
TF_SETDB1
TF_SIN3A
TF_REST
TF_NR3C1
TF_ZNF263
TF_BCLAF1
TF_ZNF274
TF_POLR3A
TF_BDP1
TF_GATA2
TF_NR2C2

#write.table(TFP, "TFP_0.5_150.txt", sep="\t")

#########################################################TRRUST#######################################################
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

#Put script and additional files in the same folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Create network
genes <- "POLG, POLG2, MGME1, SSBP1, C10orf2, TOP1MT, RNASEH1, PRIMPOL, EXOG, ATAD3A, RRM2B, TK2, DGUOK, POLRMT, TFAM, TEFM, TFB2M, MPV17" 
value <- "0.5"
l <- "150"

# load visual style from file
vizstyle.file <- file.path(getwd(), "Dynamics.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)

string.cmd = paste('string protein query query=\"',genes,'\" cutoff=\"',value,'\" species="Homo sapiens" limit="0" newNetName=\"STRING',l,'_',value,'\"', sep='')
commandsRun(string.cmd)
string2.cmd <- paste('string expand additionalNodes=\"',l,'\" network=\"CURRENT\"')
commandsRun(string2.cmd)

deleteTableColumn(column = 'stringdb::structures', table='node')
deleteTableColumn(column = 'stringdb::sequence', table='node')
deleteTableColumn(column = 'stringdb::enhancedLabel Passthrough', table='node')
deleteTableColumn(column = 'stringdb::STRING style', table='node')

# apply visual style ()
setVisualStyle('Dynamics')
#renameNetwork('STRING_MitoCarta')

#Select input nodes
clearSelection()
selectNodes(c('POLG', 'POLG2', 'MGME1', 'SSBP1', 'C10orf2', 'TOP1MT', 'RNASEH1', 'PRIMPOL', 'EXOG', 'ATAD3A', 'RRM2B', 'TK2', 'DGUOK', 'POLRMT', 'TFAM', 'TEFM', 'TFB2M', 'MPV17'), by.col='query term')

Input <- getSelectedNodes()
setNodeShapeBypass(Input,'diamond')

# heat propagation
# https://mkutmon.gitlab.io/tutorial-network-algorithms/NetworkPropagation.html
diffusionBasic()
clearSelection()
setNodeColorMapping("diffusion_output_heat", c(0,0.5), c('#FFFFFF', '#FF0000'), style.name = "Dynamics")
assign(paste("data_",l,'_',value, sep = ""), getTableColumns('node', c('display name', 'diffusion_output_heat', 'diffusion_output_rank')))

#Select one added protein
connections <- c("ACOT4", "AK2", "AK3", "AK4", "AK5", "ATAD3B", "AUH", "BTBD19", "C21orf91", "CCDC58", "CDA", "CMPK1", "CMPK2", "COLEC11", "COX4I1", "CRLF3", "CROT", "DARS2", "DCTD", "DCTPP1", "DDX28", "DHX30", "DNA2", "DTYMK", "DUT", "ENSP00000433415", "ENTPD1", "ENTPD3", "ENTPD8", "EPRS", "ERAL1", "EXO5", "FAM178A", "FASTKD5", "FEN1", "GABPA", "GABPB1", "GAL3ST3", "GARS", "GATC", "GKAP1", "GPR146", "GRSF1", "GSTK1", "GUK1", "IMMT", "INIP", "INTS3", "KIAA0020", "KLHL35", "LIG1", "LLPH", "LONP1", "LRPPRC", "LRRC37B", "MCM5", "MCM9", "MFN2", "MPV17L", "MPV17L2", "MRPL12", "MRPL21", "MT-ATP6", "MT-CO1", "MT-CO2", "MT-CYB", "MTERF1", "MTERF2", "MTERF3", "MTERF4", "MT-ND4", "MT-ND6", "MTPAP", "NDUFA5", "NDUFB8", "NDUFS3", "NIFK", "NIP7", "NKRF", "NME4", "NPTX2", "NRF1", "NSUN4", "NT5C", "NT5C1A", "NT5C1B", "NT5C2", "NT5C3A", "NT5C3B", "NT5M", "NTPCR", "NUDT19", "OXA1L", "PCNA", "PDE12", "PEAK1", "PEX10", "PEX12", "PEX13", "PEX14", "PEX2", "PIF1", "PNP", "POLDIP2", "POLE", "POLQ", "PPARGC1A", "PPRC1", "PRIM1", "PRIM2", "PTCD1", "PTCD2", "PTCD3", "RARS2", "RBM28", "REV3L", "RFWD3", "RHCG", "RNASEH2A", "RNASEH2C", "RPA1", "RPA2", "RPL30", "RRM2", "SAMHD1", "SLC25A4", "SSBP2", "SUCLA2", "SUCLG1", "SUPV3L1", "SURF1", "TACO1", "TBL3", "TFB1M", "TICRR", "TIMELESS", "TK1", "TMEM242", "TMEM70", "TOMM20", "TP53AIP1", "TRAPPC12", "TRMT10C", "TRMU", "TSSC1", "TYMP", "UCK1", "UCN", "UTP6", "ZUFSP")

for (c in connections){
  selectNodes(c, by.col='display name')
  selectFirstNeighbors(direction = "any")
  assign(paste("degree_",c, sep = ""), getSelectedNodeCount())
  n <- 0
  count <- getSelectedNodes()
  
  
  if ("9606.ENSP00000268124" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000442563" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000366939" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000419665" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000309595" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000328835" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000313350" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000313816" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000287675" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000368030" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000251810" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000299697" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000264093" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000369383" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000465759" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000420588" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000462963" %in% count) {
    n <- n+1
  }
  if ("9606.ENSP00000355471" %in% count) {
    n <- n+1
  }
  
  assign(paste("count_",c, sep = ""), n)
  
  clearSelection()
} 


Degree_0.5_150 <- rbind(degree_ACOT4, degree_AK2, degree_AK3, degree_AK4, degree_AK5, degree_ATAD3B, degree_AUH, degree_BTBD19, degree_C21orf91, degree_CCDC58, degree_CDA, degree_CMPK1, degree_CMPK2, degree_COLEC11, degree_COX4I1, degree_CRLF3, degree_CROT, degree_DARS2, degree_DCTD, degree_DCTPP1, degree_DDX28, degree_DHX30, degree_DNA2, degree_DTYMK, degree_DUT, degree_ENSP00000433415, degree_ENTPD1, degree_ENTPD3, degree_ENTPD8, degree_EPRS, degree_ERAL1, degree_EXO5, degree_FAM178A, degree_FASTKD5, degree_FEN1, degree_GABPA, degree_GABPB1, degree_GAL3ST3, degree_GARS, degree_GATC, degree_GKAP1, degree_GPR146, degree_GRSF1, degree_GSTK1, degree_GUK1, degree_IMMT, degree_INIP, degree_INTS3, degree_KIAA0020, degree_KLHL35, degree_LIG1, degree_LLPH, degree_LONP1, degree_LRPPRC, degree_LRRC37B, degree_MCM5, degree_MCM9, degree_MFN2, degree_MPV17L, degree_MPV17L2, degree_MRPL12, degree_MRPL21, degree_MTERF1, degree_MTERF2, degree_MTERF3, degree_MTERF4, degree_MTPAP, degree_NDUFA5, degree_NDUFB8, degree_NDUFS3, degree_NIFK, degree_NIP7, degree_NKRF, degree_NME4, degree_NPTX2, degree_NRF1, degree_NSUN4, degree_NT5C, degree_NT5C1A, degree_NT5C1B, degree_NT5C2, degree_NT5C3A, degree_NT5C3B, degree_NT5M, degree_NTPCR, degree_NUDT19, degree_OXA1L, degree_PCNA, degree_PDE12, degree_PEAK1, degree_PEX10, degree_PEX12, degree_PEX13, degree_PEX14, degree_PEX2, degree_PIF1, degree_PNP, degree_POLDIP2, degree_POLE, degree_POLQ, degree_PPARGC1A, degree_PPRC1, degree_PRIM1, degree_PRIM2, degree_PTCD1, degree_PTCD2, degree_PTCD3, degree_RARS2, degree_RBM28, degree_REV3L, degree_RFWD3, degree_RHCG, degree_RNASEH2A, degree_RNASEH2C, degree_RPA1, degree_RPA2, degree_RPL30, degree_RRM2, degree_SAMHD1, degree_SLC25A4, degree_SSBP2, degree_SUCLA2, degree_SUCLG1, degree_SUPV3L1, degree_SURF1, degree_TACO1, degree_TBL3, degree_TFB1M, degree_TICRR, degree_TIMELESS, degree_TK1, degree_TMEM242, degree_TMEM70, degree_TOMM20, degree_TP53AIP1, degree_TRAPPC12, degree_TRMT10C, degree_TRMU, degree_TSSC1, degree_TYMP, degree_UCK1, degree_UCN, degree_UTP6, degree_ZUFSP)
Count_0.5_150 <- rbind(count_ACOT4, count_AK2, count_AK3, count_AK4, count_AK5, count_ATAD3B, count_AUH, count_BTBD19, count_C21orf91, count_CCDC58, count_CDA, count_CMPK1, count_CMPK2, count_COLEC11, count_COX4I1, count_CRLF3, count_CROT, count_DARS2, count_DCTD, count_DCTPP1, count_DDX28, count_DHX30, count_DNA2, count_DTYMK, count_DUT, count_ENSP00000433415, count_ENTPD1, count_ENTPD3, count_ENTPD8, count_EPRS, count_ERAL1, count_EXO5, count_FAM178A, count_FASTKD5, count_FEN1, count_GABPA, count_GABPB1, count_GAL3ST3, count_GARS, count_GATC, count_GKAP1, count_GPR146, count_GRSF1, count_GSTK1, count_GUK1, count_IMMT, count_INIP, count_INTS3, count_KIAA0020, count_KLHL35, count_LIG1, count_LLPH, count_LONP1, count_LRPPRC, count_LRRC37B, count_MCM5, count_MCM9, count_MFN2, count_MPV17L, count_MPV17L2, count_MRPL12, count_MRPL21, count_MTERF1, count_MTERF2, count_MTERF3, count_MTERF4, count_MTPAP, count_NDUFA5, count_NDUFB8, count_NDUFS3, count_NIFK, count_NIP7, count_NKRF, count_NME4, count_NPTX2, count_NRF1, count_NSUN4, count_NT5C, count_NT5C1A, count_NT5C1B, count_NT5C2, count_NT5C3A, count_NT5C3B, count_NT5M, count_NTPCR, count_NUDT19, count_OXA1L, count_PCNA, count_PDE12, count_PEAK1, count_PEX10, count_PEX12, count_PEX13, count_PEX14, count_PEX2, count_PIF1, count_PNP, count_POLDIP2, count_POLE, count_POLQ, count_PPARGC1A, count_PPRC1, count_PRIM1, count_PRIM2, count_PTCD1, count_PTCD2, count_PTCD3, count_RARS2, count_RBM28, count_REV3L, count_RFWD3, count_RHCG, count_RNASEH2A, count_RNASEH2C, count_RPA1, count_RPA2, count_RPL30, count_RRM2, count_SAMHD1, count_SLC25A4, count_SSBP2, count_SUCLA2, count_SUCLG1, count_SUPV3L1, count_SURF1, count_TACO1, count_TBL3, count_TFB1M, count_TICRR, count_TIMELESS, count_TK1, count_TMEM242, count_TMEM70, count_TOMM20, count_TP53AIP1, count_TRAPPC12, count_TRMT10C, count_TRMU, count_TSSC1, count_TYMP, count_UCK1, count_UCN, count_UTP6, count_ZUFSP)
# need to perform the MT-ATP6, MT-CO2, MT-ND4 individually, also query proteins

write.table(Degree_0.5_150, "Degree_0.5_150.txt", sep = "\t")
write.table(Count_0.5_150, "Count_0.5_150.txt", sep = "\t")

#Transcription factor extention 
Repression <- file.path(getwd(), "TRRUST-hsa-repression.xgmml")
Activation <- file.path(getwd(), "TRRUST-hsa-activation.xgmml")
Unknown <- file.path(getwd(), "TRRUST-hsa-unknown.xgmml")
CTLextend.cmd = paste('cytargetlinker extend idAttribute="display name" linkSetFiles="', Repression, ',', Activation, ',', Unknown, '" network=current direction=SOURCES', sep="")
commandsRun(CTLextend.cmd)
layoutNetwork()
renameNetwork('TF_extended_TRRUST')

vizstyle.file <- file.path(getwd(), "CyTargetLinker.xml")
LoadStyle.cmd = paste('vizmap load file file="',vizstyle.file,'"', sep="")
commandsRun(LoadStyle.cmd)
setVisualStyle('CyTargetLinker')

#transcription factors
connectionsp <- c("AHR", "AR", "ATF1", "ATF3", "ATM", "CDX1", "CREB1", "CRTC1", "E2F1", "EP300", "ESR1", "EZH2", "GATA1", "HDAC1", "HMGN1", "HOXA1", "ING1", "MEF2A", "MITF", "MYB", "MYC", "MYCN", "NFKB1", "POU2AF1", "POU2F1", "POU2F2", "PRDM1", "PURA", "RB1", "RELA", "RFX1", "SP1", "SP3", "TFDP1", "TP53", "YY1")

for (c in connectionsp){
  selectNodes(c, by.col='CTL.label')
  selectFirstNeighbors(direction = "outgoing")
  assign(paste("TF_",c, sep = ""), getSelectedNodes())
  clearSelection()
  
}

#TFP <- rbind(TF_EP300, TF_AR, TF_PURA, TF_RB1, TF_RFX1, TF_HDAC1, TF_CRTC1, TF_AHR, TF_SP1, TF_MEF2A, TF_E2F1, TF_MITF, TF_MYB, TF_YY1, TF_ESR1, TF_MYC, TF_MYCN, TF_EZH2, TF_TFDP1, TF_GATA1, TF_CDX1, TF_POU2AF1, TF_POU2F2, TF_POU2F1, TF_CREB1, TF_PRDM1, TF_TP53, TF_HMGN1, TF_HOXA1, TF_ING1, TF_RELA, TF_NFKB1, TF_SP3, TF_ATF1, TF_ATF3, TF_ATM)

TF_EP300
TF_AR
TF_PURA
TF_RB1
TF_RFX1
TF_HDAC1
TF_CRTC1
TF_AHR
TF_SP1
TF_MEF2A
TF_E2F1
TF_MITF
TF_MYB
TF_YY1
TF_ESR1
TF_MYC
TF_MYCN
TF_EZH2
TF_TFDP1
TF_GATA1
TF_CDX1
TF_POU2AF1
TF_POU2F2
TF_POU2F1
TF_CREB1
TF_PRDM1
TF_TP53
TF_HMGN1
TF_HOXA1
TF_ING1
TF_RELA
TF_NFKB1
TF_SP3
TF_ATF1
TF_ATF3
TF_ATM

#write.table(TFP, "TFP_0.5_150.txt", sep="\t")