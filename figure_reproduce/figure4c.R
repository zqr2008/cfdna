library(tidyverse)
library(ggplot2)
library(sjmisc)
library(pheatmap)
library(reshape2)
library(vegan)
library(ade4)
library(WGCNA)
library(rstatix)
library(patchwork)
library(ggrepel)

Labeled_Lip_Un_Tar_total <- read.csv("~/cfDNA/Labeled_Lip_Un_Tar_total.csv")
Iden_QuantNEG <- read.delim("~/cfDNA/Iden_QuantNEG.txt")
Iden_QuantPOS <- read.delim("~/cfDNA/Iden_QuantPOS.txt")
gene50_tss <- read.delim("~/cfDNA/gene50_tss.txt", row.names=1)
data299new <- read_excel("data299new.xlsx")
group_349 <- read.csv("group_349.csv")
group_349$people <- as.character(group_349$people)

Iden_QuantNEG <- Iden_QuantNEG %>%
  dplyr::select(1,3)

Iden_QuantPOS <- Iden_QuantPOS%>%
  dplyr::select(1,3)


Iden_Quant <- rbind(Iden_QuantNEG,Iden_QuantPOS)
  
Labeled_Lip_Un_Tar_total  <- Labeled_Lip_Un_Tar_total  %>%
  column_to_rownames("ID")


labelled <- as.data.frame(t(Labeled_Lip_Un_Tar_total))

labelled <- labelled  %>%
  rownames_to_column("Metabolite.ID")

labelled$Metabolite.ID <- str_split_fixed(labelled$Metabolite.ID,'X',3)[,2]

labelled <- labelled %>%
  left_join(Iden_Quant,by = "Metabolite.ID") %>%
  drop_na("LipidIon") %>%
  #group_by(LipidIon) %>% 
  #filter(n()>1)
  distinct(LipidIon,.keep_all = T) %>%
  column_to_rownames("LipidIon") 
  
labelled <- as.data.frame(t(labelled))

heatmap <- gene50_tss %>%
  rotate_df() %>%
  mutate(across(everything(), as.numeric))

heatmap <- as.data.frame(heatmap)

xx<-str_split_fixed(colnames(heatmap),'_',6)
xx<-xx[,4]
colnames(heatmap)<-xx

#join the clinical information
heatmap <- heatmap %>%
  rownames_to_column("sample") %>%
  left_join(group_349,by = "sample") %>%
  rename(ID = people) %>%
  left_join(data299new,by= "ID") 

  
labelled <- labelled %>%
  rownames_to_column("sam") %>%
  mutate(sample = paste0("D",sam)) %>%
  dplyr::inner_join(heatmap,by="sample") 

lipid <- labelled %>%
  dplyr::select(1:PDXP,disease,group) %>%
  mutate_at(.vars = vars(1:sample), .fun = as.numeric) %>%
  dplyr::select(-c(sam))

lipid2 <- labelled %>%
  dplyr::select(1:PDXP,disease,group,time,ID,GDMgroup) %>%
  mutate_at(.vars = vars(1:sample), .fun = as.numeric) %>%
  dplyr::select(-c(sam))

logfolddata <- lipid %>%
  dplyr::select(1:sample,disease,group) %>%
  #group_by(group) %>%
  group_by(disease) %>%
  dplyr::mutate(across(where(is.numeric), sum, .names = "sum_{col}")) %>%
  dplyr:: select(starts_with("sum")) %>%
  #distinct(group,.keep_all = T) %>%
  distinct(disease,.keep_all = T) %>%
  column_to_rownames("disease") %>%
  #column_to_rownames("group") %>%
  rotate_df() %>%
  mutate(log = log2(case/ctr)) %>%
  rownames_to_column("loglipid") %>%
  mutate(lipidname = str_split_fixed(loglipid,'_',2)[,2]) 
  #mutate(log2_early = log2(case1/ctr1),
  #      log2_middle = log2(case2/ctr2),
  #      log2_late = log2(case3/ctr3)) 
 
  
filterlog <-  logfolddata %>%
    filter( abs(log) >0.4)
  # filter( abs(log2_middle) >0.6)
  

name_log <- filterlog$lipidname





lipid <- lipid %>%
  dplyr::select(-sample)




all_p_val <- list()
pv <- list()
p_values <- list()

for (col in names(lipid)) {
  if (typeof(lipid[,col]) == "double" ){
  p <- wilcox.test(get(col)~ as.factor(disease),alternative = "two.sided",
                   data = lipid)$p.value
  pv <- append(pv,p)
  all_p_val <- append(all_p_val,col)
  #print(all_p_val)

  if (p < 0.01){
  p_values <- append(p_values,col)
  }
}
}



pvalue.df <- cbind(as.data.frame(unlist(all_p_val)),as.data.frame(unlist(pv))) %>%
  dplyr::rename(lipidname =`unlist(all_p_val)`,
         p.value = `unlist(pv)`)


selected <- lipid2 %>%
  dplyr::select(any_of(name_log),RCAN3:PDXP,time,ID,group,GDMgroup) %>%
  dplyr::select(any_of(unlist(p_values)),RCAN3:PDXP,time,ID,group,GDMgroup) %>%
  drop_na() 




##################################################
#draw formal heatmap
##################################################
aqm <- melt(selected, id=c("ID","time","GDMgroup"), na.rm=TRUE) %>%
  mutate(newname = paste0(time,variable))%>%
  dplyr::select(ID,GDMgroup,newname,value)


GDM <- reshape(aqm, idvar = c("ID","GDMgroup"), timevar = "newname", direction = "wide") %>%
  drop_na() %>%
  mutate_at(.vars = vars(3:167),.funs = as.numeric) %>%
  #mutate_at(.vars = vars(3:173),.funs = as.numeric) %>%
  filter(GDMgroup == "GDM")

control <- reshape(aqm, idvar = c("ID","GDMgroup"), timevar = "newname", direction = "wide") %>%
  drop_na() %>%
  mutate_at(.vars = vars(3:167),.funs = as.numeric) %>%
  #mutate_at(.vars = vars(3:173),.funs = as.numeric) %>%
  filter(GDMgroup == "Control")

GDMcor <- WGCNA::cor(GDM[,3:167])
controlcor <- WGCNA::cor(control[,3:167])


GDMp <- WGCNA::corPvalueStudent(GDMcor,nSamples = 77)
controlp <- WGCNA::corPvalueStudent(controlcor,nSamples = 75)






################################################
#transfer to function
################################################



built_cor <- function(x,spec = "GDM"){
  if (spec == "GDM") {
    namewhole <- as.data.frame(GDMcor) %>%
      dplyr::select(ends_with(x)) %>%
      rownames_to_column("lipidion")
  }
  
  if (spec == "control") {
    namewhole <- as.data.frame(controlcor) %>%
      dplyr::select(ends_with(x)) %>%
      rownames_to_column("lipidion")
  }
  
  if (spec == "GDMp"){
    namewhole <- as.data.frame(GDMp) %>%
      dplyr::select(ends_with(x)) %>%
      rownames_to_column("lipidion")
  }
  
  if (spec == "controlp"){
    namewhole <- as.data.frame(controlp) %>%
      dplyr::select(ends_with(x)) %>%
      rownames_to_column("lipidion")
  }
  
  namefirst <- namewhole %>%
    filter(str_detect(lipidion,"first")) %>%
    dplyr::select(contains("first") |  lipidion)
  
  namesecond <- namewhole %>%
    filter(str_detect(lipidion,"second")) %>%
    dplyr::select(contains("second") |  lipidion)
  
  namethird <- namewhole %>%
    filter(str_detect(lipidion,"third")) %>%
    dplyr::select(contains("third") |  lipidion)
  
  nameend <- cbind(namefirst,namesecond)
  nameend <- cbind(nameend,namethird)
  
  nameend <- nameend[,-c(2,4)] 
  nameend$lipidion <- substr(nameend$lipidion,start = 12,stop = 50)
  nameend <- nameend %>%
    column_to_rownames("lipidion")
  print(nameend)
}


plotlipid <- function(x){
  gdmpart <- built_cor(x,spec = "GDM")
  controlpart <- built_cor(x,spec = "control")
  gdmp <- built_cor(x,spec = "GDMp")
  controlp <- built_cor(x,spec = "controlp")
  
  mapx <- t(cbind(gdmpart[1:5,],controlpart[1:5,]))
  mapp <- t(cbind(gdmp[1:5,],controlp[1:5,]))
  
  textMatrix <- paste(signif(as.matrix(mapx), 2), "\n(",
                      signif(as.matrix(mapp), 1), ")", sep = "")
  
  k <- gsub(pattern = '(rep)', x = colnames(mapx), replacement = "")
  k <- gsub("(", "", k, fixed = TRUE) 
  k <- gsub(")", "", k, fixed = TRUE) 
  
  q <- gsub(pattern = x, x = row.names(mapx), replacement = "")
  q <- gsub(pattern = 'value.first', x = q, replacement = "1st")
  q <- gsub(pattern = 'value.second', x = q, replacement = "2nd")
  q <- gsub(pattern = 'value.third', x = q, replacement = "3rd")
  
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 4, 0, 0)) 
  pic <-labeledHeatmap(mapx,
                       xLabels = k,
                       yLabels = q,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = 1,
                       zlim = c(-0.4,0.4),
                       ySymbols= row.names(mapx),
                       colorLabels= FALSE,
                       main = paste0(x," relationships with lipidomics"),
                       cex.lab = 0.8)
  return(pic)
  
}

################################################
#draw lipid-related tss score
################################################
pdf("figure4c-1.pdf",width = 8, height = 6)
plotlipid("PRSS1")
dev.off()

pdf("figure4c-2.pdf",width = 8, height = 6)
plotlipid("SPATS2L")
dev.off()

