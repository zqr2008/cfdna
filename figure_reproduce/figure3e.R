library(tidyverse)
library(ggplot2)
library(jtools)
library(forestmangr)
library(wesanderson)
library(ggsci)
library(hrbrthemes)
library(ggthemes)
library(patchwork)
library(readxl)
library(sjmisc)
library(ggpubr)
#######################################################
#load data
#######################################################
mytheme <- theme(axis.text = element_text(size = 40))+
  theme (plot.title = element_text (size = 30, face = "bold" ))+
  theme(legend.text=element_text(size=18))+
  theme(legend.title=element_text(size=20))+
  theme(plot.subtitle=element_text(size=20))+
  theme(axis.text.x=element_text(size=18))+
  theme(axis.text.y=element_text(size=16))+
  theme(axis.title.x= element_text(size = 22))+
  theme(axis.title.y = element_text(size = 22)) +
  theme(text=element_text(family="sans"))


gene50_tss <- read.delim("~/cfDNA/gene50_tss.txt", row.names=1)
data299new <- read_excel("~/cfDNA/data299new.xlsx")
group_349 <- read.csv("~/cfDNA/group_349.csv")
group_349$people <- as.character(group_349$people)



#######################################################
#set theme
#######################################################
mytheme <- theme(axis.text = element_text(size = 40))+
  theme (plot.title = element_text (size = 30, face = "bold" ))+
  theme(legend.text=element_text(size=18))+
  theme(legend.title=element_text(size=20))+
  theme(plot.subtitle=element_text(size=20))+
  theme(axis.text.x=element_text(size=16))+
  theme(axis.text.y=element_text(size=28))+
  theme(axis.title.x= element_text(size = 22))+
  theme(axis.title.y = element_text(size = 22))


heatmap <- gene50_tss %>%
  rotate_df() %>%
  mutate(across(everything(), as.numeric))

xx<-str_split_fixed(colnames(heatmap),'_',6)
xx<-xx[,4]
colnames(heatmap)<-xx




#######################################################
#join the clinical information
#######################################################
heatmap <- heatmap %>%
  rownames_to_column("sample") %>%
  left_join(group_349,by = "sample") %>%
  dplyr::rename(ID = people) %>%
  left_join(data299new,by= "ID") 


#######################################################
#select growth related genes
#######################################################
preterm.gene <- list("CBS","CDH23","CNTN4","RNF213","SMARCD1","TWSG1")
#preterm.gene <- list("CBS","CDH23")



heatmap_odds <- heatmap %>%
  plyr::mutate(pre = case_when(Gestationalage<=37  ~"preterm",
                                   Gestationalage > 37~"normal"))

heatmap_odds$pre <- factor(heatmap_odds$pre)



pic <- list()

for (i in preterm.gene){
  
mylogit <- glm(pre ~ get(i), data = heatmap_odds,family = "binomial")
modelsum <- summary(mylogit)
co <- as.data.frame(modelsum$coefficients)
pvalue <- co[2,4]


out <- as.data.frame(exp(cbind(Odds_Ratio = coef(mylogit), confint(mylogit,level = 0.95)))) 
out <- out %>%
  dplyr::rename(conf.low =`2.5 %`,
                conf.high = `97.5 %`) %>%
  rownames_to_column("term") 

out <- out[-1,]
out[1,1] <- i
out[1,5] <- pvalue
out <- out %>%
  dplyr::rename(pp = V5)

pic[[i]] <- ggplot(out, aes(x=reorder(term, Odds_Ratio), y=Odds_Ratio)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), 
                width = 0.,size  = 0.6,
                position = "dodge", color= "#508ea2") +
  geom_hline(yintercept = 1, color = "#802417", size = 1,linetype = "dotted") +
  geom_point(size = 6,shape =15,color = "#508ea2") + 
  coord_flip() +
  theme_minimal() +
  xlab("")+
  ylab("") +
  #scale_y_continuous(limits = c(0,2))+
  geom_text(aes(x = 1.5,y = conf.high*0.6,label = paste0(round(Odds_Ratio,2),"[",round(conf.low,2),"-",round(conf.high,2),"]")),size = 5.0,font = "sans")+
  geom_text(aes(x = 1.3,y = conf.high*0.6,label = paste0("P: ",round(pp,2))),size = 5,font = "sans")+
  mytheme
}


pic[["SMARCD1"]] <- pic[["SMARCD1"]]+
  ylab("Odds ratio")
  
  
pic <-(pic[["CBS"]]+pic[["CDH23"]])/(pic[["CNTN4"]]+pic[["RNF213"]])/(pic[["SMARCD1"]]+pic[["TWSG1"]])





