library(tidyverse)
library(ggplot2)
library(MetBrewer)
library(patchwork)
library(lme4)
library(lsmeans)
library(ggResidpanel)


################################################################
#set personal theme
################################################################
mytheme <- theme(axis.text = element_text(size = 40))+
  theme (plot.title = element_text (size = 30, face = "bold" ))+
  theme(legend.text=element_text(size=18))+
  theme(legend.title=element_text(size=20))+
  theme(plot.subtitle=element_text(size=20))+
  theme(axis.text.x=element_text(size=20))+
  theme(axis.text.y=element_text(size=20))+
  theme(axis.title.x= element_text(size = 22))+
  theme(axis.title.y = element_text(size = 22))

set.seed(1993)

################################################################
#load data
################################################################

df<-read.delim("feature_1794.tsv")
group_349 <- read.delim("group.tsv")
group_349$people <- as.character(group_349$people)
group_349 <- group_349 %>%
  mutate(disease= case_when(str_detect(group,"case") ~ "GDM",
                            str_detect(group,"ctr") ~ "Control"))


df <- df %>% 
  left_join(group_349,by="sample") 


################################################################
#plot MA related
################################################################

MA <- ggplot(data= df,
               aes(x =  time,
                   y =  MA))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#007e2f","#ffcd12"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("MA")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: <0.01", x = 2.5, y = 6.2), size = 5, col = "black")+
  mytheme



all <- ggplot(data= df,
       aes(x =  time,
           y =  all))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#007e2f","#ffcd12"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("all")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: 0.03", x = 2.5, y = 0.97), size = 5, col = "black")+
  mytheme




short <- ggplot(data= df,
              aes(x =  time,
                  y =  short))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#007e2f","#ffcd12"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("short")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: 0.03", x = 2.5, y = 0.97), size = 5, col = "black")+
  mytheme



peak <- ggplot(data= df,
              aes(x =  time,
                  y =  peak))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#007e2f","#ffcd12"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("peak")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: 0.08", x = 2.5, y = 0.97), size = 5, col = "black")+
  mytheme







long <- ggplot(data= df,
               aes(x =  time,
                   y =  long))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#007e2f","#ffcd12"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("long")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: 0.11", x = 2.5, y = 0.98), size = 5, col = "black")+
  mytheme+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = "bottom")



MA+all+short+peak+long
ggsave("figure2e.pdf",width = 25, height = 18,units = "cm")



################################################################
#plot motif related
################################################################
picccca <- ggplot(data= df,aes(x =  time,
                               y =  CCCA))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey",geom='errorbar')+
  scale_fill_manual(values = c("#E3242B","#17486f"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("CCCA")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: 0.04", x = 2.5, y = 0.03), size = 5, col = "black")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = "none")+
  mytheme


picactt <- ggplot(data= df,aes(x =  time,
                     y =  ACTT))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#E3242B","#17486f"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("ACTT")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: <0.01", x = 2.5, y = 0.0085), size = 5, col = "black")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = "none")+
  mytheme


picaccg <- ggplot(data= df,aes(x =  time,
                               y =  ACCG))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#E3242B","#17486f"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("ACCG")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: <0.01", x = 2.5, y = 0.00125), size = 5, col = "black")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = "none")+
  mytheme



picgcgg <- ggplot(data= df,aes(x =  time,
                               y =  GCGG))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#E3242B","#17486f"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("GCGG")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: <0.01", x = 2.5, y = 0.0007), size = 5, col = "black")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = "none")+
  mytheme



picgcgc <- ggplot(data= df,aes(x =  time,
                               y =  GCGC))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey")+
  scale_fill_manual(values = c("#E3242B","#17486f"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("GCGC")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: <0.01", x = 2.5, y = 0.0006), size = 5, col = "black")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = "bottom")+
  mytheme



pictcgg <- ggplot(data= df,aes(x =  time,
                               y =  TCGG))+
  stat_boxplot(aes(fill = disease),geom ='errorbar',width = 0.4,position = position_dodge(width = 0.75)) + 
  geom_boxplot(aes(fill = disease),outlier.size = 1,outlier.color = "grey",geom='errorbar')+
  scale_fill_manual(values = c("#E3242B","#17486f"))+ 
  theme_classic()+
  xlab("") +
  ylab("") +
  ggtitle("TCGG")+
  guides(fill="none")+
  scale_x_discrete(labels = c("first" = "1st","second" = "2nd","third" = "3rd"))+
  theme (plot.title = element_text (size = 16, face = "bold",hjust = 0.5))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))+
  geom_text(aes(label = "P value: <0.01", x = 2.5, y = 0.0006), size = 5, col = "black")+
  guides(fill=guide_legend(title=NULL))+
  theme(legend.position = "none")+
  mytheme
  
picccca+picactt+picaccg+picgcgg+picgcgc+pictcgg

ggsave("figure2d.pdf",width = 25, height = 18,units = "cm")


