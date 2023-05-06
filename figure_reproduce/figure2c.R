library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(wesanderson)
library(ggsci)
library(lme4)
library(ggResidpanel)
library(lsmeans)
library(readxl)
library(MetBrewer)
library(rstatix)
library(sjmisc)

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
seqff_1794 <- read.delim("~/cfDNA/seqff_1794.txt")
data299new <- read_excel("data299new.xlsx")
group_349 <- read.csv("group_349.csv")

group_349$people <- as.character(group_349$people)
group_349 <- group_349 %>%
  dplyr::select(people,sample)

seqff_1794 <- seqff_1794 %>%
  left_join(group_349,by="sample") %>%
  dplyr:: rename(ID = people) %>%
  left_join(data299new, by ="ID") %>%
  mutate(trimester = case_when(str_detect(group,"1")~"1st",
                               str_detect(group,"2")~"2nd",
                               str_detect(group,"3")~"3rd"),
         disease = case_when(str_detect(group,"ctr")~"Control",
                           str_detect(group,"case")~"GDM"),
         subject =  substr(number,1,7))


################################################################
#plot
################################################################
seqff_1794 %>% 
  group_by(trimester,disease) %>%
  summarise(seqff_mean = mean(seqff),
            seqff_sd = sd(seqff)) %>%
  ggplot(aes(x = trimester, y = seqff_mean,fill = disease,color = disease))+
  geom_errorbar(mapping = aes (ymax = seqff_mean + seqff_sd,
                               ymin = seqff_mean - seqff_sd),
                width = 0.1,position = position_dodge(width=0.2)) +
  geom_point(mapping = aes(x = trimester, y = seqff_mean,fill = disease),
             position = position_dodge(width=0.2))+
  geom_line(aes(x = trimester, y = seqff_mean, group = disease),
            position = position_dodge(width=0.2))+
  theme_classic() +
  scale_color_manual(values = met.brewer("Austria"))+
  xlab("trimester") +
  ylab("fetal fraction") +
  geom_text(aes(label = "disease difference: p < 0.01", x = 2.5, y = 0.3), size = 5, col = "black")+
  mytheme 

ggsave("figure2c.pdf",width = 15, height = 10,units = "cm")


