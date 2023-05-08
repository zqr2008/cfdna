library(tidyverse)
library(readr)
library(readxl)
library(RColorBrewer)
library(hrbrthemes)
library(MetBrewer)
library(cluster)    
library(factoextra) 
library(Gmisc)
library(grid)
library(easyGgplot2)


trans <- read.delim("~/cfDNA/trans.txt")
trans$people <- as.character(trans$people)
trans <- trans %>%
  dplyr::rename(ID = people) %>%
  left_join(data299new,by = "ID")%>%
  plyr::mutate(newpreterm = case_when(Gestationalage<=37 & GDMgroup =="GDM" ~"preterm with disease",
                                   Gestationalage > 37 & GDMgroup =="Control" ~"normal without disease",
                                   Gestationalage<=37 & GDMgroup =="Control"~"preterm without disease",
                                   Gestationalage >37 & GDMgroup =="GDM" ~"normal with disease"))

trans <- trans %>%
  mutate_at(.vars = vars(3:8),.funs = as.factor)

time1 <- trans %>% 
  dplyr::group_by(first_cluster) %>%
  dplyr::count(newpreterm) %>%
  dplyr::mutate(all = sum(n) ) %>%
  dplyr::mutate(percentage = n / all)  

pic1 <- ggplot2.barplot(data = time1, 
                xName="first_cluster",
                yName='percentage',
                groupName="newpreterm",
                groupColors=c("#899DA4","#C93312","#FAEFD1","#DC863B"), showLegend=TRUE,
                backgroundColor="white", color="black",
                xtitle="1st trimester",
                ytitle="Percentages of preterm birth", 
                mainTitle="",
                removePanelGrid=TRUE,removePanelBorder=TRUE,
                axisLine=c(0.5, "solid", "black"))+
  theme(legend.position = "none")

time2 <- trans %>% 
  group_by(second_cluster) %>%
  dplyr::count(newpreterm) %>%
  dplyr::mutate(all = sum(n) ) %>%
  dplyr::mutate(percentage = n / all)  

pic2 <- ggplot2.barplot(data = time2, 
                xName="second_cluster",
                yName='percentage',
                groupName="newpreterm",
                groupColors=c("#899DA4","#C93312","#FAEFD1","#DC863B"), showLegend=TRUE,
                backgroundColor="white", color="black",
                xtitle="2nd trimester",
                ytitle="", 
                mainTitle="",
                removePanelGrid=TRUE,removePanelBorder=TRUE,
                axisLine=c(0.5, "solid", "black"))+
  theme(legend.position = "none")


time3 <- trans %>% 
  group_by(third_cluster) %>%
  dplyr::count(newpreterm) %>%
  dplyr::mutate(all = sum(n) ) %>%
  dplyr::mutate(percentage = n / all)  

pic3 <- ggplot2.barplot(data = time3, 
                xName="third_cluster",
                yName='percentage',
                groupName="newpreterm",
                groupColors=c("#899DA4","#C93312","#FAEFD1","#DC863B"), showLegend=T,
                backgroundColor="white", color="black",
                xtitle="3rd trimester",
                ytitle="", 
                mainTitle="",
                removePanelGrid=TRUE,removePanelBorder=TRUE,
                axisLine=c(0.5, "solid", "black"))


pic1+pic2+pic3
ggsave("figure3c.pdf",width = 25, height = 10,units = "cm")
