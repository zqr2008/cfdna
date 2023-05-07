library(tidyverse)
library(pheatmap)
library(sjmisc)
library(readxl)
library(RColorBrewer)
library(grid)

gene50_tss <- read.delim("~/cfDNA/gene50_tss.txt", row.names=1)
data299new <- read_excel("data299new.xlsx")
group_349 <- read.csv("group_349.csv")
group_349$people <- as.character(group_349$people)



heatmap <- gene50_tss %>%
  rotate_df() %>%
  mutate(across(everything(), as.numeric))

heatmap<- base::scale(heatmap,center = TRUE, scale = TRUE)
heatmap <- as.data.frame(heatmap)


#join the clinical information
heatmap <- heatmap %>%
  rownames_to_column("sample") %>%
  left_join(group_349,by = "sample") %>%
  rename(ID = people) %>%
  left_join(data299new,by= "ID") 


heatmap_summary <- heatmap %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(across(starts_with("gene"),list(custom_mean = mean), .names = "{.col}_{.fn}")) %>%
  remove_rownames %>% 
  column_to_rownames(var="group")

xx<-str_split_fixed(colnames(heatmap_summary),'_',6)
xx<-xx[,4]
colnames(heatmap_summary)<-xx


heatmap_summary <- t(heatmap_summary) 


anno <- as.data.frame(colnames(heatmap_summary)) 

anno <- anno %>% 
  rename(group = `colnames(heatmap_summary)`) %>%
  mutate(condition = case_when(str_detect(group,"case")~"GDM",
                               str_detect(group,"ctr")~"Control")) %>%
  mutate(trimester = case_when(str_detect(group,"1")~"1st",
                               str_detect(group,"2")~"2nd",
                               str_detect(group,"3")~"3rd")) %>%
  column_to_rownames("group")



ann_colors = list(condition= c(GDM = "#E3242B",
                               Control = "#17486f"),
                  trimester = c(`1st` = "#006400",
                                `2nd` = "#f6c200",
                                `3rd` = "#d04e00"))



coul <- colorRampPalette(brewer.pal(12, "RdBu"))(50)


pdf(file = "figure3b.pdf",width = 10, height = 12)
pheatmap(heatmap_summary,
         color = rev(coul),
         clustering_distance_rows = "correlation",
         cluster_cols = F,
         cluster_rows = T,
         show_colnames  = F,
         gaps_col = 3,
         cutree_rows = 4,
         row_km = 4,
         annotation_col  = anno,
         annotation_colors =  ann_colors,
         border_color = "white",
         legend = TRUE,
         name="z-score",
         fontsize = 15
)

grid.text("z-score" , x= 0.82, y = 0.8)
dev.off()


