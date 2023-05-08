library(tidyverse)
library(readr)
library(readxl)
library(hrbrthemes)
library(factoextra) 
library(Gmisc)
library(ggraph)
library(igraph) 
library(tidygraph)
library(psych)
library(graphlayouts)
library(wesanderson)


set.seed(123)

################################################################
#load data
################################################################
gene50_tss <- read.delim("~/cfDNA/gene50_tss.txt", row.names=1)
data299new <- read_excel("data299new.xlsx")
group_349 <- read.csv("group_349.csv")
group_349$people <- as.character(group_349$people)



heatmap <- gene50_tss %>%
  rotate_df() %>%
  mutate(across(everything(), as.numeric))

heatmap<- base::scale(heatmap,center = TRUE, scale = TRUE)
heatmap <- as.data.frame(heatmap)

xx<-str_split_fixed(colnames(heatmap),'_',6)
xx<-xx[,4]
colnames(heatmap)<-xx
corr <- heatmap


################################################################
#load pathway information
################################################################
pathway <- read.csv("GO_50gene_label.csv") %>%
  select(Description,X50gene) %>%
  separate(X50gene,c("a","b","c","d","e"))

node <- read_excel("GSEA-wiki-padj0.05-50gene2pathway(1).xlsx") %>%
  dplyr:: distinct(Description,.keep_all = T) %>%
  select(Description,gene) %>%
  drop_na(gene) %>%
  separate(gene,c("a","b","c")) %>%
  pivot_longer(-Description,names_to = "term12") %>%
  drop_na(value) %>%
  select(-term12) %>%
  dplyr::rename(term = Description) %>%
  distinct(term,.keep_all = T)

hit_list <- data.frame()

for (i in (1:nrow(pathway))){
  hit_i <- as.data.frame(t(combn(pathway[i,2:6],m=2)))
  hit_i[,3] <- pathway[i,1]
  hit_list <- rbind(hit_list,hit_i)
}

hit_list <- hit_list[is.na(hit_list$V2)==FALSE,] 
hit_list$V1 <- trimws(hit_list$V1)
hit_list$V2 <- trimws(hit_list$V2)
hit_list$name1 <- paste0(hit_list$V1,hit_list$V2)
hit_list$name2 <- paste0(hit_list$V2,hit_list$V1)
hit_list <- hit_list[,3:5]


corr_matrix <- cor(corr,method = c("spearman"))
corr_matrix[abs(corr_matrix) < 0.06] <- 0
diag(corr_matrix) <- 0
corr_matrix <- abs(corr_matrix)

ig <- graph.adjacency(corr_matrix, weighted=TRUE, mode="lower") %>%
  set_edge_attr("label", index = 10,value = "A")

edeglist <-as.data.frame(get.edgelist(ig))
edeglist <- edeglist %>%
  dplyr::mutate(name=paste0(V1,V2)) %>%
  dplyr::rename(name1=name) %>%
  plyr::join(hit_list,type = 'left',match = 'first',by = "name1") %>%
  dplyr::select(-name2) %>%
  dplyr::rename(name2=name1,
                V31=V3) %>%
  plyr::join(hit_list,type = 'left',match = 'first',by = "name2") %>%
  dplyr::mutate(categroy = coalesce(V3,V31)) 


edeglist$categroy[is.na(edeglist$categroy)]<-"unmapped"

E(ig)$pathway <- edeglist$categroy
vertex.color <- read.csv("vertex.csv")
V(ig)$color <- vertex.color$class



new_graph <- graph_from_data_frame(d=as_data_frame(ig,what="edges") %>%
                                     dplyr:: arrange(desc(pathway)),
                                   vertices=as_data_frame(ig,what="vertices"))

xy <- layout_with_stress(ig)
V(ig)$x <- (xy[, 1])
V(ig)$y <- (xy[, 2])

ggraph(new_graph,"manual", x = V(ig)$x, y = V(ig)$y) +
  geom_edge_link2(aes(color=pathway,width = log(abs(weight))), 
                  start_cap = circle(2, 'mm'),
                  end_cap = circle(2, 'mm'),
                  check_overlap = TRUE) + 
  geom_node_point(aes(color = color),size = 3) + 
  scale_edge_width(range = c(1,2),name="log(weight)")+
  geom_node_label(aes(label = name),label.padding = unit(0.08, "lines"),
                  nudge_x = 0.05,nudge_y = 0.1,repel = FALSE) +
  scale_color_manual(name="type of TSS scores",values = wes_palette("BottleRocket2"))+
  scale_edge_color_manual(values = c("#a82203","#208cc0","#533d14","#cf5e4e","#3b7c70","#3b3a3e","#924099","#df9ed4","#c59349",
                                              "#637b31", "#003967","#ECECEC"))+
                                                theme_void() +
  theme(text = element_text(size=14,family="sans"),
        plot.title  = element_text(family="sans"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"))+
  guides(edge_color = guide_legend(override.aes = list(edge_width = 5))) 

ggsave("figure4d.pdf",width = 28, height = 16,units = "cm")

