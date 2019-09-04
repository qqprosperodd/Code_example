
## Drosophila promoter classification

Process the data to understand what cofactors can affect to the promoter
(Haberle et al., Nature, 2019)
<https://www.nature.com/articles/s41586-019-1210-7>.

This classification is for S2 cell.

In the article, they did STAP-seq in OSC, but STAP-seq in OSC and that
in S2 is so different (See, FigS4).  
So, I don’t know if we can apply the classification in S2 to OSC.

Anyway, if you want to improve your statistic skill, you can leran PCA,
k-means, visualization of heatmaps.

Scaling the value, k-means, main figure of heatmaps.

``` r
setwd("~/Desktop/GSE1161197_bw/")
library(tidyverse)
library(readxl)
library(gplots)
library(RColorBrewer)
library(GGally)
library(ggrepel)

haberle <- read_xlsx("41586_2019_1210_MOESM12_ESM.xlsx")
#sup-12 file contains averaged-tagcount of S2. first, scale(per COF). 
haberle2 <- haberle %>% 
  select(oligo_id, Chro, Med25, Mof, Nej_p300, p65, Lpt, trr, trx,
         CG15356_EMSY, CG7154_Brd9, Brd8, GFP, Med15, gfzf, fs1h_Brd4) %>%
  gather(-oligo_id, key = "COF", value = "nom_tag") %>% 
  group_by(oligo_id) %>% nest() %>% 
  mutate(log_z = map(data, function(x) {scale(x$nom_tag)})) %>% unnest() %>% 
  select(-nom_tag) %>% spread(COF, log_z)
#zscore can calcurate by function "scale"
haberle2 %>% select(-oligo_id) %>% gather(key = "type", value = "value") %>% 
  ggplot() + geom_density(aes(x = value, color=type)) +
  coord_cartesian(xlim = c(-2,4), ylim = c(0,4))
ggsave("~/Desktop/GSE1161197_bw/density3.png", width =6, height =6, dpi = 300)

km <- haberle2 %>% 
  select(Chro, Med25, Mof, Nej_p300, p65, Lpt, trr, trx, CG15356_EMSY,
         CG7154_Brd9, Brd8, GFP, Med15, gfzf, fs1h_Brd4) %>% 
  kmeans(centers = 5, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")
#k-means select one site randomly in each step.
#So, the start site isn't good, the solution sometimes become local solution.
#So, you must set argument "nstart" to 10 >.
compload <- km[["centers"]] %>% as_tibble() %>% mutate(cluster = paste("cluster", 1:5))
compload %>% gather(-cluster, key = "COF", value = "Mean") %>% 
  ggplot(aes(x = COF, y = Mean, fill = cluster, color = cluster)) +
  geom_bar(stat = 'identity', position = "identity", width = 0.75) +
  facet_grid(cluster ~ ., scales = 'free_y') +
  xlab("COF") + ylab("Component Loading") +
  theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1))
ggsave("~/Desktop/GSE1161197_bw/component_loading.png", width =6, height =6, dpi = 300)
#component loading
km$size %>% as_tibble() %>%
  mutate(cluster = paste("cluster", 1:5)) %>%
  ggplot(aes(x = value, y = cluster)) + geom_point() + coord_cartesian(xlim = c(0,14000)) +
  geom_text(aes(label = value), vjust = 1.5, colour = "black",size = 3)
ggsave("~/Desktop/GSE1161197_bw/cluster_component.png", width =3, height =3, dpi = 300)
#the count of component of each cluster

#heamap
#gplots::heatmap.2() takes lots of times, so sampling 2000 records randomly.
km3 <- haberle2 %>% mutate(km = km$cluster,
                                         km2 = if_else(km ==1,"red",
                                                       if_else(km ==2, "orange",
                                                               if_else(km == 3, "purple",
                                                                       if_else(km == 4, "blue", "green"))))) %>% 
  sample_n(2000) %>% arrange(km)

colfunc <- colorRampPalette(c("blue4", "white", "red4"))
#color palette for heatmaps
ppi <- 300
png("~/Desktop/GSE1161197_bw/cor190805.png", width = 15*ppi, height = 10*ppi, res = ppi)
km3 %>% select(Chro, Med25, Mof, Nej_p300, p65, Lpt, trr, trx, CG15356_EMSY,
               CG7154_Brd9, Brd8, GFP, Med15, gfzf, fs1h_Brd4) %>% 
  as.matrix() %>%  t() %>%
  heatmap.2(col=colfunc(75), ColSideColors = km3$km2,
            density.info = "none", trace = "none", dendrogram="row",
            labCol = "", Colv=NA, margin=c(5, 10))
dev.off()
#heatmap
#If you want to reduce drawing time, set "RowA = NA".

haberle2 %>% mutate(km = km$cluster,
                           km2 = if_else(km ==1,"red",
                                         if_else(km ==2, "orange",
                                                 if_else(km == 3, "purple",
                                                         if_else(km == 4, "blue", "green"))))) %>% 
  write_tsv("kmmeans_result.tsv")
km_full <- read_tsv("kmmeans_result.tsv")
#write to tsv file
#kmeans may return different result in each execution.  
#So I strongly recommend you to write the result and to use the result in below analysis.
```

Below, other figure of the article.

``` r
#######################
#when you choose the number of clusters.
#######################
df <- haberle2 %>% 
  select(Chro, Med25, Mof, Nej_p300, p65, Lpt, trr, trx, CG15356_EMSY,
         CG7154_Brd9, Brd8, GFP, Med15, gfzf, fs1h_Brd4)
pct_var <- data.frame(pct_var = 0, num_clusters = 2:10)
totalss <- kmeans(df, centers = 10, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")$totss
for(i in 2:10){
  pct_var[i-1, 'pct_var'] <- kmeans(df, centers = i, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")$betweenss/totalss
}
ggplot(pct_var, aes(x = num_clusters, y = pct_var)) + 
  geom_line() + geom_point() + theme_bw() +
  ylim(0, 1) + xlab("Number of clusters") + ylab("% Variance Explained")
ggsave("~/Desktop/GSE1161197_bw/Variance_kmeans0806.png", width =6, height =5, dpi = 300)

#################
#PCA analysis
################
PCA <- haberle %>% 
  select(oligo_id, Chro, Med25, Mof, Nej_p300, p65, Lpt, trr, trx, CG15356_EMSY,
         CG7154_Brd9, Brd8, GFP, Med15, gfzf, fs1h_Brd4) %>%
  gather(-oligo_id, key = "COF", value = "nom_tag") %>% 
  group_by(COF) %>% nest() %>% 
  mutate(log_z = map(data, function(x) {log10(x$nom_tag+1)})) %>% unnest() %>% 
  select(-nom_tag) %>% spread(COF, log_z) %>% 
  select(Chro, Med25, Mof, Nej_p300, p65, Lpt, trr, trx, CG15356_EMSY,
         CG7154_Brd9, Brd8, GFP, Med15, gfzf, fs1h_Brd4) %>%
  as.matrix() %>% t() %>% 
  prcomp()
#If you remove t(), you can get the score of each TSS region.  This score is useful when you visualize PC1 vs PC2 with the k-means result.

PCA_scrip <- PCA$x[,1:15] %>% as_tibble() %>% 
  mutate(COF = row.names(PCA$x[1:15,]))
PCA1 <- PCA$sdev %>% as_tibble() %>% 
  mutate(Comp = paste0("PC", 1:15))

PCA1 %>% mutate(perval = value^2*100/sum(value^2)) %>% 
  ggplot(aes(x = perval, y = fct_reorder(Comp, perval))) +
  geom_point() + xlab("% variance explained") + ylab("component") +
  geom_text(aes(label = round(perval,digits = 2)), vjust = 0, hjust=-0.4, colour = "black",size = 3)+
  coord_cartesian(xlim = c(0, 90))
ggsave("~/Desktop/GSE1161197_bw/Variance_PCA2.png", width =3 , height =4, dpi = 300)
#screep plot

PCA2 <- PCA_scrip %>% gather(-COF, key= "component", value = "value") 
ggpairs(data = PCA_scrip, 
        columns = 1:5,
        aes(color = COF),
        upper = "blank",
        diag = "blank",
        axisLabels = "show",
        switch = "both")
ggsave("~/Desktop/GSE1161197_bw/Variance_PCA.png", width =6, height =5, dpi = 300)
#component1~5 PCA plot

PCA2 %>% mutate(color2 = factor(component, levels=paste0("PC",1:15))) %>% 
  ggplot(aes(x = COF, y = value, color = color2, fill = color2)) +
  geom_bar(stat = "identity") + 
  facet_grid(color2 ~ ., scales = "free_y") +
  xlab("COF") + ylab("Component Loading")+
  theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1)) +
  labs(fill = "Component", color = "Component")
ggsave("~/Desktop/GSE1161197_bw/Variance_PCA3.png", width =6, height =15, dpi = 300)
#component loading

PCA_scrip %>% 
  ggplot(aes(x = PC1, y = PC2, color = COF)) +
  geom_point() +
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  xlab("PC1(76.2%)") + ylab("PC2(12%)")+
  geom_text_repel(aes(label = COF), 
                  size = 3.5,
                  box.padding = unit(0.4,"lines"),
                  segment.color = "#cccccc",
                  segment.size = 0.3,
                  fontface = "bold",
                  point.padding = unit(0.3,"lines"))
ggsave("~/Desktop/GSE1161197_bw/Variance_PCA4.png", width =7, height =5, dpi = 300)
#comp1 vs comp2
```

All R code was done on R3.6.
