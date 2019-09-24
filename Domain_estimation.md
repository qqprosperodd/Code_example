
# Multiple ChIP analysis and Genomic Domain analysis using HMMt

The goal of ChIP seq is comparing multiple datas and having insights of
Chromatin States.

Here, this is one of the example.

Filion et al., Cell, 2010 is the sourse. They analysed 55 DamID data by
PCA and HMMt. Finally they classified dm genome into 5 domains.

### Unix code

Using deeptools’ bamcompare, caluculate binalysed log2(IP/input).

``` r
#first get coverage file.
cd /Volumes/HDCZ-UT/Domain/;
for Y in 200 500 1000; do
for X in  SRR3503057 SRR3503062 SRR3503054  ; do
bamCompare -b1 ../mapping/keio_ChIPseq/${X}_sort.bam -b2 ../mapping/keio_ChIPseq/SRR3503050_sort.bam  \
-o ${Y}bp/${X}_${Y}.bedgraph -of bedgraph --binSize=${Y} --minMappingQuality 15 \
--normalizeUsing None -p 6 --scaleFactorsMethod readCount --operation log2;
done;
for X in  SRR7992744 SRR7992749 ; do
bamCompare -b1 ../mapping/Brennecke_ChIP2/${X}_sort.bam -b2 ../mapping/Brennecke_ChIP2/SRR7992739_sort.bam  \
-o ${Y}bp/${X}_${Y}.bedgraph -of bedgraph --binSize=${Y} --minMappingQuality 15 \
--normalizeUsing None -p 6 --scaleFactorsMethod readCount --operation log2;
done;
done;
#you can get binsize =200, 500, 1000 bp coverage file.
```

### R code

First, bedgraph files are read and gathered into single large table.

After loess-normalization, correlation of each samples, clustering, PCA
analysis.

``` r
library(tidyverse)
library(ggrepel)
library(GGally)
library(ggcorrplot)
library(ggdendro)
#ggdendro is hclust tools.

fig <- function(A, B){
  #This function's argument is bin-length.
  #A = "200bp/"
  #B = 200
  files3 <- dir(paste0("~/../Desktop/Domain/", A), pattern = "\\.bedgraph$", full.names = TRUE)
  files4 <- files3 %>%
    str_replace_all("^C:/Users/flust/Documents/../Desktop/Domain/|.bedgraph$", "") %>% 
    str_replace_all(A, "") %>% 
    str_replace_all(paste0("_", B), "")
  #in order to get shorter name
  
  df2 <- vector("list", length(files3))
  for (i in seq_along(files3)) {
    df2[[i]] <- read_tsv(files3[[i]], 
                         col_names = c("chr", "start", "end", "cov"))
    df2[[i]][["type"]] <- files4[[i]]
  }
  #bedgraph may have other cols, but 4 cols are necessary for the analysis.
  tx2 <- df2 %>% bind_rows() 
  txalpha <- tx2 %>%
    na_if(0) %>% 
    drop_na() %>% 
    filter(end - start == B) %>% 
    filter(chr != "chrY") %>% 
    spread(key = type, value = cov)
  #In order to delete non-coverage region, drop these region.  
  #Last, spread the dataframe. You should arrange this dataframe (Later in this example).
  tx3 <- txalpha %>% mutate(center = (start + end) /2) %>% 
    gather(-chr, -start, -end, -center, key = "type", value = "cov") %>% 
    group_by(chr, type) %>% nest() %>% 
    mutate(oberon = map(data, function(x) {stats::loess(x$cov ~ x$center, span = 0.1, na.action = na.exclude)}),
           titania = map2(data, oberon, function(x,y) {stats::predict(y, x$center)})) %>%
    dplyr::select(-oberon) %>% unnest() %>% 
    mutate(loess_nom = cov - titania) %>% 
    dplyr::select(chr, start, end, loess_nom, type)
  #loess-normalization
  
  tx4 <- tx3 %>%
    spread(key = type, value = loess_nom) %>% arrange(chr, start)
  #Don't skip arrange function.  You get wrong data.
  
  tx4 %>% write_tsv(paste0("~/../Desktop/Domain/",A, "loess_nom_", B ,".tsv"))
  
  #####
  #Correlation
  COR1 <- tx4 %>% select(-chr, -start, -end) %>% 
    drop_na() %>% 
    cor(method = "pearson")
  COR1 %>% ggcorrplot(hc.order = TRUE, method = "square",
                      colors = c("blue4", "white", "red4"), lab = TRUE)
  ggsave(paste0("~/../Desktop/Domain/",A, "pearson.png"), width =8, height =7, dpi = 300)
  #pearson
  COR2 <- tx4 %>% select(-chr, -start, -end) %>% 
    drop_na() %>% 
    cor(method = "spearman")
  COR2 %>% ggcorrplot(hc.order = TRUE, method = "square",
                      colors = c("blue4", "white", "red4"), lab = TRUE)
  ggsave(paste0("~/../Desktop/Domain/",A, "spearman.png"), width =8, height =7, dpi = 300)
  #spearman
  
  #####
  #hclust
  hcl <- tx4 %>% select(-chr, -start, -end) %>% 
    drop_na() %>% 
    as.matrix() %>% t() %>% 
    dist() %>% hclust(method = "ward.D2")
  #dist() usually caluculate eucledian distance.  method "ward" is similar to "k-means".
  ggdendrogram(hcl, rotate = TRUE, size = 2)
  ggsave(paste0("~/../Desktop/Domain/",A, "hcl_ward.png"), width =3, height =3, dpi = 300)
  
  #####
  #PCA
  PCA <- tx4 %>% select(-chr, -start, -end) %>% 
    drop_na() %>% 
    as.matrix() %>% t() %>% 
    prcomp()
  
  PCA_scrip <- PCA$x[,1:10] %>% as_tibble() %>% 
    mutate(COF = row.names(PCA$x[1:10,]))
  
  PCA1 <- PCA$sdev %>% enframe() %>% 
    mutate(Comp = paste0("PC", 1:10))
  
  PCA1 %>% mutate(perval = value^2*100/sum(value^2)) %>% 
    ggplot(aes(x = perval, y = fct_reorder(Comp, perval))) +
    geom_point() + xlab("% variance explained") + ylab("component") +
    geom_text(aes(label = round(perval,digits = 2)), vjust = 0, hjust=-0.4, colour = "black",size = 3)+
    coord_cartesian(xlim = c(0, 90))
  ggsave(paste0("~/../Desktop/Domain/",A, "Variance_PCA2.png"), width =3 , height =4, dpi = 300)
  #screep plot
  PCA2 <- PCA_scrip %>% gather(-COF, key= "component", value = "value") 
  
  ggpairs(data = PCA_scrip, 
          columns = 1:5,
          aes(color = COF),
          upper = "blank",
          diag = "blank",
          axisLabels = "show",
          switch = "both")
  ggsave(paste0("~/../Desktop/Domain/",A, "Variance_PCA.png"), width =6, height =5, dpi = 300)
  #scatter plot of component1~5
  PCA2 %>% mutate(color2 = factor(component, levels=paste0("PC",1:15))) %>% 
    ggplot(aes(x = COF, y = value, color = color2, fill = color2)) +
    geom_bar(stat = "identity") + 
    facet_grid(color2 ~ ., scales = "free_y") +
    xlab("COF") + ylab("Component Loading")+
    theme(axis.text.x = element_text(angle = 30, hjust=1, vjust=1)) +
    labs(fill = "Component", color = "Component")
  ggsave(paste0("~/../Desktop/Domain/",A, "Variance_PCA3.png"), width =6, height =15, dpi = 300)
  #component 
  
  PCA_vari <- PCA1 %>% mutate(perval = signif(value^2*100/sum(value^2), digits = 4)) %>% 
    select(perval)
  PCA_scrip %>% 
    ggplot(aes(x = PC1, y = PC2, color = COF)) +
    geom_point() +
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    xlab(paste0("PC1(", PCA_vari[1,], "%)")) + 
    ylab(paste0("PC2(", PCA_vari[2,], "%)"))+
    geom_text_repel(aes(label = COF), 
                    size = 3.5,
                    box.padding = unit(0.4,"lines"),
                    segment.color = "#cccccc",
                    segment.size = 0.3,
                    fontface = "bold",
                    point.padding = unit(0.3,"lines"))
  ggsave(paste0("~/../Desktop/Domain/",A, "Variance_PCA4.png"), width =7, height =5, dpi = 300)
  #comp1 vs comp2
}

fig("200bp/", 200)
fig("500bp/", 500)
fig("1000bp/", 1000)
```

### HMMt (R code)

In order to create PC1 vs PC2…, recomputing PCA.

``` r
fig2 <- function(A, B, C){
  #A = "200bp/"
  #B = 200
  tx4 <- read_tsv(paste0("~/../Desktop/Domain/",A, "loess_nom_", B ,".tsv"))
  PCA <- tx4 %>% select(-chr, -start, -end) %>% 
    drop_na() %>% 
    as.matrix() %>% t() %>% 
    prcomp()
  PCA1 <- PCA$sdev %>% enframe() %>% 
    mutate(Comp = paste0("PC", 1:10))
  PCA_vari <- PCA1 %>% mutate(perval = signif(value^2*100/sum(value^2), digits = 4)) %>% 
    select(perval)
  
  #HMMT
  loess_nom <- HMMt::bridge(tx4 %>% as.data.frame())
  HMM_baum <- HMMt::BaumWelchT(loess_nom[["x"]], m=C)
  #"bridge" covers non-coverage region.
  
  Result <- HMM_baum$ViterbiPath %>% 
    enframe(value = "state") %>% 
    mutate(visible = loess_nom$nonvirtuals) %>% 
    filter(visible == TRUE) %>% 
    bind_cols(tx4) 
  Result %>% write_tsv(paste0("~/../Desktop/Domain/",A, "HMMt_after_all_state", C,".tsv"))
  #save the data.
  
  PCAalpha <- Result %>% select(-chr, -start, -end, -name, -state, -visible) %>% 
    drop_na() %>% 
    as.matrix() %>% 
    prcomp()
  
  PCA_scrip <- PCAalpha$x %>% as_tibble() %>% 
    bind_cols(Result %>% drop_na() %>% select(chr, start, end, state)) %>% 
    mutate(state2 = paste0("state", state))
  
  PCA_scrip %>% sample_n(10000) %>% 
    ggplot(aes(x = PC1, y = PC2, color = state2)) +
    geom_point(alpha = 0.5)+
    xlab(paste0("PC1(", PCA_vari[1,], "%)")) + 
    ylab(paste0("PC2(", PCA_vari[2,], "%)"))
  ggsave(paste0("~/../Desktop/Domain/", A, "state", C, "PCA1.png"), width =12, height =10, dpi = 300)
  #PC1 vs PC2
  PCA_scrip %>% sample_n(10000) %>% 
    ggplot(aes(x = PC2, y = PC3, color = state2)) +
    geom_point(alpha = 0.5)+
    xlab(paste0("PC2(", PCA_vari[2,], "%)")) + 
    ylab(paste0("PC3(", PCA_vari[3,], "%)"))
  ggsave(paste0("~/../Desktop/Domain/",A, "state", C, "PCA2.png"), width =12, height =10, dpi = 300)
  #PC2 vs PC3
  PCA_scrip %>% sample_n(10000) %>% 
    ggplot(aes(x = PC1, y = PC3, color = state2)) +
    geom_point(alpha = 0.5)+
    xlab(paste0("PC1(", PCA_vari[1,], "%)")) + 
    ylab(paste0("PC3(", PCA_vari[3,], "%)"))
  ggsave(paste0("~/../Desktop/Domain/",A, "state", C, "PCA3.png"), width =12, height =10, dpi = 300)
  #PC1 vs PC3
}
fig2("200bp/", 200, 3)
fig2("200bp/", 200, 4)
fig2("200bp/", 200, 5)
fig2("500bp/", 500, 3)
fig2("500bp/", 500, 4)
fig2("500bp/", 500, 5)
fig2("1000bp/", 1000, 3)
fig2("1000bp/", 1000, 4)
fig2("1000bp/", 1000, 5)
```

### Save the domain data (R code)

Bind the bed next to each other.

``` r
Domainfun <- function(A,B,C) {
  Domain <- read_tsv(paste0("~/../Desktop/Domain/",A, "HMMt_after_all_state", C,".tsv"))
  Jurius <- Domain %>% 
    mutate(state = paste0("state", state)) %>% 
    dplyr::select(name, state, chr, start, end) %>% 
    group_by(state, chr) %>% nest() %>% 
    mutate(pack = map(data, function(x) {which(diff(x$name) != 1)}),
           end_merge = map2(data, pack, function(x,y) {x$end[y]}),
           demetrius = map(pack, function(x) {append(x+1,1,after =0) %>% head(-1)}),
           start_merge = map2(data, demetrius, function(x,y) {x$start[y]})) %>% 
    dplyr::select(chr, start_merge, end_merge, state) %>% 
    unnest() %>% arrange(chr, start_merge) %>% 
    select(chr, start_merge, end_merge,state) %>%
    filter(start_merge < end_merge) %>% 
  write_tsv(paste0("~/../Desktop/Domain/",A, "state", C, "Domain.bed"),
            col_names = FALSE)
}
Domainfun("200bp/", 200, 3)
Domainfun("200bp/", 200, 4)
Domainfun("200bp/", 200, 5)
Domainfun("500bp/", 500, 3)
Domainfun("500bp/", 500, 4)
Domainfun("500bp/", 500, 5)
Domainfun("1000bp/", 1000, 3)
Domainfun("1000bp/", 1000, 4)
Domainfun("1000bp/", 1000, 5)
```

### HMMt of each beadgraph (R code)

Each factor’s bind sites.

``` r
writeDomain <- function(A, B){
  tx4 <- read_tsv(paste0("~/../Desktop/Domain/",A, "loess_nom_", B ,".tsv"))
  test <- tx4 %>% gather(-chr, -start, -end, key = "type", value = "cov") %>% 
    group_by(type) %>% nest() %>% 
    mutate(loess_norm = map(data, function(x) {HMMt::bridge(x %>% as.data.frame())}),
           HMM_baum = map(loess_norm, function(x) {HMMt::BaumWelchT(x[["x"]], m=2)}))
  
  test2 <- test %>%
    mutate(HMM_baum2 = map(HMM_baum, function(x) {x$ViterbiPath %>%
        enframe(value = "state")}),
        visible = map(loess_norm, function(x) {x$nonvirtuals})) %>% 
    select(HMM_baum2, visible) %>% unnest() %>% 
    filter(visible == TRUE) %>% 
    bind_cols(tx4 %>% gather(-chr, -start, -end, key = "type", value = "cov")) 
  
  Jurius <- test2 %>% 
    mutate(state = paste0("state", state)) %>% 
    dplyr::select(name, state, chr, start, end, type) %>% 
    group_by(state, chr, type) %>% nest() %>% 
    mutate(pack = map(data, function(x) {which(diff(x$name) != 1)}),
           end_merge = map2(data, pack, function(x,y) {x$end[y]}),
           demetrius = map(pack, function(x) {append(x+1,1,after =0) %>% head(-1)}),
           start_merge = map2(data, demetrius, function(x,y) {x$start[y]})) %>% 
    dplyr::select(chr, start_merge, end_merge, state, type) %>% 
    unnest() %>% arrange(type, chr, start_merge) %>% 
    filter(state == "state2") %>%
    mutate(state = "bound") %>% 
    select(chr, start_merge, end_merge,state, type) %>%
    filter(start_merge < end_merge) %>%
    group_by(type) %>% 
    do(write_tsv(., paste0("~/../Desktop/Domain/",A, unique(.$type), "_Domain.bed"),
                 col_names = FALSE))
}

writeDomain("200bp/", 200)
writeDomain("500bp/", 500)
writeDomain("1000bp/", 1000)
```

All unix code was done on Anaconda3.

All R was done in R3.6 on windows10. If you use MacOS, you may have some
trouble in installing HMMt. See
<https://www.rstudio.com/products/rpackages/devtools/>
