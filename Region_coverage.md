
# Get region coverage, scatter plot and box plot.

We take so many times for analysis of Bed files. R helps us take less
time for analysis.

In this example, I get coverage of 9 bed region from 33 bam files.

### R code

You chop genome into short length fragment.

``` r
#full genome file (for CPM calculation)
library(tidyverse)
genome <- read_tsv("/Volumes/HDCZ-UT/genome/UCSC_dm6/dm6.genome", col_names = c("chr1", "end"))
genome2 <- genome %>% mutate(start = 0) %>%
  select(chr1, start, end) %>%
  filter(chr1 != "chrY")
genome2 %>% write_tsv("/Volumes/HDCZ-UT/chipmeta/dm6_full_genome.bed", col_names = FALSE)

#In this example, I've done same to 5 files.
setwd("/Volumes/HDCZ-UT/chipmeta/")
filedir <- "/Volumes/HDCZ-UT/chipmeta/"
files <- dir(paste0(filedir, "ref_summit_bed/"), pattern = "\\.bed$", full.names = TRUE)
files2 <- files %>%
  str_replace_all("^/Volumes/HDCZ-UT/chipmeta/ref_summit_bed//|.bed$", "") %>%
  str_split(pattern = "\\_", simplify = TRUE) %>%
  .[,1] %>%
  str_replace("OSC", "DHS")
#Using Regular expression, you can get short name of the file.
df <- vector("list", length(files))
for (i in seq_along(files)) {
  df[[i]] <- read_tsv(files[[i]], col_names = c("chr1", "start1", "end1"),
                      col_types = cols(chr1 = col_character(),
                                       start1 = col_double(),
                                       end1 = col_double()))
  df[[i]][["type"]] <- files2[[i]]
}
tx2 <- bind_rows(df) %>%
  filter(chr1 != "chrY") %>%
  filter(!(str_detect(chr1, pattern = "\\_"))) %>%
  distinct()
#write
writeBED <- function(z, y) {
  theseus <- z %>% group_by(type) %>% nest()
  hippolyta <- map(1:length(theseus$type), function(x) {data.frame(theseus$data[[x]])})
  names(hippolyta) <- theseus$type
  map(1:length(hippolyta), function(x) write_tsv(hippolyta[[x]],
                                                 paste0(filedir , names(hippolyta[x]), y, ".bed"),
                                                 col_names = FALSE))
}
writeBED(tx2, "")
#new peak(filter peaks)
tx2 %>% mutate(start2 = start1 - 500,
       end2 = end1 + 500) %>%
  select(chr1, start2, end2, type) %>% writeBED("_extended")
#1000bp extended  peaks
```

### Unix code

To get active/inactive enhancer region, I used bedtools intersect.

I refer to
<http://www.sciencemag.org/cgi/pmidlookup?view=long&pmid=23328393>.

I get data from <https://www.nature.com/articles/nature13994>.

``` r
cd /Volumes/HDCZ-UT/chipmeta/;
#get inactive enhancer peak.
for X in dCP hkCP; do
bedtools intersect -v -a ${X}.bed -b DHS_extended.bed > ${X}_close.bed;
done;
#enhancers that only activate hkCP/dCP.
bedtools intersect -v -a dCP.bed -b hkCP_extended.bed > dev_enhancer.bed;
bedtools intersect -v -a hkCP.bed -b dCP_extended.bed > hk_enhancer.bed;
```

### R code

you get 1000 length extended bed file of inactive enhancers and hkCP/dCP
only enhancers.

``` r
library(tidyverse)
filedir <- "/Volumes/HDCZ-UT/chipmeta/"
setwd(filedir)
hkCP <- read_tsv("hkCP_close.bed",
                 col_names = c("chr1", "start1", "end1"))
dCP <- read_tsv("dCP_close.bed",
                 col_names = c("chr1", "start1", "end1"))
dev <- read_tsv("dev_enhancer.bed",
                col_names = c("chr1", "start1", "end1"))
hk <- read_tsv("hk_enhancer.bed",
                col_names = c("chr1", "start1", "end1"))
#Using Regular expression, you get short name of the file.
writeBED2 <- function(z, y, fes) {
  z %>% mutate(start2 = start1 - 500,
               end2 = end1 + 500) %>%
    select(chr1, start2, end2) %>%
  write_tsv(paste0(filedir, fes, y, ".bed"), col_names = FALSE)
}
writeBED2(hkCP, "_extended", "hkCPclose")
writeBED2(dCP, "_extended", "dCPclose")
writeBED2(dev, "_extended", "dCPonly")
writeBED2(hk, "_extended", "hkCPonly")
```

### Unix code

you get coverage of multiple region by samtools bedcov.

``` r
cd /Volumes/HDCZ-UT/chipmeta/;
for X in [[files]]; do
samtools bedcov -Q 15 ${X}_extended4.bed Bamfile/*/*.bam > Coverage/5000bp/${X}_cov.bed;
samtools bedcov -Q 15 ${X}_extended3.bed Bamfile/*/*.bam > Coverage/2000bp/${X}_cov.bed;
samtools bedcov -Q 15 ${X}_extended2.bed Bamfile/*/*.bam > Coverage/10000bp/${X}_cov.bed;
samtools bedcov -Q 15 ${X}_extended.bed Bamfile/*/*.bam > Coverage/1000bp/${X}_cov.bed;
done;
for X in dm6_full_genome; do
samtools bedcov -Q 15 ${X}.bed Bamfile/*/*.bam > Coverage/5000bp/${X}_cov.bed;
samtools bedcov -Q 15 ${X}.bed Bamfile/*/*.bam > Coverage/2000bp/${X}_cov.bed;
samtools bedcov -Q 15 ${X}.bed Bamfile/*/*.bam > Coverage/10000bp/${X}_cov.bed;
samtools bedcov -Q 15 ${X}.bed Bamfile/*/*.bam > Coverage/1000bp/${X}_cov.bed;
done;
#samtools bedcov is not good choice, because there is bedtools coverage.
#But when you use *, bedtools doesn't work.

ls Bamfile/*/*.bam;
#output was deleted.
```

### R code

First, calculate CPM (count per million) of all experiment.

Next, get boxplot with wilcoxon rank sum test, violinplot, scatterplot
with PCC.

``` r
#get coverage file into single table.
library(tidyverse)
library(exactRankTests)
#for wilcoxon rank sum test.
library(ggpointdensity)
library(viridis)
#for pointdensity plot

setwd("/Volumes/HDCZ-UT/chipmeta/")
filedir <- "/Volumes/HDCZ-UT/chipmeta/"
files <- dir(paste0(filedir, "Coverage/"), pattern = "\\.bed$", full.names = TRUE)
files2 <- files %>%
  str_replace_all("^/Volumes/HDCZ-UT/chipmeta/Coverage//|.bed$", "") %>%
  str_replace("_extend_cov", "")
#Using Regular expression, you get short name of the file.
file_sample <- read_tsv("sample_name", col_names = "name") %>%
  pull(name) %>%
  str_split(pattern = "\\/", simplify = TRUE) %>%
  .[,3] %>%
  str_replace_all(c(".bam" = "", "_merged" = "", "_sort" = "")) %>%
  str_replace("_IN_bowtie2", "_cp190_IN") %>%
  str_replace("_bowtie2", "")
file_sample
#You have to pull the result of "ls Bamfile/*/*.bam" into text file in advance. 
df <- vector("list", length(files))
for (i in seq_along(files)) {
  df[[i]] <- read_tsv(files[[i]], col_names= c("chr", "start", "end", file_sample))
  df[[i]][["type"]] <- files2[[i]]
}
tx2 <- bind_rows(df) %>%
  unite(chr, start, end, col = "anotation", sep = "-") %>%
  distinct()

#create count file of all region for CPM calculation.
CPM_ctrl <- tx2 %>% 
  filter(type == "dm6_full_genome_cov") %>% 
  select(-anotation, -type) %>% 
  gather(key = "sample", value = "cov_ctrl") %>% 
  group_by(sample) %>% summarise(cov_ctrl = sum(cov_ctrl))

#calculate CPM.
CPM <- tx2 %>% 
  gather(-anotation, -type, key = "sample", value = "cov") %>% 
  left_join(CPM_ctrl, by = "sample")  %>% 
  group_by(sample) %>% nest() %>% 
  mutate(CPM = map(data, function(x) {log2(x$cov * 1000000/x$cov_ctrl + 1)})) %>% 
  unnest() %>% select(-cov, -cov_ctrl) %>% 
  filter(type != ("dm6_full_genome_cov"))

#wilcox_test&boxplot violin plot, scatter plot
fig <- function(x, y) {
  #this function requires: 1, annotation; 2, type; 3, sample; 4, CPM;
  #wilcoxon rank sum test
  datalabel <- x %>% group_by(type) %>%
    summarize(results = wilcox.exact(CPM ~ sample, paired = T)$p.value,
              count = n()) %>%
    mutate(result2 = paste("italic(p) ==", signif(results, digits = 2)),
           count2 = paste("italic(n) ==", count),
           star = if_else(results > 1e-2, "NS",
                          if_else(results > 1e-4, "*",
                                  if_else(results > 2.2e-16, "**","***"))))
  #violin plot
  graph2 <- x %>% ggplot(aes(x = sample, y = CPM, fill = type)) +
    geom_violin() +
    facet_wrap(~ type)
  graph2
  ggsave(paste0(filedir, y, "2.png"), width =12, height =9, dpi = 300)
  #boxplot
  graph <- x %>% ggplot(aes(x = sample, y = CPM, fill = type)) +
    geom_boxplot() + facet_wrap(~ type) +
    geom_text(x = 1.1, y = 9, data = datalabel, aes(label = result2), parse = TRUE, hjust = 0) +
    #display by 3rd decimal point using function round().
    #If you want to move the text position, change the value of geom_text().
    geom_text(x = 1.1, y = 8, data = datalabel, aes(label = count2), parse = TRUE, hjust = 0) +
    geom_text(x = 1.3, y = 10, data = datalabel, aes(label = star), hjust = 0,size = 4.5) +
    ylab(expression(log[2](CPM + 1))) +
    labs(title = paste0(y, "seq"))
  graph
  ggsave(paste0(filedir, y, ".png"), width =12, height =9, dpi = 300)
  #scatter plot
  scatter <- x %>% spread(key = sample, value = CPM) 
  scat_x <- scatter %>% colnames() %>% .[1:2] %>% append(c("siEGFP", "siPiwi"), after = 2)
  colnames(scatter) <- scat_x
  COR <- scatter %>% select(-anotation) %>%
    group_by(type) %>% summarise(pcc = cor(siEGFP, siPiwi)) %>% 
    mutate(result2 = paste("italic(PCC) ==", signif(pcc, digits = 2)))
  DIFalpha <- scatter %>% 
    mutate(threshold = if_else(siEGFP > siPiwi, "down", "up")) %>% 
    group_by(type, threshold) %>% summarise(change = n()) %>% 
    mutate(result2 = if_else(threshold == "up", paste("up =", change),
                             paste("down =", change)))
  graph3 <- scatter %>% ggplot() +
    geom_pointdensity(aes(x = siEGFP, y = siPiwi), size = 0.3, alpha = 0.6) +
    scale_color_viridis() +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_abline(intercept = c(-log2(3), log2(3)), slope = 1, linetype="dashed") +
    geom_text(x = 6, y = 3, data = COR, aes(label = result2), parse = TRUE, hjust = 0,size = 2.5) +
    geom_text(x = 3, y = 0, data = DIFalpha %>% filter(threshold =="down"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 6, data = DIFalpha %>% filter(threshold =="up"), aes(label = result2), hjust = 0,size = 2.5) +
    facet_wrap(~ type)
  graph3 
  ggsave(paste0(filedir, y, "3.png"), width =12, height =9, dpi = 300)
  #graph4
  #FC>2 or FC < 1/2 are shown 
  DIF <- scatter %>% filter(siEGFP !=0, siPiwi != 0) %>% 
    mutate(plus_thresold = log2(2 ^ (siEGFP + 1) - 1),
           minus_threshold = log2(2 ^ (siEGFP - 1) + 1/2),
           state = if_else(siPiwi > plus_thresold & (siPiwi > 2 | siEGFP >2), "plus", 
                           if_else(siPiwi <= minus_threshold & (siPiwi > 2 | siEGFP >2), "minus",
                                   if_else((siPiwi > 2 | siEGFP >2), "stay", 
                                            "NS"))))
  totalDIF <- DIF %>% group_by(type) %>% summarise(countall = n())
  DIF2 <- DIF %>% group_by(type, state) %>% summarise(count = n()) %>% 
    left_join(totalDIF, key = "type") %>% 
    mutate(result2 = paste0(state,"=", count),
           result3 = paste0(signif(count*100 / countall, digits = 2), "%"))
  graph4 <- DIF %>% ggplot() +
    geom_point(aes(x = siEGFP, y = siPiwi, color = state), size = 0.3, alpha = 0.6) + 
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_abline(intercept = c(-log2(3), log2(3)), slope = 1, linetype="dashed") +
    geom_text(x = 0, y = 8, data = DIF2 %>% filter(state =="plus"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 7.5, data = DIF2 %>% filter(state =="plus"), aes(label = result3), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 7, data = DIF2 %>% filter(state =="NS"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 6.5, data = DIF2 %>% filter(state =="NS"), aes(label = result3), hjust = 0,size = 2.5) +
    geom_text(x = 6, y = 2, data = DIF2 %>% filter(state =="stay"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 6, y = 1.5, data = DIF2 %>% filter(state =="stay"), aes(label = result3), hjust = 0,size = 2.5) +
    geom_text(x = 6, y = 1, data = DIF2 %>% filter(state =="minus"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 6, y = 0.5, data = DIF2 %>% filter(state =="minus"), aes(label = result3), hjust = 0,size = 2.5) +
    facet_wrap(~ type)
  graph4 
  ggsave(paste0(filedir, y, "4.png"), width =12, height =9, dpi = 300)
  #"plus" and "minus" region -> bed file
  write_state <- DIF %>% select(anotation, type, state) %>% 
    filter(state =="plus" | state == "minus") %>% 
    separate(anotation, into =  c("chr", "start", "end"), sep = "-")
  write_state %>% group_by(type, state) %>% 
    do(write_tsv(., paste0(filedir , "chara_region/", unique(.$type), "_", y, unique(.$state), ".bed"),
                 col_names = FALSE))
}
# * = p < 1e-2
# ** = p < 1e-4
# *** = p < 2.2e-16
# else = NS (I refer threashold to Hug et al., Cell, 2017.)
#If you want to use other way, see https://github.com/kassambara/ggpubr/issues/65.
#ATAC
ATAC <- CPM %>% filter(str_detect(sample, pattern = "ATAC"))
fig(ATAC, "ATAC")
#GRO
GRO <- CPM %>% filter(str_detect(sample, pattern = "SRR60966"))
fig(GRO, "GRO")
#RNA
poly <- CPM %>% filter(str_detect(sample, pattern = "Bre"))
fig(poly, "poly")
rRNA <- CPM %>% filter(sample == "SRR7939449" | sample == "SRR7939450")
fig(rRNA, "rRNA")

########
#ChIPseq(If you have input experiment.)
#######
CPM2 <- tx2 %>%
  gather(-anotation, -type, key = "sample", value = "cov") %>%
  left_join(CPM_ctrl, by = "sample")  %>%
  group_by(sample) %>% nest() %>%
  mutate(CPM = map(data, function(x) {x$cov * 1000000/x$cov_ctrl})) %>%
  unnest() %>% select(-cov, -cov_ctrl) %>%
  filter(type != ("dm6_full_genome_cov"))
fig2 <- function(x, y) {
  #this function requires: 1, annotation; 2, type; 3, sample; 4, CPM;
  #wilcoxon rank sum test
  datalabel <- x %>% group_by(type) %>%
    summarize(results = wilcox.exact(cov ~ sample, paired = T)$p.value,
              count = n()) %>%
    mutate(result2 = paste("italic(p) ==", signif(results, digits = 2)),
           count2 = paste("italic(n) ==", count),
           star = if_else(results > 1e-2, "NS",
                          if_else(results > 1e-4, "*",
                                  if_else(results > 2.2e-16, "**","***"))))
  #violin plot
  graph2 <- x %>% ggplot(aes(x = sample, y = cov, fill = type)) +
    geom_violin() +
    facet_wrap(~ type)
  graph2
  ggsave(paste0(filedir, y, "2.png"), width =12, height =9, dpi = 300)
  #boxplot
  graph <- x %>% ggplot(aes(x = sample, y = cov, fill = type)) +
    geom_boxplot() + facet_wrap(~ type) +
    geom_text(x = 1.1, y = 4, data = datalabel, aes(label = result2), parse = TRUE, hjust = 0) +
    geom_text(x = 1.1, y = 3, data = datalabel, aes(label = count2), parse = TRUE, hjust = 0) +
    geom_text(x = 1.3, y = 5, data = datalabel, aes(label = star), hjust = 0,size = 4.5) +
    ylab(expression(log[2](IP(CPM)/input(CPM) + 1))) +
    labs(title = paste0(y, "seq"))
  graph
  ggsave(paste0(filedir, y, ".png"), width =12, height =9, dpi = 300)
  #scatter plot
  scatter <- x %>% spread(key = sample, value = cov) 
  scat_x <- scatter %>% colnames() %>% .[1:2] %>% append(c("siEGFP", "siPiwi"), after = 2)
  colnames(scatter) <- scat_x
  COR <- scatter %>% select(-anotation) %>%
    group_by(type) %>% summarise(pcc = cor(siEGFP, siPiwi)) %>% 
    mutate(result2 = paste("italic(PCC) ==", signif(pcc, digits = 2)))
  DIFalpha <- scatter %>% 
    mutate(threshold = if_else(siEGFP > siPiwi, "down", "up")) %>% 
    group_by(type, threshold) %>% summarise(change = n()) %>% 
    mutate(result2 = if_else(threshold == "up", paste("up =", change),
                             paste("down =", change)))
  graph3 <- scatter %>% ggplot() +
    geom_pointdensity(aes(x = siEGFP, y = siPiwi), size = 0.3, alpha = 0.6) +
    scale_color_viridis() +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_abline(intercept = c(-log2(3), log2(3)), slope = 1, linetype="dashed") +
    geom_text(x = 4, y = 2, data = COR, aes(label = result2), parse = TRUE, hjust = 0,size = 2.5) +
    geom_text(x = 2, y = 0, data = DIFalpha %>% filter(threshold =="down"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 4, data = DIFalpha %>% filter(threshold =="up"), aes(label = result2), hjust = 0,size = 2.5) +
    facet_wrap(~ type)
  graph3 
  ggsave(paste0(filedir, y, "3.png"), width =12, height =9, dpi = 300)
  DIF <- scatter %>% filter(siEGFP !=0, siPiwi != 0) %>% 
    mutate(plus_thresold = log2(2 ^ (siEGFP + 1) - 1),
           minus_threshold = log2(2 ^ (siEGFP - 1) + 1/2),
           state = if_else(siPiwi > plus_thresold & (siPiwi > 1 | siEGFP >1), "plus", 
                           if_else(siPiwi <= minus_threshold & (siPiwi > 1 | siEGFP >1), "minus",
                                   if_else((siPiwi > 1 | siEGFP >1), "stay", 
                                           "unbound"))))
  totalDIF <- DIF %>% group_by(type) %>% summarise(countall = n())
  DIF2 <- DIF %>% group_by(type, state) %>% summarise(count = n()) %>% 
    left_join(totalDIF, key = "type") %>% 
    mutate(result2 = paste0(state,"=", count),
           result3 = paste0(signif(count*100 / countall, digits = 2), "%"))
  graph4 <- DIF %>% ggplot() +
    geom_point(aes(x = siEGFP, y = siPiwi, color = state), size = 0.3, alpha = 0.6) + 
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_abline(intercept = c(-log2(3), log2(3)), slope = 1, linetype="dashed") +
    geom_text(x = 0, y = 5.5, data = DIF2 %>% filter(state =="plus"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 5, data = DIF2 %>% filter(state =="plus"), aes(label = result3), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 4.5, data = DIF2 %>% filter(state =="unbound"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 0, y = 4, data = DIF2 %>% filter(state =="unbound"), aes(label = result3), hjust = 0,size = 2.5) +
    geom_text(x = 4, y = 2, data = DIF2 %>% filter(state =="stay"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 4, y = 1.5, data = DIF2 %>% filter(state =="stay"), aes(label = result3), hjust = 0,size = 2.5) +
    geom_text(x = 4, y = 1, data = DIF2 %>% filter(state =="minus"), aes(label = result2), hjust = 0,size = 2.5) +
    geom_text(x = 4, y = 0.5, data = DIF2 %>% filter(state =="minus"), aes(label = result3), hjust = 0,size = 2.5) +
    facet_wrap(~ type)
  graph4 
  ggsave(paste0(filedir, y, "4.png"), width =12, height =9, dpi = 300)
  write_state <- DIF %>% select(anotation, type, state) %>% 
    filter(state =="plus" | state == "minus") %>% 
    separate(anotation, into =  c("chr", "start", "end"), sep = "-")
  write_state %>% group_by(type, state) %>% 
    do(write_tsv(., paste0(filedir , "chara_region/", unique(.$type), "_", y, unique(.$state), ".bed"),
                 col_names = FALSE))
}
#for ChIP
#Bre ChIP
Bre <- CPM2 %>% filter(sample %in% c("SRR7992739", "SRR7992740","SRR7992744",
                                    "SRR7992745","SRR7992749", "SRR7992750")) %>%
  spread(key = sample, value = CPM) %>%
  filter(SRR7992739 != 0 & SRR7992740 != 0) %>%
  mutate(siE_H3K9me3 = log2(SRR7992744 / SRR7992739 + 1),
         siP_H3K9me3 = log2(SRR7992745 / SRR7992740 + 1),
         siE_Pol2 = log2(SRR7992749 / SRR7992739 + 1),
         siP_Pol2 = log2(SRR7992750 / SRR7992740 + 1)) %>%
  select(-contains("SRR")) %>%
  gather(-anotation, -type, key = "sample", value = "cov")
Bre_H3K9me3 <- Bre %>% filter(str_detect(sample, pattern = "H3K9me3"))
fig2(Bre_H3K9me3, "Bre_H3K9me3")
Bre_Pol2 <- Bre %>% filter(str_detect(sample, pattern = "Pol2"))
fig2(Bre_Pol2, "Bre_Pol2")
```

This code was generated on Python3.7 of Anaconda3, R3.60.
