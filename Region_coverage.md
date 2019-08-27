
# Get region coverage and, plot scatter plot and box plot.

We take many time for analysis of Bed files. R helps take less time for
analysis.

In this example, we analyse multiple Bam files on multiple BED region.

We get coverage on 9 bed region from 33 bam files.

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
#Using Regular expression, you get short name of the file.
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

I refer the method to
<http://www.sciencemag.org/cgi/pmidlookup?view=long&pmid=23328393>.

I get the data from <https://www.nature.com/articles/nature13994>.

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

``` r
cd /Volumes/HDCZ-UT/chipmeta/;
for X in SRR3503092_extended dm6_full_genome TEs-DHS-dm6_extended hkCP_extended DHS_extended dCP_extended \
dCPclose_extended hkCPclose_extended dCPonly_extended hkCPonly_extended; do
samtools bedcov -Q 15 ${X}.bed Bamfile/*/*.bam > Coverage/${X}_cov.bed;
done;
#samtools bedcov is not good choice, because there is bedtools coverage.
#But when you use *, bedtools doesn't work.

#10000bp length file.
cd /Volumes/HDCZ-UT/chipmeta/;
for X in SRR3503092_extended2 dm6_full_genome TEs-DHS-dm6_extended2 hkCP_extended2 DHS_extended2 \
dCP_extended2 dCPclose_extended2 hkCPclose_extended2 dCPonly_extended2 hkCPonly_extended2; do
samtools bedcov -Q 15 ${X}.bed Bamfile/*/*.bam > Coverage/10000bp/${X}_cov.bed;
done;

ls Bamfile/*/*.bam
Bamfile/ATACseq/siEGFP_ATAC_merged.bam
Bamfile/ATACseq/siPiwi_ATAC_merged.bam
Bamfile/ChIPseq/SRR7992739_sort.bam
Bamfile/ChIPseq/SRR7992740_sort.bam
Bamfile/ChIPseq/SRR7992744_sort.bam
Bamfile/ChIPseq/SRR7992745_sort.bam
Bamfile/ChIPseq/SRR7992749_sort.bam
Bamfile/ChIPseq/SRR7992750_sort.bam
Bamfile/ChIPseq/siE-H3K27Ac_sort.bam
Bamfile/ChIPseq/siE-H3K27me3_sort.bam
Bamfile/ChIPseq/siE-H3K4me3_sort.bam
Bamfile/ChIPseq/siE-IN_sort.bam
Bamfile/ChIPseq/siE_CP190_bowtie2_sort.bam
Bamfile/ChIPseq/siE_IN_bowtie2_sort.bam
Bamfile/ChIPseq/siP-H3K27Ac_sort.bam
Bamfile/ChIPseq/siP-H3K27me3_sort.bam
Bamfile/ChIPseq/siP-H3K4me3_sort.bam
Bamfile/ChIPseq/siP-IN_sort.bam
Bamfile/ChIPseq/siP_CP190_bowtie2_sort.bam
Bamfile/ChIPseq/siP_IN_bowtie2_sort.bam
Bamfile/DHSseq/OSC_DHS.bam
Bamfile/DHSseq/SRR569912_bowtie2_sort.bam
Bamfile/GROseq/SRR609665_sort.bam
Bamfile/GROseq/SRR609666_sort.bam
Bamfile/RIPseq/SRR3503092_bowtie2_sort.bam
Bamfile/RIPseq/SRR8791716_sort.bam
Bamfile/RNAseq/EGFPKD_Bre_poly_sort.bam
Bamfile/RNAseq/PiwiKD_Bre_poly_sort.bam
Bamfile/RNAseq/SRR7939449_sort.bam
Bamfile/RNAseq/SRR7939450_sort.bam
Bamfile/STARRseq/SRR1297299rmdup.bam
Bamfile/STARRseq/dCP.bam
Bamfile/STARRseq/hkCP.bam
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
#you can change wilcoxon rank sum test to student t-test, wilcoxon signed-rank test.
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
  ggsave(paste0(filedir, y, "2.png"), width =8, height =6, dpi = 300)
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
  ggsave(paste0(filedir, y, ".png"), width =8, height =6, dpi = 300)
  #scatter plot
  scatter <- x %>% spread(key = sample, value = CPM)
  scat_x <- scatter %>% colnames() %>% .[1:2] %>% append(c("siEGFP", "siPiwi"), after = 2)
  colnames(scatter) <- scat_x
  graph3 <- scatter %>% ggplot() +
    geom_point(aes(x = siEGFP, y = siPiwi, color = type), size = 0.3, alpha = 0.6) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_abline(intercept = c(-log2(3), log2(3)), slope = 1, linetype="dashed") +
    facet_wrap(~ type)
  graph3
  ggsave(paste0(filedir, y, "3.png"), width =8, height =6, dpi = 300)
}
# * = p < 1e-2
# ** = p < 1e-4
# *** = p < 2.2e-16
# else NS (threashold was refered to Hug et al., Cell, 2017.(<https://www.cell.com/cell/pdf/S0092-8674(17)30343-4.pdf>))
#In order to get figure by short code, I use function.
#If you want other method of display, see https://github.com/kassambara/ggpubr/issues/65.
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
  ggsave(paste0(filedir, y, "2.png"), width =8, height =6, dpi = 300)
  #boxplot
  graph <- x %>% ggplot(aes(x = sample, y = cov, fill = type)) +
    geom_boxplot() + facet_wrap(~ type) +
    geom_text(x = 1.1, y = 4, data = datalabel, aes(label = result2), parse = TRUE, hjust = 0) +
    geom_text(x = 1.1, y = 3, data = datalabel, aes(label = count2), parse = TRUE, hjust = 0) +
    geom_text(x = 1.3, y = 5, data = datalabel, aes(label = star), hjust = 0,size = 4.5) +
    ylab(expression(log[2](IP(CPM)/input(CPM) + 1))) +
    labs(title = paste0(y, "seq"))
  graph
  ggsave(paste0(filedir, y, ".png"), width =8, height =6, dpi = 300)
  #scatter plot
  scatter <- x %>% spread(key = sample, value = cov)
  scat_x <- scatter %>% colnames() %>% .[1:2] %>% append(c("siEGFP", "siPiwi"), after = 2)
  colnames(scatter) <- scat_x
  graph3 <- scatter %>% ggplot() +
    geom_point(aes(x = siEGFP, y = siPiwi, color = type), size = 0.3, alpha = 0.6) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_abline(intercept = c(-log2(3), log2(3)), slope = 1, linetype="dashed") +
    facet_wrap(~ type)
  graph3
  ggsave(paste0(filedir, y, "3.png"), width =8, height =6, dpi = 300)
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
#keio chip seq
KEIO <- CPM2 %>% filter(str_detect(sample, pattern = "si")) %>%
  filter(!(str_detect(sample, pattern = "ATAC"))) %>%
  spread(key = sample, value = CPM) %>%
  filter(`siE-IN` != 0 & `siP-IN` != 0) %>%
  mutate(keio_E_H3K4me3 = log2(`siE-H3K4me3` / `siE-IN` + 1),
         keio_P_H3K4me3 = log2(`siP-H3K4me3` / `siP-IN` + 1),
         keio_E_H3K27me3 = log2(`siE-H3K27me3` / `siE-IN` + 1),
         keio_P_H3K27me3 = log2(`siP-H3K27me3` / `siP-IN` + 1),
         keio_E_H3K27Ac = log2(`siE-H3K27Ac` / `siE-IN` + 1),
         keio_P_H3K27Ac = log2(`siP-H3K27Ac` / `siP-IN` + 1)) %>%
  select(-contains("si")) %>%
  gather(-anotation, -type, key = "sample", value = "cov")
keio_H3K4me3 <- KEIO %>% filter(str_detect(sample, pattern = "H3K4me3"))
fig2(keio_H3K4me3, "keio_H3K4me3")
keio_H3K27me3 <- KEIO %>% filter(str_detect(sample, pattern = "H3K27me3"))
fig2(keio_H3K27me3, "keio_H3K27me3")
keio_H3K27Ac <- KEIO %>% filter(str_detect(sample, pattern = "H3K27Ac"))
fig2(keio_H3K27Ac, "keio_H3K27Ac")
```

This code was generated on Python3.7, R3.60.