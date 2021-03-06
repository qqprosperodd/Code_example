
# Mappable region bed file

To make mapping quality bedgraph, we chop genome into short length and
map those sequences to the genome and count the fraction of mapped reads
in each binsize.

Mappable region bed file is also useful in Hi-C filtering.

This code procedure refers to
<https://bitbucket.org/tanaylab/schic2/src/default/map3c/>.

The article of schic2 is <https://www.nature.com/articles/nature23001>
(Nagano et al., Nature, 2017).

### R code

You chop genome into short length fragment.

``` r
library(tidyverse)
genome <- read_tsv("/Volumes/HDCZ-UT/genome/UCSC_dm6/dm6.genome", col_names = c("chr", "size"))
filedir <- "/Volumes/HDCZ-UT/genome/mapping_quality/"

#50bp, 100bp, 150bp length bed file make
fastq_make <- function(res, length){
  #res = shift size, length = fragment length
  delete <- length / res
  ref_res <- genome %>%
    group_by(chr) %>% nest() %>%
    mutate(binres = map(data, function(x) {seq(from=0, to=(x$size %/% res - delete)*res, by = res)}),
           binres2 = map(binres, function(x) {x + length})) %>%
    select(-data) %>%
    unnest()
  ref_res %>%
    write_tsv(paste0(filedir, "dm6_genome_", res, "bin_", length, "len.bed"), 
              col_names = FALSE)
}
fastq_make(10, 150)
fastq_make(10, 100)
fastq_make(10, 50)
```

To generate rest-sites bed with enogth mapping quality, I deside to get
coverage of rest-site.

For bedtools coverage, you get 50-lenth rest-site fragment.

``` r
filedir <- "/Volumes/HDCZ-UT/genome/mapping_quality/rest_bed/"
setwd(filedir)
library(tidyverse)
DpnII <- read_tsv("/Volumes/HDCZ-UT/genome/UCSC_dm6/DpnIIrestsites.bed",
                  col_names = c("chr", "start", "end"))
HindIII <- read_tsv("/Volumes/HDCZ-UT/genome/UCSC_dm6/HindIIIrestsites.bed",
                  col_names = c("chr", "start", "end"))
DpnII %>% mutate(start2 = start-48, end2 = end + 48) %>%
  select(chr, start2, end2) %>%
  filter(start2 > 0) %>%
  write_tsv("DpnIIrestsites_100.bed")
HindIII %>% mutate(start2 = start-47, end2 = end + 47) %>%
  select(chr, start2, end2) %>%
  write_tsv("HindIIIrestsites_100.bed")
```

### Unix code

First, you get fasta sequence of short length fragment.

Next, map them to the genome.

You get bedgraph coverage of full genome and coverage of 100
length-extended rest-site.

If you haven’t get rest-site bed file yet, please create it by
HiCExplorer
findRestSites.

<https://hicexplorer.readthedocs.io/en/latest/content/tools/findRestSite.html#findrestsite>

``` r
#Mapping
cd /Volumes/HDCZ-UT/genome/mapping_quality/;
mkdir Discas/;
for X in dm6_genome_10bin_150len dm6_genome_10bin_100len dm6_genome_10bin_50len; do
bedtools getfasta -fi /Volumes/HDCZ-UT/genome/UCSC_dm6/dm6_chromosome.fa -bed ${X}.bed > ${X}.fa;
bwa mem -A 1 -B 4 -E 50 -L 0 -t 4 /Volumes/HDCZ-UT/genome/UCSC_dm6/bwaIndex/bwa_dm6_genome \
 ${X}.fa | samtools view -@ 6 -bS - > Discas/bwa_${X}.bam;
samtools sort -@ 6 -m 4G Discas/bwa_${X}.bam > bwa_${X}_sort.bam;
samtools index bwa_${X}_sort.bam;
bowtie2 -p 4 -f -x /Volumes/HDCZ-UT/genome/UCSC_dm6/Bowtie2Index/genome \
 -U ${X}.fa | samtools view -@ 6 -bS - > Discas/bowtie2_${X}.bam;
samtools sort -@ 6 -m 4G Discas/bowtie2_${X}.bam > bowtie2_${X}_sort.bam;
samtools index bowtie2_${X}_sort.bam;
done;
#Max mapping quality of bowtie2 is 42 and Also of bwa is 60.

#Generate mapping quality bedgraph
mkdir Coverage/;
for Y in bowtie2 bwa ;do
for X in dm6_genome_10bin_150len dm6_genome_10bin_100len dm6_genome_10bin_50len; do
bamCoverage -b ${Y}_${X}_sort.bam -o Coverage/${Y}_${X}.bedgraph \
-of bedgraph --binSize=10 --minMappingQuality 15 -p 6 ;
done;
done;

#Get coverage of 100 length-etended rest-site with mapping quality > 15.
cd /Volumes/HDCZ-UT/genome/mapping_quality/;
for Y in bowtie2 bwa ;do
for X in HindIII DpnII ; do
for Z in dm6_genome_10bin_150len ; do
samtools view -@ 6 -q 15 -bS ${Y}_${Z}_sort.bam | bedtools coverage -b stdin \
 -a rest_bed/${X}restsites_100.bed > Coverage/${X}_${Y}_${Z}_bedcov.bed;
done;
done;
done;
```

### R code

Filter Coverage \>= 20 and you get final rest-site.

``` r
filedir <- "/Volumes/HDCZ-UT/genome/mapping_quality/Coverage/"
setwd(filedir)
files <- dir(filedir, pattern = "\\.bed$", full.names = TRUE)
files2 <- files %>%
  str_replace_all("^/Volumes/HDCZ-UT/genome/mapping_quality/Coverage//|.bed$", "") %>%
  str_replace("_bedcov", "")

df <- vector("list", length(files))
for (i in seq_along(files)) {
  df[[i]] <- read_tsv(files[[i]], col_names= c("chr", "start", "end", "cov"))
  df[[i]][["type"]] <- files2[[i]]
}
tx2 <- bind_rows(df) %>%
  filter(cov >= 20) %>%
  mutate(start2 = if_else(str_detect(type, pattern = "Hind") , start+47, start+48),
         end2 = if_else(str_detect(type, pattern = "Hind") , end-47, end-48))

cov_e <- bind_rows(df) %>% spread(type, cov)
Den <- cov_e %>%
  mutate(order1 = seq(from=1, to=length(cov_e$start),by=1)) %>%
  gather(-chr, -start, -end, -order1, key = "type", value ="coverage") %>%
  group_by(type) %>% nest() %>%
  mutate(threshold = map(data, function(x) {if_else(x$coverage >= 20, "valid", "discard")})) %>% 
  unnest() %>%
  filter(threshold =="discard") %>% 
  select(type, order1) %>%
  group_by(type) %>% nest()
Den4 <- flatten(Den$data)
names(Den4) <- Den$type
library(UpSetR)
ppi <- 300
png("upset.png", width = 10*ppi, height = 6*ppi, res = ppi)
UpSetR::fromList(Den4) %>%
  UpSetR::upset(nsets = length(files),sets.bar.color = "#56B4E9", order.by="freq", nintersects = NA)
dev.off()

Den5 <- tx2 %>% select(type, chr ,start2, end2) %>% 
  group_by(type) %>% nest()
Den6 <- map(1:length(Den5$type), function(x) data.frame(Den5$data[[x]]))
names(Den6) <- Den5$type
map(1:length(Den6), 
    function(x) write_tsv(Den6[[x]],
                          paste0("/Volumes/HDCZ-UT/genome/mapping_quality/filterbed/", names(Den6[x]), ".bed"),
                                          col_names = FALSE))
```

This code was generated on Python3.7, R3.60.
