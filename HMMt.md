
# DamIDseq analysis and broad peak detection using HMMt

DamID method was first developed by Bas Van Steensel (Van Steensel et
al., Nat Gen, 2001).

If you want to know the method, see (Vogel et al., Nature Method, 2007).

The analysis of DamID-chip was developed in (Filion et al., Cell, 2010).

And this example code uses HMMt package Filion used in the article.

In Van Steensel’s later article, they used chromHMM for DamIDseq, but
HMMt gives fine result.

### Unix code

Read filtering and Mapping.

Merge bam file is not used.

``` r
cd /Volumes/HDCZ-UT/mapping/Piwi_DamIDseq/;
mkdir Discas;
for X in SRR4296242 SRR4296243 SRR4296244 SRR4296245 ; do
seqkit grep -j 3 -sirp ^GATC ${X}.fastq -o Discas/${X}fill.fastq;
#DamIDseq read starts ^GATC because of sample procedure.
cutadapt -j 3 -a GATCCTCGGCCGCGACC -m 50 -q 20 -o ${X}_trimed.fastq Discas/${X}fill.fastq;
#PCR primer::GGTCGCGGCCGAGGATC
bowtie2 -p 6 -x ../../genome/UCSC_dm6/Bowtie2Index/genome -U ${X}_trimed.fastq -S ../Discas/${X}.sam;
grep -v "XS:" ../Discas/${X}.sam > ../Discas/${X}uniq.sam;
#only take unique reads into account.
samtools view -@ 6 -bS ../Discas/${X}uniq.sam > ../Discas/${X}.bam;
samtools sort -@ 6 -m 4G ../Discas/${X}.bam > ${X}_sort.bam;
samtools index ${X}_sort.bam;
done;
samtools merge -f Dam_merged.bam SRR4296242_sort.bam SRR4296243_sort.bam;
samtools merge -f Dam_Piwi_merged.bam SRR4296244_sort.bam SRR4296245_sort.bam;
for X in Dam Dam_Piwi; do
samtools index ${X}_merged.bam;
done;
for X in $(ls . | grep .bam$); do
tex2=`echo ${X} | sed 's/_sort.bam//'`;
bedtools bamtobed -i ${X} > ${tex2}_count.bed;
done;
```

### R code

First, filter mapping quality \>25, and get coverage on GATC site. Write
the result to bedgraph.

Test if the input Data fits the student’s t model.

Next, caluculate log2(DamID / input), and normalize it by loess method.
If you want to know thee loess method, see
<https://rafalab.github.io/dsbook/smoothing.html>. Write the normalized
result to bedgraph.

Last, use HMMt to detect binding site. Write the merged binding site to
bed file.

HMMt can get from gui11aume/HMMt <https://github.com/gui11aume/HMMt>.
Bas Van Steensel and his PhD student, Tom introduced me to the package.

``` r
#If HMMt isn't installed,
library(devtools)
install_github("gui11aume/HMMt")
library(HMMt)

#To know how HMMt works.
x <- c(rt(1000, df = 3), rt(1000, df = 3) + 1)
#t distribution (mean = 0) in 1 to 1000, and t distribution (mean = 1) in 1001 to 2000
plot(x, type ="l")
y <- BaumWelchT(x)
lines(BaumWelchT(x)$ViterbiPath - 1, col = 2)
#details and code :: github/gui11aume

#about HMM -> speech and lamguage processing (Dan Jurafsky and james H Martin)
#about HMM -> https://shop.ohmsha.co.jp/shopdetail/000000000574/
#about HMM of t-distribution -> (Filion et al., Cell, 2010)

###HMMt needs log2-normalized and loess-fitted data
###Procedure :: filter (^GATC), map, samtools bedcov(covert to bedgraph), load R, log2(Dam-Dam)
###Procedure2 :: loess-normalization, bridge for HMMt, finally HMMt
library(tidyverse)
library(data.table)
files3 <- dir("~/../Desktop/190702/", pattern = "\\_count.bed$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("^C:/Users/flust/Documents/../Desktop/190702/|_count.bed$", "")
df2 <- vector("list", length(files3))
DpnII <- read_tsv("~/../Desktop/190702/DpnIIrestsites.bed",
                  col_names = c("chr", "start", "end")) %>%
  mutate(center = start +2) %>% dplyr::select(chr, center)
for (i in seq_along(files3)) {
  df2[[i]] <- fread(files3[[i]],
                    col.names = c("chr", "start", "end", "dis", "quality", "strand"),
                    sep = ("\t")) %>%
    as_tibble() %>%
    filter(quality >= 25) %>%
    filter(chr !="chrY") %>%
    mutate(center = if_else(strand =="+", start +2, end - 2)) %>%
    group_by(chr, center) %>% summarize(coverage = n())
  df2[[i]][["type"]] <- files4[[i]]
}
#bedtools bamtobed result bed contains mapping quality in 5th col.
#http://www.acgt.me/blog/2014/12/16/understanding-mapq-scores-in-sam-files-does-37-42
Coverage1 <- DpnII %>%
  left_join(df2 %>%
              bind_rows() %>%
              spread(key = type, value = coverage),
            by = c("chr", "center")) %>%
  filter(!(SRR4296242 =="NA" & SRR4296243 =="NA" & SRR4296244 =="NA" & SRR4296245 =="NA"))
#First, gather GATC site using left_join with GATC restsites.
Coverage1 %>% dplyr::select(-SRR4296244,-SRR4296245) %>%
  drop_na() %>%
  mutate(Dam = SRR4296242 + SRR4296243,
         start = center -2,
         end = center +2) %>%
  dplyr::select(chr, start, end, Dam) %>%
  write_tsv("~/../Desktop/190702/Dam.bedgraph",
            col_names = FALSE)
Coverage1 %>% dplyr::select(-SRR4296242,-SRR4296243) %>%
  drop_na() %>%
  mutate(Dam = SRR4296244 + SRR4296245,
         start = center -2,
         end = center +2) %>%
  dplyr::select(chr, start, end, Dam) %>%
  write_tsv("~/../Desktop/190702/Dam_Piwi.bedgraph",
            col_names = FALSE)
############################################
#Confirm if student's t estimation can fit to the distribution.
test_t <- Coverage1 %>% dplyr::select(-SRR4296244,-SRR4296245) %>% filter(chr !="chrY") %>%
  filter(!(SRR4296242 =="NA"|SRR4296243 =="NA")) %>%
  mutate(dis_cov = log2(SRR4296242/ SRR4296243)) %>%
  dplyr::select(chr, center, dis_cov) %>%
  group_by(chr) %>% nest() %>%
  mutate(oberon = map(data, function(x) {stats::loess(x$dis_cov ~ x$center, span = 0.4)}),
         titania = map2(data, oberon, function(x,y) {stats::predict(y, x$center)})) %>%
  dplyr::select(-oberon) %>% unnest() %>%
  mutate(loess_nom = dis_cov - titania) %>%
  dplyr::select(chr, center, loess_nom)
library(MASS)
fit_norm <- fitdistr(test_t$loess_nom, "normal")
fit_t <- fitdistr(test_t$loess_nom, "t",
                  start = list(m=mean(test_t$loess_nom),
                               s=sd(test_t$loess_nom),
                               df=3),
                  lower=c(-1, 0.001,1))
test_t %>%
  ggplot() + geom_density(aes(x = loess_nom)) +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  stat_function(fun = dnorm, args = list(mean = fit_norm[["estimate"]][["mean"]],
                                         sd = fit_norm[["estimate"]][["sd"]])) +
  stat_function(fun = dnorm, args = list(mean = fit_t[["estimate"]][["m"]],
                                      sd = fit_t[["estimate"]][["s"]]))
ggsave("~/../Desktop/190702/t_distri.png", width =3, height = 4, dpi = 300)
#You can get the figure S5A in (Filion et al., Cell, 2010).
Coverage1 %>% filter(chr !="chrY") %>%
  filter(!(SRR4296242 =="NA"|SRR4296243 =="NA" | SRR4296244 =="NA"|SRR4296245 =="NA")) %>%
  mutate(dis_cov1 = log2(SRR4296245/ SRR4296243),
         dis_cov2 = log2(SRR4296244/ SRR4296242),
         dis_cov = (dis_cov1 + dis_cov2)/2) %>%
  dplyr::select(chr, center, dis_cov) %>%
  filter(chr =="chr2L") %>%
  ggplot(aes(x = center, y = dis_cov)) + geom_line()+
  geom_smooth(method = "loess", span = 0.4, se = FALSE)
ggsave("~/../Desktop/190702/distri_chr2L.png", width =10, height = 4, dpi = 300)
###################################################
#replicate merge and LOESS normalization
Merge_LOESS <- Coverage1 %>% filter(chr !="chrY") %>%
  filter(!(SRR4296242 =="NA"|SRR4296243 =="NA" | SRR4296244 =="NA"|SRR4296245 =="NA")) %>%
  mutate(dis_cov1 = log2(SRR4296245/ SRR4296243),
         dis_cov2 = log2(SRR4296244/ SRR4296242),
         dis_cov = (dis_cov1 + dis_cov2)/2) %>%
  dplyr::select(chr, center, dis_cov) %>%
  group_by(chr) %>% nest() %>%
  mutate(oberon = map(data, function(x) {stats::loess(x$dis_cov ~ x$center, span = 0.4)}),
         titania = map2(data, oberon, function(x,y) {stats::predict(y, x$center)})) %>%
  dplyr::select(-oberon) %>% unnest() %>%
  mutate(loess_nom = dis_cov - titania,
         start = center -2,
         end = center +2) %>%
  dplyr::select(chr, start, end, loess_nom)
#distribution after normalization
Merge_LOESS %>% mutate(center = start + 2) %>%
  filter(chr =="chr2L") %>%
  ggplot(aes(x = center, y = loess_nom)) + geom_line()
ggsave("~/../Desktop/190702/distri_chr2L_after.png", width =10, height = 4, dpi = 300)
###################################################

loess_nom <- HMMt::bridge(Merge_LOESS %>% as.data.frame())
HMM_baum <- HMMt::BaumWelchT(loess_nom[["x"]])
#HMMt bridge uniformize sampling rate

Result <- HMM_baum$ViterbiPath %>%
  enframe(value = "state") %>%
  mutate(cov_norm = loess_nom$x,
         visible = loess_nom$nonvirtuals) %>%
  drop_na() %>%
  bind_cols(Merge_LOESS) %>%
  rownames_to_column(var = "order_num") %>%
  dplyr::select(chr, start,end, order_num, state,  loess_nom)
#sum the result
Result %>% filter(state ==2) %>%
  dplyr::select(chr, start,end, loess_nom) %>%
  write_tsv("~/../Desktop/190702/Piwi_bound_site_binless.bedgraph",
            col_names = FALSE)
Result %>%
  dplyr::select(chr, start,end, loess_nom) %>%
  write_tsv("~/../Desktop/190702/Piwi_coverage.bedgraph",
            col_names = FALSE)
Jurius <- Result %>% filter(state ==2) %>% mutate(order_num2 = as.double(order_num)) %>%
  dplyr::select(-state, -order_num) %>%
  group_by(chr) %>% nest() %>%
  mutate(pack = map(data, function(x) {which(diff(x$order_num2) != 1)}),
         end_merge = map2(data, pack, function(x,y) {x$end[y]}),
         demetrius = map(pack, function(x) {append(x+1,1,after =0) %>% head(-1)}),
         start_merge = map2(data, demetrius, function(x,y) {x$start[y]})) %>%
  dplyr::select(chr, start_merge, end_merge) %>%
  unnest() %>%
  write_tsv("~/../Desktop/190702/Piwi_bound_site.bed",
            col_names = FALSE)
```

# Sample code of Piwi CLIP

HMMt can develop other experiments if the sample distribution can be
estimated by student’s t distribution.

Here, an example of CLIP experiment.

Also, there are other experiment tools that detect broad peak of ChIP or
CLIP by Hidden Markov Model.

See chromHMM or HitoneHMM.

chromHMM : (Ernest and Kellis, Nature Method, 2012)

HistoneHMM : (Heinig et al., BMC, 2015)

### Unix code

``` r
mkdir Discas/;
for X in PiwiCLIP_IN PiwiCLIP_IP_rep2 PiwiCLIP_IP_rep1; do
cutadapt -j 3 -m 27 --max-n 0 -a TGGAATTCTCGGGTGCCAAGG -o Discas/${X}_trimed.fastq.gz ${X}.fastq.gz;
seqkit rmdup -s Discas/${X}_trimed.fastq.gz -o Discas/${X}_trimed2.fastq.gz;
cutadapt -j 3 -u 9 -m 18 -q 20 -o ${X}_f.fastq.gz Discas/${X}_trimed2.fastq.gz;
bowtie2 -p 6 -N 1 -x ../../genome/UCSC_dm6/Bowtie2Index/genome -U ${X}_f.fastq.gz -S Discas/${X}.sam;
samtools view -@ 6 -bS Discas/${X}.sam > Discas/${X}.bam;
samtools sort -@ 6 -m 4G Discas/${X}.bam > ${X}_sort.bam;
samtools index ${X}_sort.bam;
done;
for X in $(ls . | grep .bam$); do
tex2=`echo ${X} | sed 's/_sort.bam//'`;
bedtools bamtobed -i ${X} > ${tex2}_count.bed;
done;
```

### R code

``` r
genome <- read_tsv("~/../Desktop/genome/UCSC_dm6/dm6.genome",
                   col_names = c("chr", "size"))

files3 <- dir("~/../Desktop/190709/genome", pattern = "\\_count.bed$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("^C:/Users/flust/Documents/../Desktop/190709/genome/|_count.bed$", "")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[i]] <- read_tsv(files3[[i]],
                    col_names = c("chr", "start", "end", "dis", "quality", "strand")) %>%
    filter(quality >= 15) %>%
    filter(chr !="chrY")
  df2[[i]][["type"]] <- files4[[i]]
}
Coverage1 <- df2 %>% bind_rows() %>%
  mutate(length = end - start,
         center = (start + end) %/% 2) %>%
  dplyr::select(-dis) %>% filter(type != "PiwiCLIP_IP_rep2")
chck1 <- Coverage1 %>% group_by(length) %>% summarise(n())

Resolution_make <- function(res){
  ref_res <- genome %>% group_by(chr) %>% nest() %>%
    mutate(binres = map(data, function(x) {seq(from=0, to=(x$size %/% res)*res, by = res)}),
           binres2 = map2(binres,data, function(x,y) {append(tail(x, -1), y$size)})) %>%
    dplyr::select(-data) %>%
    unnest()
  ref_ano <- Coverage1 %>% mutate(binres = (center %/% res) * res) %>%
    dplyr::select(chr, binres, type) %>% group_by(chr, binres, type) %>%
    summarise(coverage = n()) %>% spread(key = type, value = coverage) %>%
    rename(ctrl = PiwiCLIP_IN, Piwi = PiwiCLIP_IP_rep1)
  Cov_bin <- left_join(ref_res, ref_ano, by = c("chr", "binres")) %>%
    filter(chr !="chrY")
  Cov_bin
}
result <- Resolution_make(1000)
#1kb resolution
#chop genome in 1000 length fragment.
test_make <- function(res){
  ref_res <- genome %>% group_by(chr) %>% nest() %>%
    mutate(binres = map(data, function(x) {seq(from=0, to=(x$size %/% res)*res, by = res)}),
           binres2 = map2(binres,data, function(x,y) {append(tail(x, -1), y$size)})) %>%
    dplyr::select(-data) %>% unnest()
  ref_ano1 <- df2 %>% bind_rows() %>%
    mutate(center = (start + end) %/% 2,
           binres = (center %/% res) * res) %>%
    filter(type == "PiwiCLIP_IN")
  ref_ano2 <- ref_ano1$dis %>% sample(150000, replace = FALSE) %>% enframe(value = "dis")
  ref_ano <- ref_ano2 %>% left_join(ref_ano1, by ="dis") %>%
    dplyr::select(chr, binres, type) %>% group_by(chr, binres, type) %>%
    summarise(coverage = n()) %>% spread(key = type, value = coverage) %>%
    rename(ctrl = PiwiCLIP_IN)
  Cov_bin <- left_join(ref_res, ref_ano, by = c("chr", "binres")) %>%
    filter(chr !="chrY")
  Cov_bin
}
result_t1 <- test_make(1000)
result_t1 %>% group_by(ctrl) %>% summarise(n())
result_t2 <- test_make(1000) %>% rename(ctrl1 = ctrl)
t_test <- left_join(result_t1, result_t2, by = c("chr", "binres")) %>%
  drop_na() %>% mutate(dis_cov = log2(ctrl1/ctrl),
                                      center = binres + 500) %>%
  dplyr::select(chr, center, dis_cov) %>%
  group_by(chr) %>% nest() %>%
  mutate(oberon = map(data, function(x) {stats::loess(x$dis_cov ~ x$center, span = 0.4)}),
         titania = map2(data, oberon, function(x,y) {stats::predict(y, x$center)})) %>%
  dplyr::select(-oberon) %>% unnest() %>%
  mutate(loess_nom = dis_cov - titania,
         start = center -500,
         end = center +500) %>%
  dplyr::select(chr, start, end, loess_nom)
fit_norm <- fitdistr(t_test$loess_nom, "normal")
fit_t <- fitdistr(t_test$loess_nom, "t",
                  start = list(m=mean(t_test$loess_nom),
                               s=sd(t_test$loess_nom),
                               df=3),
                  lower=c(-1, 0.001,1))
t_test %>%
  ggplot() + geom_density(aes(x = loess_nom)) +
  coord_cartesian(xlim = c(-0.1, 0.1)) +
  stat_function(fun = dnorm, args = list(mean = fit_norm[["estimate"]][["mean"]],
                                         sd = fit_norm[["estimate"]][["sd"]])) +
  stat_function(fun = dnorm, args = list(mean = fit_t[["estimate"]][["m"]],
                                         sd = fit_t[["estimate"]][["s"]]))
ggsave("~/../Desktop/190709/t_distri.png", width =3, height = 4, dpi = 300)

result %>% filter(chr !="chrY") %>%
  drop_na %>%
  mutate(dis_cov = log2(Piwi/ctrl),
         center = binres + 500) %>%
  dplyr::select(chr, center, dis_cov) %>%
  filter(chr =="chr2L") %>%
  ggplot(aes(x = center, y = dis_cov)) + geom_line()+
  geom_smooth(method = "loess", span = 0.4, se = FALSE)
ggsave("~/../Desktop/190709/distri_chr2L.png", width =10, height = 4, dpi = 300)
########################################################

Merge_LOESS <- result %>%
  drop_na() %>%
  mutate(dis_cov = log2(Piwi/ctrl),
         center = binres + 500) %>%
  dplyr::select(chr, center, dis_cov) %>%
  group_by(chr) %>% nest() %>%
  mutate(oberon = map(data, function(x) {stats::loess(x$dis_cov ~ x$center, span = 0.4)}),
         titania = map2(data, oberon, function(x,y) {stats::predict(y, x$center)})) %>%
  dplyr::select(-oberon) %>% unnest() %>%
  mutate(loess_nom = dis_cov - titania,
         start = center -500,
         end = center +500) %>%
  dplyr::select(chr, start, end, loess_nom)
loess_nom <- HMMt::bridge(Merge_LOESS %>% as.data.frame())
HMM_baum <- HMMt::BaumWelchT(loess_nom[["x"]])

Result <- HMM_baum$ViterbiPath %>%
  enframe(value = "state") %>%
  mutate(cov_norm = loess_nom$x,
         visible = loess_nom$nonvirtuals) %>%
  drop_na() %>%
  bind_cols(Merge_LOESS) %>%
  rownames_to_column(var = "order_num") %>%
  dplyr::select(chr, start,end, order_num, state,  loess_nom)

Result %>% filter(state ==2) %>%
  dplyr::select(chr, start,end, loess_nom) %>%
  write_tsv("~/../Desktop/190709/Piwi_bound_site_CLIP.bedgraph",
            col_names = FALSE)
Result %>%
  dplyr::select(chr, start,end, loess_nom) %>%
  write_tsv("~/../Desktop/190709/Piwi_CLIP_coverage.bedgraph",
            col_names = FALSE)
Jurius <- Result %>% filter(state ==2) %>% mutate(order_num2 = as.double(order_num)) %>%
  dplyr::select(-state, -order_num) %>%
  group_by(chr) %>% nest() %>%
  mutate(pack = map(data, function(x) {which(diff(x$order_num2) != 1)}),
         end_merge = map2(data, pack, function(x,y) {x$end[y]}),
         demetrius = map(pack, function(x) {append(x+1,1,after =0) %>% head(-1)}),
         start_merge = map2(data, demetrius, function(x,y) {x$start[y]})) %>%
  dplyr::select(chr, start_merge, end_merge) %>%
  unnest() %>%
  write_tsv("~/../Desktop/190709/Piwi_bound_site_CLIP.bed",
            col_names = FALSE)
```

All R was done in R3.6 on windows10. If you use MacOS, you may have some
trouble in installing HMMt. See
<https://www.rstudio.com/products/rpackages/devtools/>
