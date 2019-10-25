
setwd("/media/meg/Titania/Piwi_CLASH/")
filedir <- "/media/meg/Titania/Piwi_CLASH/"
library(tidyverse)
library(gplots)
library(RColorBrewer)
library(gtools)
library(gridExtra)
piRNA_length <- read_tsv("~/hyb/data/fastq/piRNA_all.fasta.fai", 
                         col_names = c("piRNA", "pi_length")) %>% 
  separate(piRNA, sep = "\\_", into = c("gene_name"))
Vienna <- read_csv("Piwi_CLASH_comp_dm6_25_2_hybrids_ua.vienna", 
                   col_names = "data") %>% 
  mutate(index = rep(1:(length(.$data) / 3), each = 3),
         type = if_else(row_number() %% 3 == 1,"name",
                        if_else(row_number() %% 3 == 2, "seq", "score"))) %>% 
  spread(key = "type", value = "data") 
#Vienna fileは最初の行にリードの名前、２行目にリードの配列、３行目に結合箇所.

Vienna2 <- Vienna %>% 
  mutate(name = str_replace_all(name, "transposable_element", "TE"),
         name = str_replace(name, "pre_miRNA", "premiRNA")) %>% 
  separate(score, sep = "\\t", into= c("contact", "deltaG")) %>% 
  separate(name, sep = "\\_", into= c("dis11", "count", "gene_name", "dis", "dis2", "type2", "start", "end", "dis4", 
                                      "target", "dis5", "dis6", "type3", "start2", "end2")) %>% 
  separate(end, sep = "-", into = c("end", "dis10")) %>% 
  select(-contains("dis")) %>% 
  mutate(deltaG = str_replace(deltaG, "\\(", "") %>% str_replace("\\)", "") %>% as.double(),
         count = count %>% as.integer(),
         start = start %>% as.integer(),
         end = end %>% as.integer())
Vienna2 %>% group_by(type2, type3) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(y =type3, x=type2))+
  geom_tile(aes(fill = count))+
  geom_text(aes(label = count))+
  scale_fill_gradient(low = "white", high = "green")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  xlab("3' RNA type") + ylab("5' RNA type")
ggsave(paste0(filedir,"vienna/", "RNAtype_test_collapse.png"), width =7, height =6, dpi = 300)
Vienna2 %>% group_by(type2, type3) %>% 
  summarise(count = sum(count)) %>% 
  ggplot(aes(y =type3, x=type2))+
  geom_tile(aes(fill = count))+
  geom_text(aes(label = count))+
  scale_fill_gradient(low = "white", high = "green")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  xlab("3' RNA type") + ylab("5' RNA type")
ggsave(paste0(filedir,"vienna/", "RNAtype_test_total.png"), width =7, height =6, dpi = 300)
Vienna3 <- Vienna2 %>% filter(type2 == "piRNA") %>% 
  filter(start==1) %>% 
  mutate(contact2 = str_sub(contact, 1, (end - start + 1)),
         length = nchar(seq),
         piRNAlen = end-start+1)%>% 
  filter(str_detect(contact2, "\\(")) %>% 
  left_join(piRNA_length, by = "gene_name") %>% 
  filter(pi_length == piRNAlen)
Vienna3 %>% mutate(deltaG2 = if_else(deltaG %% 1 >=0.5, deltaG %/% 1+1, deltaG %/% 1)) %>% 
  group_by(deltaG2) %>% summarise(total = sum(count)) %>% 
  ggplot(aes(x = deltaG2, y = total)) + geom_line(color= "blue")+
  xlab("dG (kcal/mol)") + ylab("Number of chimeric reads")
ggsave(paste0(filedir,"vienna/", "distri.png"), width =7, height =6, dpi = 300)
Vienna3 %>% ggplot(aes(x = deltaG)) + geom_freqpoly(binwidth=1)
ggsave(paste0(filedir,"vienna/", "distri.png"), width =7, height =6, dpi = 300)
Vienna3 %>% group_by(type3) %>% 
  summarise(count_collapse = n()) %>% 
  ggplot(aes(x = fct_reorder(type3, count_collapse) %>% fct_rev(), y = count_collapse)) +
  geom_col(fill = "orange")+
  xlab("RNA type") + ylab("count of collapsed reads") +
  geom_text(aes(label = count_collapse), vjust = -0.5, hjust=0.5, colour = "black",size = 3)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
ggsave(paste0(filedir,"vienna/", "piRNA_target_RNA_type_collapse.png"), width =5, height =4, dpi = 300)
Vienna3 %>% group_by(type3) %>% 
  summarise(count = sum(count)) %>% 
  ggplot(aes(x = fct_reorder(type3, count) %>% fct_rev(), y = count)) +
  geom_col(fill = "orange")+
  xlab("RNA type") + ylab("count of reads") +
  geom_text(aes(label = count), vjust = -0.5, hjust=0.5, colour = "black",size = 3)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
ggsave(paste0(filedir,"vienna/", "piRNA_target_RNA_type.png"), width =5, height =4, dpi = 300)

Vienna3 %>% group_by(piRNAlen) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(x = piRNAlen, y = count)) + geom_col(fill = "orange")+
  xlab("piRNA length") + ylab("Number of chimeric reads") +
  geom_text(aes(label = count), vjust = -0.5, hjust=0.5, colour = "black",size = 3)+
  scale_x_continuous(breaks = c(24,25,26,27,28,29,30))
ggsave(paste0(filedir,"vienna/", "piRNAlen_distri_collapse.png"), width =5, height =4, dpi = 300)
Vienna3 %>% group_by(piRNAlen) %>% 
  summarise(count = sum(count)) %>% 
  ggplot(aes(x = piRNAlen, y = count)) + geom_col(fill = "orange")+
  xlab("piRNA length") + ylab("Number of chimeric reads") +
  geom_text(aes(label = count), vjust = -0.5, hjust=0.5, colour = "black",size = 3)+
  scale_x_continuous(breaks = c(24,25,26,27,28,29,30))
ggsave(paste0(filedir,"vienna/", "piRNAlen_distri.png"), width =5, height =4, dpi = 300)

Vienna3 %>% write_tsv(paste0(filedir,"vienna/", "Piwi_CLASH.tsv"))

Vienna4 <- Vienna3 %>% 
  select(index, gene_name,target, piRNAlen,seq, contact2, deltaG) %>% 
  mutate(heat = contact2 %>% str_replace_all("\\.","0-") %>% str_replace_all("\\(","1-")) 

Vienna5 <- Vienna4 %>%
  mutate(score = map(heat, function(x) {x %>% str_sub(end = -2) %>% str_split(pattern = "-")})) %>% 
  unnest(cols = score) %>% 
  mutate(base = map(score, function(x) {seq(1, length(x), by = 1)})) %>% 
  select(index, deltaG, piRNAlen, gene_name, score, base) %>% 
  unnest(cols = c(score, base)) %>% 
  mutate(score_G = as.integer(score) * deltaG,
         gene_name = paste0(gene_name, "_", index))
Vienna5 %>%  
  ggplot() +
  geom_tile(aes(y =gene_name, fill = score, x=base))+
  facet_wrap(~piRNAlen, scales = "free_y") + 
  scale_fill_manual(values = c("0"="white","1"= "black"))
Vienna5 %>%  filter(piRNAlen %in% c(25,26,27,28)) %>%
  mutate(piRNAlen2 = paste0("piRNA (", piRNAlen, "nt)")) %>% 
  ggplot() +
  geom_tile(aes(y =gene_name, fill = score_G, x=base))+
  facet_wrap(~piRNAlen2, strip.position ="top", scales = "free_y", ncol=4)+
  scale_fill_gradient(low = "black", high = "white")+
  labs(fill = "dG (kcal/mol)", x = "bases (nt)", y= "ordered by k-means (center = 5)")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        strip.background=element_blank())
ggsave(paste0(filedir,"vienna/", "test.png"),  width =6, height =4, dpi = 300)

hcl <- Vienna5 %>% select(-score, -deltaG) %>% 
  group_by(piRNAlen) %>% nest() %>% 
  mutate(miranda = map(data, function(x) {x %>% spread(base, score_G)}),
         hclust = map(miranda, function(x) {x %>% select(-gene_name, -index) %>%
             as.matrix() %>% dist() %>% hclust(method = "ward.D2")}),
         order = map(hclust, function(x) {x[["order"]]})) %>% 
  select(-hclust, -data) %>% 
  unnest(cols = c(miranda, order))
hcl2 <- hcl %>% gather(-piRNAlen, -index, -order, -gene_name, key = "base", value = "score_G") %>% 
  drop_na() %>% 
  mutate(base = base %>% as.integer()) %>% 
  arrange(piRNAlen, order, base)
hcl2 %>% filter(piRNAlen %in% c(25,26,27,28)) %>%
  mutate(piRNAlen2 = paste0("piRNA (", piRNAlen, "nt)")) %>% 
  ggplot() +
  geom_tile(aes(y =fct_reorder(gene_name, order), fill = score_G, x=base))+
  facet_wrap(~piRNAlen2, strip.position ="top", scales = "free_y", ncol=4)+
  scale_fill_gradient(low = "black", high = "white")+
  labs(fill = "dG (kcal/mol)", x = "bases (nt)", y= "ordered by hcl")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        strip.background=element_blank())
ggsave(paste0(filedir, "vienna/","test_hcl.png"), width =6, height =4, dpi = 300)

#k-means
kmeans <- Vienna5 %>% select(-score) %>% 
  filter(piRNAlen %in% c(25, 26, 27, 28)) %>% 
  group_by(piRNAlen) %>% nest() %>% 
  mutate(miranda = map(data, function(x) {x %>% spread(base, score_G)}),
         km = map(miranda, function(x) {x %>% select(-index, -gene_name, -deltaG) %>% 
             kmeans(centers = 3, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")}),
         cluster = map(km, function(x) {x[["cluster"]]})) %>% 
  select(-km, -data) %>% 
  unnest(cols = c(miranda, cluster))
kmeans2 <- kmeans %>% gather(-piRNAlen, -index, -cluster, -gene_name, -deltaG, key = "base", value = "score_G") %>% 
  drop_na() %>% 
  mutate(base = base %>% as.integer()) %>% 
  group_by(cluster) %>% nest() %>% 
  mutate(center_dG = map(data, function(x) {x$deltaG %>% mean()})) %>% 
  unnest(cols = c(data, center_dG)) %>% 
  arrange(piRNAlen, cluster, gene_name, base) %>%
  mutate(piRNAlen2 = paste0("piRNA (", piRNAlen, "nt)"))
kmeans2 %>% ggplot() +
  geom_tile(aes(y =fct_reorder(gene_name, center_dG) %>% fct_rev(), fill = score_G, x=base))+
  facet_wrap(~piRNAlen2, strip.position ="top", scales = "free_y", ncol=4)+
  scale_fill_gradient(low = "black", high = "white")+
  labs(fill = "dG (kcal/mol)", x = "bases (nt)", y= "ordered by k-means (center = 3)")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        strip.background=element_blank())
ggsave(paste0(filedir, "vienna/","test_kmeans.png"), width =6, height =4, dpi = 300) 

#elbow plot
kmeans_elbow <- Vienna5 %>% select(-score, -deltaG) %>% 
  filter(piRNAlen == 26) %>% 
  spread(base, score_G)
kmeans_elbow2 <- kmeans_elbow %>% select(-index, -piRNAlen, -gene_name)
pct_var <- data.frame(pct_var = 0, num_clusters = 2:20)
totalss <- kmeans(kmeans_elbow2, centers = 20, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")$totss
for(i in 2:20){
  pct_var[i-1, 'pct_var'] <- kmeans(kmeans_elbow2, centers = i, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")$betweenss/totalss
}
ggplot(pct_var, aes(x = num_clusters, y = pct_var)) + 
  geom_line() + geom_point() + theme_bw() +
  ylim(0, 1) + xlab("Number of clusters") + ylab("% Variance Explained")
ggsave(paste0(filedir,"vienna/", "Variance_kmeanselbow_piRNA.png"),width =6, height =5, dpi = 300)

######################################
#shuffle
######################################
filedir2 <- "/media/meg/Titania/Piwi_CLASH/shuffle"
files1 <- dir(filedir2, pattern = "\\.vienna$", full.names = TRUE)
files2 <- files1 %>%
  str_replace_all("/media/meg/Titania/Piwi_CLASH/shuffle/Piwi_CLASH_comp_dm6_25_2_hybrids_ua_|.vienna$", "")
df2 <- vector("list", length(files1))
for (i in seq_along(files1)) {
  df2[[i]] <- read_csv(files1[[i]], 
                       col_names = "data") %>% 
    mutate(index = rep(1:(length(.$data) / 3), each = 3),
           type = if_else(row_number() %% 3 == 1,"name",
                          if_else(row_number() %% 3 == 2, "seq", "score"))) %>% 
    spread(key = "type", value = "data") 
  df2[[i]][["sample"]] <- files2[[i]]
}
Vienna_shuffle <- df2 %>% bind_rows() %>% 
  mutate(name = str_replace_all(name, "transposable_element", "TE"),
         name = str_replace(name, "pre_miRNA", "premiRNA"),
         name = str_replace(name, "_shuf", "")) %>% 
  separate(score, sep = "\\t", into= c("contact", "deltaG")) %>% 
  separate(name, sep = "\\_", into= c("dis11", "count", "gene_name", "dis", "dis2", "type2", "start", "end", "dis4", 
                                      "target", "dis5", "dis6", "type3", "start2", "end2")) %>% 
  separate(end, sep = "-", into = c("end", "dis10")) %>% 
  select(-contains("dis")) %>% 
  mutate(deltaG = str_replace(deltaG, "\\(", "") %>% str_replace("\\)", "") %>% as.double(),
         count = count %>% as.integer(),
         start = start %>% as.integer(),
         end = end %>% as.integer()) %>%
  filter(type2 == "piRNA") %>% 
  filter(start==1) %>% 
  mutate(contact2 = str_sub(contact, 1, (end - start + 1)),
         length = nchar(seq),
         piRNAlen = end-start+1)%>% 
  filter(str_detect(contact2, "\\(")) %>% 
  left_join(piRNA_length, by = "gene_name") %>% 
  filter(pi_length == piRNAlen)

Vienna3 %>% mutate(sample = "chimera") %>% 
  bind_rows(Vienna_shuffle) %>% 
  mutate(sample2 = if_else(sample =="chimera", sample,
                           if_else(sample == "5shuffle", "5' were shuffled",
                                   if_else(sample == "3shuffle", "3' were shuffled",
                                           "both were shuffled")))) %>% 
  mutate(deltaG2 = if_else(deltaG %% 1 >=0.5, deltaG %/% 1+1, deltaG %/% 1)) %>% 
  group_by(sample2, deltaG2) %>% summarise(total = sum(count)) %>% 
  ggplot(aes(x = deltaG2, y = total)) + geom_line(aes(color = sample2), size =1)+
  geom_point(aes(color = sample2), size = 2)+
  xlab("dG (kcal/mol)") + ylab("Number of chimeric reads")+
  guides(color = guide_legend(title =NULL))+
  theme(legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave(paste0(filedir, "vienna/", "shuffle_dG.png"),width =5, height =4, dpi = 300)

Vienna3 %>% mutate(sample = "chimera") %>% 
  bind_rows(Vienna_shuffle) %>% 
  mutate(sample2 = if_else(sample =="chimera", sample,
                           if_else(sample == "5shuffle", "5' were shuffled",
                                   if_else(sample == "3shuffle", "3' were shuffled",
                                           "both were shuffled")))) %>% 
  filter(type3 =="mRNA") %>% 
  mutate(deltaG2 = if_else(deltaG %% 1 >=0.5, deltaG %/% 1+1, deltaG %/% 1)) %>% 
  group_by(sample2, deltaG2) %>% summarise(total = sum(count)) %>% 
  ggplot(aes(x = deltaG2, y = total)) + geom_line(aes(color = sample2), size =1)+
  geom_point(aes(color = sample2), size = 2)+
  xlab("dG (kcal/mol)") + ylab("Number of chimeric reads")+
  guides(color = guide_legend(title =NULL))+
  theme(legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave(paste0(filedir, "vienna/", "shuffle_dG_mRNA.png"),width =5, height =4, dpi = 300)
#no more use code. for reference
#colfunc <- colorRampPalette(c("white", "black"))
#Vienna6 <- Vienna5 %>% select(-score) %>% 
#  spread(base, score_G) %>% select(-gene_name) %>% 
# column_to_rownames(var = "index") %>% 
# as.matrix()
#Vienna6 %>%  heatmap.2(col=colfunc(50), 
#         density.info = "none", trace = "none", dendrogram="row",
#         labCol = "", Colv =NA, margin=c(1, 1))
#test4 <- Vienna3 %>% filter(type3 %in%c("TE", "emsenbleTE"))

################################################################
###############################################################
###################################################################

setwd("/media/meg/Titania/mapping/AGO1_CLASH/")
library(tidyverse)
filedir <- "/media/meg/Titania/mapping/AGO1_CLASH/"
files3 <- dir(filedir, pattern = "\\_comp_hOH7_hybrids_ua.vienna$", full.names = TRUE)
files4 <- files3 %>%
  str_replace_all("/media/meg/Titania/mapping/AGO1_CLASH//|_comp_hOH7_hybrids_ua.vienna$", "") %>% 
  str_replace_all("_trim", "")
df2 <- vector("list", length(files3))
for (i in seq_along(files3)) {
  df2[[i]] <- read_csv(files3[[i]], 
                       col_names = "data") %>% 
    mutate(index = rep(1:(length(.$data) / 3), each = 3),
           type = if_else(row_number() %% 3 == 1,"name",
                          if_else(row_number() %% 3 == 2, "seq", "score"))) %>% 
    spread(key = "type", value = "data") 
  df2[[i]][["sample"]] <- files4[[i]]
}


Vienna_ago <- df2 %>% bind_rows()
miRNA_length <- read_tsv("~/hyb/data/fastq/hOH7-microRNA.fasta.fai", 
                         col_names = c("miRNA", "mi_length")) %>% 
  separate(miRNA, sep = "\\_", into = c("gene_name"))
Vienna2_ago <- Vienna_ago %>% 
  mutate(name = str_replace_all(name, "transposable_element", "TE"),
         name = str_replace(name, "pre_miRNA", "premiRNA")) %>% 
  separate(score, sep = "\\t", into= c("contact", "deltaG")) %>% 
  separate(name, sep = "\\_", into= c("dis11", "count", "gene_name", "dis", "dis2", "type2", "start", "end", "dis4", 
                                      "target", "dis5", "dis6", "type3", "start2", "end2")) %>% 
  separate(end, sep = "-", into = c("end", "dis10")) %>% 
  select(-contains("dis")) %>% 
  mutate(deltaG = str_replace(deltaG, "\\(", "") %>% str_replace("\\)", "") %>% as.double(),
         count = count %>% as.integer(),
         start = start %>% as.integer(),
         end = end %>% as.integer())
Vienna2_ago %>% group_by(type2, type3) %>% 
  summarise(count = n()) %>% 
  ggplot(aes(y =type3, x=type2))+
  geom_tile(aes(fill = count))+
  geom_text(aes(label = count), size = 3)+
  scale_fill_gradient(low = "white", high = "green")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  xlab("3' RNA type") + ylab("5' RNA type")
ggsave(paste0(filedir,"vienna/", "RNAtype_test_collapse.png"), width =7, height =6, dpi = 300)
Vienna2_ago %>% group_by(type2, type3) %>% 
  summarise(count = sum(count)) %>% 
  ggplot(aes(y =type3, x=type2))+
  geom_tile(aes(fill = count))+
  geom_text(aes(label = count), size = 3)+
  scale_fill_gradient(low = "white", high = "green")+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))+
  xlab("3' RNA type") + ylab("5' RNA type")
ggsave(paste0(filedir,"vienna/", "RNAtype_test_total.png"), width =7, height =6, dpi = 300)
Vienna3_ago <- Vienna2_ago %>% filter(type2 == "microRNA") %>% 
  mutate(contact2 = str_sub(contact, 1, (end - start + 1)),
         length = nchar(seq),
         miRNAlen = end-start+1) %>% 
  filter(str_detect(contact2, "\\(")) %>% 
  left_join(miRNA_length, by = "gene_name") %>% 
  filter(mi_length == miRNAlen)
Vienna3_ago %>% mutate(deltaG2 = if_else(deltaG %% 1 >=0.5, deltaG %/% 1+1, deltaG %/% 1)) %>% 
  group_by(deltaG2) %>% summarise(total = sum(count)) %>% 
  ggplot(aes(x = deltaG2, y = total)) + geom_line(color= "blue")+
  xlab("dG (kcal/mol)") + ylab("Number of chimeric reads")
ggsave(paste0(filedir,"vienna/", "distri.png"), width =7, height =6, dpi = 300)
Vienna3_ago %>% ggplot(aes(x = deltaG)) + geom_freqpoly(binwidth=1,color= "blue")
ggsave(paste0(filedir,"vienna/", "distri_collapse.png"), width =7, height =6, dpi = 300)

Vienna3_ago %>% write_tsv(paste0(filedir,"vienna/", "AGO1_CLASH.tsv"))

Vienna3_ago %>% group_by(type3) %>% 
  summarise(count_collapse = n()) %>% 
  ggplot(aes(x = fct_reorder(type3, count_collapse) %>% fct_rev(), y = count_collapse)) +
  geom_col(fill = "orange")+
  xlab("RNA type") + ylab("count of collapsed reads") +
  geom_text(aes(label = count_collapse), vjust = -0.5, hjust=0.5, colour = "black",size = 3)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
ggsave(paste0(filedir,"vienna/", "miRNA_target_RNA_type_collapse.png"), width =5, height =4, dpi = 300)
Vienna3_ago %>% group_by(type3) %>% 
  summarise(count = sum(count)) %>% 
  ggplot(aes(x = fct_reorder(type3, count) %>% fct_rev(), y = count)) +
  geom_col(fill = "orange")+
  xlab("RNA type") + ylab("count of reads") +
  geom_text(aes(label = count), vjust = -0.5, hjust=0.5, colour = "black",size = 3)+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
ggsave(paste0(filedir,"vienna/", "miRNA_target_RNA_type.png"), width =5, height =4, dpi = 300)

Vienna4_ago <- Vienna3_ago %>% filter(type3 == "mRNA") %>% 
  select(index, gene_name,target, miRNAlen,seq, contact2, deltaG, sample, count) %>% 
  mutate(heat = contact2 %>% str_replace_all("\\.","0-") %>% str_replace_all("\\(","1-"))
Vienna5_ago <- Vienna4_ago %>% mutate(score = map(heat, function(x) {x %>% str_sub(end = -2) %>% str_split(pattern = "-")})) %>% 
  unnest(cols = score) %>% 
  mutate(base = map(score, function(x) {seq(1, length(x), by = 1)})) %>% 
  select(sample, index, deltaG, miRNAlen, gene_name, score, base, count) %>% 
  unnest(cols = c(score, base)) %>% 
  mutate(score_G = as.integer(score) * deltaG,
         gene_name = paste0(sample, "_", gene_name, "_", index))
Vienna5_ago %>%  
  ggplot() +
  geom_tile(aes(y =gene_name, fill = score, x=base))+
  facet_wrap(~miRNAlen, scales = "free_y") + 
  scale_fill_manual(values = c("0"="white","1"= "black"))
Vienna5_ago %>%  filter(miRNAlen %in% c(21,22,23)) %>%
  mutate(miRNAlen2 = paste0("miRNA (", miRNAlen, "nt)")) %>% 
  ggplot() +
  geom_tile(aes(y =gene_name, fill = score_G, x=base))+
  facet_wrap(~miRNAlen, strip.position ="top", scales = "free_y")+
  scale_fill_gradient(low = "black", high = "white")+
  labs(fill = "dG (kcal/mol)", x = "bases (nt)", y= "ordered by k-means (center = 5)")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top",
        strip.background=element_blank())
ggsave(paste0(filedir,"vienna/", "miRNA_test.png"), width =12, height =20, dpi = 300)

hcl_ago <- Vienna5_ago %>% select(-score, -deltaG) %>%
  filter(miRNAlen %in% c(21,22,23)) %>% 
  group_by(miRNAlen) %>% nest()  %>% 
  mutate(miranda = map(data, function(x) {x %>% spread(base, score_G)}),
         hclust = map(miranda, function(x) {x %>% select(-gene_name, -index) %>%
             as.matrix() %>% dist() %>% hclust(method = "ward.D2")}),
         order = map(hclust, function(x) {x[["order"]]})) %>% 
  select(-hclust, -data) %>% 
  unnest(cols = c(miranda, order))
hcl2_ago <- hcl_ago %>% gather(-miRNAlen, -index, -order, -gene_name, -sample, key = "base", value = "score_G") %>% 
  drop_na() %>% 
  mutate(base = base %>% as.integer()) %>% 
  arrange(miRNAlen, order, base)
hcl2_ago %>% ggplot() +
  geom_tile(aes(y =fct_reorder(gene_name, order), fill = score_G, x=base))+
  facet_wrap(~miRNAlen, scales = "free_y")+
  scale_fill_gradient(low = "black", high = "white")
ggsave("~/Desktop/test_hcl_ago.png", width =12, height =20, dpi = 300)     

#k-means
kmeans_ago <- Vienna5_ago %>% select(-score) %>% 
  filter(miRNAlen ==22) %>% 
  group_by(miRNAlen) %>% nest() %>% 
  mutate(miranda = map(data, function(x) {x %>% spread(base, score_G)}),
         km = map(miranda, function(x) {x %>% select(-index, -gene_name,-sample, -deltaG, -count) %>% 
             kmeans(centers = 5, nstart = 1000, iter.max = 1000, algorithm = "MacQueen")}),
         cluster = map(km, function(x) {x[["cluster"]]})) %>% 
  select(-km, -data) %>% 
  unnest(cols = c(miranda, cluster))
kmeans2_ago <- kmeans_ago %>% 
  gather(-miRNAlen, -index, -cluster, -gene_name, -sample, -deltaG, -count, key = "base", value = "score_G") %>% 
  drop_na() %>% 
  mutate(base = base %>% as.integer()) %>% 
  group_by(cluster) %>% nest() %>% 
  mutate(center_dG = map(data, function(x) {x$deltaG %>% mean()})) %>% 
  unnest(cols = c(data, center_dG)) %>% 
  arrange(miRNAlen, cluster, gene_name, base) %>%
  mutate(miRNAlen2 = paste0("miRNA (", miRNAlen, "nt)"),
         cluster2 = paste0("cluster", cluster))
p1 <- kmeans2_ago %>% ggplot() +
  geom_tile(aes(y =fct_reorder(gene_name, center_dG) %>% fct_rev(), fill = score_G, x=base))+
  labs(fill = "dG (kcal/mol)", x = "bases (nt)", y= "ordered by k-means (center = 5)")+
  scale_fill_gradient(low = "black", high = "white")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top") +
  scale_x_continuous(breaks = c(1,5,10,15,20,22))
p2 <- kmeans2_ago  %>% ggplot() +
  geom_tile(aes(y =fct_reorder(gene_name,center_dG) %>% fct_rev(), fill = cluster2, x = base))+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(color = "white"),
        axis.text.x = element_text(color = "white"),
        axis.ticks.x = element_line(color = "white"))
g_legend <- function(A) {
  tmp <- ggplot_gtable(ggplot_build(A))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
mylegend <- g_legend(p1)
ppi <- 300
png(paste0(filedir,"vienna/", "kmeans_test.png"), width = 5*ppi, height = 6*ppi, res = ppi)
grid.arrange(p1 + theme(legend.position = "none"),
             p2 + theme(legend.position = "none"),
             ncol = 2, widths = c(15,1),
             mylegend, nrow = 2, heights = c(10,1))
dev.off()

###################################
filedir2 <- "/media/meg/Titania/mapping/AGO1_CLASH/shuffle"
files1 <- dir(filedir2, pattern = "\\.vienna$", full.names = TRUE)
files2 <- files1 %>%
  str_replace_all("/media/meg/Titania/mapping/AGO1_CLASH/shuffle/|.vienna$", "") %>% 
  str_replace_all("_comp_hOH7_hybrids_ua", "") %>% 
  str_replace_all("_trim", "")
df2 <- vector("list", length(files1))
for (i in seq_along(files1)) {
  df2[[i]] <- read_csv(files1[[i]], 
                       col_names = "data") %>% 
    mutate(index = rep(1:(length(.$data) / 3), each = 3),
           type = if_else(row_number() %% 3 == 1,"name",
                          if_else(row_number() %% 3 == 2, "seq", "score"))) %>% 
    spread(key = "type", value = "data") 
  df2[[i]][["sample"]] <- files2[[i]]
}
Vienna_shuffle_ago <- df2 %>% bind_rows() %>% 
  mutate(name = str_replace_all(name, "transposable_element", "TE"),
         name = str_replace(name, "pre_miRNA", "premiRNA"),
         name = str_replace(name, "_shuf", "")) %>% 
  separate(score, sep = "\\t", into= c("contact", "deltaG")) %>% 
  separate(name, sep = "\\_", into= c("dis11", "count", "gene_name", "dis", "dis2", "type2", "start", "end", "dis4", 
                                      "target", "dis5", "dis6", "type3", "start2", "end2")) %>% 
  separate(end, sep = "-", into = c("end", "dis10")) %>% 
  select(-contains("dis")) %>% 
  mutate(deltaG = str_replace(deltaG, "\\(", "") %>% str_replace("\\)", "") %>% as.double(),
         count = count %>% as.integer(),
         start = start %>% as.integer(),
         end = end %>% as.integer()) %>%
  filter(type2 == "microRNA") %>% 
  filter(start==1) %>% 
  mutate(contact2 = str_sub(contact, 1, (end - start + 1)),
         length = nchar(seq),
         miRNAlen = end-start+1)%>% 
  filter(str_detect(contact2, "\\(")) %>% 
  left_join(miRNA_length, by = "gene_name") %>% 
  filter(mi_length == miRNAlen)

Vienna3_ago %>% mutate(sample = "chimera") %>% 
  bind_rows(Vienna_shuffle_ago) %>% 
  mutate(sample2 = if_else(sample =="chimera", sample,
                           if_else(str_detect(sample, "5shuffle_3shuffle"),  "both were shuffled",
                                   if_else(str_detect(sample, "3shuffle"), "3' were shuffled",
                                           "5' were shuffled")))) %>% 
  mutate(deltaG2 = if_else(deltaG %% 1 >=0.5, deltaG %/% 1+1, deltaG %/% 1)) %>% 
  group_by(sample2, deltaG2) %>% summarise(total = sum(count)) %>% 
  ggplot(aes(x = deltaG2, y = total)) + geom_line(aes(color = sample2), size =1)+
  geom_point(aes(color = sample2), size = 2)+
  xlab("dG (kcal/mol)") + ylab("Number of chimeric reads")+
  guides(color = guide_legend(title =NULL))+
  theme(legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave(paste0(filedir, "vienna/", "shuffle_dG.png"),width =5, height =4, dpi = 300)

Vienna3_ago %>% mutate(sample = "chimera") %>% 
  bind_rows(Vienna_shuffle_ago) %>% 
  mutate(sample2 = if_else(sample =="chimera", sample,
                           if_else(str_detect(sample, "5shuffle_3shuffle"),  "both were shuffled",
                                   if_else(str_detect(sample, "3shuffle"), "3' were shuffled",
                                           "5' were shuffled")))) %>% 
  filter(type3 =="mRNA") %>% 
  mutate(deltaG2 = if_else(deltaG %% 1 >=0.5, deltaG %/% 1+1, deltaG %/% 1)) %>% 
  group_by(sample2, deltaG2) %>% summarise(total = sum(count)) %>% 
  ggplot(aes(x = deltaG2, y = total)) + geom_line(aes(color = sample2), size =1)+
  geom_point(aes(color = sample2), size = 2)+
  xlab("dG (kcal/mol)") + ylab("Number of chimeric reads")+
  guides(color = guide_legend(title =NULL))+
  theme(legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave(paste0(filedir, "vienna/", "shuffle_dG_mRNA.png"),width =5, height =4, dpi = 300)
