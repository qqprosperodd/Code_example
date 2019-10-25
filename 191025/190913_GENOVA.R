library(GENOVA)
library(tidyverse)
library(bindrcpp)
library(BiocStyle)
library(bigwrig)
library(RColorBrewer)
#color.bar関数はパッケージのない関数らしい。そこでここで予め定義する。
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}

insulation.callTAD <- function(experiment,  BEDcolor = "127,201,127", verbose = F){
  if(is.null(experiment[["INSULATION"]])){ stop("Call insulation score first and store in experiment$INSULATION")}
  res <- experiment[["RES"]]
  scooch <- floor(100e3 / res)
  entries <- list()
  df = NULL
  CHROMS <- unique(experiment[["INSULATION"]][,1])
  
  experiment[["INSULATION"]]$V2 <- experiment[["INSULATION"]][,2] + experiment[["RES"]]
  
  for(CCC in CHROMS){
    if(verbose){   message("Starting chromosome",CCC, "\n")  }
    
    INSU <- experiment[["INSULATION"]][experiment[["INSULATION"]][,1] == CCC ,]
    insCol <- INSU[,4]
    
    #determine minimum and maximum values
    min_value <- min(insCol[insCol != -Inf])
    max_value <- max(insCol[insCol != Inf])
    #check for every insulatoin score that is -Inf or Inf if
    #    the value before or after it is the same
    #to prevent "sudden" peaks to be filtered out
    for (i in (2:length(insCol))){
      insCol[i==-Inf & i - 1 != -Inf & i + 1 != -Inf] <-  min_value *2
      insCol[i==Inf & i - 1 != Inf & i + 1 != Inf] <-  max_value * 2
    }
    # set values of Inf and -Inf that are left to the minimal and maximal value of the chromosome
    insCol[insCol == -Inf] <- min_value
    insCol[insCol == Inf] <- max_value
    
    #set any remaining values that ar not finite (e.g. NA)
    # to NaN (only for chrY as this only has -Inf and Inf)
    insCol[!(is.finite(unlist(insCol)))] <- NaN
    INSU <- cbind(INSU[,1:3],scale(insCol, center = TRUE, scale = TRUE))
    colnames(INSU)[4] <- "V4"
    if( length(INSU$V4) == length(which(is.nan(INSU$V4)))){break}
    
    if(verbose){   message("Computing the delta-vector...\n")  }
    add <- 1:scooch
    i <- (scooch+1):(nrow(INSU)-scooch)
    i.rep <- rep(i, each=length(add))
    right <- matrix(INSU$V4[i.rep+add], ncol=scooch, byrow=T)
    left <- matrix(INSU$V4[i.rep-add], ncol=scooch, byrow=T)
    right <- apply(right, 1, mean)
    left <- apply(left, 1, mean)
    delta = left - right
    deltaDF <- data.frame(INSU[i,1:4], delta)
    colnames(deltaDF)[1:4] <- paste0("V", 1:4)
    colnames(deltaDF)[5] <- "delta"
    deltaDF <- dplyr::arrange(deltaDF,  V1,V2)
    deltaDF$ID <- 1:nrow(deltaDF)
    
    ####
    # First find peaks
    ###
    VALLEYS <- diff(c(.Machine$integer.max, deltaDF$V4)) > 0L
    VALLEYS <- cumsum(rle(VALLEYS)$lengths)
    VALLEYS <- VALLEYS[seq.int(1L, length(VALLEYS), 2L)]
    if (deltaDF$V4[[1]] == deltaDF$V4[[2]]) {
      VALLEYS <- VALLEYS[-1]
    }
    if(verbose){message("Calling borders...\n")}
    
    boundaryCalls <- NULL
    VALLEYS <- sort(unique(c(VALLEYS+1, VALLEYS, VALLEYS-1)))
    VALLEYS <- VALLEYS[VALLEYS > 2]
    for(j in VALLEYS[2:length(VALLEYS)]){
      dat <- deltaDF[(j-1):(j+1),]
      normalOrder <- order(dat$delta)
      dat0 <- dat
      dat0$delta[2] <- 0
      zeroOrder <- order(dat0$delta)
      if(! all(complete.cases(dat))){next}
      if(! all(normalOrder == c(3,2,1))){next}
      if(!all(zeroOrder == normalOrder)){next}
      UP <- max(deltaDF[(j-1),5],na.rm = T)
      DOWN <- min(deltaDF[(j+1),5],na.rm = T)
      if((UP - DOWN) < 0.1){next}
      if(is.null(boundaryCalls)) {
        boundaryCalls <- rbind(boundaryCalls,dat[2,])
      } else{
        nearbyPoints <- dplyr::filter(boundaryCalls, V2  >= (dat[2,2] - (res*scooch)) , V2 <= (dat[2,3] + (res*scooch) ))  # res*scooch
        if(nrow(nearbyPoints) == 0){
          boundaryCalls <- rbind(boundaryCalls,dat[2,])
        } else {
          closeToZeroDelta <- min(abs(c(nearbyPoints$delta, dat[2,]$delta)))
          winner <- rbind(nearbyPoints, dat[2,])[which(abs(c(nearbyPoints$delta, dat[2,]$delta)) == closeToZeroDelta),]
          losers <- rbind(nearbyPoints, dat[2,])[-which(abs(c(nearbyPoints$delta, dat[2,]$delta)) == closeToZeroDelta),]
          boundaryCalls <- rbind(boundaryCalls,dat[2,])
          boundaryCalls <- dplyr::anti_join(boundaryCalls, losers, by = c("V1", "V2", "V3", "V4", "delta"))
        }
      }
    }
    
    for(i in 1:nrow(boundaryCalls)){
      if( boundaryCalls[i,3] < boundaryCalls[i,2]  ){
        tmp3 <- boundaryCalls[i,3]
        tmp2 <- boundaryCalls[i,2]
        boundaryCalls[i,3] <- tmp2
        boundaryCalls[i,2] <- tmp3
      }
    }
    
    boundaryCalls <- boundaryCalls[with(boundaryCalls, order(V1, V2)), ]
    
    if(verbose){message("Generating bedgraph...\n")}
    for(i in 2:nrow(boundaryCalls)){
      if(boundaryCalls[i-1,1] != boundaryCalls[i,1]){next}
      prev <- boundaryCalls[i-1,3]
      now <- boundaryCalls[i,1:2]
      now[,2] <- now[,2]
      prev <- prev
      ddd <-cbind(now, prev,now,prev, BEDcolor)[c(1,3,2,1,3,2,7)]
      colnames(ddd) <- c("a", 'b', 'c', 'd', 'e', 'f', 'g')
      df <- rbind(df, ddd)
    }
  }
  
  return(bedgraph = df)
}





setwd("/media/meg/Titania/HiC_Pro/HiC_Out/hic_results/matrix/")
filedir <- "/media/meg/Titania/HiC_Pro/HiC_Out/"
EGFP_10kb <- construct.experiment(ignore.checks = T, # time-saver for vignette.
                                     signalPath = 'OSC1_EGFP_rep1/iced/10000/OSC1_EGFP_rep1_10000_iced.matrix',
                                     indicesPath = 'OSC1_EGFP_rep1/raw/10000/OSC1_EGFP_rep1_10000_abs.bed',
                                     name = "OSC1_EGFP",
                                     color = "black")
Piwi_10kb <- construct.experiment(ignore.checks = T, # time-saver for vignette.
                                     signalPath = 'OSC1_Piwi_rep1/iced/10000/OSC1_Piwi_rep1_10000_iced.matrix',
                                     indicesPath = 'OSC1_Piwi_rep1/raw/10000/OSC1_Piwi_rep1_10000_abs.bed',
                                     name = "Piwi_EGFP",
                                     color = "blue")

#cis-quantation
cisChrom_out <- cisTotal.perChrom(EGFP_10kb)
ppi <- 300
png(paste0(filedir, "cis_quantation.png"), width = 5*ppi, height = 5*ppi, res = ppi)
plot(cisChrom_out$perChrom, las = 2)
abline(h = cisChrom_out$genomeWide, col = "red")
dev.off()
#chromosomee plot
png(paste0(filedir, "chromosomme_contact.png"), width = 5*ppi, height = 5*ppi, res = ppi)
EGFP_10kb %>% chromosomeMatrix()
dev.off()
#RCP
RCP1 <- RCP(experimentList = list(EGFP_10kb, Piwi_10kb),
            chromsToUse = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"))
RCP1 %>% visualise.RCP.ggplot(smooth = T, combine = F)
ggsave(paste0(filedir, "RSP_each.png"), width =8, height =6, dpi = 300)
#each chromosome
RCP1 %>% visualise.RCP.ggplot(smooth = T, combine = T)
ggsave(paste0(filedir, "RSP_all.png"), width =5, height =4, dpi = 300)
#all chromosome

##################################
#insulation
##################################
cols <- c("#f03b20", "#ffeda0", "white", "#31a354")
#heatmap color
png(paste0(filedir, "dendro_test.png"), width = 10*ppi, height = 5*ppi, res = ppi)
layout(matrix(c(1,2), nrow = 1, ncol = 2), widths = c(1,5))
color.bar(colorRampPalette(cols)(100), -1, nticks = 5)
insulation.domainogram(EGFP_10kb,
                       'chr2L',
                       14e6,
                       18e6,
                       window.size1 = 1,
                       window.size2 = 101,
                       step = 2)
dev.off()

#genome-wide insulation score
#insulation scoreのために、chr4をのぞく
filter_chrom <- function(data, chr) {
  data1 <- data
  ICE <- data1[["ICE"]]
  ABS <- data1[["ABS"]]
  filter_chr <- ABS %>% filter(V1 == chr)
  ICE_fill <- ICE %>% filter(!(V1 %in% filter_chr$V4)) %>% 
    filter(!(V2 %in% filter_chr$V4))
  ABS_fill <- ABS %>% filter(V1 != chr)
  CHRS_fill <- data1[["CHRS"]] %>% subset(. != chr)
  new_name <- ABS_fill %>% rownames_to_column(var = "V5") %>% 
    mutate(V5 = V5 %>% as.integer()) %>% 
    select(V4, V5)
  data1[["ICE"]] <- ICE_fill %>% left_join(new_name %>% rename(V1 = V4), by = "V1") %>% 
    mutate(V1 = V5) %>% select(-V5) %>% left_join(new_name %>% rename(V2 = V4), by = "V2") %>% 
    mutate(V2 = V5) %>% select(V1,V2,V3) %>% data.table::as.data.table()
  data1[["ABS"]] <- ABS_fill %>% left_join(new_name, by = "V4") %>% 
    mutate(V4 = V5) %>% select(V1,V2,V3,V4) %>% as.data.frame()
  data1[["CHRS"]] <- CHRS_fill
  return(data1)
}

EGFP_10kb2 <- EGFP_10kb %>% filter_chrom('chr4') %>% filter_chrom('chrY')
Piwi_10kb2 <- Piwi_10kb %>% filter_chrom('chr4') %>% filter_chrom('chrY')


######################################################
EGFP_10kb_insulation <- genome.wide.insulation(hic = EGFP_10kb2,
                                               window.size = 25)
EGFP_10kb_insulation %>% write_tsv(paste0(filedir, "EGFP_10kb_insulation.tsv"))
Piwi_10kb_insulation <- genome.wide.insulation(hic = Piwi_10kb2,
                                               window.size = 25)
Piwi_10kb_insulation %>% write_tsv(paste0(filedir, "Piwi_10kb_insulation.tsv"))
EGFP_10kb_insulation <- read_tsv(paste0(filedir, "EGFP_10kb_insulation.tsv"),
                                 col_types = cols(start = col_integer(),
                                                  end = col_integer()))
Piwi_10kb_insulation <- read_tsv(paste0(filedir, "Piwi_10kb_insulation.tsv"),
                                 col_types = cols(start = col_integer(),
                                                  end = col_integer()))
#call TAD
EGFP_10kb[["INSULATION"]] <- EGFP_10kb_insulation
EGFP_10kb
TADcalls <- insulation.callTAD(EGFP_10kb, BEDcolor = "#7ec0ee", verbose = T)

  #get insulation heatmaps
insulation.heatmap_out = insulation.heatmap(
  insulationList = list(WT = Hap1_WT_10kb_insulation,
                        SCC4 = Hap1_SCC4_10kb_insulation ),
  bed = WT_TADs,
  zlim = c(-.5,0.25), # zlim.
  profileZlim = c(-.75,-.1) # zlim for the profile
)
  

