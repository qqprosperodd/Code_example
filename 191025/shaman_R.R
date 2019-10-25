  filedir <- "/media/meg/Titania/HiCExplorer/Ulianov/"      
  setwd(filedir)
  library(tidyverse)
  files <- dir(filedir, pattern = "\\.tsv$", full.names = TRUE)
  files2 <- files %>%
    str_replace_all("^/media/meg/Titania/HiCExplorer/hug//|.tsv$", "")
  
  filedir <- "/media/meg/Titania/HiCExplorer/cool/"      
  setwd(filedir)
  files <- dir(filedir, pattern = "\\.tsv$", full.names = TRUE)
  files2 <- files %>%
    str_replace_all("^/media/meg/Titania/HiCExplorer/cool//|.tsv$", "")
  
  df <- vector("list", length(files))
  for (i in seq_along(files)) {
    df[[i]] <- read_tsv(files[[i]], col_names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "cov"),
                        col_types = cols(chrom1 = col_character(),
                                         start1 = col_double(),
                                         end1 = col_double(),
                                         chrom2 = col_character(),
                                         start2 = col_double(),
                                          end2 = col_double(),
                                         cov = col_double())) %>%
      mutate(start1 = (start1 + end1)%/%2,
             end1 = start1 + 1,
             start2 = (start2 + end2)%/%2,
             end2 = start2 + 1)
    df[[i]][["type"]] <- files2[[i]]
  }
  tx2 <- bind_rows(df) %>%
    filter(!(chrom1 == "chrY" | chrom2 == "chrY")) %>% 
    group_by(type) %>% nest()
  rm(df)
  Den6 <- map(1:length(tx2$type), function(x) data.frame(tx2$data[[x]]))
  names(Den6) <- tx2$type
  rm(tx2)
  map(1:length(Den6), function(x) write_tsv(Den6[[x]], 
                                            paste0(filedir, "shaman/", names(Den6[x]), ".txt"),
                                            col_names = TRUE))
#################
################

setwd("~/Desktop/pra/")
library(misha)
#dm genome can be obtein from flybase, but the fa isnt separated.
#so I manually make fa and chromosizes.
#

#genome load
gdb.init('dm6_25')

gsetroot("dm6_25")
gtrack.2d.import(track="hic_obs_test1",
                 description="observed hic data", 
                 file="dm6_25/tracks/OSCtest10000merge.txt")
gtrack.2d.import(track="hic_obs_test_siEGFP",
                 description="observed hic data", 
                 file="dm6_25/tracks/siEGFPmerge2.txt")
gtrack.2d.import(track="hic_obs_test_siPiwi",
                 description="observed hic data", 
                 file="dm6_25/tracks/siPiwimerge2.txt")
gdb.reload()
library(shaman)
options(shaman.mc_support=1)
Ret <- shaman_shuffle_hic_track(track_db="~/Desktop/pra/dm6_25",
                         obs_track_nm="hic_obs_test_siEGFP", 
                         work_dir=tempdir())

#plot expected
gsetroot("~/Desktop/pra/dm6_25")
point_score <-  gextract("hic_obs_testhug_shuffle", 
                       gintervals.2d("2L", 19e06, 22e06, "2L", 12e06, 2e06),
                       colnames="score")

shaman_plot_map_score_with_annotations("dm6_25",
                                       point_score,
                                       gintervals("2L", 19e06, 22e06),
                                       a_colors=c("#4572A7", "#AA4643"))
ggsave(filename="ulianov_test.png",width=10,height=6,dpi=300)


  options(shaman.mc_support=1)    #configuring multi-core mode
options(gmax.data.size=1e+09)
if (gtrack.exists("hic_score_new")) {
  gtrack.rm("hic_score_new", force=TRUE)
  gdb.reload()
} 
Ret2 <- shaman_score_hic_track(track_db="~/Desktop/pra/dm6_25", 
                       work_dir=tempdir(),
                       score_track_nm="hic_score_new", 
                       obs_track_nms="hic_obs_test1",
                       exp_track_nms="hic_obs_test1_shuffle")
  
    

##############################
filedir <- "/media/meg/Titania/HiCExplorer/hug/"   
files <- files[3]
files2 <- files2[3]


