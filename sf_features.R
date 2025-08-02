library(tidyverse)
library(caret)
source("rf_script.R")
codon <- read_csv("codon.csv")
codon_common <- codon %>% 
  drop_na(AA_3l) %>% 
  mutate(AA_3l=factor(AA_3l, unique(AA_3l))) %>% 
  group_by(AA_3l) %>% 
  slice_max(ecoli, n=1) %>% 
  ungroup()
aa_all <- codon_common$AA_1l
dna_all <- codon_common$codon_dna

# Features -----

zscale <- read_csv("zscale.csv") %>% 
  rename(AA_3l=abbrev) %>% 
  select(AA_3l, z1:z5) %>% 
  left_join(codon) %>% 
  rename(aa=AA_1l) %>% 
  distinct(aa, .keep_all=TRUE) %>% 
  select(aa, z1:z5) %>% 
  filter(!is.na(aa))

tscale <- read_csv("tscale.csv") %>% 
  select(abv, starts_with("T")) %>% 
  filter(!is.na(abv)) %>% 
  rename(aa=abv)

aaindex <- read_csv("aaindex.csv") %>% 
  pivot_longer(-AAIndex, names_to="aa", values_to="aai") %>% 
  pivot_wider(names_from=AAIndex, values_from=aai)

stscale <- read_csv("stscale.csv") %>% 
  select(-aminoacid)

vhsescale  <- read_csv("vhsescale.csv") %>% 
  select(-aminoacid)

# lum_2.8ms data -----

# sf_mean <- read_csv("sf_mean.csv") %>% 
#   select(-sd_lum)
# 
# sf_mean_z <- sf_mean %>% 
#   add_1mer() %>% 
#   rep_z(zscale) %>% 
#   select(-aa, -(aa_0:aa_4)) %>% 
#   rowwise %>% 
#   mutate(z1=mean(c_across(starts_with("z1")))) %>% 
#   mutate(z2=mean(c_across(starts_with("z2")))) %>% 
#   mutate(z3=mean(c_across(starts_with("z3")))) %>% 
#   mutate(z4=mean(c_across(starts_with("z4")))) %>% 
#   mutate(z5=mean(c_across(starts_with("z5"))))
# 
# sf_mean_z_3mer <- sf_mean %>% 
#   mutate(aa=paste0(aa,"FS")) %>% 
#   add_3mer2() %>% 
#   pivot_longer(-c(aa,mean_lum),names_to="mer",values_to="seq") %>% 
#   left_join(scale_mer(aa1=zscale,aa2=zscale,aa3=zscale) %>% 
#               rename(seq=aa)) %>% 
#   mutate(aa=substr(aa,1,5)) %>% 
#   pivot_wider(id_cols=aa,names_from=mer,values_from=z1:z5)
#   
# sf_mean_t <- sf_mean %>% 
#   add_1mer() %>% 
#   rep_z(tscale) %>% 
#   select(-aa, -(aa_0:aa_4)) %>% 
#   rowwise %>% 
#   mutate(T1=mean(c_across(starts_with("T1")))) %>% 
#   mutate(T2=mean(c_across(starts_with("T2")))) %>% 
#   mutate(T3=mean(c_across(starts_with("T3")))) %>% 
#   mutate(T4=mean(c_across(starts_with("T4")))) %>% 
#   mutate(T5=mean(c_across(starts_with("T5"))))
# 
# sf_mean_t_3mer <- sf_mean %>% 
#   mutate(aa=paste0(aa,"FS")) %>% 
#   add_3mer2() %>% 
#   pivot_longer(-c(aa,mean_lum),names_to="mer",values_to="seq") %>% 
#   left_join(scale_mer(aa1=tscale,aa2=tscale,aa3=tscale) %>% 
#               rename(seq=aa)) %>% 
#   mutate(aa=substr(aa,1,5)) %>% 
#   pivot_wider(id_cols=aa,names_from=mer,values_from=T1:T5)
# 
# sf_mean_aai <- sf_mean %>% 
#   add_1mer() %>% 
#   rep_z(aaindex) %>% 
#   t_addmean(10) %>% 
#   select(-c(ends_with("_0"), 
#             ends_with("_1"),
#             ends_with("_2"),
#             ends_with("_3"),
#             ends_with("_4"))) %>% 
#   select(-aa)
# 
# sf_mean_st <- sf_mean %>% 
#   add_1mer() %>% 
#   rep_z(stscale) %>% 
#   t_addmean(3) %>% 
#   select(-starts_with("aa"), -mea)
# 
# sf_mean_st_3mer <- sf_mean %>% 
#   mutate(aa=paste0(aa,"FS")) %>% 
#   add_3mer2() %>% 
#   pivot_longer(-c(aa,mean_lum),names_to="mer",values_to="seq") %>% 
#   left_join(scale_mer(aa1=stscale,aa2=stscale,aa3=stscale) %>% 
#               rename(seq=aa)) %>% 
#   mutate(aa=substr(aa,1,5)) %>% 
#   pivot_wider(id_cols=aa,names_from=mer,values_from=ST1:ST8)
# 
# sf_mean_vhse <- sf_mean %>% 
#   add_1mer() %>% 
#   rep_z(vhsescale) %>% 
#   t_addmean(5) %>% 
#   select(-starts_with("aa"), -mean_)
# 
# sf_mean_vhse_3mer <- sf_mean %>% 
#   mutate(aa=paste0(aa,"FS")) %>% 
#   add_3mer2() %>% 
#   pivot_longer(-c(aa,mean_lum),names_to="mer",values_to="seq") %>% 
#   left_join(scale_mer(aa1=vhsescale,aa2=vhsescale,aa3=vhsescale) %>% 
#               rename(seq=aa)) %>% 
#   mutate(aa=substr(aa,1,5)) %>% 
#   pivot_wider(id_cols=aa,names_from=mer,values_from=VHSE1:VHSE8)
# 
# sf_mean_all <- sf_mean_z %>% 
#   bind_cols(sf_mean_t %>% select(-mean_lum))
# 
# sf_mean_all2 <- sf_mean_z %>% 
#   bind_cols(sf_mean_t %>% select(-mean_lum)) %>% 
#   bind_cols(sf_mean_st %>% select(-mean_lum)) %>% 
#   bind_cols(sf_mean_vhse %>% select(-mean_lum))
# 
# write_csv(sf_mean_all2, "sf_mean_all2.csv")
# 
# sf_mean_allres <- sf_mean_z %>% select(!z1:z5) %>% 
#   bind_cols(sf_mean_t %>% select(-mean_lum) %>% select(!T1:T5)) %>% 
#   bind_cols(sf_mean_st %>% select(-mean_lum) %>% select(!ST1:ST8)) %>% 
#   bind_cols(sf_mean_vhse %>% select(-mean_lum) %>% select(!VHSE1:VHSE8))
# 
# write_csv(sf_mean_allres, "sf_mean_allres.csv")
# 
# sf_mean_a2m <- sf_mean_z %>% 
#   bind_cols(sf_mean_t %>% select(-mean_lum)) %>% 
#   bind_cols(sf_mean_st %>% select(-mean_lum)) %>% 
#   bind_cols(sf_mean_vhse %>% select(-mean_lum)) %>% 
#   bind_cols(sf_mean_z_3mer %>% select(-aa)) %>% 
#   bind_cols(sf_mean_t_3mer %>% select(-aa)) %>% 
#   bind_cols(sf_mean_st_3mer %>% select(-aa)) %>% 
#   bind_cols(sf_mean_vhse_3mer %>% select(-aa))
# 
# write_csv(sf_mean_a2m, "sf_mean_a2m.csv")

# cf_aadna -----

bind_all_aadna <- function(cf_aadna){
  # cf_aadna_z <- cf_aadna %>% 
  #   add_1mer() %>% 
  #   rep_z(zscale) %>% 
  #   select(-aa, -dna, -(aa_0:aa_4)) %>% 
  #   rowwise %>% 
  #   mutate(z1=mean(c_across(starts_with("z1")))) %>% 
  #   mutate(z2=mean(c_across(starts_with("z2")))) %>% 
  #   mutate(z3=mean(c_across(starts_with("z3")))) %>% 
  #   mutate(z4=mean(c_across(starts_with("z4")))) %>% 
  #   mutate(z5=mean(c_across(starts_with("z5"))))
  
  cf_aadna_z <- cf_aadna %>% 
    add_1mer() %>% 
    rep_z(zscale) %>% 
    select(-dna) %>% 
    t_addmean(2, "z") %>% 
    select(-starts_with("aa"))
  
  cf_aadna_t <- cf_aadna %>% 
    add_1mer() %>% 
    rep_z(tscale) %>% 
    select(-dna) %>% 
    t_addmean(2, "T") %>% 
    select(-starts_with("aa"))
  
  # cf_aadna_t <- cf_aadna %>% 
  #   add_1mer() %>% 
  #   rep_z(tscale) %>% 
  #   select(-aa, -dna, -(aa_0:aa_4)) %>% 
  #   rowwise %>% 
  #   mutate(T1=mean(c_across(starts_with("T1")))) %>% 
  #   mutate(T2=mean(c_across(starts_with("T2")))) %>% 
  #   mutate(T3=mean(c_across(starts_with("T3")))) %>% 
  #   mutate(T4=mean(c_across(starts_with("T4")))) %>% 
  #   mutate(T5=mean(c_across(starts_with("T5"))))
  
  cf_aadna_st <- cf_aadna %>% 
    add_1mer() %>% 
    rep_z(stscale) %>% 
    select(-dna) %>% 
    t_addmean(3, "ST") %>% 
    select(-starts_with("aa"))
  
  cf_aadna_vhse <- cf_aadna %>% 
    add_1mer() %>% 
    rep_z(vhsescale) %>% 
    select(-dna) %>% 
    t_addmean(5, "VHSE") %>% 
    select(-starts_with("aa"))
  
  cf_aadna_all2 <- cf_aadna_z %>% 
    bind_cols(cf_aadna_t %>% select(-energy)) %>% 
    bind_cols(cf_aadna_st %>% select(-energy)) %>% 
    bind_cols(cf_aadna_vhse %>% select(-energy))
  
  cf_aadna_all2
}

# cf_aadna <- read_csv("cf_aadna.csv") %>% 
#   filter(nchar(aa)==5) %>% 
#   select(-sd_lum)
# 
# cf_aadna_new <- cf_aadna %>% 
#   bind_rows(read_csv("energy_240524.csv") %>% 
#               select(dna, aa, energy, mean_lum) %>% 
#               na.omit)
# 
# cf_aadna_new_all <- bind_all_aadna(cf_aadna_new)
# 
# write_csv(cf_aadna_all2, "cf_aadna_all2.csv")
# write_csv(cf_aadna_new_all, "cf_aadna_new_all.csv")
# 
# energy_all <- read_csv("energy2_clean.csv")
# 
# cf_aadna_2409 <- read_csv("lum_240911.csv") %>% 
#   select(aa, dna, average) %>% 
#   rename(mean_lum=average) %>% 
#   left_join(energy_all)
# 
# cf_aadna_2409_bind <- cf_aadna_new %>% 
#   rows_upsert(cf_aadna_2409) %>% 
#   rows_upsert(cf_aadna_new %>% select(aa,energy))
# 
# cf_aadna_2409_all <- bind_all_aadna(cf_aadna_2409_bind)
# 
# write_csv(cf_aadna_2409_bind, "cf_aadna_2409_bind.csv")
# write_csv(cf_aadna_2409_all, "cf_aadna_2409_all.csv")
