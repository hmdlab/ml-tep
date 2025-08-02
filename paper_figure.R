library(tidyverse)
library(ggseqlogo)
library(patchwork)
library(ggpubr)
library(caret)
source("rf_script.R")
library(ggpmisc)
library(ggsankey)
source("sankey_script.R")
source("sf_features.R")
library(forcats)

windowsFonts(fig_family = windowsFont("Arial"))
fig_family <- "Arial"

energy2_clean <- read_csv("energy2_clean.csv")
aadna_all <- read_csv("aadna_all.csv")

# fig 2 -----

pep_neg_seq <- read_csv("peptides_neg.csv") %>% 
  mutate(peptides=substr(peptides,2,5))

pep_neg_seq_all <- pep_neg_seq %>% 
  uncount(count)

g2bb <- ggseqlogo(pep_neg_seq$peptides) + 
  theme(text = element_text(family = fig_family)) + 
  scale_x_continuous(name="Position")

lum_2.8_seq <- read_csv("lum_2.8ms_merge.csv") %>% 
  filter(nchar(aa)==5) %>% 
  mutate(aa=substr(aa,2,5))
lum_pos_2.8_skik <- lum_2.8_seq %>% 
  filter(lum_2.8ms>=86) %>% 
  group_by(aa) %>% 
  summarize(lum_2.8ms=mean(lum_2.8ms), .groups="keep") %>% 
  ungroup

g2ba <- ggseqlogo(lum_pos_2.8_skik$aa) + 
  theme(text = element_text(family = fig_family)) + 
  scale_x_continuous(name="Position")

g2ba + g2bb +
  plot_annotation(tag_levels = list("a","b"))

ggsave("svg/2b.svg", width=15, height=5)

# fig 5 -----

# 1回目
cf_aadna <- read_csv("cf_aadna.csv")
cf_aadna_all2 <- read_csv("cf_aadna_all2.csv")
# 2回目
cf_seqbind <- read_csv("sf_seqbind.csv")
cf_aadna_new <- read_csv("cf_aadna_new_all.csv")
# 3回目
cf_seqbind2 <- read_csv("sf_seqbind2.csv")
cf_aadna_new2 <- read_csv("cf_aadna_2409_all.csv")

load("rfcv_cfall2.Rdata")

g5a <- rfcv_cfall2 %>% 
  cor_rfcv_plot()
rfcv_cfall2 %>% 
  cor_rfcv_table
rfcv_cfall2$pred %>% 
  inner_join(rfcv_cfall2$bestTune) %>% 
  summarise(cor=cor(pred,obs), RMSE=RMSE(pred,obs))

# g5b <- rfcv_cfall2 %>% 
#   importance_plot("Feature Importance", group, 1,1) +
#   scale_fill_hue(labels = c("Energy", "ST-scale", "T-scale", "VHSE-scale", "Z-scale"))

g5b <- rfcv_cfall2 %>% 
  importance_plot_paper("Feature Importance", group, 1,1) +
  scale_fill_hue(labels = c("Energy", "ST-scale", "T-scale", "VHSE-scale", "Z-scale"))

lib_all_aadna <- energy2_clean %>% 
  left_join(aadna_all) %>% 
  anti_join(cf_aadna %>% 
              distinct(aa))

lib_all_energy <- lib_all_aadna %>% 
  bind_all_aadna()

rfcv_cfall2_pred <- lib_all_aadna %>% 
  select(-energy) %>% 
  mutate(pred=predict(rfcv_cfall2, lib_all_energy))

rfcv_cfall2_pred_top <- rfcv_cfall2_pred %>% 
  mutate(aa=substr(aa,2,5)) %>% 
  arrange(-pred) %>% 
  head(100)

g5c <- rfcv_cfall2_pred_top %>% 
  pull(aa) %>% 
  ggseqlogo::ggseqlogo() + 
  scale_x_continuous(name="Position")

g5d <- seq_sankey(rfcv_cfall2_pred_top) + 
  scale_x_discrete(name="Position")

g5a + g5b + g5c + g5d +
  plot_annotation(tag_levels = list("A","B","C","D"))

ggsave("svg/5.svg", width=15, height=9.5)
  
# Fig 6 -----

load("rfcv_cfd_new.Rdata")

g6a <- rfcv_cfd %>% 
  cor_rfcv_plot()
rfcv_cfd %>% 
  cor_rfcv_table
rfcv_cfd$pred %>% 
  inner_join(rfcv_cfd$bestTune) %>% 
  summarise(cor=cor(pred,obs), RMSE=RMSE(pred,obs))

# g6b <- rfcv_cfd %>% 
#   importance_plot("Feature Importance", group, 1,1) +
#   scale_fill_hue(labels = c("Energy", "ST-scale", "T-scale", "VHSE-scale", "Z-scale"))
g6b <- rfcv_cfd %>% 
  importance_plot_paper("Feature Importance", group, 1,1)

rfcv_cfd_pred <- lib_all_aadna %>% 
  select(-energy) %>% 
  mutate(pred=predict(rfcv_cfd, lib_all_energy))

rfcv_cfd_pred_top <- rfcv_cfd_pred %>% 
  mutate(aa=substr(aa,2,5)) %>% 
  arrange(-pred) %>% 
  head(100)

g6c <- rfcv_cfd_pred_top %>% 
  pull(aa) %>% 
  ggseqlogo::ggseqlogo() + 
  scale_x_continuous(name="Position")

g6d <- seq_sankey(rfcv_cfd_pred_top) + 
  scale_x_discrete(name="Position")

g6a + g6b + g6c + g6d +
  plot_annotation(tag_levels = list("A","B","C","D"))

ggsave("svg/6.svg", width=15, height=9.5)

# Fig 7 -----

lum_240911 <- read_csv("lum_240911.csv") %>% 
  mutate(group=factor(グループ,unique(グループ))) %>% 
  mutate(group=fct_recode(group, "toppred_rf"="RF予測上位", "toppred_XGB"="XGB予測上位",
                          "toppred_common"="共通予測上位", "random"="予測ランダム")) %>% 
  rename(predict_rf=RF予測値, measured=average)

g7b <- lum_240911 %>% 
  drop_na(predict_rf) %>% 
  ggplot(aes(x=predict_rf, y=measured, color=group)) + 
  geom_point(size=2) +
  coord_fixed() + 
  theme_pubr()

lum_240911 %>% 
  drop_na(predict_rf) %>% 
  summarize(cor=cor(predict_rf, measured))

g7b

ggsave("svg/7b.svg", width=7, height=5)

load("rfcv_cfd_new2.Rdata")

g7ca <- rfcv_cfd %>% 
  cor_rfcv_plot()
rfcv_cfd %>% 
  cor_rfcv_table
rfcv_cfd$pred %>% 
  inner_join(rfcv_cfd$bestTune) %>% 
  summarise(cor=cor(pred,obs), RMSE=RMSE(pred,obs))

g7cb <- rfcv_cfd %>% 
  importance_plot_paper("Feature Importance", group, 1,1)# +
  # scale_fill_hue(labels = c("Energy", "ST-scale", "T-scale", "VHSE-scale", "Z-scale"))

g7ca + g7cb +
  plot_annotation(tag_levels = list("a","b"))

ggsave("svg/7c.svg", width=15, height=5)


lib_all_aadna <- energy2_clean %>% 
  left_join(aadna_all) %>% 
  anti_join(cf_seqbind2 %>% 
              distinct(aa))

lib_all_energy <- lib_all_aadna %>% 
  bind_all_aadna()

rfcv_cfd_pred <- lib_all_aadna %>% 
  select(-energy) %>% 
  mutate(pred=predict(rfcv_cfd, lib_all_energy))

rfcv_cfd_pred_top <- rfcv_cfd_pred %>% 
  mutate(aa=substr(aa,2,5)) %>% 
  arrange(-pred) %>% 
  head(100)

g7e <- rfcv_cfd_pred_top %>% 
  pull(aa) %>% 
  ggseqlogo::ggseqlogo() + 
  scale_x_continuous(name="Position")

g7f <- seq_sankey(rfcv_cfd_pred_top) + 
  scale_x_discrete(name="Position")

g7ca + g7cb + g7e + g7f +
  plot_annotation(tag_levels = list("C","D","E","F"))

ggsave("svg/7new.svg", width=15, height=9.5)

# Fig S2 -----

load("xgbcv_cfall2.Rdata")

gs2a <- xgbcv_cfall2 %>% 
  cor_rfcv_plot()
xgbcv_cfall2 %>% 
  cor_rfcv_table
xgbcv_cfall2$pred %>% 
  inner_join(xgbcv_cfall2$bestTune) %>% 
  summarise(cor=cor(pred,obs), RMSE=RMSE(pred,obs))

gs2b <- xgbcv_cfall2 %>% 
  importance_plot_paper("Feature Importance", group, 1,1)# +
  #scale_fill_hue(labels = c("Energy", "ST-scale", "T-scale", "VHSE-scale", "Z-scale"))

xgbcv_cfall2_pred <- lib_all_aadna %>% 
  select(-energy) %>% 
  mutate(pred=predict(xgbcv_cfall2, lib_all_energy))

xgbcv_cfall2_pred_top <- xgbcv_cfall2_pred %>% 
  mutate(aa=substr(aa,2,5)) %>% 
  arrange(-pred) %>% 
  head(100)

gs2c <- xgbcv_cfall2_pred_top %>% 
  pull(aa) %>% 
  ggseqlogo::ggseqlogo() + 
  scale_x_continuous(name="Position")

gs2d <- seq_sankey(xgbcv_cfall2_pred_top) + 
  scale_x_discrete(name="Position")

gs2a + gs2b + gs2c + gs2d + 
  plot_annotation(tag_levels = list("A", "B", "C", "D"))

ggsave("svg/s2.svg", width=15, height=9.5)

# Fig S3 -----

load("xgbcv_cfd_new.Rdata")

gs3a <- xgbcv_cfd %>% 
  cor_rfcv_plot()
xgbcv_cfd %>% 
  cor_rfcv_table
xgbcv_cfd$pred %>% 
  inner_join(xgbcv_cfd$bestTune) %>% 
  summarise(cor=cor(pred,obs), RMSE=RMSE(pred,obs))

gs3b <- xgbcv_cfd %>% 
  importance_plot_paper("Feature Importance", group, 1,1)# +
  #scale_fill_hue(labels = c("Energy", "ST-scale", "T-scale", "VHSE-scale", "Z-scale"))

xgbcv_cfd_pred <- lib_all_aadna %>% 
  select(-energy) %>% 
  mutate(pred=predict(xgbcv_cfd, lib_all_energy))

xgbcv_cfd_pred_top <- xgbcv_cfd_pred %>% 
  mutate(aa=substr(aa,2,5)) %>% 
  arrange(-pred) %>% 
  head(100)

gs3c <- xgbcv_cfd_pred_top %>% 
  pull(aa) %>% 
  ggseqlogo::ggseqlogo() + 
  scale_x_continuous(name="Position")

gs3d <- seq_sankey(xgbcv_cfd_pred_top) + 
  scale_x_discrete(name="Position")

gs3a + gs3b + gs3c + gs3d + 
  plot_annotation(tag_levels = list("A", "B", "C", "D"))

ggsave("svg/s3.svg", width=15, height=9.5)

# Fig S4 -----

load("xgbcv_cfd_new2.Rdata")

gs4a <- xgbcv_cfd %>% 
  cor_rfcv_plot()
xgbcv_cfd %>% 
  cor_rfcv_table
xgbcv_cfd$pred %>% 
  inner_join(xgbcv_cfd$bestTune) %>% 
  summarise(cor=cor(pred,obs), RMSE=RMSE(pred,obs))

gs4b <- xgbcv_cfd %>% 
  importance_plot_paper("Feature Importance", group, 1,1)# +
  #scale_fill_hue(labels = c("Energy", "ST-scale", "T-scale", "VHSE-scale", "Z-scale"))

gs4a + gs4b + 
  plot_annotation(tag_levels = list("A", "B")) + 
  plot_layout(ncol=1)

ggsave("svg/s4.svg", width=8, height=9.5)


# Fig S5 -----

all_pred_bind <- read_csv("all_pred_bind.csv")

all_pred_bind %>% 
  # left_join(seq_candidate_bind) %>% 
  # slice_max(pred_mean, n=10000) %>%
  ggplot(aes(pred_rf, pred_xgb)) + 
  geom_density_2d_filled() + 
  theme_minimal()

ggsave("svg/s5.svg", width=6.5, height=5)
