library(tidyverse)
library(caret)
source("rf_script.R")
# source("sf_features.R")
library(ggpmisc)
library(ggsankey)
source("sankey_script.R")

# 1回目
cf_aadna <- read_csv("cf_aadna.csv")
cf_aadna_all2 <- read_csv("cf_aadna_all2.csv")
# 2回目
cf_aadna_new <- read_csv("cf_aadna_new_all.csv")
# 3回目
cf_aadna_new2 <- read_csv("cf_aadna_2409_all.csv")

ntree <- 1000
seed <- 1

# rf_all -----

tune_grid_rf <- expand.grid(.mtry=seq(5,75,5))

print("all")
set.seed(seed)

rfcv_cfall2 <- train(
  mean_lum ~ .,
  data = cf_aadna_all2,
  method = "rf",
  ntree=ntree,
  tuneGrid = tune_grid_rf,
  trControl = trainControl(method = "cv", 
                           number = 5,
                           savePredictions = T)
)
save(rfcv_cfall2, file="rfcv_cfall2.Rdata")

rfcv_cfall2 %>% cor_rfcv_table() %>% 
  ggplot(aes(Resample, RMSE)) + 
  geom_point() + 
  ggtitle(paste("seed =", seed))

rfcv_cfall2$bestTune
plot(rfcv_cfall2)
plot(rfcv_cfall2$finalModel)

rfcv_cfall2 %>% 
  importance_plot("cf1回目_RF(回帰)変数重要度", group, 1,1)

rfcv_cfall2 %>% 
  cor_rf_plot()
rfcv_cfall2 %>% 
  cor_rfcv_plot()
rfcv_cfall2 %>% 
  cor_rfcv_plot2("cfnew_RF(回帰)予測")
rfcv_cfall2 %>% 
  cor_rfcv_table
rfcv_cfall2$pred %>% 
  inner_join(rfcv_cfall2$bestTune) %>% 
  summarise(rmse=RMSE(pred,obs))

# pred -----

# lib_all_all2 <- read_csv("lib_all_all2.csv")
# lib_all_aa <- read_csv("lib_all_aa.csv")
energy2_clean <- read_csv("energy2_clean.csv")
aadna_all <- read_csv("aadna_all.csv")
# sf_seqbind <- read_csv("sf_seqbind2.csv")

# lib_all_aadna <- lib_all_aa %>% 
#   left_join(energy2_clean %>% 
#               left_join(aadna_all))

lib_all_aadna <- energy2_clean %>% 
  left_join(aadna_all) %>% 
  anti_join(cf_aadna %>% 
              distinct(aa))

# lib_all_energy <- lib_all_aadna %>% 
#   select(energy) %>% 
#   bind_cols(lib_all_all2)

source("sf_features.R")

lib_all_energy <- lib_all_aadna %>% 
  bind_all_aadna()

rfcv_cfall2_pred <- lib_all_aadna %>% 
  select(-energy) %>% 
  mutate(pred=predict(rfcv_cfall2, lib_all_energy))

rfcv_cfall2_pred_top <- rfcv_cfall2_pred %>% 
  mutate(aa=substr(aa,2,5)) %>% 
  arrange(-pred) %>% 
  head(100)

rfcv_cfall2_pred_top %>% 
  pull(aa) %>% 
  ggseqlogo::ggseqlogo() + 
  ggtitle("候補上位100配列")

seq_sankey(rfcv_cfall2_pred_top)

set.seed(3)
top_cluster_rf <- rfcv_cfd_all_pred %>% 
  arrange(-pred) %>% 
  head(5000) %>% 
  select(aa) %>% 
  separate(aa, c("blank","zero","one","two","three","four"), sep="") %>% 
  select(-"blank") %>% 
  fastDummies::dummy_cols(remove_selected_columns=TRUE) %>% 
  kmeans(5)

seq_candidate_rfnew <- rfcv_cfd_all_pred %>% 
  arrange(-pred) %>% 
  head(5000) %>% 
  mutate(clust=top_cluster_rf$cluster) %>% 
  group_by(clust) %>% 
  slice_max(pred, n=10) %>% 
  ungroup %>% 
  mutate(pred=round(pred,3)) %>% 
  rename(pred_lum=pred)

seq_candidate_rfnew %>% 
  ggplot(aes(clust, pred_lum, group=clust)) + 
  geom_boxplot() +
  labs(y="pred_rm")

write_csv(seq_candidate_rfnew, "seq_candidate_rfnew.csv")

seq_sankey(seq_candidate_rfnew %>% 
             mutate(aa=substr(aa,2,5)))
seq_sankey(seq_candidate_rfnew %>% 
             mutate(aa=substr(aa,2,5)) %>% 
             filter(clust==5))



# xgb_all -----

tune_grid_xgb <- expand.grid(
  nrounds = seq(from = 50, to = 250, by = 25),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(1, 2, 3, 4, 5, 6),
  gamma = 0,
  colsample_bytree = c(0.4, 0.7, 1),
  min_child_weight = 1,
  subsample = 1
)
tuneplot <- function(x, probs = .90) {
  ggplot(x) +
    coord_cartesian(ylim = c(quantile(x$results$RMSE, probs = probs), min(x$results$RMSE))) +
    theme_bw()
}


print("all")
set.seed(seed)

xgbcv_cfall2 <- train(
  mean_lum ~ .,
  data = cf_aadna_all2,
  method = "xgbTree",
  tuneGrid = tune_grid_xgb,
  trControl = trainControl(method = "cv", 
                           number = 5,
                           savePredictions = T)
)
save(xgbcv_cfall2, file="xgbcv_cfall2.Rdata")

xgbcv_cfall2 %>% cor_rfcv_table() %>% 
  ggplot(aes(Resample, RMSE)) + 
  geom_point() + 
  ggtitle(paste("seed =", seed))

tuneplot(xgbcv_cfall2)
xgbcv_cfall2$bestTune
plot(xgbcv_cfall2)

xgbcv_cfall2 %>% 
  importance_plot("cf1回目_xgb(回帰)変数重要度", group, 1,1)

xgbcv_cfall2 %>% 
  cor_rf_plot()
xgbcv_cfall2 %>% 
  cor_rfcv_plot()
xgbcv_cfall2 %>% 
  cor_rfcv_plot2("cfnew_xgb(回帰)予測")
xgbcv_cfall2 %>% 
  cor_rfcv_table
xgbcv_cfall2$pred %>% 
  inner_join(xgbcv_cfd$bestTune) %>% 
  summarise(rmse=RMSE(pred,obs))

# pred -----

# lib_all_all2 <- read_csv("lib_all_all2.csv")
# lib_all_aa <- read_csv("lib_all_aa.csv")
energy2_clean <- read_csv("energy2_clean.csv")
aadna_all <- read_csv("aadna_all.csv")
# sf_seqbind <- read_csv("sf_seqbind.csv")

# lib_all_aadna <- lib_all_aa %>% 
#   left_join(energy2_clean %>% 
#               left_join(aadna_all))

lib_all_aadna <- energy2_clean %>% 
  left_join(aadna_all) %>% 
  anti_join(cf_aadna %>% 
              distinct(aa))

# lib_all_energy <- lib_all_aadna %>% 
#   select(energy) %>% 
#   bind_cols(lib_all_all2)

lib_all_energy <- lib_all_aadna %>% 
  bind_all_aadna()# %>% 
  # bind_cols(lib_all_aadna %>% select(energy))

xgbcv_cfall2_pred <- lib_all_aadna %>% 
  select(-energy) %>% 
  mutate(pred=predict(xgbcv_cfall2, lib_all_energy))

xgbcv_cfall2_pred_top <- xgbcv_cfall2_pred %>% 
  mutate(aa=substr(aa,2,5)) %>% 
  arrange(-pred) %>% 
  head(100)

xgbcv_cfall2_pred_top %>% 
  pull(aa) %>% 
  ggseqlogo::ggseqlogo() + 
  ggtitle("候補上位100配列")

seq_sankey(xgbcv_cfall2_pred_top)

set.seed(3)
top_cluster_xgb <- xgbcv_cfd_all_pred %>% 
  arrange(-pred) %>% 
  head(5000) %>% 
  select(aa) %>% 
  separate(aa, c("blank","zero","one","two","three","four"), sep="") %>% 
  select(-"blank") %>% 
  fastDummies::dummy_cols(remove_selected_columns=TRUE) %>% 
  kmeans(5)

seq_candidate_xgbnew <- xgbcv_cfd_all_pred %>% 
  arrange(-pred) %>% 
  head(5000) %>% 
  mutate(clust=top_cluster_xgb$cluster) %>% 
  group_by(clust) %>% 
  slice_max(pred, n=10) %>% 
  ungroup %>% 
  mutate(pred=round(pred,3)) %>% 
  rename(pred_lum=pred)

seq_candidate_xgbnew %>% 
  ggplot(aes(clust, pred_lum, group=clust)) + 
  geom_boxplot() +
  labs(y="pred_xgb")

write_csv(seq_candidate_xgbnew, "seq_candidate_xgbnew.csv")

seq_sankey(seq_candidate_xgbnew %>% 
             mutate(aa=substr(aa,2,5)))
seq_sankey(seq_candidate_xgbnew %>% 
             mutate(aa=substr(aa,2,5)) %>% 
             filter(clust==5))

# Bind -----

all_pred_bind <- rfcv_cfd_all_pred %>% 
  rename(pred_rf=pred) %>% 
  left_join(xgbcv_cfd_all_pred %>% 
              rename(pred_xgb=pred)) %>% 
  mutate(pred_mean=(pred_rf+pred_xgb)/2)

write_csv(all_pred_bind, "all_pred_bind.csv")

all_pred_bind %>% 
  slice_max(pred_mean, n=1000) %>%
  ggplot(aes(pred_rf, pred_xgb)) + 
  geom_point()

seq_candidate_bind <- seq_candidate_rfnew %>% 
  rename(pred_rf=pred_lum, clust_rf=clust) %>% 
  full_join(seq_candidate_xgbnew %>% 
              rename(pred_xgb=pred_lum, clust_xgb=clust)) %>% 
  mutate(clust=ifelse(!is.na(clust_rf),ifelse(!is.na(clust_xgb),"rf_xgb","rf"),ifelse(!is.na(clust_xgb),"xgb",NA))) %>% 
  select(dna, aa, clust)

all_pred_bind %>% 
  left_join(seq_candidate_bind) %>% 
  slice_max(pred_mean, n=5000) %>%
  ggplot(aes(pred_rf, pred_xgb, color=clust)) + 
  geom_point()

all_pred_bind %>% 
  left_join(seq_candidate_bind) %>% 
  # slice_max(pred_mean, n=10000) %>%
  ggplot(aes(pred_rf, pred_xgb)) + 
  geom_point()

all_pred_bind %>% 
  left_join(seq_candidate_bind) %>% 
  # slice_max(pred_mean, n=10000) %>%
  ggplot(aes(pred_rf, pred_xgb)) + 
  geom_density_2d_filled()

all_pred_bind %>% 
  inner_join(seq_candidate_bind) %>% 
  filter(clust=="rf_xgb") %>% 
  left_join(seq_candidate_rfnew %>% 
              select(-pred_lum) %>% 
              rename(clust_rf=clust)) %>% 
  left_join(seq_candidate_xgbnew %>% 
              select(-pred_lum) %>% 
              rename(clust_xgb=clust)) %>% 
  select(-pred_mean, -clust)

seq_candidate_rfnew %>% 
  inner_join(all_pred_bind %>% 
               inner_join(seq_candidate_bind) %>% 
               filter(clust=="rf_xgb") %>% 
               select(-clust))
seq_candidate_xgbnew %>% 
  inner_join(all_pred_bind %>% 
               inner_join(seq_candidate_bind) %>% 
               filter(clust=="rf_xgb") %>% 
               select(-clust))
all_pred_bind %>% 
  left_join(seq_candidate_rfnew %>% 
              mutate(clust=factor(clust))) %>% 
  slice_max(pred_mean, n=5000) %>%
  ggplot(aes(pred_rf, pred_xgb, color=clust)) + 
  geom_point()
all_pred_bind %>% 
  left_join(seq_candidate_xgbnew %>% 
              mutate(clust=factor(clust))) %>% 
  slice_max(pred_mean, n=5000) %>%
  ggplot(aes(pred_rf, pred_xgb, color=clust)) + 
  geom_point()
