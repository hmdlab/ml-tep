library(tidyverse)
library(caret)
library(ggpubr)

rep_z <- function(df, zscale){
  df %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_0")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_1")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_2")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_3")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_4"))))
}
rep_z2 <- function(df, zscale){
  df %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_01")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_12")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_23")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_34"))))
}
rep_z3 <- function(df, zscale){
  df %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_02")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_13")))) %>% 
    left_join(zscale %>% 
                rename_with(~(zscale %>% colnames() %>% paste0("_24"))))
}

add_1mer <- function(df){
  df %>% 
    mutate(aa_0=substr(aa,1,1)) %>% 
    mutate(aa_1=substr(aa,2,2)) %>% 
    mutate(aa_2=substr(aa,3,3)) %>% 
    mutate(aa_3=substr(aa,4,4)) %>% 
    mutate(aa_4=substr(aa,5,5))
}

add_2mer <- function(df){
  df %>% 
    mutate(aa_01=substr(aa,1,2)) %>% 
    mutate(aa_12=substr(aa,2,3)) %>% 
    mutate(aa_23=substr(aa,3,4)) %>% 
    mutate(aa_34=substr(aa,4,5))
}

add_3mer <- function(df){
  df %>% 
    mutate(aa_02=substr(aa,1,3)) %>% 
    mutate(aa_13=substr(aa,2,4)) %>% 
    mutate(aa_24=substr(aa,3,5))
}

add_3mer2 <- function(df){
  df %>% 
    mutate(aa_02=substr(aa,1,3)) %>% 
    mutate(aa_13=substr(aa,2,4)) %>% 
    mutate(aa_24=substr(aa,3,5)) %>% 
    mutate(aa_35=substr(aa,4,6)) %>% 
    mutate(aa_46=substr(aa,5,7))
}

t_addmean <- function(df, nlength, stchar){
  df %>% 
    left_join(df %>% 
                select(starts_with(stchar),starts_with("aa")) %>% 
                pivot_longer(-starts_with("aa")) %>% 
                mutate(nmean=substr(name,1,nlength)) %>% 
                group_by(aa,nmean) %>% 
                summarize(value=mean(value)) %>% 
                ungroup %>% 
                pivot_wider(everything(), names_from="nmean"))
}

bind_all_aadna <- function(aadna){
  aadna <- aadna %>% 
    select(aa)
  
  lib_all_z <- aadna %>% 
    add_1mer() %>% 
    rep_z(zscale) %>% 
    t_addmean(2, "z") %>% 
    select(-starts_with("aa"))
  
  lib_all_t <- aadna %>% 
    add_1mer() %>% 
    rep_z(tscale) %>% 
    t_addmean(2, "T") %>% 
    select(-starts_with("aa"))
  
  lib_all_st <- aadna %>% 
    add_1mer() %>% 
    rep_z(stscale) %>% 
    t_addmean(3, "ST") %>% 
    select(-starts_with("aa"))
  
  lib_all_vhse <- aadna %>% 
    add_1mer() %>% 
    rep_z(vhsescale) %>% 
    t_addmean(5, "VHSE") %>% 
    select(-starts_with("aa"))
  
  lib_all_all2 <- lib_all_z %>% 
    bind_cols(lib_all_t) %>% 
    bind_cols(lib_all_st) %>% 
    bind_cols(lib_all_vhse)
  
  lib_all_all2
}

scale_mer <- function(...){
  expand_grid(...) %>% 
    unnest(everything(),names_sep="_") %>% 
    unite(aa, ends_with("_aa"), sep="") %>% 
    pivot_longer(-aa,names_to=c("n","z"),names_sep="_") %>% 
    group_by(aa,z) %>% 
    summarize(value=mean(value)) %>% 
    pivot_wider(names_from=z, values_from=value) %>% 
    ungroup
}

auc_rf <- function(pred, cv=FALSE){
  pred$pred %>% 
    inner_join(pred$bestTune) %>% 
    {if(cv) group_by(.,Resample)
      else .} %>% 
    summarize(AUC=round(as.numeric(
      pROC::roc(obs, positive, levels=c("positive","negative"), direction=">")$auc),3)) %>% 
    {if(cv) unite(.,AUC, Resample, AUC, sep="; AUC=")
      else .} %>% 
    pull(AUC)
}

auc_plot <- function(pred, guide, title="AUC ="){
  auc_ans <- pred %>% 
    auc_rf(TRUE)
  aucmean_ans <- pred %>% 
    auc_rf(FALSE)
  pred$pred %>% 
    inner_join(pred$bestTune) %>% 
    group_by(Resample) %>% 
    group_split() %>% 
    map(~lift(obs ~ positive, data = ., class = "positive")$data) %>% 
    map2(unique(pred$pred$Resample), ~mutate(.x, fold=.y)) %>% 
    bind_rows() %>% 
    ggplot(aes(1 - Sp, Sn, color = fold)) + 
    geom_line(size=1) +
    coord_fixed() + 
    scale_color_discrete(guide = guide_legend(title = guide), labels=auc_ans) + 
    ggtitle(paste(title, aucmean_ans))
}

importance_plot <- function(pred, title, fill, sub1=1, sub2=2){
  fill <- enquo(fill)
  varImp(pred)$importance %>% 
    rownames_to_column("Feature") %>% 
    arrange(-Overall) %>% 
    slice_head(n=20) %>% 
    mutate(Feature=factor(Feature, rev(Feature))) %>% 
    dplyr::mutate(!!fill:=substr(Feature, sub1,sub2)) %>% 
    ggplot(aes(Overall, Feature, fill=!!fill)) +
    geom_bar(stat="identity") + 
    ggtitle(title) + 
    theme_gray(base_family = "Arial")
}

importance_plot_paper <- function(pred, title, fill, sub1=1, sub2=2){
  fill <- enquo(fill)
  varImp(pred)$importance %>% 
    rownames_to_column("Feature") %>% 
    arrange(-Overall) %>% 
    slice_head(n=10) %>% 
    mutate(Feature=factor(Feature, rev(Feature))) %>% 
    dplyr::mutate(!!fill:=substr(Feature, sub1,sub2)) %>% 
    ggplot(aes(Overall, Feature)) +
    geom_bar(stat="identity") + 
    ggtitle(title) + 
    theme_pubr(base_family = "Arial")
}

importance_plot_var <- function(pred_list, title, fill, sub1=1, sub2=2){
  fill <- enquo(fill)
  pred_list %>% 
    map(~(varImp(.)$importance)) %>% 
    map(rownames_to_column, "Feature") %>% 
    bind_rows() %>% 
    group_by(Feature) %>% 
    summarize(imp_mean=mean(Overall), imp_sd=sd(Overall)) %>% 
    arrange(-imp_mean) %>% 
    slice_head(n=20) %>% 
    mutate(Feature=factor(Feature, rev(Feature))) %>% 
    mutate(!!fill:=substr(Feature, sub1,sub2)) %>% 
    ggplot(aes(imp_mean, Feature, fill=!!fill)) +
    geom_bar(stat="identity") + 
    geom_errorbar(aes(xmin=imp_mean-imp_sd, xmax=imp_mean+imp_sd)) + 
    ggtitle(title) + 
    theme_pubr(base_family = "HiraKakuPro-W3")
}

cor_rf_plot <- function(pred){
  # corall <- pred$pred %>% 
  #   inner_join(pred$bestTune) %>% 
  #   summarize(cor= round(cor(pred, obs),3)) %>% 
  #   pull(cor) %>% 
  #   paste0("cor=",.)
  pred$pred %>% 
    inner_join(pred$bestTune) %>% 
    ggplot(aes(obs, pred)) +
    geom_point() +
    coord_fixed() +
    theme_pubr() +
    ggtitle("Overall")
}

predict_cv <- function(pred){
  pred$pred %>% 
    inner_join(pred$bestTune) %>% 
    arrange(rowIndex)
}

cor_rf_table <- function(pred){
  pred %>% 
    predict_cv() %>% 
    summarize(cor= round(cor(pred, obs),3)) %>% 
    bind_cols(
      pred$results %>% 
        inner_join(pred$bestTune)
    )
}

cor_rfcv_table <- function(pred){
  pred %>% 
    predict_cv %>% 
    group_by(Resample) %>% 
    summarize(cor= round(cor(pred, obs),3)) %>% 
    left_join(pred$resample)
}

cor_rfcv_plot <- function(pred, grid=TRUE){
  corcv <- pred %>% 
    cor_rfcv_table %>% 
    mutate(title=Resample) %>% 
    pull(title)
  xyscale <- pred$pred %>% 
    inner_join(pred$bestTune) %>% 
    select(obs, pred) %>% 
    summarise(across(everything(), list(min = min, max = max)))
  g <- pred$pred %>% 
    inner_join(pred$bestTune) %>% 
    group_by(Resample) %>% 
    group_split() %>% 
    map2(corcv,~{ggplot(.x, aes(obs, pred)) +
        geom_point() +
        coord_fixed() +
        theme_pubr() + 
        ggtitle(.y)+
        scale_x_continuous(limits=c(xyscale$obs_min, xyscale$obs_max))+
        scale_y_continuous(limits=c(xyscale$pred_min, xyscale$pred_max))})
  g <- c(g, list(cor_rf_plot(pred)))
  if(grid) cowplot::plot_grid(plotlist=g, ncol=3)
  else g
}

cor_rfcv_paperplot <- function(pred){
  predtable <- pred$pred %>% 
    inner_join(pred$bestTune)
  xyscale <- pred$pred %>% 
    inner_join(xgbcv_cfd$bestTune) %>% 
    select(obs, pred) %>% 
    summarise(across(everything(), list(min = min, max = max)))
  g1 <- predtable %>% 
    filter(Resample=="Fold1") %>% 
    ggplot(aes(obs,pred))+geom_point()+coord_fixed()+
    scale_x_continuous(limits=c(xyscale$obs_min, xyscale$obs_max))+
    scale_y_continuous(limits=c(xyscale$pred_min, xyscale$pred_max))
  g2 <- predtable %>% 
    filter(Resample=="Fold2") %>% 
    ggplot(aes(obs,pred))+geom_point()+coord_fixed()+
    scale_x_continuous(limits=c(xyscale$obs_min, xyscale$obs_max))+
    scale_y_continuous(limits=c(xyscale$pred_min, xyscale$pred_max))
  g3 <- predtable %>% 
    filter(Resample=="Fold3") %>% 
    ggplot(aes(obs,pred))+geom_point()+coord_fixed()+
    scale_x_continuous(limits=c(xyscale$obs_min, xyscale$obs_max))+
    scale_y_continuous(limits=c(xyscale$pred_min, xyscale$pred_max))
  g4 <- predtable %>% 
    filter(Resample=="Fold4") %>% 
    ggplot(aes(obs,pred))+geom_point()+coord_fixed()+
    scale_x_continuous(limits=c(xyscale$obs_min, xyscale$obs_max))+
    scale_y_continuous(limits=c(xyscale$pred_min, xyscale$pred_max))
  g5 <- predtable %>% 
    filter(Resample=="Fold5") %>% 
    ggplot(aes(obs,pred))+geom_point()+coord_fixed()+
    scale_x_continuous(limits=c(xyscale$obs_min, xyscale$obs_max))+
    scale_y_continuous(limits=c(xyscale$pred_min, xyscale$pred_max))
  g6 <- predtable %>% 
    ggplot(aes(obs,pred))+geom_point()+coord_fixed()+
    scale_x_continuous(limits=c(xyscale$obs_min, xyscale$obs_max))+
    scale_y_continuous(limits=c(xyscale$pred_min, xyscale$pred_max))
  g1 + g2 + g3 + g4 + g5 + g6 +
    plot_annotation(tag_levels = list("Fold1","Fold2","Fold3","Fold4","Fold5","Overall"))
}

tplot <- function(table, size=5){
  table %>% 
    ggplot() + 
    annotate(geom = 'table',
             x=0,
             y=0,
             size=size,
             label=list(table %>% 
                          mutate_if(is.numeric, round, 3))) + 
    theme_void()
}
cor_rfcv_plot2 <- function(pred, title=NULL){
  pred %>% 
    predict_cv %>% 
    ggplot(aes(obs, pred, color=Resample)) +
    geom_point() +
    coord_fixed() +
    ggtitle(title) + 
    theme_pubr(base_family = "HiraKakuPro-W3")
}
