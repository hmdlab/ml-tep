library(ggsankey)

seq_sankey <- function(seq){
  top_sep <- seq %>% 
    select(aa) %>% 
    separate(aa, c("0","1","2","3","4"), sep="") %>% 
    select(-`0`) %>% 
    mutate(`1`=paste0(`1`,"_1"),
           `2`=paste0(`2`,"_2"),
           `3`=paste0(`3`,"_3"),
           `4`=paste0(`4`,"_4")) %>% 
    mutate(id=row_number())
  
  top_list <- top_sep %>% 
    select(-id) %>% 
    pivot_longer(everything(), names_to="position", values_to="node") %>% 
    pull(node) %>% 
    table() %>% 
    sort(decreasing = T) %>% 
    names() %>% 
    factor(.,.)
  
  top_table <- top_sep %>% make_long(`1`,`2`,`3`,`4`) %>% 
    mutate(node=factor(node, top_list), next_node=factor(next_node, top_list))
  
  ggplot(top_table, aes(x = x, 
                        next_x = next_x, 
                        node = node, 
                        next_node = next_node,
                        fill = node,
                        label = substr(node,1,1))) +
    geom_sankey(flow.alpha = 0.5, node.color = 1) +
    geom_sankey_label(size = 3.5, color = 1, fill = "white") +
    scale_fill_viridis_d(option = "A", alpha = 0.95) +
    # scale_fill_continuous(type="viridis") + 
    theme_sankey(base_size = 16) +
    theme(legend.position = "none")
}
