f_module_overlap <- function(df,
                             mod_set1, 
                             mod_set2, 
                             mod_no1, 
                             mod_no2) {
  
  # find intersect between our gene lists associated with each module set
  gene_list1 <- df %>% filter(mod_set == mod_set1 & !is.na(module)) %>% pull(ensembl_gene_id)
  gene_list2 <- df %>% filter(mod_set == mod_set2 & !is.na(module)) %>% pull(ensembl_gene_id)
  
  gene_universe_intersect <- intersect(gene_list1, gene_list2) 
  
  # pull genes in mo_no1 of mod_set1
  module1 <- df %>% 
    filter(mod_set == mod_set1) %>% 
    filter(module == mod_no1) %>% 
    pull(ensembl_gene_id)
  
  # pull genes in mo_no2 of mod_set2
  module2 <- df %>% 
    filter(mod_set == mod_set2) %>% 
    filter(module == mod_no2) %>% 
    pull(ensembl_gene_id)
  
  # find intersect of module1 with gene_universe_intersect genes
  module1_intersect <- intersect(gene_universe_intersect, module1)
  
  # find intersect of module2 with gene_universe_intersect genes
  module2_intersect <- intersect(gene_universe_intersect, module2)
  
  # find intersect of module intersects
  module1_module2_intersect <- intersect(module1_intersect, module2_intersect)
  
  # set variables for hypergeometric function
  q <- (module1_module2_intersect %>% length) - 1 # overlap between two modules - 1
  m <- module1_intersect %>% length # module 1 
  n <- (gene_universe_intersect %>% length) - m # pop size - module 1
  k <- module2_intersect %>% length # module 2
  
  # return tibble
  tibble(mod_set1 = mod_set1,
         mod_set2 = mod_set2,
         mod1 = mod_no1,
         mod2 = mod_no2,
         overlap = q + 1,
         q = q,
         m = m,
         n = n,
         k = k,
         p_val = phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE))
  
}

