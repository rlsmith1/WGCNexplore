

# libraries ---------------------------------------------------------------

    library(tidyverse)
    library(GO.db)
    library(topGO)
    library(biomaRt)



# function ----------------------------------------------------------------

#### function to run topGO on modules from WGCNA ####
## inputs: gene_module (character vector of ensembl_gene_id for genes in module of interest)
##        gene_universe (all genes in experiment after filtration)
##        df_gene_go (matrix of ensembl attribute go_id mapped to ensembl gene IDs, external gene name, and GO term (BP, CC, MF))
##        l_gene_2_GO (list of all ensembl gene IDs associated with each GO ID, where the name of each object in the list is the GO ID)
##        go_ontology (the GO ontology namespace of interest, on of "BP", "MF", or "CC")
## ouputs: a tibble containing the GO terms for that module with weightFisher < 0.05


f_top_GO_modules <- function(gene_module, 
                             gene_universe, 
                             df_gene_go, 
                             l_gene_2_GO, 
                             go_ontology = c("BP", "MF", "CC")) {
  
  # remove genes in list of interest (i.e., module) without GO annotation
  keep <- gene_module %in% (df_gene_go %>% pull(ensembl_gene_id) %>% unique)
  keep <- which(keep == TRUE)
  gene_module <- gene_module[keep]
  
  # make named factor indicating genes that are in module out of all the genes in the universe (lol not that universe)
  gene_module_universe <- as.integer(gene_universe %in% gene_module) %>% factor()
  names(gene_module_universe) <- gene_universe
  
  # create a topGO data object
  GOdata <- new("topGOdata",
                ontology = go_ontology,
                allGenes = gene_module_universe,
                annotationFun = annFUN.GO2genes,
                GO2genes = l_gene_2_GO,
                description = "GO analysis of module")
  
  # test for significance
  weight_fisher_result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  # generate results table
  all_GO <- usedGO(GOdata)
  all_res <- GenTable(GOdata,
                      weightFisher = weight_fisher_result,
                      orderBy = "weightFisher",
                      topNodes = length(all_GO)) %>% 
    as_tibble() %>% 
    mutate(weightFisher = as.numeric(weightFisher),
           weightFisher = ifelse(is.na(weightFisher), 0, weightFisher)) %>% 
    filter(weightFisher < 0.05) %>% 
    mutate(p_adj = p.adjust(weightFisher, method = "BH"))
  
  # return the final tibble
  return(all_res)
  
  
}

