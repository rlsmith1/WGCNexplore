

library(tidyverse)

setwd("/Users/rachelsmith/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/WGCNA/wgcna_app")

##
## load in and rename objects based on dx
##

shiny_objects <- paste0("objects/", list.files("objects"))

for (object in shiny_objects) {
        
        # gene expression matrix with residuals
        if (grepl("resids", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_resids_scz <- df_resids %>% mutate(dx = "schizo", .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_resids_bd <- df_resids %>% mutate(dx = "bipolar", .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_resids_mdd <- df_resids %>% mutate(dx = "mdd", .before = 1) 
                        
                }
                
        }
        
        # module assignments
        if (grepl("module_assignments", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_modules_scz <- df_modules %>% mutate(mod_set = paste0("schizo_", mod_set), .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_modules_bd <- df_modules %>% mutate(mod_set = paste0("bipolar_", mod_set), .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_modules_mdd <- df_modules %>% mutate(mod_set = paste0("mdd_", mod_set), .before = 1) 
                        
                }
                
        }
        
        # cell type enrichment
        if (grepl("cell_type", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_cell_type_scz <- df_modules_cell_types %>% mutate(mod_set = paste0("schizo_", mod_set), .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_cell_type_bd <- df_modules_cell_types %>% mutate(mod_set = paste0("bipolar_", mod_set), .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_cell_type_mdd <- df_modules_cell_types %>% mutate(mod_set = paste0("mdd_", mod_set), .before = 1) 
                        
                }
                
        }
        
        # hypergeometric overlap
        if (grepl("module_overlap", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_overlap_scz <- df_mods_overlap %>% mutate(mod_set1 = paste0("schizo_", mod_set1),
                                                                     mod_set2 = paste0("schizo_", mod_set2),
                                                                     .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_overlap_bd <- df_mods_overlap %>% mutate(mod_set1 = paste0("bipolar_", mod_set1),
                                                                     mod_set2 = paste0("bipolar_", mod_set2),
                                                                     .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_overlap_mdd <- df_mods_overlap %>% mutate(mod_set1 = paste0("mdd_", mod_set1),
                                                                     mod_set2 = paste0("mdd_", mod_set2),
                                                                     .before = 1) 
                        
                }
                
        }
        
        # module enrichment
        if (grepl("module_go_results", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_mods_go_scz <- df_mods_go %>% mutate(mod_set = paste0("schizo_", mod_set), .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_mods_go_bd <- df_mods_go %>% mutate(mod_set = paste0("bipolar_", mod_set), .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_mods_go_mdd <- df_mods_go %>% mutate(mod_set = paste0("mdd_", mod_set), .before = 1) 
                        
                }
                
        }
        
        # module enrichment bigrams
        if (grepl("module_go_bigrams", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_mods_go_bigrams_scz <- df_mods_go_bigrams %>% mutate(mod_set = paste0("schizo_", mod_set), .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_mods_go_bigrams_bd<- df_mods_go_bigrams %>% mutate(mod_set = paste0("bipolar_", mod_set), .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_mods_go_bigrams_mdd <- df_mods_go_bigrams %>% mutate(mod_set = paste0("mdd_", mod_set), .before = 1) 
                        
                }
                
        }
        
        # overlap enrichment
        if (grepl("overlap_go_results", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_mods_overlap_go_scz <- df_mods_overlap_go %>% mutate(mod_set1 = paste0("schizo_", mod_set1),
                                                                                mod_set2 = paste0("schizo_", mod_set2),
                                                                                .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_mods_overlap_go_bd <- df_mods_overlap_go %>% mutate(mod_set1 = paste0("bipolar_", mod_set1),
                                                                                mod_set2 = paste0("bipolar_", mod_set2),
                                                                                .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_mods_overlap_go_mdd <- df_mods_overlap_go %>% mutate(mod_set1 = paste0("mdd_", mod_set1),
                                                                                mod_set2 = paste0("mdd_", mod_set2),
                                                                                .before = 1) 
                        
                }
                
        }
        
        # overlap enrichment bigrams
        if (grepl("overlap_go_bigrams", object)) {
                
                if (grepl("schizo", object)) {
                        
                        load(object)
                        df_mods_overlap_go_bigrams_scz <- df_mods_overlap_go_bigrams %>% mutate(mod_set1 = paste0("schizo_", mod_set1),
                                                                                                mod_set2 = paste0("schizo_", mod_set2),
                                                                                                .before = 1) 
                        
                } else if (grepl("bipolar", object)) {
                        
                        load(object)
                        df_mods_overlap_go_bigrams_bd <- df_mods_overlap_go_bigrams %>% mutate(mod_set1 = paste0("bipolar_", mod_set1),
                                                                                                mod_set2 = paste0("bipolar_", mod_set2),
                                                                                                .before = 1) 
                        
                } else if (grepl("mdd", object)) {
                        
                        load(object)
                        df_mods_overlap_go_bigrams_mdd <- df_mods_overlap_go_bigrams %>% mutate(mod_set1 = paste0("mdd_", mod_set1),
                                                                                                mod_set2 = paste0("mdd_", mod_set2),
                                                                                                .before = 1) 
                        
                }
                
        }
        
        
}

##
## combine
##

df_resids <- df_resids_scz %>% bind_rows(df_resids_bd, df_resids_mdd) %>% distinct()
df_modules <- df_modules_scz %>% bind_rows(df_modules_bd, df_modules_mdd) %>% distinct() %>% 
        dplyr::rename("module" = "merged_module")
df_cell_type <- df_cell_type_scz %>% bind_rows(df_cell_type_bd, df_cell_type_mdd) %>% distinct()
df_overlap <- df_overlap_scz %>% bind_rows(df_overlap_bd, df_overlap_mdd) %>% distinct()
df_mods_go <- df_mods_go_scz %>% bind_rows(df_mods_go_bd, df_mods_go_mdd) %>% distinct()
df_mods_go_bigrams <- df_mods_go_bigrams_scz %>% bind_rows(df_mods_go_bigrams_bd, df_mods_go_bigrams_mdd) %>% distinct()
df_overlap_go <- df_mods_overlap_go_scz %>% bind_rows(df_mods_overlap_go_bd, df_mods_overlap_go_mdd) %>% distinct()
df_overlap_go_bigrams <- df_mods_overlap_go_bigrams_scz %>% bind_rows(df_mods_overlap_go_bigrams_bd, df_mods_overlap_go_bigrams_mdd) %>% distinct()

## 
## save objects for shiny
##

save(df_resids, 
     df_modules, df_cell_type, df_overlap,
     df_mods_go, df_mods_go_bigrams,
     df_overlap_go, df_overlap_go_bigrams,
     file = "shiny/wgcna_app/objects/2022-07-05_all_shiny_objects.RDS")



##
## export df_modules to .csv for shiny
##
load("shiny/wgcna_app/objects/2022-07-05_all_shiny_objects.RDS")
df_modules %>% 
        pivot_wider(id_cols = ensembl_gene_id, names_from = mod_set, values_from = module) %>% 
        write.csv("shiny/wgcna_module_set_comparison/wgcna_module_set_comparison/data/modules.csv",
                  row.names = FALSE)



