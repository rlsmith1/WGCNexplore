
#### RUN THIS SCRIPT FOR THE MODULE SETS OF INTEREST TO LOAD INTO SHINY!!!! ####


# libraries ---------------------------------------------------------------

        library(tidyverse)
        library(purrr)
        library(janitor)
        library(RColorBrewer)
        library(ggalluvial)



# set theme for plots -----------------------------------------------------


        theme_set(theme_bw() +
                          theme(plot.title = element_text(size = 20),
                                axis.title = element_text(size = 18),
                                axis.text = element_text(size = 15),
                                strip.text = element_text(size = 18),
                                legend.title = element_text(size = 18),
                                legend.text = element_text(size = 15)))


# ### SET DX HERE!!! ------------------------------------------------------

        dx <- "bipolar"
        

# load module data --------------------------------------------------------


        load(paste0("objects/2022-07-01_", dx, "_module_assignments.RDS")) # df_modules
        # note: 2022-06-30 for schizo



# select module sets of interest based on sft, min_size, and cut_height results --------


        df_modules <- df_modules %>% 
                filter(mod_set %in% c("sft8_minSize30_cutHeight0.2", "sft9_minSize30_cutHeight0.15",
                                      "sft9_minSize30_cutHeight0.2")) %>% 
                dplyr::rename("module" = "merged_module")


# plot module sizes -------------------------------------------------------

        
        # df_modules %>%
        #         group_by(mod_set) %>% 
        #         dplyr::count(module) %>% 
        #         
        #         ggplot(aes(x = module, y = n)) +
        #         geom_point(aes(color = mod_set, size = n)) +
        #         labs(y = "n_genes") +
        #         scale_x_discrete(labels = labels,
        #                          breaks = seq(from = 0, to = n_mods, by = 1)) +
        #         facet_wrap(~mod_set, scales = "free_x") +
        #         guides(color = "none") +
        #         ggtitle("Number of genes in each module per set")
        


# alluvial overlap plot --------------------------------------------------------

# 
#         # COLORS
#         qual_col_pals <- brewer.pal.info %>% 
#                 rownames_to_column("palette") %>% 
#                 filter(category == "qual", !grepl("Pastel", palette))
#         
#         col_vector <- c("#636363", unlist(mapply(brewer.pal, 
#                                                  (qual_col_pals$maxcolors - 1), 
#                                                  qual_col_pals$palette)))
#         df_colors <- tibble(module = levels(df_modules$module),
#                color = col_vector[1:length(levels(df_modules$module))]) %>% 
#                 mutate(module = factor(module, levels = 0:max(as.numeric(module))))
#         
#         # DF
#         df_alluvial <- df_modules %>% 
#                 dplyr::select(ensembl_gene_id, mod_set, module) %>% 
#                 left_join(df_colors, by = "module")
#         
#         # PLOT
#         df_alluvial %>% 
#                 
#                 ggplot(aes(x = mod_set, stratum = module, alluvium = ensembl_gene_id)) +
#                 stat_stratum() +
#                 geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#                 stat_flow(aes(fill = I(color))) +
#                 labs(x = "module set")
#         


# hypergeometric overlap test ---------------------------------------------

        
        # SOURCE FUNCTION
        source("functions/module_set_hypergeometric_overlap.R")

        # CALCULATE HYPERGEOMETRIC OVERLAP BETWEEN ALL MODULES IN EACH MODULE SET
        mod_sets <- df_modules %>% pull(mod_set) %>% unique
        mod_set_combos <- cross2(mod_sets, mod_sets, .filter = `==`)
        l_mod_set_combos <- seq_len(length(mod_set_combos)) %>% purrr::map(~mod_set_combos[[.x]] %>% unlist)
        
        df_mods_overlap <- tibble()
        for (sets in l_mod_set_combos)  {
                
                print(paste0("calculating overlap between ", sets[1], " and ", sets[2]))
                
                df_mod_sets <- df_modules %>% filter(mod_set %in% sets)
                
                modules1 <- df_mod_sets %>%
                        filter(mod_set == sets[1]) %>%
                        arrange(module) %>%
                        pull(module) %>%
                        unique
                
                modules2 <- df_mod_sets %>%
                        filter(mod_set == sets[2]) %>%
                        arrange(module) %>%
                        pull(module) %>%
                        unique
                
                df_tmp2 <- tibble()
                for (i in modules1) {
                        
                        print(paste0("calculating overlap for ", sets[1], " module ", i))
                        
                        df_tmp1 <- modules2 %>% map_dfr(~f_module_overlap(df = df_mod_sets, # makes sure merged_module is named "module"
                                                                          mod_set1 = sets[1],
                                                                          mod_set2 = sets[2],
                                                                          mod_no1 = i,
                                                                          mod_no2 = .x)) %>%
                                mutate(p_adj = p.adjust(p_val, method = "BH")) %>%
                                filter(p_adj < 0.05)
                        
                        df_tmp2 <- df_tmp2 %>% bind_rows(df_tmp1)
                        
                }
                
                df_mods_overlap <- df_mods_overlap %>% bind_rows(df_tmp2)
                
        }
        
        
        save(df_mods_overlap,
             file = paste0("objects/", Sys.Date(), "_", dx, "_module_overlap.RDS"))
        
                
        
# hypergeometric overlap matrix plot --------------------------------------

        # # SOURCE PLOT FUNCTION
        # source("functions/f_module_overlap_matrix_plot.R")
        # 
        # # CREATE DF FOR PLOT
        # df_overlap_plot <-
        #         df_mods_overlap %>%
        #         mutate(is_big = ifelse(q >= 50, 1 , 0) %>% factor()) %>%
        #         # mutate(p_adj = ifelse(p_adj > 0.999, NA, p_adj)) %>% # set p_adj = 1 to NA
        #         mutate(p_adj = ifelse(p_adj == 0, 1e-300, p_adj)) %>%
        #         mutate(mod_set1 = factor(mod_set1, levels = unique(.$mod_set1))) %>%
        #         mutate(mod1 = factor(mod1, levels = levels(df_mods_overlap$mod2)),
        #                mod2 = factor(mod2, levels = levels(df_mods_overlap$mod2))) %>%
        #         dplyr::select(-c(q, m, n, k))
        # 
        # # PLOT
        # f_plot_module_matrix(df_overlap_plot, "sft10_minSize30_cutHeight0.15", "sft10_minSize30_cutHeight0.2")
        


# cell-type enrichment ----------------------------------------------------

        # LOAD CELL-TYPE DATA
        load("objects/cell_type_data.RDS") # df_cell_type

        # CREATE DF
        df_modules_cell_types <- df_modules %>%
                dplyr::select(mod_set, ensembl_gene_id, module) %>%
                left_join(df_cell_type) %>%
                filter(!is.na(type) & grepl("lake", type)) %>% 
                mutate(type = str_remove(type, "_lake"))
        
        save(df_modules_cell_types,
             file = paste0("objects/", Sys.Date(), "_", dx, "_module_cell_type.RDS"))
        
        # # PLOT Z-SCORE OF CELL-TYPE ENRICHMENT
        # df_modules_cell_types %>% 
        #         group_by(mod_set, module) %>% 
        #         dplyr::count(type) %>% 
        #         mutate(z_score = (n - mean(n, na.rm = TRUE))/sd(n, na.rm = TRUE)) %>% 
        #         
        #         # filter(!grepl("synapses", type)) %>% 
        #         
        #         ggplot(aes(x = type, y = module,
        #                    fill = z_score)) +
        #         geom_tile(aes(width = 0.95, height = 0.95)) +
        #         scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-3, 3)) +
        #         facet_wrap(~mod_set, scales = "free", nrow = 2) +
        #         labs(x = "Lake cell type") + 
        #         coord_flip() +
        #         theme_classic()
        

        