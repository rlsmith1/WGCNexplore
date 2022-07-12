

#### RUN THIS SCRIPT FOR THE MODULE SETS OF INTEREST TO LOAD INTO SHINY!!!! ####


# libraries ---------------------------------------------------------------

        library(tidyverse)
        library(purrr)
        library(janitor)
        library(RColorBrewer)
        library(biomaRt)
        library(tidytext)
        library(ggwordcloud)



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
        
        
        
# select module sets of interest based on sft, min_size, and cut_height results --------
        
        
        df_modules <- df_modules %>% 
                filter(mod_set %in% c("sft8_minSize30_cutHeight0.2", "sft9_minSize30_cutHeight0.15",
                                      "sft9_minSize30_cutHeight0.2")) %>% 
                dplyr::rename("module" = "merged_module")

        
        
# load GO data from biomaRt -----------------------------------------------

     
        # PULL ALL WGCNA GENES FOR BACKGROUND
        gene_universe <- df_modules %>% pull(ensembl_gene_id) %>% unique
        
        # GET GENE INFO FROM ENSEMBL DATABASE
        hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                       dataset = "hsapiens_gene_ensembl",
                                       version = 88) # GRCh38.87

        # MAP EACH GENE TO GO TERM
        df_gene_go <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "namespace_1003", "name_1006"),
                            filters = "ensembl_gene_id",
                            values = gene_universe,
                            mart = hsapiens_ensembl) %>%
                as.data.frame() %>%
                as_tibble() %>%
                filter(go_id != "") %>%
                dplyr::rename("gene_ontology" = "namespace_1003",
                              "term" = "name_1006")

        # save(df_gene_go, file = "objects/gene_ontology_data.RDS")
        
        # BUILD ANNOTATION LIST
        load("objects/gene_ontology_data.RDS")
        
        l_gene_2_GO <- df_gene_go %>% 
                dplyr::select(ensembl_gene_id, go_id) %>% 
                unstack()


# functional enrichment on modules ----------------------------------------

        
        
        # RUN ENRICHMENT FUNCTION ON ALL MODULES IN EACH MODULE SET
        
                doParallel::registerDoParallel()
                
                mod_sets <- df_modules$mod_set %>% unique
                df_mods_go <- tibble()
                for (set in mod_sets) {
                        
                        print(paste0("determining functional enrichment for module set ", set))
                        
                        df_mod_set <- df_modules %>% filter(mod_set == set)
                        modules <- df_mod_set %>% filter(!is.na(module)) %>% arrange(module) %>% pull(module) %>% unique
                        n_mods <- length(modules) - 1
                        
                        df_tmp <- 0:n_mods %>%
                                map_dfr(~f_top_GO_modules(gene_module = df_mod_set %>%
                                                                  filter(module == .x) %>%
                                                                  pull(ensembl_gene_id),
                                                          gene_universe = gene_universe,
                                                          df_gene_go = df_gene_go,
                                                          l_gene_2_GO = l_gene_2_GO,
                                                          go_ontology = "BP") %>%
                                                mutate(module = .x, .before = 1)
                                ) %>%
                                clean_names() %>%
                                dplyr::select(-term) %>%
                                left_join(df_gene_go, by = "go_id") %>%
                                dplyr::select(module, go_id, term, everything()) %>%
                                dplyr::select(-gene_ontology) %>% 
                                distinct() %>%
                                mutate(mod_set = set, .before = 1)
                        
                        df_mods_go <- df_mods_go %>% bind_rows(df_tmp)
                        
                }
                
                # save
                save(df_mods_go,
                     file = paste0("objects/", Sys.Date(), "_", dx, "_module_go_results.RDS"))
                
        # WORD CLOUD DATA
                
                # remove stop words
                data("stop_words")
                
                # remove stop words
                c_stop_words <- stop_words$word
                collapsed_stop_words <- paste0("\\b", c_stop_words, "\\b", collapse = "|")
                
                df_mods_go_no_stop_words <- df_mods_go %>%
                        group_by(mod_set, module) %>%
                        mutate(text = str_remove_all(term, collapsed_stop_words) %>% trimws())
                
                # bigrams
                df_mods_go_bigrams <-  df_mods_go_no_stop_words %>%
                        unnest_tokens(ngram, text, token = "ngrams", n = 2) %>%
                        dplyr::count(ngram, sort = TRUE)
                
                save(df_mods_go_bigrams,
                     file = paste0("objects/", Sys.Date(), "_", dx, "_module_go_bigrams.RDS"))
                
                

# module GO enrichment plots -----------------------------------------------------


        # STRENGTH OF ENRICHMENT

                df_mods_go %>% filter(mod_set == "sft10_minSize30_cutHeight0.2") %>%

                        ggplot(aes(x = -log10(p_adj))) +
                        geom_histogram() +
                        facet_wrap(~module, scales = "free") +

                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
                        theme(strip.text = element_text(size = 8))

        # ENRICHMENT TERMS
                df_mods_go %>% filter(mod_set == "sft10_minSize30_cutHeight0.2") %>%

                        group_by(module) %>% arrange(module, p_adj) %>% top_n(n = 5) %>% # filter(module == 19) %>%
                        ggplot(aes(x = -log10(p_adj), y = reorder(str_wrap(term, width = 15), -p_adj))) +
                        geom_point(aes(size = significant/annotated, color = p_adj)) +
                        labs(y = "") +
                        facet_wrap(~module, scales = "free") +

                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +

                        theme(axis.text = element_text(size = 6),
                              strip.text = element_text(size = 6))

        # WORD CLOUDS

                df_mods_go_bigrams %>%
                        ungroup %>%
                        filter(!grepl("regulation|pathway|process|signal", ngram)) %>%
                        filter(mod_set == "sft10_minSize30_cutHeight0.2") %>%
                        filter(module != 0) %>%

                        group_by(module) %>%
                        mutate(size = n/max(n)) %>%
                        top_n(n = 10) %>% # filter(module == 12) %>%

                        mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(70, 30)),
                               angle = ifelse(row_number() <= 3, 0, angle)) %>%
                        mutate(module = factor(module, levels = 0:max(module))) %>%

                        ggplot(aes(label = ngram,
                                   size = size,
                                   angle = angle,
                                   x = module,
                                   color = size)) +
                        geom_text_wordcloud_area(rm_outside = TRUE) +
                        scale_size_area(max_size = 6) +
                        facet_wrap(~module, scales = "free") +
                        theme_minimal() +
                        theme(strip.text = element_blank())

                
                                

# functional enrichment on overlaps ---------------------------------------

                
        # RUN ENRICHMENT FUNCTION ON ALL OVERLAPS BETWEEN MODULES     
                
                df_mods_overlap_filt <- df_mods_overlap %>% filter(q > 50)
                df_mods_overlap_go <- tibble()
                
                # for (i in seq_len(nrow(df_mods_overlap_filt))) {
                for (i in 1:(nrow(df_mods_overlap_filt))) {
                        
                        print(paste0("determining functional enrichment for module set overlap ", i, " of ", nrow(df_mods_overlap_filt)))
                        
                        df_current_overlap <- df_mods_overlap_filt[i, ]
                        
                        mod_set1 <- df_current_overlap %>% pull(mod_set1)
                        mod_set2 <- df_current_overlap %>% pull(mod_set2)
                        mod1 <- df_current_overlap %>% pull(mod1)
                        mod2 <- df_current_overlap %>% pull(mod2)
                        
                        genes1 <- df_modules %>% filter(mod_set == mod_set1 & module == mod1) %>% pull(ensembl_gene_id)
                        genes2 <- df_modules %>% filter(mod_set == mod_set2 & module == mod2) %>% pull(ensembl_gene_id)
                        genes <- intersect(genes1, genes2)
                        
                        df_tmp <- f_top_GO_modules(gene_module = genes,
                                                   gene_universe = gene_universe,
                                                   df_gene_go = df_gene_go,
                                                   l_gene_2_GO = l_gene_2_GO,
                                                   go_ontology = "BP") %>%
                                clean_names() %>%
                                dplyr::select(-term) %>%
                                left_join(df_gene_go, by = "go_id") %>%
                                dplyr::select(-gene_ontology) %>% 
                                distinct() %>%
                                mutate(mod_set1 = mod_set1,
                                       mod_set2 = mod_set2,
                                       mod1 = mod1,
                                       mod2 = mod2,
                                       .before = 1)
                        
                        df_mods_overlap_go <- df_mods_overlap_go %>% bind_rows(df_tmp)
                        
                }
                
                # save
                save(df_mods_overlap_go,
                     file = paste0("objects/", Sys.Date(), "_", dx, "_overlap_go_results.RDS"))
                
                
        # WORD CLOUD DATA
                
                # remove stop words
                df_mods_overlap_go_no_stop_words <- df_mods_overlap_go %>%
                        group_by(mod_set1, mod_set2, mod1, mod2) %>%
                        mutate(text = str_remove_all(term, collapsed_stop_words) %>% trimws())
                
                # bigrams
                df_mods_overlap_go_bigrams <-  df_mods_overlap_go_no_stop_words %>%
                        unnest_tokens(ngram, text, token = "ngrams", n = 2) %>%
                        dplyr::count(ngram, sort = TRUE)
                
                save(df_mods_overlap_go_bigrams,
                     file = paste0("objects/", Sys.Date(), "_", dx, "_overlap_go_bigrams.RDS"))
                
        
                
# overlap GO enrichment plots -----------------------------------------------------
                
                
        # FORMAT DF FOR PLOTTING
                df_mods_overlap_go_grouped <- df_mods_overlap_go %>%
                        mutate(mod_set1 = factor(mod_set1, levels = unique(.$mod_set1)),
                               mod_set2 = factor(mod_set2, levels = unique(.$mod_set2)),
                               mod1 = factor(mod1, levels = levels(.$mod2))) %>%
                        group_by(mod_set1, mod_set2) %>%
                        mutate(group = cur_group_id(), .before = 1) %>%
                        ungroup

        # STRENGTH OF ENRICHMENT

                df_mods_overlap_go_grouped %>%
                        filter(group == 1 & mod1 == 1 & mod2 == 2) %>%

                        ggplot(aes(x = -log10(p_adj))) +
                        geom_histogram() +

                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
                        theme(strip.text = element_text(size = 8))

                # ENRICHMENT TERMS

                df_mods_overlap_go_grouped %>%

                        # select mod_sets and modules
                        # filter(mod_set1 == "sft10_minSize30_cutHeight0.15" & mod_set2 == "sft10_minSize30_cutHeight0.2") %>%
                        filter(mod1 == 1) %>%

                        # arrange to plot
                        group_by(mod2) %>%
                        arrange(mod1, mod2, p_adj) %>%
                        top_n(n = 10) %>%

                        ggplot(aes(x = -log10(p_adj), y = reorder(str_wrap(term, width = 15), -p_adj))) +
                        geom_point(aes(size = significant/annotated, color = p_adj)) +
                        labs(y = "") +
                        facet_wrap(~mod2, scales = "free") +

                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +

                        theme(axis.text = element_text(size = 6),
                              strip.text = element_text(size = 6))

                # WORD CLOUDS

                df_mods_overlap_go_bigrams %>%
                        ungroup %>%
                        filter(!grepl("regulation|pathway|process|signal", ngram)) %>%
                        # filter(mod_set1 == "sft10_minSize30_cutHeight0.15" & mod_set2 == "sft10_minSize30_cutHeight0.2") %>%
                        filter(mod1 == 1) %>%

                        group_by(mod2) %>%
                        mutate(size = n/max(n)) %>%
                        top_n(n = 10) %>%

                        mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(70, 30)),
                               angle = ifelse(row_number() <= 3, 0, angle)) %>%
                        # mutate(module = factor(module, levels = 0:max(module))) %>%

                        ggplot(aes(label = ngram,
                                   size = size,
                                   angle = angle,
                                   x = mod2,
                                   color = mod2)) +
                        geom_text_wordcloud_area(rm_outside = TRUE) +
                        scale_size_area(max_size = 6) +
                        facet_wrap(~mod2, scales = "free") +
                        theme_minimal() +
                        theme(strip.text = element_blank())

        
        
        