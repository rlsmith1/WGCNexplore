##
## Describe app here
##

# libraries ---------------------------------------------------------------

        library(shiny)
        library(tidyverse)
        library(purrr)
        library(janitor)
        library(ggalluvial)
        library(RColorBrewer)
        library(ggwordcloud)

# set theme for plots -----------------------------------------------------

        theme_set(theme_bw() +
                          theme(plot.title = element_text(size = 20),
                                axis.title = element_text(size = 18),
                                axis.text = element_text(size = 15),
                                strip.text = element_text(size = 18),
                                legend.title = element_text(size = 18),
                                legend.text = element_text(size = 15)))


# data --------------------------------------------------------------------

        # setwd("/Users/rachelsmith/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/WGCNA/wgcna_app/shiny/wgcna_app/")

        load("objects/2022-07-05_all_shiny_objects.RDS")
        # df_resids, 
        # df_modules, df_cell_type, df_overlap,
        # df_mods_go, df_mods_go_bigrams,
        # df_overlap_go, df_overlap_go_bigrams



# Shiny server --------------------------------------------------------

        
shinyServer(function(input, output) {
        
        
        ##
        ## SET UP REACTIVES
        ##
        
        # select module sets
        
        mod_set1 <- reactive({input$mod_set1})
        mod_set2 <- reactive({input$mod_set2})
        
        # create df_modules using these sets
        df_modules_reactive <- reactive({
                
                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                
                df_modules %>% 
                        filter(mod_set %in% c(mod_set1, mod_set2)) %>% 
                        mutate(mod_set = factor(mod_set, levels = c(mod_set1, mod_set2)))
                
        })
        
        # get cell-type data
        df_cell_type_reactive <- reactive({
                
                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                
                df_cell_type %>% 
                        filter(mod_set %in% c(mod_set1, mod_set2)) %>% 
                        mutate(mod_set = factor(mod_set, levels = c(mod_set1, mod_set2)))
                
        })
        
        # get hypergeometric overlap data
        df_overlap_reactive <- reactive({
                
                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                
                df_overlap %>% filter(mod_set1 == mod_set1 & mod_set2 == mod_set2)
                
        })
        
        # find maximum number of modules in either set
        n_mods <- reactive({
                
                df_modules <- df_modules_reactive()
                
                df_modules %>% 
                        pull(module) %>% 
                        as.character %>% 
                        as.numeric %>% 
                        max
                
        })
        
        # create color df for alluvial plots based on number of modules
        df_colors_reactive <- reactive({
                
                df_modules <- df_modules_reactive()
                
                qual_col_pals <- brewer.pal.info %>% 
                        rownames_to_column("palette") %>% 
                        filter(category == "qual", !grepl("Pastel", palette))
                
                col_vector <- c("#636363", unlist(mapply(brewer.pal, 
                                                         (qual_col_pals$maxcolors - 1), 
                                                         qual_col_pals$palette)))
                
                tibble(module = levels(df_modules$module),
                       color = col_vector[1:length(levels(df_modules$module))]) %>% 
                        mutate(module = factor(module, levels = 0:max(as.numeric(module))))
                
        })
        
        
        ##
        ## SIDEBAR 1 PLOTS
        ##
        
        
        # plot module sizes
        output$module_sizes <- renderPlot({
                
                df_modules <- df_modules_reactive()
                n_mods <- n_mods()
                labels <- 0:n_mods %>% ifelse(. %% 5 == 0, ., "")
                
                df_modules %>%
                        group_by(mod_set) %>% 
                        dplyr::count(module) %>% 
                        
                        ggplot(aes(x = module, y = n)) +
                        geom_point(aes(color = mod_set, size = n)) +
                        labs(y = "n_genes") +
                        scale_x_discrete(labels = labels,
                                         breaks = seq(from = 0, to = n_mods, by = 1)) +
                        facet_wrap(~mod_set, scales = "free_x") +
                        guides(color = "none") +
                        ggtitle("Number of genes in each module per set")
                
                
        }, height = 500, width = 900) %>% bindEvent(input$go1)
        
        
        # main alluvial plot
        output$main_alluvial <- renderPlot({
                
                df_modules <- df_modules_reactive()
                df_colors <- df_colors_reactive()
                
                df_alluvial <- df_modules %>% 
                        dplyr::select(ensembl_gene_id, mod_set, module) %>% 
                        left_join(df_colors, by = "module")
                
                # remove gray module?
                if (input$grey_mod == FALSE) {
                        
                        df_alluvial <- df_alluvial %>% filter(module != 0)
                        
                } else if (input$grey_mod == TRUE) {
                        
                        df_alluvial <- df_alluvial
                        
                }
                
                # plot
                df_alluvial %>% 
                        
                        ggplot(aes(x = mod_set, stratum = module, alluvium = ensembl_gene_id)) +
                        stat_stratum() +
                        geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                        stat_flow(aes(fill = I(color))) +
                        labs(x = "module set") +
                        ggtitle(paste0(mod_set1(), " intersects with ", mod_set2()))
                
                
        }, height = 600, width = 900) %>% bindEvent(input$go1)
        
        
        # overlap matrix plot
        output$overlap_matrix <- renderPlot({
                
                df_overlap <- df_overlap_reactive()
                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                
                df_overlap_plot <-
                        df_overlap %>%
                        # filter(mod_set1 == mod_set1 & mod_set2 == mod_set2) %>% 
                        mutate(is_big = ifelse(q >= 50, 1 , 0) %>% factor()) %>%
                        # mutate(p_adj = ifelse(p_adj > 0.999, NA, p_adj)) %>% # set p_adj = 1 to NA
                        mutate(p_adj = ifelse(p_adj == 0, 1e-300, p_adj)) %>%
                        mutate(mod_set1 = factor(mod_set1, levels = unique(.$mod_set1))) %>%
                        mutate(mod1 = factor(mod1, levels = levels(df_overlap$mod2)),
                               mod2 = factor(mod2, levels = levels(df_overlap$mod2))) %>%
                        dplyr::select(-c(q, m, n, k))
                
                # plot
                x_axis_scale <- df_overlap_plot %>% 
                        filter(mod_set1 == mod_set1) %>% 
                        arrange(mod1) %>% 
                        pull(mod1) %>% unique
                
                y_axis_scale <- df_overlap_plot %>% 
                        filter(mod_set2 == mod_set2) %>%
                        arrange(mod2) %>% 
                        pull(mod2) %>% unique
                
                df_overlap_plot %>% 
                        filter(mod_set1 == mod_set1 & mod_set2 == mod_set2) %>% 
                        
                        ggplot(aes(x = mod1, y = mod2)) +
                        geom_tile(aes(fill = -log10(p_adj))) +
                        geom_text(aes(label = overlap)) +
                        
                        scale_fill_gradientn(colors = brewer.pal(5, "Reds"), na.value = "white") +
                        scale_x_discrete(limits = x_axis_scale) +
                        scale_y_discrete(limits = y_axis_scale) +
                        labs(x = mod_set1, y = mod_set2) +
                        theme_classic() +
                        theme(plot.title = element_text(size = 18),
                              axis.title = element_text(size = 15),
                              axis.text = element_text(size = 12),
                              strip.text = element_text(size = 15),
                              legend.title = element_text(size = 15),
                              legend.text = element_text(size = 12))
                
                
        }, height = 600, width = 900) %>% bindEvent(input$go1)
        
        
        # module cell type enrichment
        output$cell_type <- renderPlot({
                
                n_mods <- n_mods()
                df_cell_type <- df_cell_type_reactive()
                
                df_cell_type %>% 
                        group_by(mod_set, module) %>%
                        dplyr::count(type) %>%
                        mutate(z_score = (n - mean(n, na.rm = TRUE))/sd(n, na.rm = TRUE)) %>%

                        # filter(!grepl("synapses", type)) %>%

                        ggplot(aes(x = type, y = module,
                                   fill = z_score)) +
                        geom_tile(aes(width = 0.95, height = 0.95)) +
                        scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-3, 3)) +
                        facet_wrap(~mod_set, scales = "free", nrow = 2) +
                        labs(x = "Lake cell type") +
                        coord_flip() +
                        ggtitle("Z-score of cell-type enrichment") +
                        theme_classic() +
                        theme(axis.title = element_text(size = 15),
                              axis.text = element_text(size = 12),
                              strip.text = element_text(size = 15),
                              legend.title = element_text(size = 15),
                              legend.text = element_text(size = 12))
                              
        }, height = 800, width = 900) %>% bindEvent(input$go1)
        
        ##
        ## SIDEBAR 2 PLOTS
        ##
        
        # alluvial zoom plot
        output$zoom_alluvial <- renderPlot({


                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                df_overlap <- df_overlap_reactive()

                df_colors <- df_colors_reactive()
                df_modules <- df_modules_reactive()

                # q > 50?
                if (input$q50 == FALSE) {

                        df_overlap <- df_overlap

                } else if (input$q50 == TRUE) {

                        df_overlap <- df_overlap %>% filter(q >= 50)

                }

                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", mod_set1, mod_set2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", mod_set2, mod_set1)

                m <- df_overlap %>%
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)

                df_alluvial_zoom <- df_modules %>%
                        filter((mod_set == mod_set_zoom & module == mod_zoom) |
                                       (mod_set == mod_set_other & module %in% m))

                # plot
                df_alluvial_zoom %>%
                        left_join(df_colors) %>%

                        ggplot(aes(x = mod_set, stratum = module, alluvium = ensembl_gene_id)) +
                        stat_stratum() +
                        geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
                        stat_flow(aes(fill = I(color))) +
                        labs(x = "module set") +
                        ggtitle(paste0(mod_set1, " intersects with ", mod_set2))


        }, height = 600, width = 900) %>% bindEvent(input$go2)

        
        # word cloud 1: enrichment of zoomed module in mod set
        output$wordcloud1 <- renderPlot({

                module_set <- ifelse(input$mod_set_zoom == "Module set 1",
                                     mod_set1(),
                                     mod_set2())

                mod_zoom <- input$mod_zoom

                df_wordcloud1 <- df_mods_go_bigrams %>% filter(mod_set == module_set)

                df_wordcloud1 %>%
                        ungroup %>%
                        filter(!grepl("regulation|pathway|process|signal", ngram)) %>%
                        filter(module == mod_zoom) %>%

                        group_by(module) %>%
                        mutate(size = n/max(n)) %>%
                        top_n(n = 10) %>%

                        mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(70, 30)),
                               angle = ifelse(row_number() <= 3, 0, angle)) %>%
                        mutate(module = factor(module, levels = 0:max(module))) %>%

                        ggplot(aes(label = ngram,
                                   size = size,
                                   angle = angle,
                                   x = module,
                                   color = size)) +
                        geom_text_wordcloud_area(rm_outside = TRUE) +
                        scale_size_area(max_size = 10) +
                        facet_wrap(~module, scales = "free") +
                        labs(x = "") +
                        ggtitle(paste0(module_set, " module ", mod_zoom)) +
                        theme_minimal() +
                        theme(strip.text = element_blank(),
                              plot.title = element_text(size = 18),
                              axis.text = element_blank())

        }, height = 500, width = 900) %>% bindEvent(input$go2)

        
        # histogram 1: strength of enrichment of zoomed module in mod set
        output$histogram1 <- renderPlot({

                module_set <- ifelse(input$mod_set_zoom == "Module set 1",
                                     mod_set1(),
                                     mod_set2())

                mod_zoom <- input$mod_zoom

                df_histogram1 <- df_mods_go %>% filter(mod_set == module_set & module == mod_zoom)

                df_histogram1 %>%

                        ggplot(aes(x = -log10(p_adj))) +
                        geom_histogram() +
                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
                        ggtitle(paste0("strength of ", module_set, " module ", mod_zoom, " enrichment "))

        }, height = 500, width = 900) %>% bindEvent(input$go2)

        
        # word cloud 2: enrichment of significant overlapping modules in other mod set
        output$wordcloud2 <- renderPlot({

                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                df_overlap <- df_overlap_reactive()

                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", mod_set1, mod_set2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", mod_set2, mod_set1)

                m <- df_overlap %>%
                        filter(q > 50) %>%
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)

                df_wordcloud2 <- df_mods_go_bigrams %>%
                        filter(mod_set == mod_set_other & module %in% m)

                # plot
                df_wordcloud2 %>%
                        ungroup %>%
                        filter(!grepl("regulation|pathway|process|signal", ngram)) %>%

                        group_by(module) %>%
                        mutate(size = n/max(n)) %>%
                        top_n(n = 10) %>%

                        mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(70, 30)),
                               angle = ifelse(row_number() <= 3, 0, angle)) %>%
                        mutate(module = factor(module, levels = 0:max(module))) %>%

                        ggplot(aes(label = ngram,
                                   size = size,
                                   angle = angle,
                                   x = module,
                                   color = module)) +
                        geom_text_wordcloud_area(rm_outside = TRUE) +
                        scale_size_area(max_size = 10) +
                        # facet_wrap(~module, scales = "free") +
                        ggtitle(paste0(ifelse(input$mod_set_zoom == "Module set 1",
                                              mod_set1,
                                              mod_set2), " intersecting modules enrichment")) +
                        theme_minimal() +
                        theme(plot.title = element_text(size = 20),
                              axis.text = element_text(size = 12),
                              axis.title = element_text(size = 15))

        }, height = 500, width = 900) %>% bindEvent(input$go2)

        # histogram 2: strength of enrichment of significant overlapping modules in other mod set
        output$histogram2 <- renderPlot({

                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                df_overlap <- df_overlap_reactive()

                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", mod_set1, mod_set2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", mod_set2, mod_set1)

                m <- df_overlap %>%
                        filter(q > 50) %>%
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)

                df_histogram2 <- df_mods_go %>%
                        filter(mod_set == mod_set_other & module %in% m)

                # plot
                df_histogram2 %>%

                        ggplot(aes(x = -log10(p_adj), fill = factor(module))) +
                        geom_histogram() +

                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +

                        guides(fill = "none") +
                        facet_wrap(~module) +
                        ggtitle(paste0("strength of ", ifelse(input$mod_set_zoom == "Module set 1",
                                                              mod_set1,
                                                              mod_set2), " intersecting modules enrichment strength"))

        }, height = 500, width = 900) %>% bindEvent(input$go2)

        # word cloud 3: enrichment of significant overlaps
        output$wordcloud3 <- renderPlot({

                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                df_overlap <- df_overlap_reactive()

                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", mod_set1, mod_set2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", mod_set2, mod_set1)

                m <- df_overlap %>%
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)

                df_wordcloud3 <- df_overlap_go_bigrams %>%
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other & mod2 %in% m)

                # plot
                df_wordcloud3 %>%
                        ungroup %>%
                        filter(!grepl("regulation|pathway|process|signal", ngram)) %>%
                        dplyr::rename("module" = "mod2") %>%

                        group_by(module) %>%
                        mutate(size = n/max(n)) %>%
                        top_n(n = 10) %>%

                        mutate(angle = 90 * sample(c(0, 1), n(), replace = TRUE, prob = c(70, 30)),
                               angle = ifelse(row_number() <= 3, 0, angle)) %>%

                        ggplot(aes(label = ngram,
                                   size = size,
                                   angle = angle,
                                   x = module,
                                   color = module)) +
                        geom_text_wordcloud_area(rm_outside = TRUE) +
                        scale_size_area(max_size = 10) +
                        # facet_wrap(~module, scales = "free") +
                        ggtitle(paste0(ifelse(input$mod_set_zoom == "Module set 1",
                                              mod_set1,
                                              mod_set2), " module OVERLAPS enrichment")) +
                        theme_minimal() +
                        theme(plot.title = element_text(size = 20),
                              axis.text = element_text(size = 12),
                              axis.title = element_text(size = 15))

        }, height = 500, width = 900) %>% bindEvent(input$go2)

        # histogram 3: strength of enrichment of significant overlaps
        output$histogram3 <- renderPlot({

                mod_set1 <- mod_set1()
                mod_set2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                df_overlap <- df_overlap_reactive()

                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", mod_set1, mod_set2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", mod_set2, mod_set1)

                m <- df_overlap %>%
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)

                df_histogram3 <- df_overlap_go %>%
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other & mod2 %in% m)

                # plot
                df_histogram3 %>%
                        dplyr::rename("module" = "mod2") %>%

                        ggplot(aes(x = -log10(p_adj), fill = factor(module))) +
                        geom_histogram() +

                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +

                        guides(fill = "none") +
                        facet_wrap(~module) +
                        ggtitle(paste0("strength of ", ifelse(input$mod_set_zoom == "Module set 1",
                                                              mod_set1,
                                                              mod_set2), " module OVERLAPS enrichment strength"))

        }, height = 500, width = 900) %>% bindEvent(input$go2)
        

})        
        


        
        
        
        
        
        
        
        

        





