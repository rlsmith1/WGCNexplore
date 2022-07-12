##
## describe application here
##


# required libraries ------------------------------------------------------

require(shiny)
require(tidyverse)
require(janitor)
require(RColorBrewer)
require(ggalluvial)
require(purrr)
require(biomaRt)
require(tidytext)


# source functions --------------------------------------------------------

source("functions/module_set_hypergeometric_overlap.R")
source("functions/f_top_GO_modules.R")

# data --------------------------------------------------------------------


# setwd("/Users/rachelsmith/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/WGCNA/wgcna_app/shiny/wgcna_module_set_comparison/wgcna_module_set_comparison/")

## read in data/modules.csv
df_modules <- read.csv("data/modules.csv") %>% 
        as_tibble() %>% 
        mutate_if(is.numeric, ~factor(.x, levels = 0:max(.x)))

## read in cell-type data (data/cell_types.csv)
df_cell_type <- read.csv("data/cell_type.csv") %>% 
        as_tibble()

## if data is in wide format, pivot to long with column (module set) names going to column "mod_set"
if (ncol(df_modules) > 3) {
        
        df_modules <- df_modules %>%
                pivot_longer(2:ncol(.), names_to = "mod_set", values_to = "module")
        
} else {df_modules <- df_modules}

# set theme for plots -----------------------------------------------------

theme_set(theme_bw() +
                  theme(plot.title = element_text(size = 20),
                        axis.title = element_text(size = 18),
                        axis.text = element_text(size = 15),
                        strip.text = element_text(size = 18),
                        legend.title = element_text(size = 18),
                        legend.text = element_text(size = 15)))

# UI ----------------------------------------------------------------------


ui <- fluidPage(
        
        ##
        ## APPLICATION TITLE
        ##
        
        titlePanel("WGCNA module set comparison"),
        
        ##
        ## SIDEBAR LAYOUT
        ##
        tabsetPanel(
                
                ## TAB 1
                tabPanel("Module sets overview", fluid = TRUE,
                         
                         ## sidebar 1
                         sidebarLayout(
                                 
                                 sidebarPanel(
                                         numericInput(inputId = "n_mod_sets",
                                                      label = "Number of module sets to compare:",
                                                      value = 2,
                                                      min = 2,
                                                      step = 0.5
                                         ),
                                         actionButton("go0", "GO"),
                                         width = 3
                                 ),
                                 
                                 mainPanel = NULL
                                 
                         ),
                         
                         ## sidebar 2
                         sidebarLayout(
                                 
                                 sidebarPanel(
                                         selectInput(inputId = "mod_set1",
                                                     label = "Module set 1",
                                                     choices = df_modules %>% pull(mod_set) %>% unique(),
                                                     selected = NULL),
                                         selectInput(inputId = "mod_set2",
                                                     label = "Module set 2",
                                                     choices = df_modules %>% pull(mod_set) %>% unique(),
                                                     selected = NULL),
                                         checkboxInput(inputId = "grey_mod",
                                                       label = "Show grey?",
                                                       value = TRUE),
                                         actionButton(inputId = "go1",
                                                      label = "GO"),
                                         width = 3
                                 ),
                                 
                                 mainPanel(
                                         
                                         tabsetPanel(
                                                 id = "tabetPanelID",
                                                 type = "tabs",
                                                 
                                                 tabPanel("Module number & size", plotOutput("module_sizes")),
                                                 tabPanel("Module overlap", 
                                                          plotOutput("main_alluvial", height = 600, width = 900),
                                                          plotOutput("overlap_matrix", height = 600, width = 900),
                                                 ),
                                                 tabPanel("Cell-type enrichment",
                                                          plotOutput("cell_type", height = 800, width = 900))
                                                 
                                         )
                                         
                                 )
                                 
                         )
                ),
                
                ## TAB 2
                tabPanel("Module zoom", fluid = TRUE,
                         
                         ## sidebar 3
                         sidebarLayout(
                                 
                                 sidebarPanel(
                                         
                                         selectInput(inputId = "mod_set_zoom",
                                                     label = "Module set to zoom in on",
                                                     choices = c("Module set 1", "Module set 2"),
                                                     selected = NULL,
                                                     multiple = FALSE),
                                         selectInput(inputId = "mod_zoom",
                                                     label = "Module to zoom in on",
                                                     choices = levels(df_modules$module)),
                                         checkboxInput(inputId = "q50",
                                                       label = "overlap >= 50 genes",
                                                       value = FALSE),
                                         actionButton(inputId = "go2",
                                                      label = "GO"),
                                         width = 3
                                 ),
                                 
                                 mainPanel(
                                         
                                         tabsetPanel(
                                                 id = "tabetPanelID",
                                                 type = "tabs",
                                                 
                                                 tabPanel("Zoomed module overlaps", 
                                                          plotOutput("zoom_alluvial", height = 600, width = 900)),
                                                 tabPanel("Zoomed module functional enrichment",
                                                          plotOutput("zoom_module_go", height = 800, width = 900)),
                                                 tabPanel("Overlapping modules functional enrichment",
                                                          plotOutput("zoom_module_overlaps_go", height = 800, width = 900)),
                                                 tabPanel("Overlaps only functional enrichment",
                                                          plotOutput("intersects_go", height = 800, width = 900))  
                                                 
                                         )
                                         
                                 )
                         )
                         
                )
                
        )
        
)




# Server ------------------------------------------------------------------


server <- function(input, output) {
        
        ##
        ## UPDATE SIDEBAR 2 BASED ON SIDEBAR 1 (NUMBER OF MOD_SETS) INPUT
        ##
        
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
        
        # pull gene universe (all genes in all module sets)
        gene_universe_reactive <- reactive({
                
                df_modules <- df_modules_reactive()
                df_modules %>% pull(ensembl_gene_id) %>% unique 
                
        })
        
        
        
        ##
        ## CHECK IF MODULE SET HAS BEEN RUN FOR HYPERGEOMETRIC OVERLAP, CELL-TYPE ENRICHMENT, AND FUNCTIONAL ENRICHMENT
        ##
        
        data <- reactiveValues()
        
        data$df_overlap <- tibble(mod_set1 = as.character(),
                                  mod_set2 = as.character(),
                                  mod1 = factor(),
                                  mod2 = factor(),
                                  overlap = as.numeric(),
                                  q = as.numeric(),
                                  m = as.numeric(),
                                  n = as.numeric(),
                                  k = as.numeric(),
                                  p_val = as.numeric(),
                                  p_adj = as.numeric())
        
        data$df_gene_go <- tibble(ensembl_gene_id = as.character(),
                                  external_gene_name = as.character(),
                                  go_id = as.character(),
                                  term = as.character())
        
        data$l_gene_2_GO <- list()
        
        data$df_mods_go <- tibble(mod_set = as.character(),
                                  module = as.numeric(),
                                  go_id = as.character(),
                                  term = as.character(),
                                  annotated = as.numeric(),
                                  significant = as.numeric(),
                                  expected = as.numeric(),
                                  weight_fisher = as.numeric(),
                                  p_adj = as.numeric())
        
        data$df_overlaps_go <- tibble(mod_set1 = as.character(),
                                      mod_set2 = as.character(),
                                      mod1 = as.numeric(),
                                      mod2 = as.numeric(),
                                      go_id = as.character(),
                                      term = as.character(),
                                      annotated = as.numeric(),
                                      significant = as.numeric(),
                                      expected = as.numeric(),
                                      weight_fisher = as.numeric(),
                                      p_adj = as.numeric())
        
        # If module overlap has not been calculated, run hypergeometric tests here
        observeEvent(input$go1, {
                
                df_modules <- df_modules_reactive()
                
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                
                modules1 <- df_modules %>%
                        filter(mod_set == modSet1) %>%
                        arrange(module) %>%
                        pull(module) %>%
                        unique
                
                modules2 <- df_modules %>%
                        filter(mod_set == modSet2) %>%
                        arrange(module) %>%
                        pull(module) %>%
                        unique
                
                if (
                        (nrow(data$df_overlap) == 0) |
                        (nrow(filter(data$df_overlap, mod_set1 == modSet1, mod_set2 == modSet2)) == 0)
                ) {
                        
                        print("calculating new overlap")
                        
                        for (i in modules1) {
                                
                                print(paste0("module ", i, " of ", max(as.numeric(modules1)) - 1))
                                
                                df_tmp <- modules2 %>% map_dfr(~f_module_overlap(df = df_modules,
                                                                                 mod_set1 = modSet1,
                                                                                 mod_set2 = modSet2,
                                                                                 mod_no1 = i,
                                                                                 mod_no2 = .x)) %>%
                                        mutate(p_adj = p.adjust(p_val, method = "BH"))
                                
                                data$df_overlap <- bind_rows(data$df_overlap, df_tmp)
                                
                        }
                        
                        for (i in modules2) {
                                
                                print(paste0("module ", i, " of ", max(as.numeric(modules2)) - 1))
                                
                                df_tmp <- modules1 %>% map_dfr(~f_module_overlap(df = df_modules,
                                                                                 mod_set1 = modSet2,
                                                                                 mod_set2 = modSet1,
                                                                                 mod_no1 = i,
                                                                                 mod_no2 = .x)) %>%
                                        mutate(p_adj = p.adjust(p_val, method = "BH"))
                                
                                data$df_overlap <- bind_rows(data$df_overlap, df_tmp)
                                
                        }
                        
                        
                } 
                
        })
        
        # If no functional enrichment has been run, create mart and get GO IDs here
        observeEvent(input$go2, {
                
                df_modules <- df_modules_reactive()
                gene_universe <- gene_universe_reactive()
                
                if(nrow(data$df_gene_go) == 0) {
                        
                        doParallel::registerDoParallel()
                        
                        print("creating mart")
                        
                        # get gene info from ensembl database
                        hsapiens_ensembl <- useEnsembl(biomart = "genes",
                                                       dataset = "hsapiens_gene_ensembl",
                                                       version = 88) # GRCh38.87
                        
                        print("mapping genes to GO terms")
                        
                        # map each gene to go term
                        data$df_gene_go <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006"),
                                                 filters = "ensembl_gene_id",
                                                 values = gene_universe,
                                                 mart = hsapiens_ensembl) %>%
                                as.data.frame() %>%
                                as_tibble() %>%
                                filter(go_id != "") %>%
                                dplyr::rename("term" = "name_1006")
                        
                        # unstack
                        data$l_gene_2_GO <- data$df_gene_go %>% 
                                dplyr::select(ensembl_gene_id, go_id) %>% 
                                unstack()
                        
                        
                }
                
        })
        
        # If functional enrichment on selected module and its overlaps hasn't been run, run topGO here
        observeEvent(input$go2, {
                
                
                # bring in reactives & inputs
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                df_modules <- df_modules_reactive()
                gene_universe <- gene_universe_reactive()
                mod_zoom <- input$mod_zoom

                # bring in data
                df_gene_go <- data$df_gene_go
                l_gene_2_GO <- data$l_gene_2_GO
                df_overlap <- data$df_overlap
                
                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", modSet1, modSet2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", modSet2, modSet1)
                
                overlap_mods <- df_overlap %>%
                        filter(p_adj < 0.05) %>% 
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)
                
                # run topGO on zoomed module
                if ((nrow(filter(data$df_mods_go, mod_set == mod_set_zoom, module == mod_zoom)) == 0)) {
                        
                        doParallel::registerDoParallel()
                        
                        print(paste0("calculating functional enrichment for ", mod_set_zoom, " module ", mod_zoom))
                        
                        module <- df_modules %>% 
                                filter(mod_set == mod_set_zoom & module == mod_zoom) %>% 
                                pull(ensembl_gene_id) %>% 
                                unique()
                        
                        df_tmp <- f_top_GO_modules(gene_module = module,
                                                   gene_universe = gene_universe,
                                                   df_gene_go = df_gene_go,
                                                   l_gene_2_GO = l_gene_2_GO,
                                                   go_ontology = "BP") %>%
                                mutate(module = mod_zoom, .before = 1) %>%
                                clean_names() %>%
                                dplyr::select(-term) %>%
                                left_join(df_gene_go %>% dplyr::select(go_id, term), 
                                          by = "go_id") %>%
                                dplyr::select(module, go_id, term, everything()) %>%
                                distinct() %>%
                                mutate(mod_set = mod_set_zoom, .before = 1) %>% 
                                mutate(module = as.numeric(module))
                        
                        data$df_mods_go <- bind_rows(data$df_mods_go, df_tmp)
                        
                } 
                
                        
                for (m in overlap_mods) {
                        
                        # run topGO on overlapping modules
                        if ((nrow(filter(data$df_mods_go, mod_set == mod_set_other, module == m)) == 0)) {
                                
                                doParallel::registerDoParallel()
                                
                                print(paste0("calculating functional enrichment for ", mod_set_other, " module ", m))
                                
                                module <- df_modules %>% 
                                        filter(mod_set == mod_set_other & module == m) %>% 
                                        pull(ensembl_gene_id) %>% 
                                        unique()
                                
                                df_tmp <- f_top_GO_modules(gene_module = module,
                                                           gene_universe = gene_universe,
                                                           df_gene_go = df_gene_go,
                                                           l_gene_2_GO = l_gene_2_GO,
                                                           go_ontology = "BP") %>%
                                        mutate(module = m, .before = 1) %>%
                                        clean_names() %>%
                                        dplyr::select(-term) %>%
                                        left_join(df_gene_go %>% dplyr::select(go_id, term), 
                                                  by = "go_id") %>%
                                        dplyr::select(module, go_id, term, everything()) %>%
                                        distinct() %>%
                                        mutate(mod_set = mod_set_other, .before = 1) %>% 
                                        mutate(module = as.numeric(module))
                                
                                data$df_mods_go <- bind_rows(data$df_mods_go, df_tmp)
                                
                        } 
                        
                        # run topGO on gene intersects
                        if ((nrow(filter(data$df_overlaps_go, 
                                         mod_set1 == mod_set_zoom, mod1 == mod_zoom,
                                         mod_set2 == mod_set_other, mod2 == m)) == 0)) {
                                
                                doParallel::registerDoParallel()
                                
                                print(paste0("calculating functional enrichment for ", mod_set_zoom, " module ", mod_zoom,
                                             " and ", mod_set_other, " module ", m, " intersect"))
                                
                                module1 <- df_modules %>% 
                                        filter(mod_set == mod_set_zoom & module == mod_zoom) %>% 
                                        pull(ensembl_gene_id) %>% 
                                        unique()
                                
                                module2 <- df_modules %>% 
                                        filter(mod_set == mod_set_other & module == m) %>% 
                                        pull(ensembl_gene_id) %>% 
                                        unique()
                                
                                module_intersect <- intersect(module1, module2)
                                
                                df_tmp <- f_top_GO_modules(gene_module = module_intersect,
                                                           gene_universe = gene_universe,
                                                           df_gene_go = df_gene_go,
                                                           l_gene_2_GO = l_gene_2_GO,
                                                           go_ontology = "BP") %>%
                                        clean_names() %>%
                                        mutate(mod_set1 = mod_set_zoom,
                                               mod_set2 = mod_set_other,
                                               mod1 = mod_zoom,
                                               mod2 = m) %>% 
                                        dplyr::select(-term) %>%
                                        left_join(df_gene_go %>% dplyr::select(go_id, term), 
                                                  by = "go_id") %>%
                                        # dplyr::select(module, go_id, term, everything()) %>%
                                        distinct() %>%
                                        mutate(mod1 = as.numeric(mod1),
                                               mod2 = as.numeric(mod2))
                                
                                data$df_overlaps_go <- bind_rows(data$df_overlaps_go, df_tmp)
                                
                        } 
                        
                }
                
        })
        
        ##
        ## TAB 1 PLOTS
        ##
        
        # number & size of modules
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
                        ggtitle("Alluvial plot")
                
                
        }, height = 600, width = 900) %>% bindEvent(input$go1)
        
        # overlap matrix plot
        output$overlap_matrix <- renderPlot({
                
                df_overlap <- data$df_overlap
                df_modules <- df_modules_reactive()
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                
                modules1 <- df_modules %>%
                        filter(mod_set == modSet1) %>%
                        arrange(module) %>%
                        pull(module) %>%
                        unique
                
                modules2 <- df_modules %>%
                        filter(mod_set == modSet2) %>%
                        arrange(module) %>%
                        pull(module) %>%
                        unique
                
                # format for plotting
                df_overlap_plot <-
                        df_overlap %>%
                        filter(mod_set1 == modSet1 & mod_set2 == modSet2) %>% 
                        mutate(is_big = ifelse(q >= 50, 1 , 0) %>% factor()) %>%
                        mutate(p_adj = ifelse(p_adj > 0.05, NA, p_adj)) %>% # set p_adj = 1 to NA
                        mutate(p_adj = ifelse(p_adj == 0, 1e-300, p_adj)) %>%
                        mutate(mod_set1 = factor(mod_set1, levels = unique(.$mod_set1))) %>%
                        mutate(mod1 = factor(mod1, levels = levels(df_overlap$mod2)),
                               mod2 = factor(mod2, levels = levels(df_overlap$mod2))) %>%
                        dplyr::select(-c(q, m, n, k))
                
                # plot
                df_overlap_plot %>% 
                        
                        ggplot(aes(x = mod1, y = mod2)) +
                        geom_tile(aes(fill = -log10(p_adj))) +
                        geom_text(data = df_overlap_plot %>% filter(p_adj < 0.05), 
                                  mapping = aes(label = overlap)) +
                        
                        scale_fill_gradientn(colors = brewer.pal(5, "Reds"), na.value = "white") +
                        scale_x_discrete(limits = modules1) +
                        scale_y_discrete(limits = modules2) +
                        labs(x = modSet1, y = modSet2) +
                        ggtitle("Hypergeometric overlaps")
                
                
        }, height = 600, width = 900) %>% bindEvent(input$go1)
        
        # cell-type enrichment plot
        output$cell_type <- renderPlot({
                
                df_modules <- df_modules_reactive()
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                
                df_modules %>%
                        filter(mod_set %in% c(modSet1, modSet2)) %>% 
                        
                        dplyr::select(mod_set, ensembl_gene_id, module) %>%
                        left_join(df_cell_type, by = "ensembl_gene_id") %>%
                        group_by(mod_set, module) %>%
                        dplyr::count(type) %>%
                        mutate(z_score = (n - mean(n, na.rm = TRUE))/sd(n, na.rm = TRUE)) %>%
                        na.omit() %>% 
                        
                        # filter(!grepl("synapses", type)) %>%
                        
                        ggplot(aes(x = type, y = module,
                                   fill = z_score)) +
                        geom_tile(aes(width = 0.95, height = 0.95)) +
                        scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu")), limits = c(-3, 3)) +
                        facet_wrap(~mod_set, scales = "free", nrow = 2) +
                        labs(x = "Lake cell type") +
                        coord_flip()
                
                
        }, height = 800, width = 900) %>% bindEvent(input$go1)
        
        ##
        ## TAB 2 PLOTS
        ##
        
        # alluvial zoom plot
        output$zoom_alluvial <- renderPlot({
                
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                
                df_colors <- df_colors_reactive()
                df_modules <- df_modules_reactive()
                df_overlap <- data$df_overlap
                
                # q > 50?
                if (input$q50 == FALSE) {
                        
                        df_overlap <- df_overlap
                        
                } else if (input$q50 == TRUE) {
                        
                        df_overlap <- df_overlap %>% filter(q >= 50)
                        
                }
                
                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", modSet1, modSet2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", modSet2, modSet1)
                
                m <- df_overlap %>%
                        filter(p_adj < 0.05) %>% 
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
                        ggtitle("Zoomed module overlaps")
                
                
        }, height = 600, width = 900) %>% bindEvent(input$go2)
        
        # zoomed module enrichment
        output$zoom_module_go <- renderPlot({
                
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                df_mods_go <- data$df_mods_go
                
                # select module set
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", modSet1, modSet2)
                
                # plot
                print(df_mods_go %>% dplyr::count(mod_set, module))
                
                df_mods_go %>% 
                        filter(mod_set == mod_set_zoom & module == mod_zoom) %>%
                        arrange(p_adj) %>% 
                        dplyr::slice(1:10) %>% 
                        
                        ggplot(aes(x = -log10(p_adj), y = reorder(str_wrap(term, width = 35), -p_adj))) +
                        geom_point(aes(size = significant/annotated, color = -log10(p_adj))) +
                        labs(y = "") +

                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
                        scale_color_gradientn(colors = brewer.pal(9, "Reds")) +
                        
                        ggtitle(paste0(mod_set_zoom, " module ", mod_zoom, " enrichment"))

                
        }, height = 800, width = 900) %>% bindEvent(input$go2)
        
        # zoomed module overlapping modules (full) enrichment
        output$zoom_module_overlaps_go <- renderPlot({
                
                # bring in reactives & inputs
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                
                # bring in data
                df_overlap <- data$df_overlap
                df_mods_go <- data$df_mods_go
                
                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", modSet1, modSet2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", modSet2, modSet1)
                
                overlap_mods <- df_overlap %>%
                        filter(p_adj < 0.05) %>% 
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)
                
                # plot
                df_mods_go %>% 
                        
                        filter(mod_set == mod_set_other & module %in% overlap_mods) %>%
                        group_by(module) %>% 
                        arrange(module, p_adj) %>% 
                        dplyr::top_n(n = 10) %>% 
                        
                        ggplot(aes(x = -log10(p_adj), y = reorder_within(str_wrap(term, width = 35), -p_adj, module))) +
                        geom_point(aes(size = significant/annotated, color = -log10(p_adj))) +
                        labs(y = "") +
                        facet_wrap(~module, scales = "free", ncol = 2) +
                        scale_y_discrete(labels = function(x) gsub("*_.", "", x)) +
                        
                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
                        scale_color_gradientn(colors = brewer.pal(9, "Reds")[3:9]) +
                        
                        ggtitle(paste0(mod_set_other, " module ", paste(overlap_mods, sep = ", "), " enrichment"))
                
                
                
        }, height = 800, width = 900) %>% bindEvent(input$go2)
        
        # intersecting gene sets enrichment
        output$intersects_go <- renderPlot({
                
                # bring in reactives & inputs
                modSet1 <- mod_set1()
                modSet2 <- mod_set2()
                mod_zoom <- input$mod_zoom
                
                # bring in data
                df_overlap <- data$df_overlap
                df_mods_go <- data$df_mods_go
                df_overlaps_go <- data$df_overlaps_go
                
                # select modules to look at
                mod_set_zoom <- ifelse(input$mod_set_zoom == "Module set 1", modSet1, modSet2)
                mod_set_other <- ifelse(input$mod_set_zoom == "Module set 1", modSet2, modSet1)
                
                overlap_mods <- df_overlap %>%
                        filter(p_adj < 0.05) %>% 
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom) %>%
                        filter(mod_set2 == mod_set_other) %>%
                        pull(mod2)
                
                # plot
                df_overlaps_go %>% 
                        
                        filter(mod_set1 == mod_set_zoom & mod1 == mod_zoom,
                               mod_set2 == mod_set_other & mod2 %in% overlap_mods) %>%
                        group_by(mod2) %>% 
                        arrange(mod2, p_adj) %>% 
                        dplyr::top_n(n = 10) %>% 
                        
                        ggplot(aes(x = -log10(p_adj), y = reorder_within(str_wrap(term, width = 35), -p_adj, mod2))) +
                        geom_point(aes(size = significant/annotated, color = -log10(p_adj))) +
                        labs(y = "") +
                        facet_wrap(~mod2 , scales = "free", ncol = 2) +
                        scale_y_discrete(labels = function(x) gsub("*_.", "", x)) +
                        
                        geom_vline(aes(xintercept = -log10(0.05)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.01)), lty = 2, color = "black") +
                        geom_vline(aes(xintercept = -log10(0.001)), lty = 2, color = "black") +
                        scale_color_gradientn(colors = brewer.pal(9, "Reds")[3:9]) +
                        
                        ggtitle("Intersecting gene enrichment")
                
                
                
        }, height = 800, width = 900) %>% bindEvent(input$go2)
        
}


# Run app -----------------------------------------------------------------


shinyApp(ui = ui, server = server)
