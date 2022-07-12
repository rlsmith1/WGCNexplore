##
## Describe app here
##

library(shiny)

load("objects/2022-07-05_all_shiny_objects.RDS")
# df_resids, 
# df_modules, df_cell_type, df_overlap,
# df_mods_go, df_mods_go_bigrams,
# df_overlap_go, df_overlap_go_bigrams

shinyUI(fluidPage(
        
        
        ## APPLICATION TITLE
        titlePanel("WGCNA module set comparisons"),
        
        ## SIDEBAR 1: SPECIFY MODULE SETS
        sidebarLayout(
                
                sidebarPanel(
                        
                        selectInput(inputId = "mod_set1",
                                    label = "Module set 1",
                                    choices = df_cell_type %>% pull(mod_set) %>% unique(),
                                    selected = NULL),
                        selectInput(inputId = "mod_set2",
                                    label = "Module set 2",
                                    choices = df_cell_type %>% pull(mod_set) %>% unique(),
                                    selected = NULL),
                        checkboxInput(inputId = "grey_mod",
                                      label = "Show grey?",
                                      value = TRUE),
                        
                        actionButton(inputId = "go1",
                                     label = "GO"),
                        
                        width = 4
                        
                ),
                
                # Main Panel showing plots based on sidebar inputs
                mainPanel(plotOutput("module_sizes", width = "100%", height = "100%"),
                          plotOutput("main_alluvial", width = "100%", height = "100%"),
                          plotOutput("overlap_matrix", width = "100%", height = "100%"),
                          plotOutput("cell_type", width = "100%", height = "100%")
                          
                )
                
        ),
        
        # SIDEBAR 2: ZOOM IN ON MODULE FROM ONE SET AND LOOK AT OVERLAP ENRICHMENT
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
                                     label = "GO")
                ),
                
                # Show alluvial plot of the selected modules
                mainPanel(plotOutput("zoom_alluvial", width = "100%", height = "100%"),
                          plotOutput("wordcloud1", width = "100%", height = "100%"),
                          plotOutput("histogram1", width = "100%", height = "100%"),
                          plotOutput("wordcloud2", width = "100%", height = "100%"),
                          plotOutput("histogram2", width = "100%", height = "100%"),
                          plotOutput("wordcloud3", width = "100%", height = "100%"),
                          plotOutput("histogram3", width = "100%", height = "100%")
                          #
                )
                
        )
        
        
)

)




