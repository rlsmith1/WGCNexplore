

library(shiny)
library(tidyverse)

# UI
ui <- fluidPage(
        
        # Application title
        titlePanel("Practice app"),
        
        # Sidebar with a slider input for number of bins 
        sidebarLayout(
                sidebarPanel(
                        numericInput(inputId = "x",
                                     label = "X:",
                                     value = 0,
                                     step = 0.5
                        ),
                        numericInput(inputId = "y",
                                     label = "Y:",
                                     value = 0,
                                     step = 0.5
                        ),
                        actionButton("go", "GO")
                        
                ),
                
                # Show a plot of the generated distribution
                mainPanel(
                        plotOutput("plot")
                )
        )
)

##
## SERVER
##
server <- function(input, output) {
        
        data <- reactiveValues()
        
        data$df <- tibble(x = as.numeric(),
                          y = as.numeric(),
                          sum = as.numeric())
        
        observeEvent(input$go, {
                
                x1 <- input$x
                y1 <- input$y
                
                if (
                        (nrow(data$df) == 0) |
                        (nrow(filter(data$df, x == x1, y == y1)) == 0)
                ) {
                        
                        df_tmp <- tibble(x = x1,
                                         y = y1,
                                         sum = sum(x1, y1))
                        
                        data$df <- bind_rows(data$df, df_tmp)
                        
                        
                } 
                
        })
        
        # plot output
        output$plot <- renderPlot({
                
                print("df_plot")
                print(data$df)
                
                data$df %>%
                        ggplot(aes(x = x, y = y, color = sum)) +
                        geom_point(size = 3)
                
        }) %>% bindEvent(input$go)
        
}


# Run the application 
shinyApp(ui = ui, server = server)
