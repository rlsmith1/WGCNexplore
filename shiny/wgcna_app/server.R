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


