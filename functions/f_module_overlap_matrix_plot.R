

# libraries ---------------------------------------------------------------

  library(tidyverse)
  library(RColorBrewer)


# function ----------------------------------------------------------------


  f_plot_module_matrix <- function(df_overlap_plot, module_set1, module_set2) {
    
          x_axis_scale <- df_overlap_plot %>% 
                  filter(mod_set1 == module_set1) %>% 
                  arrange(mod1) %>% 
                  pull(mod1) %>% unique
          
          y_axis_scale <- df_overlap_plot %>% 
                  filter(mod_set2 == module_set2) %>%
                  arrange(mod2) %>% 
                  pull(mod2) %>% unique
          
          df_overlap_plot %>% 
                  filter(mod_set1 == module_set1 & mod_set2 == module_set2) %>% 
                  
          ggplot(aes(x = mod1, y = mod2)) +
                  geom_tile(aes(fill = -log10(p_adj))) +
                  geom_text(aes(label = overlap)) +
                         
                  scale_fill_gradientn(colors = brewer.pal(5, "Reds"), na.value = "white") +
                  scale_x_discrete(limits = x_axis_scale) +
                  scale_y_discrete(limits = y_axis_scale) +
                  labs(x = module_set1, y = module_set2) +
                  theme_classic() +
                  theme(plot.title = element_text(size = 18),
                        axis.title = element_text(size = 15),
                        axis.text = element_text(size = 12),
                        strip.text = element_text(size = 15),
                        legend.title = element_text(size = 15),
                        legend.text = element_text(size = 12))
          
  }     

