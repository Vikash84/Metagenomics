library(tidyverse)
library(ggtree)
library(ape)
library(phangorn)

input_astro_tree <- "data/process/trees/astrovirus/trees/astro_tree.treefile"

astro_tree <- read.tree(file = input_astro_tree)
midpoint_astro_tree <- midpoint(astro_tree)
astro_tree_data <- ggtree(midpoint_astro_tree)
astro_tree_table <- astro_tree_data$data

astro_tip_labels <- astro_tree_data %>% filter(isTip == TRUE) %>%
  mutate(label = str_replace_all(label, pattern = "_", replacement = " ")) %>% 
  mutate(label = str_replace(label, pattern = "Bat astrovirus Tm Guangxi LD71 2007", "Bat astrovirus Tm Guangxi LD71")) %>% 
  mutate(label = str_replace(label, pattern = "Bat astrovirus Tm Guangxi LD77 2007", "Bat astrovirus Tm Guangxi LD77")) %>%
  mutate(label = str_replace(label, pattern = "Bat astrovirus Tm Guangxi LD38 2007", "Bat astrovirus Tm Guangxi LD38"))
  
astro_node_labels <- astro_tree_data %>% filter(isTip == FALSE) %>%
  mutate(label = as.numeric(label)) %>%
  filter(label >= 50)

astro_tree_figure <- ggtree(midpoint_astro_tree) +
  geom_text(data = astro_node_labels, aes(label = label), hjust = 1.3, vjust = -0.6, size = 1) +
  geom_tiplab(data = astro_tip_labels, aes(label = label, fill = grepl("Murine astrovirus JS1", label)),
              label.padding = unit(0.05, "lines"), label.size = 0, label.shape = 4,
              align = TRUE, size = 1.8, geom = "label") +
  scale_fill_manual(values = c(NA, "lightpink")) +
  theme(legend.position="none") +
  geom_treescale(x = 0.05, y = 17, color = "red", fontsize=2.5) +
  xlim_tree(3)

astro_tree_figure

ggsave(filename = "data/process/trees/mitovirus/trees/astro_tree.png", plot = astro_tree_figure, width = 8, height = 8)
