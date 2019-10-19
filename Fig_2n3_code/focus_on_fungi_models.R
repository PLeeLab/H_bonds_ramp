library("tidyverse")
theme_new <- function() {
  theme_bw() %+replace%
    theme(
      plot.title = element_text(color = "black", size = 21, face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "black", size = 21, face = "bold", hjust = 0),
      axis.title = element_text(color = "black", size = 7),
      axis.text = element_text(color = "black", size = 6),
      legend.text = element_text(color = "black", size = 5),
      legend.title = element_text(color = "black", size = 5),
      strip.text = element_text(color = "black", size = 6),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      complete = TRUE
    )
}

models_df_Fungi <- read.delim(file = "inputs/model_fitness/OUT_Fungi_Model_fitness_5end_ALL_1_to_293_v0.0.2.txt")

models_df_Fungi %>% 
  select(Modeling) %>% 
  table() %>% 
  View()

rbind(
  models_df_Fungi %>% 
  mutate(Best_AIC = apply(.[,7:9], 1, function(x) names(x)[which.min(x)])) %>% 
  mutate(Best_BIC = apply(.[,10:12], 1, function(x) names(x)[which.min(x)])) %>% 
  select(Best_AIC) %>%
  table() %>% 
  as_tibble(),
  models_df_Fungi %>% 
    mutate(Best_AIC = apply(.[,7:9], 1, function(x) names(x)[which.min(x)])) %>% 
    mutate(Best_BIC = apply(.[,10:12], 1, function(x) names(x)[which.min(x)])) %>% 
    select(Best_BIC) %>%
    table() %>% 
    as_tibble()
  )


models_df_Fungi %>% 
  mutate(Best_AIC = apply(.[,7:9], 1, function(x) names(x)[which.min(x)]), 
         Best_BIC = apply(.[,10:12], 1, function(x) names(x)[which.min(x)])) %>% 
  mutate(Best_AIC = str_replace(Best_AIC, "AIC_", ""),
         Best_BIC = str_replace(Best_BIC, "BIC_", "")) %>% 
  select(Best_AIC) %>%
  table() %>% as_tibble()


Fig_S3 <- cbind(
  models_df_Fungi %>% 
    mutate(Best_AIC = apply(.[,7:9], 1, function(x) names(x)[which.min(x)])) %>% 
    mutate(Best_AIC = str_replace(Best_AIC, "AIC_", "")) %>%
    rename(AIC = Best_AIC) %>% 
    select(AIC) %>% 
    table() %>% 
    as.data.frame(row.names = 1) %>%
    rename(AIC = Freq),
  
  models_df_Fungi %>% 
    mutate(Best_BIC = apply(.[,10:12], 1, function(x) names(x)[which.min(x)])) %>% 
    mutate(Best_BIC = str_replace(Best_BIC, "BIC_", "")) %>%
    rename(BIC = Best_BIC) %>% 
    select(BIC) %>% 
    table() %>% 
    as.data.frame(row.names = 1) %>%
    rename(BIC = Freq)
  ) %>%
  as_tibble(rownames = "Model") %>%
  pivot_longer(cols = c(AIC, BIC), names_to = "Criterion", values_to = "Freq") %>% 
  group_by(Criterion) %>%
  mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
  mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
  ggplot(aes(x = Criterion, y = Freq, fill = Model)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1)+
  theme_new() +
  theme(legend.position = "bottom") +
  NULL

ggsave(
  plot = Fig_S3,
  filename = "outputs/Fig_S3.pdf",
  device = "pdf",
  width = 89,
  height = 89,
  units = "mm",
  useDingbats = FALSE
)





  