#-------- LIBRARIES -------#
# install.packages("tidyverse")
# install.packages("EnvStats")
# install.packages("ggpubr")
# install.packages("cowplot")
# install.packages("ggcorrplot")
# install.packages("scales")
library("tidyverse")
library("EnvStats")
library("ggpubr")
library("cowplot")
library("ggcorrplot")
library("scales") # to access break formatting functions
#-------------------------#
rm(list = ls()) # clear environment

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

#-------- INPUTS -------#
# source(file = "get_model_parameters.R")
Elements_df_boot <- read.delim(file = "inputs/OUT_Elements_df_boot.txt")
#Elements_df_boot <- read_csv(file = "~/Downloads/data.csv")

tax_df <- rbind(
  read.delim(file = "inputs/taxonomy/integrated_tax_df_Bacteria_genomes.txt"),
  read.delim(file = "inputs/taxonomy/integrated_tax_df_Archaea_genomes.txt"),
  read.delim(file = "inputs/taxonomy/integrated_tax_df_Fungi_genomes.txt")
)

models_df_Bacteria <- rbind(
  read.delim(file = "inputs/model_fitness/OUT_Model_fitness_5end_ALL_1_to_3000.txt"),
  read.delim(file = "inputs/model_fitness/OUT_Model_fitness_5end_ALL_3001_to_7000.txt"),
  read.delim(file = "inputs/model_fitness/OUT_Model_fitness_5end_ALL_7001_to_12001.txt"),
  read.delim(file = "inputs/model_fitness/OUT_Model_fitness_5end_ALL_12002_to_13000.txt"),
  read.delim(file = "inputs/model_fitness/OUT_Model_fitness_5end_ALL_13001_to_13500.txt"),
  read.delim(file = "inputs/model_fitness/OUT_Model_fitness_5end_ALL_13501_to_13900.txt"),
  read.delim(file = "inputs/model_fitness/OUT_Model_fitness_5end_ALL_13901_to_13921.txt")
)
models_df_Archaea <- read.delim(file = "inputs/model_fitness/OUT_Archaea_Model_fitness_5end_ALL_1_to_297.txt")
models_df_Fungi <- read.delim(file = "inputs/model_fitness/OUT_Fungi_Model_fitness_5end_ALL_1_to_293.txt")


ORFeome_features <- rbind(
  read.delim(file = "inputs/ORFeome_features/ORFeome_features_Archaea.txt"),
  read.delim(file = "inputs/ORFeome_features/ORFeome_features_Bacteria.txt"),
  read.delim(file = "inputs/ORFeome_features/ORFeome_features_Fungi.txt")
)
#-------------------------#
source(file = "get_model_parameters.R")
model_fitness_5end <- get_model_parameters(ELEMENT = "Hydrogen bonds", START_AFTER_CODON_POSITION = 1, END_AT_CODON_POSITION = 100, RETURN = "data")
model_fitness_5end
Fig_2A <- get_model_parameters(ELEMENT = "Hydrogen bonds", START_AFTER_CODON_POSITION = 1, END_AT_CODON_POSITION = 100, RETURN = "plot")

dir.create(path = "outputs", showWarnings = TRUE)
write.table(x = model_fitness_5end, file = "outputs/OUT_Model_fitness_5end.txt", sep = "\t", row.names = FALSE, quote = FALSE)

models_df <- rbind(models_df_Bacteria, models_df_Archaea, models_df_Fungi)
models_df$Genome <- gsub(x = gsub(x = gsub(x = models_df$Genome, pattern = "GCF_", replacement = "GCF-"), pattern = "_.*", replacement = ""), pattern = "GCF-", replacement = "GCF_")

integrated_df <- merge(x = models_df, y = tax_df, by.x = "Genome", by.y = "AssemblyAccession")

integrated_df$complete_tax <- paste(integrated_df$phylum, integrated_df$class, integrated_df$order, integrated_df$family, integrated_df$genus, integrated_df$species, sep = "__")

ORFeome_features$Genome <- gsub(x = gsub(x = gsub(x = ORFeome_features$Genome, pattern = "GCF_", replacement = "GCF-"), pattern = "_.*", replacement = ""), pattern = "GCF-", replacement = "GCF_")

integrated_df <- merge(x = integrated_df, y = ORFeome_features, by.x = "Genome", by.y = "Genome")
# write.table(x = integrated_df, file = "integrated_df.txt", sep = "\t", row.names = FALSE)

Fig_S2A <- integrated_df %>%
  ggplot(mapping = aes(
    x = superkingdom,
    y = Length_total,
    fill = superkingdom
  )) +
  geom_hline(yintercept = c(3500000), lty = 2) +
  geom_boxplot(show.legend = FALSE, width = 0.5, colour = "black", notch = TRUE) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_brewer(palette = "Set1") +
  xlab(label = "") +
  ylab(label = "Total ORFeome length (bp)") +
  stat_n_text(y.expand.factor = 0.1, geom = "text", fontface = "italic", vjust = 0, size = 4, color = "black") +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

Fig_S2B <- integrated_df %>%
  ggplot(mapping = aes(
    x = superkingdom,
    y = Length_mean,
    fill = superkingdom
  )) +
  geom_hline(yintercept = 1000, lty = 2) +
  geom_boxplot(show.legend = FALSE, width = 0.5, colour = "black", notch = TRUE) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_brewer(palette = "Set1") +
  xlab(label = "") +
  ylab(label = "Mean length of CDSs per genome (bp)") +
  stat_n_text(y.expand.factor = 0.1, geom = "text", fontface = "italic", vjust = 0, size = 4, color = "black") +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

Fig_S2C <- integrated_df %>%
  ggplot(mapping = aes(
    x = superkingdom,
    y = GC3_mean / GC_mean,
    fill = superkingdom
  )) +
  geom_hline(yintercept = c(1), lty = 2) +
  geom_boxplot(show.legend = FALSE, width = 0.5, colour = "black", notch = TRUE) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_brewer(palette = "Set1") +
  xlab(label = "") +
  ylab(label = "GC3/GC") +
  stat_n_text(y.expand.factor = 0.1, geom = "text", fontface = "italic", vjust = 0, size = 3, color = "black") +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

Fig_S2D <- integrated_df %>%
  filter(phylum != "NA") %>%
  ggplot(mapping = aes(
    x = reorder(phylum, GC_mean, FUN = mean),
    y = GC_mean * 100, # Value
    colour = superkingdom
  )) +
  geom_hline(yintercept = 50, lty = 2) +
  stat_summary(fun.y = mean, fun.ymax = function(x) {
    mean(x) + (1 * sd(x))
  }, fun.ymin = function(x) {
    mean(x) - (1 * sd(x))
  }, geom = "errorbar") +
  stat_summary(fun.y = mean, geom = "point", size = 3, shape = 15) +
  scale_color_brewer(palette = "Set1") +
  xlab(label = "Phylum") +
  ylab(label = "GC per genome (%)") +
  stat_n_text(y.expand.factor = 0.05, geom = "text", fontface = "italic", angle = 90, hjust = 0, size = 3, color = "grey") +
  theme_new() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NULL

Fig_S2E <- integrated_df %>%
  filter(phylum != "NA") %>%
  ggplot(mapping = aes(
    x = reorder(phylum, GC_mean, FUN = mean),
    y = (GC3_mean * 100) / (GC_mean * 100), # Value
    colour = superkingdom
  )) +
  geom_hline(yintercept = 1, lty = 2) +
  stat_summary(fun.y = mean, fun.ymax = function(x) {
    mean(x) + (1 * sd(x))
  }, fun.ymin = function(x) {
    mean(x) - (1 * sd(x))
  }, geom = "errorbar") +
  stat_summary(fun.y = mean, geom = "point", size = 3, shape = 15) +
  scale_color_brewer(palette = "Set1") +
  xlab(label = "Phylum") +
  ylab(label = "GC3/GC per genome") +
  stat_n_text(y.expand.factor = 0.05, geom = "text", fontface = "italic", angle = 90, hjust = 0, size = 3, color = "grey") +
  theme_new() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NULL

Fig_S2ABC <- plot_grid(Fig_S2A, Fig_S2B, Fig_S2C, labels = c("A", "B", "C"), label_size = 22, ncol = 3, nrow = 1)
Fig_S2 <- plot_grid(Fig_S2ABC, Fig_S2D, Fig_S2E, labels = c("", "D", "E"), label_size = 22, ncol = 1, nrow = 3, rel_heights = c(1, 1.75, 1.75))

ggsave(
  plot = Fig_S2,
  filename = "outputs/Fig_S2.pdf",
  device = "pdf",
  width = 10,
  height = 12,
  useDingbats = FALSE
)

Fig_S4A <- integrated_df %>%
  ggplot(mapping = aes(
    x = superkingdom,
    y = Length_total,
    fill = Modeling,
    colour = superkingdom
  )) +
  geom_hline(yintercept = c(3500000), lty = 2) +
  geom_boxplot(width = 0.5, notch = TRUE) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  # scale_fill_brewer(palette = "Set3") +
  scale_fill_manual(values = c("white", "grey")) +
  scale_color_brewer(palette = "Set1") +
  stat_compare_means(label = "p.format", method = "t.test") + # "p.format"
  xlab(label = "") +
  ylab(label = "Total ORFeome length (bp)") +
  stat_n_text(y.expand.factor = 0.1, geom = "text", fontface = "italic", vjust = 0, size = 4, color = "black") +
  theme_new() +
  theme(legend.position = "bottom") +
  NULL

Fig_S4B <- integrated_df %>%
  ggplot(mapping = aes(
    x = superkingdom,
    y = Length_mean,
    fill = Modeling,
    colour = superkingdom
  )) +
  geom_hline(yintercept = 1000, lty = 2) +
  geom_boxplot(width = 0.5, notch = TRUE) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = c("white", "grey")) +
  scale_color_brewer(palette = "Set1") +
  stat_compare_means(label = "p.format", method = "t.test") + # "p.format"
  xlab(label = "") +
  ylab(label = "Mean length of CDSs per genome (bp)") +
  stat_n_text(y.expand.factor = 0.1, geom = "text", fontface = "italic", vjust = 0, size = 4, color = "black") +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

Fig_S4C <- integrated_df %>%
  ggplot(mapping = aes(
    x = superkingdom,
    y = GC3_mean / GC_mean,
    fill = Modeling,
    colour = superkingdom
  )) +
  geom_hline(yintercept = c(1), lty = 2) +
  geom_boxplot(width = 0.5, notch = TRUE) +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_manual(values = c("white", "grey")) +
  scale_color_brewer(palette = "Set1") +
  stat_compare_means(label = "p.format", method = "t.test") + # "p.format"
  xlab(label = "") +
  ylab(label = "GC3/GC") +
  stat_n_text(y.expand.factor = 0.1, geom = "text", fontface = "italic", vjust = 0, size = 3, color = "black") +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

Fig_S4D <- integrated_df %>%
  filter(phylum != "NA") %>%
  ggplot(mapping = aes(
    x = reorder(phylum, GC_mean, FUN = mean),
    y = GC_mean * 100,
    colour = superkingdom,
    shape = Modeling,
    lty = Modeling
  )) +
  geom_hline(yintercept = 50, lty = 2) +
  stat_summary(fun.y = mean, fun.ymax = function(x) {
    mean(x) + (1 * sd(x))
  }, fun.ymin = function(x) {
    mean(x) - (1 * sd(x))
  }, geom = "errorbar", position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = mean, geom = "point", size = 2, position = position_dodge(width = 0.75)) +
  scale_color_brewer(palette = "Set1") +
  xlab(label = "Phylum") +
  ylab(label = "GC per genome (%)") +
  stat_n_text(y.expand.factor = 0.05, geom = "text", fontface = "italic", angle = 90, hjust = 0, size = 3, color = "grey") +
  theme_new() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NULL

Fig_S4E <- integrated_df %>%
  filter(phylum != "NA") %>%
  ggplot(mapping = aes(
    x = reorder(phylum, GC_mean, FUN = mean),
    y = (GC3_mean * 100) / (GC_mean * 100), # Value
    colour = superkingdom,
    shape = Modeling,
    lty = Modeling
  )) +
  geom_hline(yintercept = 1, lty = 2) +
  stat_summary(fun.y = mean, fun.ymax = function(x) {
    mean(x) + (1 * sd(x))
  }, fun.ymin = function(x) {
    mean(x) - (1 * sd(x))
  }, geom = "errorbar", position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = mean, geom = "point", size = 2, position = position_dodge(width = 0.75)) +
  scale_color_brewer(palette = "Set1") +
  xlab(label = "Phylum") +
  ylab(label = "GC3/GC per genome") +
  stat_n_text(y.expand.factor = 0.05, geom = "text", fontface = "italic", angle = 90, hjust = 0, size = 3, color = "grey") +
  theme_new() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  NULL

Fig_S4ABC <- plot_grid(Fig_S4A, Fig_S4B, Fig_S4C, labels = c("A", "B", "C"), label_size = 22, ncol = 3, nrow = 1)
Fig_S4 <- plot_grid(Fig_S4ABC, Fig_S4D, Fig_S4E, labels = c("", "D", "E"), label_size = 22, ncol = 1, nrow = 3, rel_heights = c(1, 1.75, 1.75))

ggsave(
  plot = Fig_S4,
  filename = "outputs/Fig_S4.pdf",
  device = "pdf",
  width = 10,
  height = 14,
  useDingbats = FALSE
)


Fig_2B <- integrated_df %>%
  ggplot(mapping = aes(
    x = superkingdom,
    y = No_of_CDSs,
    fill = superkingdom
  )) +
  geom_hline(yintercept = c(3500), lty = 2, lwd = 0.25) +
  geom_boxplot(width = 0.5, colour = "black", notch = TRUE, size = 0.25) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_brewer(palette = "Set1") +
  xlab(label = "") +
  ylab(label = "CDSs per ORFeome") +
  stat_n_text(y.expand.factor = 0.1, geom = "text", fontface = "italic", vjust = 0, size = 2, color = "black") +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

Fig_2C <- integrated_df %>%
  select(Modeling, superkingdom) %>%
  table() %>%
  as.data.frame() %>%
  group_by(superkingdom) %>%
  mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
  mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
  ggplot(aes(x = superkingdom, y = Percentage, fill = Modeling)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(x = superkingdom, y = Percentage, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
  # geom_text(aes(x = superkingdom, y = Freq, label = Labels), colour = "white", position = position_stack(vjust=0.5), size=4) +
  xlab(label = "") +
  # labs(subtitle = "All data") +
  scale_fill_brewer(palette = "Dark2") +
  theme_new() +
  theme(legend.position = "bottom") +
  NULL

integrated_df$P_Rate <- as.numeric(integrated_df$P_Rate)
integrated_df$P_Initial_cost <- as.numeric(integrated_df$P_Initial_cost)
P_VAL_THRESH <- 0.001
Fig_2D <- integrated_df %>%
  mutate(Significant = ifelse(test = as.numeric(P_Rate) < P_VAL_THRESH, yes = "Yes", no = "No")) %>%
  select(Significant, superkingdom) %>%
  table() %>%
  as.data.frame() %>%
  group_by(superkingdom) %>%
  mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
  mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
  ggplot(aes(x = superkingdom, y = Percentage, fill = Significant)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(x = superkingdom, y = Percentage, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
  xlab(label = "") +
  # labs(subtitle = "Subset of succesfull modeling", fill = paste0("Signif (", 100-(P_VAL_THRESH*100), "%)")) + # paste0("Significant at\nconfidence\nlevel = ", 100-(P_VAL_THRESH*100), "%")
  labs(fill = paste0("Signif (", 100 - (P_VAL_THRESH * 100), "%)")) + # paste0("Significant at\nconfidence\nlevel = ", 100-(P_VAL_THRESH*100), "%")
  scale_fill_brewer(palette = "Accent") +
  theme_new() +
  theme(legend.position = "bottom") +
  NULL

# colnames(integrated_df)
Fig_2E <- cbind(
  cor(
    x = integrated_df %>%
      filter(superkingdom == "Bacteria") %>% # Eukaryota Bacteria
      select(No_of_CDSs, Length_total, Length_mean, GC_mean, GC1_mean, GC2_mean, GC3_mean),
    y = integrated_df %>%
      filter(superkingdom == "Bacteria") %>%
      select(P_Rate) %>%
      mutate(Bacteria = -log10(P_Rate)) %>%
      select(Bacteria),
    use = "complete.obs"
  ),

  cor(
    x = integrated_df %>%
      filter(superkingdom == "Archaea") %>% # Bacteria  Archaea
      select(No_of_CDSs, Length_total, Length_mean, GC_mean, GC1_mean, GC2_mean, GC3_mean),
    y = integrated_df %>%
      filter(superkingdom == "Archaea") %>%
      select(P_Rate) %>%
      mutate(Archaea = -log10(P_Rate)) %>%
      select(Archaea),
    use = "complete.obs"
  ),

  cor(
    x = integrated_df %>%
      filter(superkingdom == "Eukaryota") %>% # Eukaryota Bacteria  Archaea
      select(No_of_CDSs, Length_total, Length_mean, GC_mean, GC1_mean, GC2_mean, GC3_mean),
    y = integrated_df %>%
      filter(superkingdom == "Eukaryota") %>%
      select(P_Rate) %>%
      mutate(Eukaryota = -log10(P_Rate)) %>%
      select(Eukaryota),
    use = "complete.obs"
  )
) %>% ggcorrplot(outline.col = "black", lab = TRUE, lab_size = 2) +
  xlab(label = "") +
  ylab(label = "Significance of Rate parameter") +
  coord_flip() +
  theme_new() +
  theme(
    legend.position = "bottom",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  ) +
  NULL

linear_model_bact <- lm(formula = -log10(P_Rate) ~ No_of_CDSs, data = integrated_df, subset = superkingdom == "Bacteria")
linear_model_arch <- lm(formula = -log10(P_Rate) ~ No_of_CDSs, data = integrated_df, subset = superkingdom == "Archaea")
linear_model_fung <- lm(formula = -log10(P_Rate) ~ No_of_CDSs, data = integrated_df, subset = superkingdom == "Eukaryota")

summary(linear_model_bact)
summary(linear_model_arch)
summary(linear_model_fung)

summary(-log10(integrated_df$P_Rate))

Fig_2F <- integrated_df %>%
  ggplot(aes(
    x = No_of_CDSs, # No_of_CDSs , Length_total
    y = -log10(P_Rate),
    # y = Rate,
    colour = superkingdom
  )) +
  geom_point(alpha=0.1) +
  geom_hline(yintercept = (0.01), lty = 2, colour = "black", lwd = 0.25) +
  #geom_smooth(method = "lm", lwd = 0.25) +
  stat_smooth(method = "lm", lwd = 0.25, fullrange = FALSE) +
  annotate("text", x = 10^2.5, y = 30, label = paste0("Adj R^2 = ", round(summary(linear_model_bact)$adj.r.squared, 3)), color = "red", size = 1.5) +
  annotate("text", x = 10^2.5, y = 25, label = paste0("Adj R^2 = ", round(summary(linear_model_arch)$adj.r.squared, 3)), color = "blue", size = 1.5) +
  annotate("text", x = 10^2.5, y = 20, label = paste0("Adj R^2 = ", round(summary(linear_model_fung)$adj.r.squared, 3)), color = "green", size = 1.5) +
  annotate("text", x = 10^3, y = 30, label = paste0("Pval = ", round(summary(linear_model_bact)$coefficients[2, 4], 4)), color = "red", size = 1.5) +
  annotate("text", x = 10^3, y = 25, label = paste0("Pval = ", round(summary(linear_model_arch)$coefficients[2, 4], 4)), color = "blue", size = 1.5) +
  annotate("text", x = 10^3, y = 20, label = paste0("Pval = ", round(summary(linear_model_fung)$coefficients[2, 4], 4)), color = "green", size = 1.5) +
  # annotate("text", x = 10^3, y = 30, label = paste0("Fstat = ", round(summary(linear_model_bact)$fstatistic[1],3)), color="red", size=2 ) +
  # annotate("text", x = 10^3, y = 25, label = paste0("Fstat = ", round(summary(linear_model_arch)$fstatistic[1],3)), color="blue" , size=2) +
  # annotate("text", x = 10^3, y = 20, label = paste0("Fstat = ", round(summary(linear_model_fung)$fstatistic[1],3)), color="green", size=2 ) +
  scale_color_brewer(palette = "Set1") +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  #ylim(c(-5, 35)) +
  # facet_wrap(facets = ~superkingdom, scales = "free_x") +
  theme_new() +
  theme(legend.position = "bottom") +
  NULL

ggsave(
  plot = plot_grid(Fig_2A, Fig_2B, Fig_2C, Fig_2D, Fig_2E, Fig_2F, labels = c("A", "B", "C", "D", "E", "F"), label_size = 8, ncol = 3, nrow = 2),
  filename = "outputs/Fig_2.pdf",
  device = "pdf",
  width = 170,
  height = 120,
  units = "mm",
  useDingbats = FALSE
)

# Rate, Carrying_capacity, Initial_cost
my_comparisons <- list(c("Archaea", "Eukaryota"), c("Bacteria", "Archaea"), c("Bacteria", "Eukaryota"))

Fig_3A <- ggplot(
  data = integrated_df,
  mapping = aes(x = superkingdom, y = Carrying_capacity)
) + #  "Carrying_capacity" "Rate" "Initial_cost"
  # geom_hline(yintercept = 0, lty=2, colour="black") +
  geom_violin(aes(colour = superkingdom, fill = superkingdom), alpha = 0.5, size = 0.25) +
  stat_n_text(y.expand.factor = 0.2, geom = "text", fontface = "italic", vjust = 0, size = 2, color = "grey") +
  stat_summary(fun.y = median, fun.ymin = function(x) {
    median(x) - sd(x)
  }, fun.ymax = function(x) {
    median(x) + sd(x)
  }, colour = "black", geom = "pointrange", size = 0.25) +
  stat_compare_means(label = "p.format", comparisons = my_comparisons, method = "wilcox", size = 2) + # "p.format" "p.signif"
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  xlab(label = "") +
  ylab(label = "Carrying capacity") +
  # labs(title = "Unsuccessful modeling as Rate=0") +
  theme_new() +
  theme(legend.position = "none") +
  NULL

Fig_3B <- ggplot(
  data = integrated_df,
  mapping = aes(x = superkingdom, y = Initial_cost)
) + #  "Carrying_capacity" "Rate" "Initial_cost"
  # geom_hline(yintercept = 0, lty=2, colour="black") +
  geom_violin(aes(colour = superkingdom, fill = superkingdom), alpha = 0.5, size = 0.25) +
  stat_n_text(y.expand.factor = 0.2, geom = "text", fontface = "italic", vjust = 0, size = 2, color = "grey") +
  stat_summary(fun.y = median, fun.ymin = function(x) {
    median(x) - sd(x)
  }, fun.ymax = function(x) {
    median(x) + sd(x)
  }, colour = "black", geom = "pointrange", size = 0.25) +
  stat_compare_means(label = "p.format", comparisons = my_comparisons, method = "wilcox", size = 2) + # "p.format" "p.signif"
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  xlab(label = "") +
  ylab(label = "Initial cost") +
  # labs(title = "Unsuccessful modeling as Rate=0") +
  theme_new() +
  theme(legend.position = "none") +
  NULL

Fig_3C <- ggplot(
  data = integrated_df %>%
    filter(Modeling == "Successful"),
  mapping = aes(x = superkingdom, y = Rate)
) + # superkingdom         phylum          class           order            family           genus
  geom_hline(yintercept = 0, lty = 2, colour = "black", lwd = 0.25) +
  geom_violin(aes(colour = superkingdom, fill = superkingdom), alpha = 0.5, size = 0.25) +
  stat_n_text(y.expand.factor = 0.2, geom = "text", fontface = "italic", vjust = 0, size = 2, color = "grey") +
  stat_summary(fun.y = median, fun.ymin = function(x) {
    median(x) - sd(x)
  }, fun.ymax = function(x) {
    median(x) + sd(x)
  }, colour = "black", geom = "pointrange", size = 0.25) +
  stat_compare_means(label = "p.format", comparisons = my_comparisons, method = "wilcox", size = 2) + # "p.format" "p.signif"
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  xlab(label = "") +
  # labs(title = "Unsuccessful modeling as Rate=0") +
  theme_new() +
  theme(legend.position = "none") +
  NULL

Fig_3D <- ggplot(
  data = integrated_df, # %>%
  # filter(Modeling == "Successful"),
  mapping = aes(x = Initial_cost, y = Carrying_capacity, colour = superkingdom)
) +
  # geom_point(alpha=0.25) +
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+
  geom_smooth(method = "lm", show.legend = FALSE, level = 0.99, lwd = 0.25) +
  geom_abline(slope = 1, lty = 2, lwd = 0.25) +
  # geom_rug(alpha=0.1) +
  xlim(c(5.5, 9)) +
  ylim(c(5.5, 9)) +
  xlab(label = "Initial cost (H bonds per codon)") +
  ylab(label = "Capacity (H bonds per codon)") +
  scale_color_brewer(palette = "Set1") +
  # scale_fill_viridis_c(option = "plasma") +
  # facet_wrap(~superkingdom) +
  theme_new() +
  # theme_bw() +
  theme(legend.position = "none") +
  NULL

ggsave(
  plot = plot_grid(Fig_3A, Fig_3B, Fig_3C, Fig_3D, labels = c("A", "B", "C", "D"), label_size = 8, ncol = 2, nrow = 2),
  filename = "outputs/Fig_3.pdf",
  device = "pdf",
  width = 89,
  height = 89,
  units = "mm",
  useDingbats = FALSE
)

p.adjust(p = c(0.04, 0.0082, 0.00025), method = "BY")


# Fig_S5A <- integrated_df %>%
#   filter(superkingdom == "Eukaryota") %>% 
#   select(Modeling, phylum) %>%
#   table() %>%
#   as.data.frame() %>%
#   group_by(phylum) %>%
#   filter(Freq > 0) %>% 
#   mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
#   #mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
#   mutate(Labels = Freq) %>%
#   #ggplot(aes(x = phylum, y = Percentage, fill = Modeling)) +
#   ggplot(aes(x = phylum, y = Freq, fill = Modeling)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   #geom_text(aes(x = phylum, y = Percentage, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   geom_text(aes(x = phylum, y = Freq, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   xlab(label = "Phylum") +
#   coord_flip() +
#   # labs(subtitle = "All data") +
#   #scale_fill_brewer(palette = "Dark2") +
#   theme_new() +
#   theme(legend.position = "bottom",
#         axis.text.y = element_text(face = "italic")) +
#   NULL
# 
# Fig_S5B <- integrated_df %>%
#   filter(superkingdom == "Eukaryota", phylum == "Ascomycota") %>% 
#   select(Modeling, class) %>%
#   table() %>%
#   as.data.frame() %>%
#   group_by(class) %>%
#   filter(Freq > 0) %>% 
#   mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
#   #mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
#   mutate(Labels = Freq) %>%
#   #ggplot(aes(x = class, y = Percentage, fill = Modeling)) +
#   ggplot(aes(x = class, y = Freq, fill = Modeling)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   #geom_text(aes(x = class, y = Percentage, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   geom_text(aes(x = class, y = Freq, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   xlab(label = "Class") +
#   coord_flip() +
#   # labs(subtitle = "All data") +
#   #scale_fill_brewer(palette = "Dark2") +
#   theme_new() +
#   theme(legend.position = "bottom",
#         axis.text.y = element_text(face = "italic")) +
#   NULL
# 
# Fig_S5C <- integrated_df %>%
#   filter(superkingdom == "Eukaryota", phylum == "Ascomycota", class == "Saccharomycetes") %>% 
#   select(Modeling, order) %>%
#   table() %>%
#   as.data.frame() %>%
#   group_by(order) %>%
#   filter(Freq > 0) %>% 
#   mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
#   #mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
#   mutate(Labels = Freq) %>%
#   #ggplot(aes(x = order, y = Percentage, fill = Modeling)) +
#   ggplot(aes(x = order, y = Freq, fill = Modeling)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   #geom_text(aes(x = order, y = Percentage, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   geom_text(aes(x = order, y = Freq, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   xlab(label = "Order") +
#   coord_flip() +
#   # labs(subtitle = "All data") +
#   #scale_fill_brewer(palette = "Dark2") +
#   theme_new() +
#   theme(legend.position = "bottom",
#         axis.text.y = element_text(face = "italic")) +
#   NULL

# Fig_S5D <- integrated_df %>%
#   filter(superkingdom == "Eukaryota", phylum == "Ascomycota", class == "Saccharomycetes", order == "Saccharomycetales") %>% 
#   select(Modeling, family) %>%
#   table() %>%
#   as.data.frame() %>%
#   group_by(family) %>%
#   filter(Freq > 0) %>% 
#   mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
#   #mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
#   mutate(Labels = Freq) %>%
#   #ggplot(aes(x = family, y = Percentage, fill = Modeling)) +
#   ggplot(aes(x = family, y = Freq, fill = Modeling)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   #geom_text(aes(x = family, y = Percentage, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   geom_text(aes(x = family, y = Freq, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
#   xlab(label = "Family") +
#   coord_flip() +
#   # labs(subtitle = "All data") +
#   #scale_fill_brewer(palette = "Dark2") +
#   theme_new() +
#   theme(legend.position = "bottom",
#         axis.text.y = element_text(face = "italic")) +
#   NULL

#Fig_S5 <- plot_grid(Fig_S5A, Fig_S5B, Fig_S5C, Fig_S5D, labels = c("A", "B", "C", "D"), label_size = 8, ncol = 2, nrow = 2, rel_heights = c(1, 1, 1, 1))

# ggsave(
#   plot = Fig_S5,
#   filename = "outputs/Fig_S5.pdf",
#   device = "pdf",
#   width = 183,
#   height = 183,
#   units = "mm",
#   useDingbats = FALSE
# )

Fig_S5 <- integrated_df %>%
  filter(superkingdom == "Eukaryota") %>% 
  select(Modeling, family) %>%
  table() %>%
  as.data.frame() %>%
  group_by(family) %>%
  filter(Freq > 0) %>% 
  mutate(Percentage = 100 * (Freq / sum(Freq))) %>%
  #mutate(Labels = paste(Freq, "(", round(Percentage, digits = 1), "%)")) %>%
  mutate(Labels = Freq) %>%
  #ggplot(aes(x = family, y = Percentage, fill = Modeling)) +
  ggplot(aes(x = family, y = Freq, fill = Modeling)) +
  geom_bar(stat = "identity", width = 0.8) +
  #geom_text(aes(x = family, y = Percentage, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
  #geom_text(aes(x = family, y = Freq, label = Labels), colour = "white", position = position_stack(vjust = 0.5), size = 1) +
  xlab(label = "Family") +
  coord_flip() +
  # labs(subtitle = "All data") +
  #scale_fill_brewer(palette = "Dark2") +
  theme_new() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(face = "italic")) +
  NULL

ggsave(
  plot = Fig_S5,
  filename = "outputs/Fig_S5.pdf",
  device = "pdf",
  width = 89,
  height = 189,
  units = "mm",
  useDingbats = FALSE
)

