#install.packages('tidyverse')
#library("corrr")
library("tidyverse")
library("cowplot")
library("ggridges")
library("EnvStats")
library("ggpubr")
#-------------------------#

rm(list=ls()); # clear environment

theme_new <- function(){
  theme_bw() %+replace%
    theme(
      plot.title = element_text(   color = "black", size = 7, face = "plain", hjust = 0),
      axis.title = element_text(   color = "black", size = 7),
      axis.text = element_text(    color = "black", size = 6),
      legend.text = element_text(  color = "black", size = 5),
      legend.title = element_text( color = "black", size = 5),
      strip.text = element_text(   color = "black", size = 6), 
      #strip.background = element_blank(),
      #panel.grid.minor = element_blank(),
      #panel.grid.major = element_blank(),
      complete = TRUE
    )
}
#-------------------------#

pd_folding_energy_BACT <- read_csv(file = "inputs/mRNA/complete_fold_prob_dataset_BACT.csv") %>% 
  mutate(Domain = "Bacteria") 
pd_folding_energy_ARCH <- read_csv(file = "inputs/mRNA/complete_fold_prob_dataset_ARCH.csv") %>% 
  mutate(Domain = "Archaea")
pd_folding_energy_FUNG <- read_csv(file = "inputs/mRNA/complete_fold_prob_dataset_FUNG.csv") %>% 
  mutate(Domain = "Eukaryote")

pd_folding_energy_dataset <- bind_rows(pd_folding_energy_BACT,
                                       pd_folding_energy_ARCH,
                                       pd_folding_energy_FUNG) %>% 
  mutate(Metric = "(A) P(mRNA folding)") %>% 
  mutate(Mean  = 1 - Mean) %>% 
  mutate(Upper = 1 - Upper) %>% 
  mutate(Lower = 1 - Lower) %>% 
  rename(ORFeome = Organism) %>% 
  pivot_wider(names_from = Metric, values_from = c(Mean, Lower, Upper))


pd_hydrogen_bonding_BACT <- read_tsv(file = "inputs/hb/OUT_Bacteria_Elements_df_boot_ALL_1_to_13922.txt") %>% 
  mutate(Domain = "Bacteria") %>% 
  filter(Element == "Hydrogen bonds")
pd_hydrogen_bonding_ARCH <- read_tsv(file = "inputs/hb/OUT_Archaea_Elements_df_boot_ALL_1_to_297.txt") %>% 
  mutate(Domain = "Archaea") %>% 
  filter(Element == "Hydrogen bonds")
pd_hydrogen_bonding_FUNG <- read_tsv(file = "inputs/hb/OUT_Fungi_Elements_df_boot_ALL_1_to_293.txt") %>% 
  mutate(Domain = "Eukaryote") %>% 
  filter(Element == "Hydrogen bonds")

pd_hydrogen_bonding_dataset <- bind_rows(pd_hydrogen_bonding_BACT,
                                         pd_hydrogen_bonding_ARCH,
                                         pd_hydrogen_bonding_FUNG) %>% 
  mutate(Metric = "(B) H bonding") %>% 
  select(-Element) %>% 
  rename(ORFeome = Genome) %>% 
  pivot_wider(names_from = Metric, values_from = c(Mean, Lower, Upper))


pd_metrics_dataset <- left_join(pd_folding_energy_dataset, pd_hydrogen_bonding_dataset)
#length(unique(pd_metrics_dataset$ORFeome))
glimpse(pd_metrics_dataset)

tax_df <- bind_rows(read_tsv(file = "inputs/taxonomy/integrated_tax_df_Bacteria_genomes.txt"),
                    read_tsv(file = "inputs/taxonomy/integrated_tax_df_Archaea_genomes.txt"),
                    read_tsv(file = "inputs/taxonomy/integrated_tax_df_Fungi_genomes.txt"))
tax_df$AssemblyAccession <- gsub(x = gsub(x = gsub(x = tax_df$AssemblyAccession, pattern = "GCF_", replacement = "GCF-"), pattern = "\\..", replacement = ""), pattern = "GCF-", replacement = "GCF_")

pd_metrics_dataset$ORFeome <- gsub(x = gsub(x = gsub(x = pd_metrics_dataset$ORFeome, pattern = "GCF_", replacement = "GCF-"), pattern = ".._.*", replacement = ""), pattern = "GCF-", replacement = "GCF_")

pd_complete_dataset <- left_join(pd_metrics_dataset, tax_df, by = c("ORFeome" = "AssemblyAccession") )
#length(unique(pd_complete_dataset$ORFeome))
#glimpse(pd_complete_dataset)
#colnames(pd_complete_dataset)

ORFeome_features <- bind_rows(
  read_tsv(file = "inputs/ORFeome_features/ORFeome_features_Archaea.txt"),
  read_tsv(file = "inputs/ORFeome_features/ORFeome_features_Bacteria.txt"),
  read_tsv(file = "inputs/ORFeome_features/ORFeome_features_Fungi.txt")
) %>% rename(ORFeome = Genome)
ORFeome_features$ORFeome <- gsub(x = gsub(x = gsub(x = ORFeome_features$ORFeome, pattern = "GCF_", replacement = "GCF-"), pattern = ".._.*", replacement = ""), pattern = "GCF-", replacement = "GCF_")

pd_complete_dataset <- left_join(pd_complete_dataset, ORFeome_features, by = c("ORFeome") )
#length(unique(pd_complete_dataset$ORFeome))
#glimpse(pd_complete_dataset)
#colnames(pd_complete_dataset)

correlations_dataframe <- bind_rows(pd_complete_dataset %>%
  filter(Position > 1 & Position <= 120) %>% 
  group_by(ORFeome) %>% 
  summarise(Pearson = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "pearson"),
            Spearman = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "spearman")) %>% 
  mutate(Range = "AP2-P120"),
pd_complete_dataset %>%
  filter(Position > 1 & Position <= 60) %>% 
  group_by(ORFeome) %>% 
  summarise(Pearson = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "pearson"),
            Spearman = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "spearman")) %>% 
  mutate(Range = "BP2-P60"),
pd_complete_dataset %>%
  filter(Position > 1 & Position <= 30) %>% 
  group_by(ORFeome) %>% 
  summarise(Pearson = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "pearson"),
            Spearman = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "spearman")) %>% 
  mutate(Range = "CP2-P30"),
pd_complete_dataset %>%
  filter(Position > 1 & Position <= 15) %>% 
  group_by(ORFeome) %>% 
  summarise(Pearson = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "pearson"),
            Spearman = cor(`Mean_(A) P(mRNA folding)`, `Mean_(B) H bonding`, use = "pairwise.complete.obs", method = "spearman")) %>% 
  mutate(Range = "DP2-P15"))

correlations_dataframe <- left_join(correlations_dataframe, tax_df, by = c("ORFeome" = "AssemblyAccession") )
correlations_dataframe <- left_join(correlations_dataframe, ORFeome_features )
#colnames(correlations_dataframe)

n_of_samples <- correlations_dataframe %>% 
  drop_na(superkingdom) %>% 
  mutate(Nucleus = ifelse(test = superkingdom=="Bacteria" | superkingdom=="Archaea", yes = "Prokaryote", no = "Eukaryote")) %>% 
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Cor_index", values_to = "Value") %>% 
  filter(Cor_index == "Pearson") %>% 
  filter(Range == "AP2-P120") %>% 
  group_by(superkingdom) %>% 
  rename(Domain = superkingdom) %>% 
  summarise(n = n())

Fig_3B <-
  correlations_dataframe %>% 
  drop_na(superkingdom) %>% 
  mutate(Nucleus = ifelse(test = superkingdom=="Bacteria" | superkingdom=="Archaea", yes = "Prokaryote", no = "Eukaryote")) %>% 
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Cor_index", values_to = "Value") %>% 
  #filter(Cor_index == "Pearson") %>% 
  #filter(Range == "AP2-P120") %>% 
  ggplot(mapping = aes(x = Value, fill=superkingdom, colour=superkingdom)) +
  #geom_vline(xintercept = c(-0.5, 0.0, 0.5), lty=2, colour="black") +
  geom_vline(xintercept = c(0.0), lty=2, colour="black") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  geom_density(alpha=1/2) +
  #geom_histogram(alpha=1/2) +
  #geom_freqpoly(alpha=1) +
  annotate("text", x = 0.6, y = 2.5, label = paste0("n = ", n_of_samples %>% filter(Domain == "Archaea") %>% select(n) ), colour = "red", size=3) +
  annotate("text", x = 0.7, y = 8, label = paste0("n = ", n_of_samples %>% filter(Domain == "Bacteria") %>% select(n) ), colour = "blue", size=3 ) +
  annotate("text", x = -0.2, y = 2.0, label = paste0("n = ", n_of_samples %>% filter(Domain == "Eukaryota") %>% select(n) ), colour = "darkgreen", size=3 ) +
  xlim(-1, 1) +
  xlab(label = "Correlation between H bonding and P(mRNA folding)") +
  #theme_bw() +
  #xlab(label = "Codon position (5'-end)") +
  #ylab(label = "Mean + 95% CI") +
  #facet_wrap(~Cor_index) +
  facet_grid(Range~Cor_index) +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

my_comparisons <- list( c("Eukaryote", "Prokaryote") )
Fig_3C <- correlations_dataframe %>% 
  drop_na(superkingdom) %>% 
  mutate(Nucleus = ifelse(test = superkingdom=="Bacteria" | superkingdom=="Archaea", yes = "Prokaryote", no = "Eukaryote")) %>% 
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Cor_index", values_to = "Value") %>% 
  ggplot(mapping = aes(x = Nucleus, y = Value, fill=Nucleus)) +
  #geom_hline(yintercept = c(-0.5, 0.0, 0.5), lty=1, colour="grey") +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2") +
  #geom_point() +
  geom_boxplot() +
  stat_n_text(y.expand.factor = 0.05, geom = "text", fontface = "italic", angle = 0, hjust = 0.5, size = 2, color = "grey") +
  stat_compare_means(comparisons = my_comparisons, angle=90, label = "p.signif", method = "t.test")+ # Add pairwise comparisons p-value
  #stat_compare_means(label = "p.format", method = "t.test", label.x = 1.5) + # label =  "p.signif"
  coord_flip() +
  ylab(label = "Correlation between H bonding and P(mRNA folding)") +
  xlab(label = "") +
  ylim(c(-1.3, 1.3)) +
  #theme_bw() +
  facet_grid(Range~Cor_index) +
  theme_new() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  NULL

Fig_3D1 <- correlations_dataframe %>% 
  drop_na(superkingdom) %>% 
  drop_na(class) %>% 
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Cor_index", values_to = "Value") %>% 
  filter(Range == "AP2-P120") %>% 
  filter(Cor_index == "Pearson") %>% 
  group_by(class, Range) %>% 
  filter(n() >= 3) %>% 
  ggplot(mapping = aes(y = reorder(x = class, FUN = median, X = Value), x = Value, fill = superkingdom )) +
  geom_vline(xintercept = c(-0.5, 0, 0.5), lty=2, colour="grey") +
  #geom_boxplot() +
  geom_density_ridges(alpha=1/2) +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(Range~Cor_index, nrow=1) +
  ylab(label = "Class") +
  xlab(label = "Pearson correlation") +
  #theme_bw() +
  theme_new() +
  theme(legend.position = "bottom") +
  NULL

Fig_3D2 <- correlations_dataframe %>% 
  drop_na(superkingdom) %>% 
  drop_na(class) %>% 
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Cor_index", values_to = "Value") %>% 
  filter(Range == "AP2-P120") %>% 
  filter(Cor_index == "Pearson") %>% 
  group_by(class, Range) %>% 
  filter(n() >= 3) %>% 
  ggplot(mapping = aes(x = reorder(x = class, FUN = median, X = Value), y = GC_mean, fill = superkingdom)) +
  geom_boxplot() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(Range~Cor_index, nrow=1) +
  coord_flip() +
  #theme_bw() +
  theme_new() +
  theme(legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  NULL

Fig_3D3 <- correlations_dataframe %>% 
  drop_na(superkingdom) %>% 
  drop_na(class) %>% 
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Cor_index", values_to = "Value") %>% 
  filter(Range == "AP2-P120") %>% 
  filter(Cor_index == "Pearson") %>% 
  group_by(class, Range) %>% 
  filter(n() >= 3) %>% 
  ggplot(mapping = aes(x = reorder(x = class, FUN = median, X = Value), y = GC3_mean, fill = superkingdom)) +
  geom_boxplot() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(Range~Cor_index, nrow=1) +
  coord_flip() +
  #theme_bw() +
  theme_new() +
  theme(legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  NULL

Fig_3D4 <- correlations_dataframe %>% 
  drop_na(superkingdom) %>% 
  drop_na(class) %>% 
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Cor_index", values_to = "Value") %>% 
  filter(Range == "AP2-P120") %>% 
  filter(Cor_index == "Pearson") %>% 
  group_by(class, Range) %>% 
  filter(n() >= 3) %>% 
  ggplot(mapping = aes(x = reorder(x = class, FUN = median, X = Value), y = GC3_mean/GC_mean, fill = superkingdom)) +
  geom_boxplot() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  facet_wrap(Range~Cor_index, nrow=1) +
  coord_flip() +
  #theme_bw() +
  theme_new() +
  theme(legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  NULL

Fig_3_TOP <- plot_grid(Fig_3B, Fig_3B, Fig_3C, labels = c("A", "B", "C"), label_size = 8, ncol = 3, nrow = 1)

Fig_3_BOTTOM <- plot_grid(Fig_3D1, Fig_3D2, Fig_3D3, Fig_3D4, labels = c("C", "", "", ""), label_size = 8, ncol = 4, nrow = 1, rel_widths = c(0.4, 0.2, 0.2, 0.2))

Fig_3 <- plot_grid(Fig_3_TOP, Fig_3_BOTTOM, labels = c("", ""), ncol = 1, nrow = 2, rel_heights = c(0.8, 1.2))

dir.create(path = "outputs", showWarnings = TRUE)

ggsave(plot = Fig_3,
       filename = "outputs/Fig_3.pdf",
       device = "pdf",
       width = 180,
       height = 200,
       units = "mm",
       useDingbats = FALSE)


