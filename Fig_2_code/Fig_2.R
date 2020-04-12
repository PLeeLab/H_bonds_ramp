#install.packages('tidyverse')
#install.packages("Hmisc")
#install.packages("RcmdrMisc")
library("tidyverse")
library("corrr")
library("cowplot")
library("Hmisc")
library("RcmdrMisc")
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
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      complete = TRUE
    )
}
#-------------------------#

pd_RSCU <- read_csv(file = "inputs/complete_rscu_dataset.csv")
pd_RSCU$Metric <- "(A) RSCU"

pd_W <- read_csv(file = "inputs/complete_w_dataset.csv")
pd_W$Metric <- "(C) CAI"

pd_folding_energy_dataset <- read_csv(file = "inputs/complete_fold_prob_dataset.csv")
pd_folding_energy_dataset$Mean  <- 1 - pd_folding_energy_dataset$Mean
pd_folding_energy_dataset$Upper <- 1 - pd_folding_energy_dataset$Upper
pd_folding_energy_dataset$Lower <- 1 - pd_folding_energy_dataset$Lower
pd_folding_energy_dataset$Metric <- "(B) P(mRNA folding)"

pd_tAI_bact_ecoli  <- read_csv(file = "inputs/wtai_tibble_Escherichia_coli.csv")
pd_tAI_bact_verru  <- read_csv(file = "inputs/wtai_tibble_Methylacidiphilum_kamchatkense_Kam1.csv")
pd_tAI_arch_metha  <- read_csv(file = "inputs/wtai_tibble_Methanosarcina_acetivorans_c2a.ASM734v1.csv")
pd_tAI_arch_halof  <- read_csv(file = "inputs/wtai_tibble_Haloferax_volcanii.csv")
pd_tAI_fung_sacch  <- read_csv(file = "inputs/wtai_tibble_Saccharomyces_cerevisiae.csv")
pd_tAI_fung_schiz  <- read_csv(file = "inputs/wtai_tibble_Schizosaccharomyces_pombe.csv")
pd_tAI_fung_neuro  <- read_csv(file = "inputs/wtai_tibble_Neurospora_crassa_73.csv")
pd_tAI_fung_asper  <- read_csv(file = "inputs/wtai_tibble_Aspergillus_nidulans.csv")
pd_tAI_dataset <- bind_rows(pd_tAI_bact_ecoli, pd_tAI_bact_verru, pd_tAI_arch_metha, pd_tAI_arch_halof, pd_tAI_fung_sacch, pd_tAI_fung_schiz, pd_tAI_fung_neuro, pd_tAI_fung_asper)
pd_tAI_dataset$Metric <- "(D) tAI"

unique(pd_folding_energy_dataset$Organism)
pd_HBs_bact_ecoli <- read_csv(file = "inputs/hbs_Escherichia_coli.csv")
pd_HBs_bact_ecoli$Organism <- "Escherichia_coli"
pd_HBs_bact_verru <- read_csv(file = "inputs/hbs_Methylacidiphilum_kamchatkense_Kam1.csv")
pd_HBs_bact_verru$Organism <- "Methylacidiphilum_kamchatkense_Kam1"
pd_HBs_arch_metha <- read_csv(file = "inputs/hbs_Methanosarcina_acetivorans_c2a.csv")
pd_HBs_arch_metha$Organism <- "Methanosarcina_acetivorans_c2a.ASM734v1"
pd_HBs_arch_halof <- read_csv(file = "inputs/hbs_Haloferax_volcanii.csv")
pd_HBs_arch_halof$Organism <- "Haloferax_volcanii"
pd_HBs_fung_sacch <- read_csv(file = "inputs/hbs_Saccharomyces_cerevisiae.csv")
pd_HBs_fung_sacch$Organism <- "Saccharomyces_cerevisiae" 
pd_HBs_fung_schiz <- read_csv(file = "inputs/hbs_Schizosaccharomyces_pombe.csv")
pd_HBs_fung_schiz$Organism <- "Schizosaccharomyces_pombe"
pd_HBs_fung_neuro <- read_csv(file = "inputs/hbs_Neurospora_crassa_73.csv")
pd_HBs_fung_neuro$Organism <- "Neurospora_crassa_73"
pd_HBs_fung_asper <- read_csv(file = "inputs/hbs_Aspergillus_nidulans.csv")
pd_HBs_fung_asper$Organism <- "Aspergillus_nidulans"

pd_HBs_dataset <- bind_rows(pd_HBs_bact_ecoli, pd_HBs_bact_verru, pd_HBs_arch_metha, pd_HBs_arch_halof, pd_HBs_fung_sacch, pd_HBs_fung_schiz, pd_HBs_fung_neuro, pd_HBs_fung_asper)
pd_HBs_dataset$Element <- NULL
pd_HBs_dataset$Metric <- "(E) H bonds"
  
integrated_metrics_df <- bind_rows(pd_RSCU,
                                   pd_W,
                                   pd_folding_energy_dataset,
                                   pd_tAI_dataset,
                                   pd_HBs_dataset)

integrated_metrics_df <- integrated_metrics_df %>% 
  mutate(sk = ifelse(test = Organism == "Escherichia_coli" | Organism == "Methylacidiphilum_kamchatkense_Kam1" | Organism == "Haloferax_volcanii" | Organism == "Methanosarcina_acetivorans_c2a.ASM734v1",
                     yes = "Prok", no = "Euk"))

Fig_2A <- ggplot(data = integrated_metrics_df %>%
         filter(Position > 1 & Position < 100), # 0 = no codon is excluded, 1 = only start codon is excluded,
       mapping = aes(x = Position, colour=Organism)) +
  #geom_vline(xintercept = 1, lty=2, colour="grey") +
  #geom_hline(yintercept = 1, lty=2, colour="grey") +
  geom_errorbar(mapping = aes(ymin=Lower, ymax=Upper), size=0.1) +
  geom_point(mapping = aes(y=Mean), size=0.15) +
  geom_line(mapping = aes(y=Mean), size=0.5) +
  scale_color_brewer(palette="Dark2") +
  #geom_smooth(aes(y=Mean))+
  xlab(label = "Codon position (5'-end)") +
  ylab(label = "Mean + 95% CI") +
  facet_wrap(Metric~Organism, scales = "free_y", ncol = 8) +
  #facet_wrap(Organism~Metric, scales = "free_y", ncol = 5) +
  #facet_grid(Metric~Organism, scales = "free_y") +
  theme_new() +
  theme(legend.position = "none") +
  NULL

ggsave(plot = Fig_2A,
       filename = "outputs/Fig_2A.pdf",
       device = "pdf",
       width = 180,
       height = 90,
       units = "mm",
       useDingbats = FALSE)

#-------------------------- Correlation analysis -------------------------#
plot_corr_net_fun <- function(CURRENT_ORGANISM = 1) {
  integrated_metrics_df %>% 
    dplyr::filter(Organism == unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]) %>% 
    dplyr::filter(Position > 1 & Position < 100) %>% 
    dplyr::select(-Lower, -Upper, -sk) %>% 
    tidyr::pivot_wider(names_from = Metric, values_from = Mean) %>%
    #group_by(Organism) %>%
    dplyr::select(-Organism, -Position) %>% 
    corrr::correlate(method = "spearman", use ="pairwise.complete.obs") %>% # spearman, kendall, pearson
    corrr::network_plot(curved = TRUE, repel = TRUE, min_cor = 0) +
    ggplot2::ggtitle(label = unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM])
}

Fig_2B <- plot_grid(plot_corr_net_fun(CURRENT_ORGANISM = 1),
                    plot_corr_net_fun(CURRENT_ORGANISM = 2),
                    plot_corr_net_fun(CURRENT_ORGANISM = 3),
                    plot_corr_net_fun(CURRENT_ORGANISM = 4),
                    plot_corr_net_fun(CURRENT_ORGANISM = 5),
                    plot_corr_net_fun(CURRENT_ORGANISM = 6),
                    plot_corr_net_fun(CURRENT_ORGANISM = 7),
                    plot_corr_net_fun(CURRENT_ORGANISM = 8), nrow = 4, col = 2, labels = "B", label_size = 8)

ggsave(plot = Fig_2B,
       filename = "outputs/Fig_2B.pdf",
       device = "pdf",
       width = 300,
       height = 400,
       units = "mm",
       useDingbats = FALSE)

#-------------- Figure 2C1 --- Spearman's --------- focus on correlation between the ramp of Hbonds and other ramps
plot_corr_barplot_fun <- function(CURRENT_ORGANISM = 1) {
  integrated_metrics_df %>% 
    dplyr::filter(Organism == unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]) %>% 
    dplyr::filter(Position > 1 & Position < 100) %>% 
    dplyr::select(-Lower, -Upper, -sk) %>% 
    tidyr::pivot_wider(names_from = Metric, values_from = Mean) %>%
    #group_by(Organism) %>%
    dplyr::select(-Organism, -Position) %>% 
    corrr::correlate(method = "spearman", use ="pairwise.complete.obs") %>% # spearman, kendall, pearson
    corrr::focus(`(E) H bonds`) %>% 
    #dplyr::mutate(rowname = reorder(rowname, `(E) H bonds`)) %>%
    ggplot(aes(x = rowname, y = `(E) H bonds`)) +
    geom_hline(yintercept = 0, lty = 1, colour = "black") +
    geom_col(aes(fill = `(E) H bonds` >= 0)) + 
    ylim(c(-1,1))+
    geom_text(aes(label = round(`(E) H bonds`, 2)), size=2) +
    scale_fill_brewer(palette="Set1") +
    xlab(label = "") +
    coord_flip() +
    ggtitle(label = unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]) +
    theme_bw() +
    theme_new() +
    theme(legend.position = "none")
}

Fig_2C_spearman <- plot_grid(plot_corr_barplot_fun(CURRENT_ORGANISM = 1),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 2),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 3),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 4),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 5),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 6),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 7),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 8),
                    nrow = 2, ncol = 4)

ggsave(plot = Fig_2C_spearman,
       filename = "outputs/Fig_2C1.pdf",
       device = "pdf",
       width = 180,
       height = 60,
       units = "mm",
       useDingbats = FALSE)

#---- Get the adjusted P value from the multiple correlation analyses ----#
unique(integrated_metrics_df$Organism)
CURRENT_ORGANISM = 1 # Need to manually iterate over each organism number
unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]
integrated_metrics_df %>% 
  dplyr::filter(Organism == unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]) %>% 
  dplyr::filter(Position > 1 & Position < 100) %>% 
  dplyr::select(-Lower, -Upper, -sk) %>% 
  tidyr::pivot_wider(names_from = Metric, values_from = Mean) %>%
  #group_by(Organism) %>%
  dplyr::select(-Organism, -Position) %>% 
  as.matrix() %>%
  #Hmisc::rcorr(type = "pearson") %>% 
  RcmdrMisc::rcorr.adjust(type = "spearman") # spearman, pearson

#-------------- Figure 2C2 - Pearson's --------- focus on correlation between the ramp of Hbonds and other ramps
plot_corr_barplot_fun <- function(CURRENT_ORGANISM = 1) {
  integrated_metrics_df %>% 
    dplyr::filter(Organism == unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]) %>% 
    dplyr::filter(Position > 1 & Position < 100) %>% 
    dplyr::select(-Lower, -Upper, -sk) %>% 
    tidyr::pivot_wider(names_from = Metric, values_from = Mean) %>%
    #group_by(Organism) %>%
    dplyr::select(-Organism, -Position) %>% 
    corrr::correlate(method = "pearson", use ="pairwise.complete.obs") %>% # spearman, kendall, pearson
    corrr::focus(`(E) H bonds`) %>% 
    #dplyr::mutate(rowname = reorder(rowname, `(E) H bonds`)) %>%
    ggplot(aes(x = rowname, y = `(E) H bonds`)) +
    geom_hline(yintercept = 0, lty = 1, colour = "black") +
    geom_col(aes(fill = `(E) H bonds` >= 0)) + 
    ylim(c(-1,1))+
    geom_text(aes(label = round(`(E) H bonds`, 2)), size=2) +
    scale_fill_brewer(palette="Set1") +
    xlab(label = "") +
    coord_flip() +
    ggtitle(label = unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]) +
    theme_bw() +
    theme_new() +
    theme(legend.position = "none")
}

Fig_2C_pearson <- plot_grid(plot_corr_barplot_fun(CURRENT_ORGANISM = 1),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 2),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 3),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 4),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 5),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 6),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 7),
                    plot_corr_barplot_fun(CURRENT_ORGANISM = 8),
                    nrow = 2, ncol = 4)


ggsave(plot = Fig_2C_pearson,
       filename = "outputs/Fig_2C2.pdf",
       device = "pdf",
       width = 180,
       height = 60,
       units = "mm",
       useDingbats = FALSE)

#---- Get the adjusted P value from the multiple correlation analyses ----#
unique(integrated_metrics_df$Organism)
CURRENT_ORGANISM = 1
unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]
integrated_metrics_df %>% 
  dplyr::filter(Organism == unique(integrated_metrics_df$Organism)[CURRENT_ORGANISM]) %>% 
  dplyr::filter(Position > 1 & Position < 100) %>% 
  dplyr::select(-Lower, -Upper) %>% 
  tidyr::pivot_wider(names_from = Metric, values_from = Mean) %>%
  #group_by(Organism) %>%
  dplyr::select(-Organism, -Position) %>% 
  as.matrix() %>%
  #Hmisc::rcorr(type = "pearson") %>% 
  RcmdrMisc::rcorr.adjust(type = "pearson") # spearman, pearson

#------ GAM regression analysis ----#
Fig_2D <- integrated_metrics_df %>%
  filter(Position > 1 & Position < 100) %>% 
  filter(Metric == "(B) P(mRNA folding)" | Metric == "(E) H bonds") %>% 
  #mutate(Organism = paste0(sk, "_", Organism)) %>% 
  pivot_wider(names_from = Metric, values_from = c(Mean, Lower, Upper)) %>% 
  #ggplot(mapping = aes(x = `Mean_(E) H bonds`, y = `Mean_(B) P(mRNA folding)`, colour = Organism, shape = sk)) +
  #ggplot(mapping = aes(x = `Mean_(E) H bonds`, y = `Mean_(B) P(mRNA folding)`, colour = sk)) +
  ggplot(mapping = aes(x = `Mean_(E) H bonds`, y = `Mean_(B) P(mRNA folding)`, colour = Organism)) +
  geom_point(size=0.15) +
  geom_errorbar(mapping = aes(ymin  = `Lower_(B) P(mRNA folding)`,  ymax = `Upper_(B) P(mRNA folding)`), size=0.1) +
  geom_errorbarh(mapping = aes(xmin = `Lower_(E) H bonds`,          xmax = `Upper_(E) H bonds`), size=0.1) +
  geom_smooth(method = "gam") +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(~Organism, nrow = 1, scales = "free") +
  #facet_wrap(~Organism, scales = "free") +
  theme_bw() +
  theme_new() +
  theme(legend.position = "none") +
  NULL

ggsave(plot = Fig_2D,
       filename = "outputs/Fig_2D.pdf",
       device = "pdf",
       width = 180,
       height = 30,
       units = "mm",
       useDingbats = FALSE)
