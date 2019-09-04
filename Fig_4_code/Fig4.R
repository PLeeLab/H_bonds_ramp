#-------- LIBRARIES -------#
library("tidyverse")
library("cowplot")

rm(list=ls()); # clear environment

theme_new <- function(){
  theme_bw() %+replace%
    theme(
      plot.title = element_text(   color = "black", size = 6, hjust = 0),
      axis.title = element_text(   color = "black", size = 6),
      axis.text = element_text(    color = "black", size = 6),
      legend.text = element_text(  color = "black", size = 6),
      legend.title = element_text( color = "black", size = 6),
      strip.text = element_text(   color = "black", size = 6),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      complete = TRUE
    )
}

bacteria_df <- read.delim("inputs/Integrated_all_df_Bacteria.txt")
length(unique(bacteria_df$Organism))
bacteria_df$superkingdom <- "Bacteria"
archaea_df <- read.delim("inputs/Integrated_all_df_Archaea.txt")
length(unique(archaea_df$Organism))
archaea_df$superkingdom <- "Archaea"

Integrated_all_df <- rbind(bacteria_df, archaea_df)

normalit <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

# GCF_000005845.2	511145	GCF_000005845.2	ASM584v2	Escherichia coli str. K-12 substr. MG1655 (E. coli)	562	Escherichia coli	GCA_000005845.2GCF_000005845.2identical	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia	Escherichia coli
Fig_4A <- ggplot(data = Integrated_all_df %>% 
                   filter(Element == "Hydrogen bonds" & Position != 1 & Organism == "GCF_000005845.2_ASM584v2"),
                 mapping = aes(x=Position, y=z_square)) +
  #geom_col(show.legend = FALSE) +
  geom_linerange(aes(ymin=0, ymax=z_square), lwd=0.25) +
  #geom_point() +
  xlab(label = "Codon position (5'-end)")+
  ylab(expression(z^2)) +
  labs(title = "Escherichia coli str. K-12 substr. MG1655") +
  theme_new() +
  NULL

Fig_S4A <- Integrated_all_df %>%
  #group_by(Organism, Element) %>%
  filter(Element == "Hydrogen bonds" & Position != 1) %>% 
  group_by(Organism) %>%
  mutate(norm_Z = normalit(z_square)) %>% 
  na.omit() %>% 
  ggplot(aes(x=Position, y=Organism)) +
  geom_tile(aes(fill = norm_Z)) +
  scale_fill_viridis_c() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_wrap(~superkingdom, ncol=1, scales = "free_y") + #scales = "free_y", 
  ylab("Organism") +
  xlab(label = "Codon position (5'-end)")+
  labs(fill=expression(Normalized~z^2)) +
  theme_new() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom", legend.key.size = unit(4, "mm")) +
  NULL

Fig_4B <- Integrated_all_df %>%
  #group_by(Organism, Element) %>%
  filter(Element == "Hydrogen bonds" & Position != 1) %>% 
  group_by(Organism) %>%
  mutate(norm_Z = normalit(z_square)) %>% 
  #na.omit() %>% 
  ggplot(aes(x=Position, y=norm_Z, colour=superkingdom, shape=superkingdom)) +
  stat_summary(fun.data = "mean_cl_boot", size=0.005) +
  #scale_colour_brewer(palette = "Set1") +
  #facet_wrap(~superkingdom, ncol=1, scales = "free_y") + #scales = "free_y", 
  ylab(expression(Normalized~z^2)) +
  xlab(label = "Codon position (5'-end)")+
  theme_new() +
  theme(legend.position = c(0.5, 0.75),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  NULL

Fig_4C <- ggplot(data = Integrated_all_df %>% 
                   filter(Element == "Hydrogen bonds" & Position != 1 & Organism == "GCF_000005845.2_ASM584v2"),
                 mapping = aes(x=Position, y=chi_gram)) +
  geom_hline(yintercept = c(0) , lty=2, lwd=0.25) +
  geom_linerange(aes(ymin=0, ymax=chi_gram), lwd=0.25) +
  xlab(label = "Codon position (5'-end)")+
  ylab(expression( over(Observed-Expected, sqrt(Expected) ))) + #"Observed - Expected / sqrt(Expected)"
  labs(title = "Escherichia coli str. K-12 substr. MG1655") +
  #facet_wrap(~Element, scales = "free_x", ncol = 6) +
  theme_new() +
  NULL

Integrated_all_df2 <- Integrated_all_df %>%
  #group_by(Organism, Element) %>%
  filter(Element == "Hydrogen bonds" & Position != 1) %>% 
  group_by(Organism) %>%
  mutate(scaled_Chigram = scale(chi_gram)) %>% 
  na.omit()

Fig_S4B <- ggplot(data = Integrated_all_df2, aes(x=Position, y=Organism)) +
  geom_tile(aes(fill = scaled_Chigram)) +
  #scale_fill_distiller(palette = "Spectral") +
  #scale_fill_gradientn(colours = terrain.colors(10)) +
  #scale_fill_gradientn(colours=c("red","white", "green")) +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(-1,1)*max(abs(Integrated_all_df2$scaled_Chigram), na.rm = TRUE),
                       name="Scaled chi gram") +
  facet_wrap(~superkingdom, ncol=1, scales = "free_y") + #scales = "free_y", 
  ylab("Organism") +
  xlab(label = "Codon position (5'-end)")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_new()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom", legend.key.size = unit(4, "mm")) +
  NULL

Fig_4D <- ggplot(data = Integrated_all_df2, aes(x=Position, y=scaled_Chigram, colour=superkingdom, shape=superkingdom)) +
  geom_hline(yintercept = c(0) , lty=2, lwd=0.25) +
  stat_summary(fun.data = "mean_cl_boot", size=0.005) +
  #scale_colour_brewer(palette = "Set1") +
  #facet_wrap(~superkingdom, ncol=1, scales = "free_y") + #scales = "free_y", 
  ylab(expression( over(Observed-Expected, sqrt(Expected) ))) + # Scaled "Observed - Expected / sqrt(Expected)"
  xlab(label = "Codon position (5'-end)")+
  theme_new() +
  theme(legend.position = c(0.5, 0.5),
        legend.background = element_blank(),
        legend.title = element_blank()) +
  NULL

Fig_S4 <- plot_grid(Fig_S4A, Fig_S4B, labels = c("A", "B"), nrow = 1, rel_widths = c(1, 1), label_size = 8)
ggsave(plot = Fig_S4,
       filename = "outputs/Fig_S4.pdf",
       device = "pdf",
       width = 89,
       height = 120,
       units = "mm",
       useDingbats = FALSE)

Fig_4 <- plot_grid(Fig_4A, Fig_4B, Fig_4C, Fig_4D, labels = c("A", "B", "C", "D"), ncol = 2, rel_widths = c(0.5, 0.5, 1, 1), label_size = 8, align = "hv", axis = "tblr")

ggsave(plot = Fig_4,
       filename = "outputs/Fig_4.pdf",
       device = "pdf",
       width = 89,
       height = 70,
       units = "mm",
       useDingbats = FALSE)