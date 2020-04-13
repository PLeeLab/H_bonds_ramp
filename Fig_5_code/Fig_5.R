#-------- LIBRARIES -------#
# BiocManager::install("ggtree")
# install.packages("tidyverse")
# install.packages("EnvStats")
# install.packages("ggpubr")
# install.packages("cowplot")
# install.packages("ggcorrplot")
# install.packages("scales")
library("EnvStats")
library("ggpubr")
library("scales") # to access break formatting functions
library("ggtree")
library("tidyverse")
library("ggstance")
library("ggcorrplot")
library("cowplot")
library("ggridges")
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

tree_ref_orfeomes <- read.newick(file = "inputs/tree_Bact_Arch_Fung_refDataSet/RAxML_bestTree.refDataSet_bacteria_archaea_fungi_refined_MOD.tre")
tree_ref_orfeomes$tip.label <- gsub(x = gsub(x = gsub(x = tree_ref_orfeomes$tip.label, pattern = "GCF_", replacement = "GCF-"), pattern = ".._.*", replacement = ""), pattern = "GCF-", replacement = "GCF_")
#plot(tree_ref_orfeomes)


#-------- INPUTS -------#
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

models_df <- rbind(models_df_Bacteria,
                   models_df_Archaea,
                   models_df_Fungi)

models_df$Genome[1]
models_df$Genome <- gsub(x = gsub(x = gsub(x = models_df$Genome, pattern = "GCF_", replacement = "GCF-"), pattern = ".._.*", replacement = ""), pattern = "GCF-", replacement = "GCF_")
# prev: "_.*"
# new: ".._.*"

tax_df$AssemblyAccession <- gsub(x = gsub(x = gsub(x = tax_df$AssemblyAccession, pattern = "GCF_", replacement = "GCF-"), pattern = "\\..", replacement = ""), pattern = "GCF-", replacement = "GCF_")

integrated_df <- merge(x = models_df, y = tax_df, by.x = "Genome", by.y = "AssemblyAccession")

integrated_df$complete_tax <- paste(integrated_df$phylum, integrated_df$class, integrated_df$order, integrated_df$family, integrated_df$genus, integrated_df$species, sep = "__")

ORFeome_features$Genome <- gsub(x = gsub(x = gsub(x = ORFeome_features$Genome, pattern = "GCF_", replacement = "GCF-"), pattern = ".._.*", replacement = ""), pattern = "GCF-", replacement = "GCF_")
# prev: "_.*"

integrated_df <- merge(x = integrated_df, y = ORFeome_features, by.x = "Genome", by.y = "Genome")

#integrated_df <- as_tibble(integrated_df)

#--------- Ramp characteristics ---------#
# Rate, Carrying_capacity, Initial_cost
my_comparisons <- list(c("Archaea", "Eukaryota"), c("Bacteria", "Archaea"), c("Bacteria", "Eukaryota"))

Fig_5B <- ggplot(
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

Fig_5C <- ggplot(
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

Fig_5D <- ggplot(
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

Fig_5E <- ggplot(
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

Fig_5ABCDE <- plot_grid(Fig_5B, Fig_5B, Fig_5C, Fig_5D, Fig_5E, labels = c("A", "B", "C", "D", "E"), nrow = 1, ncol = 5, label_size = 8)



#---------------  Phylogenomics section -----------------#
rownames(integrated_df) <- integrated_df$Genome

groupInfo <- list(Bacteria = integrated_df$Genome[integrated_df$superkingdom=="Bacteria"],
                  Archaea  = integrated_df$Genome[integrated_df$superkingdom=="Archaea"],
                  Fungi    = integrated_df$Genome[integrated_df$superkingdom=="Eukaryota"])

tree_grouped <- groupOTU(.data = tree_ref_orfeomes,
                         .node = groupInfo)

p_tree <- 
  ggtree(tr = tree_grouped,
       mapping = aes(color = group),
       layout = "rectangular") + # layout="equal_angle" | layout = "fan", open.angle=180 | layout = "rectangular"
  geom_tiplab(size=0) +
  #geom_tippoint(aes(color = group)) +
  geom_treescale(x = 0, y = 1000, width = 0.1, offset = 10) +
  guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) +
  theme_tree2(legend.position = c(0.1, 0.88)) +
  #geom_text(aes(label=node)) +
  NULL

p_tree <- flip(tree_view = p_tree, node1 = 2319, node2 = 2952) %>% flip(node1 = 2978, node2 = 3283) %>% flip(node1 = 2979, node2 = 3275)
#rotate_tree(p_tree, 180)
#integrated_df$Rate[is.na(integrated_df$Rate)] <- 0
p1 <- facet_plot(p     = p_tree,
                 panel = "Rate",
                 data  = integrated_df,
                 geom  = geom_point, 
                 mapping = aes(x = Rate, group = superkingdom))

p2 <- facet_plot(p     = p1,
                 panel = "Rate box",
                 data  = integrated_df,
                 geom  = geom_boxploth, 
                 mapping = aes(x = Rate, group = superkingdom)) #class

#-------------------------------------
i2 <- integrated_df %>% dplyr::select(Rate)
p_tree2 <- 
  ggtree(tr = tree_grouped,
         #mapping = aes(color = group),
         layout = "circular") + # layout="equal_angle" | layout = "fan", open.angle=180 | layout = "rectangular"
  geom_tiplab(size = 0.0, align=TRUE, linetype = "solid", linesize = 0.05) + #0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
  geom_tippoint(aes(color = group)) +
  geom_treescale(x = -0.5, y = 0, width = 0.2, offset = 10) +
  guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) +
  theme_tree2(legend.position = c(0.1, 0.88)) +
  #geom_text(aes(label=node)) +
  NULL

p_tree2 <- flip(tree_view = p_tree2, node1 = 2319, node2 = 2952) %>%
  flip(node1 = 2978, node2 = 3283) %>%
  flip(node1 = 2979, node2 = 3275)

Fig_5FG <- plot_grid(gheatmap(p        = p_tree2,
                              data     = i2,
                              offset   = -0.1,
                              width    = 0.2, 
                              colnames = FALSE) +
            scale_fill_viridis_c(option="C", name="Ramp rate", limits = c(min(integrated_df$Rate, na.rm = TRUE), max(integrated_df$Rate, na.rm = TRUE))) +
              theme_new() +
            theme(legend.position = "bottom"),
          
          integrated_df %>%
            dplyr::filter(!is.na(Rate)) %>%
            dplyr::filter(!is.na(phylum)) %>%
            dplyr::mutate(superkingdom2 = str_extract(superkingdom, "^.{1}")) %>% 
            dplyr::mutate(phylum = paste0(superkingdom2, "_", phylum)) %>% 
            dplyr::group_by(phylum) %>% 
            dplyr::filter(n() >= 3) %>%
            ggplot(mapping = aes(y = reorder(x = phylum, FUN = median, X = Rate), x = Rate, fill = stat(x))) +
            geom_density_ridges_gradient() +
            scale_fill_viridis_c(option = "C", limits = c(min(integrated_df$Rate, na.rm = TRUE), max(integrated_df$Rate, na.rm = TRUE))) +
            ylab(label = "Phylum") +
            #coord_flip() +
            #geom_vline(xintercept = 0, lty=2) +
            #theme_bw() +
            theme_new() +
            theme(legend.position = "bottom") +
            NULL,
          labels = c("F", "G"), label_size = 8, ncol = 2, rel_widths = c(1.25, 0.75))

ggsave(plot = plot_grid(Fig_5ABCDE, Fig_5FG, nrow = 2, ncol = 1, rel_heights = c(1/5, 4/5)),
       filename = "outputs/Fig_5.pdf",
       device = "pdf",
       width = 180,
       height = 200,
       units = "mm",
       useDingbats = FALSE)


#----------------------- Tree with namex and taxonomy names for reference
tree_grouped_with_tax <- tree_grouped

#i = 1
for (i in 1:length(tree_grouped_with_tax$tip.label)) {
  if (length(integrated_df$complete_tax[which(tree_grouped_with_tax$tip.label[i]==integrated_df$Genome)])==1) {
    tree_grouped_with_tax$tip.label[i] <- integrated_df$complete_tax[which(tree_grouped_with_tax$tip.label[i]==integrated_df$Genome)]
  } else {
    tree_grouped_with_tax$tip.label[i] <- tree_grouped_with_tax$tip.label[i]
  }
}

p_tree3 <- 
  ggtree(tr = tree_grouped_with_tax,
         #mapping = aes(color = group),
         layout = "circular") + # layout="equal_angle" | layout = "fan", open.angle=180 | layout = "rectangular"
  geom_tiplab(aes(angle=angle), size=1.0, align=TRUE, linetype = "solid", linesize = 0.05) + #0 = blank, 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
  geom_tippoint(aes(color = group)) +
  geom_treescale(x = -0.5, y = 0, width = 0.1, offset = 10) +
  #guides(color = guide_legend(override.aes = list(size = 4, shape = 15))) +
  #theme_tree2(legend.position = c(0.1, 0.88)) +
  #geom_text(aes(label=node)) +
  NULL

p_tree3 <- flip(tree_view = p_tree3, node1 = 2319, node2 = 2952) %>%
  flip(node1 = 2978, node2 = 3283) %>%
  flip(node1 = 2979, node2 = 3275)

ggsave(plot = p_tree3, filename = "outputs/Fig_treeTAX.pdf", device = "pdf", width = 20, height = 20, useDingbats=FALSE)