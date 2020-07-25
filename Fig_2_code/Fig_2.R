library("tidyverse")
library("Biostrings")
library("seqinr")
library("ggpubr")
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
      #panel.grid.minor = element_blank(),
      #panel.grid.major = element_blank(),
      complete = TRUE
    )
}

ORFs_file <- readDNAStringSet(filepath = "inputs/predicted_orfs_7929")
operon_data <- read_delim(file = "inputs/list_of_operons_7929", delim = "\t")
ORGANISM_NAME = "Escherichia_coli_str_k_12_substr_mg1655"

new_operon_id <- c()
for (i in operon_data$Operon) {
  if (!is.na(i)) { operon_id <- i }
  #print(i)
  #print(operon_id)
  new_operon_id <- c(new_operon_id, operon_id)
}

NUMBER_OF_CODONS_TO_ANALYZE = 100 # 100
MINIMUM_NUMBER_OF_ORFs_PER_OPERON = 2

operons_tibble <- operon_data %>% 
  mutate(Operon_ID = new_operon_id) %>% # bind the column with the operon numbers
  drop_na(IdGene) %>%                           # remove the rows that were empty because of the numbering of the operons
  select(-Operon) %>%                           # remove the old numbering of the operons
  select(Operon_ID, IdGene:Function) %>%    # reorder the tibble so the new numbering appears upfront
  filter( Type == "CDS" ) %>%               # make sure to subset only CDSs
  group_by(Operon_ID) %>%                   # count the number of ORFs in each operon
  filter( n() >= MINIMUM_NUMBER_OF_ORFs_PER_OPERON) %>% # keep only operons with more thant X ORFs
  mutate(IdGene = gsub(x = IdGene,
                       pattern = "CDS:",
                       replacement = "")) %>%  # simplify the name of ORFs
  mutate(ORF_length = postRight - PosLeft) %>% # calculate the lengt of every ORF
  group_by(Operon_ID) %>%                  # group by operon id
  filter(all(ORF_length >= NUMBER_OF_CODONS_TO_ANALYZE*3)) # keep only operons that contain all ORF with lengths equal or greater than the NUMBER_OF_CODONS_TO_ANALYZE
#min(operons_tibble$ORF_length)

# simplify the name of ORFs in the FASTA file to match the operons tibble names
names(ORFs_file) <- gsub(x = names(ORFs_file), pattern = "CDS:", replacement = "")
# Fix names of ORFs in the FASTA file
names(ORFs_file) <- gsub(x = names(ORFs_file), pattern = "\t", replacement = "")

#----------------------------------------------------------------------#
#----------------------------------------------------------------------#
#--------------- looping over every detected operon -------------------#
#----------------------------------------------------------------------#
#----------------------------------------------------------------------#
operons_HBs_tibble <- tibble()
HBs_per_codon_tibble <- tibble()
#OPERON = 623
for (OPERON in unique(operons_tibble$Operon_ID) ) {
  print(paste("Analizing operon: ", OPERON) )
  
  ORFs_in_each_operon <- operons_tibble %>% 
    filter(Operon_ID == OPERON)
  ORFs_in_each_operon <- ORFs_in_each_operon$IdGene
  
  # get only the nucleotide sequence out of the seqinR list object
  subset_ORFs_in_each_operon <- subseq(x = ORFs_file[ORFs_in_each_operon])
  HBs_per_codon_tibble <- bind_rows(HBs_per_codon_tibble,
                                    as_tibble(
                                      c(letterFrequency(subset_ORFs_in_each_operon, "A", as.prob=FALSE) * 2 +
                                          letterFrequency(subset_ORFs_in_each_operon, "T", as.prob=FALSE) * 2 +
                                          letterFrequency(subset_ORFs_in_each_operon, "G", as.prob=FALSE) * 3 +
                                          letterFrequency(subset_ORFs_in_each_operon, "C", as.prob=FALSE) * 3) / (width(subset_ORFs_in_each_operon)/3)
                                    ) %>% 
                                      mutate(Organism = ORGANISM_NAME) %>% 
                                      mutate(HBs_per_codon = value) %>% 
                                      mutate(Operon_ID = paste0("Operon_", OPERON) ) %>% 
                                      mutate(ORF_position_within_operon = 1:n() ) %>% # name each ORF within an operon with consequtive numbers 
                                      select(Organism, Operon_ID, ORF_position_within_operon, HBs_per_codon)
                                    )
  
  
  # Subset now to the number of codons to analyze
  subset_ORFs_in_each_operon <- subseq(x = ORFs_file[ORFs_in_each_operon], start = 1, width = NUMBER_OF_CODONS_TO_ANALYZE * 3)
  
  #----------------- PROCESS -----------------#
  # Create a matrix in which we will count the occurence of nucleotides (and hence codons) and their HB cost
  xx <- enframe(as.character(subset_ORFs_in_each_operon))$value
  xx <- lapply(X=xx, FUN=seqinr::s2c)
  xx <- data.frame(matrix(unlist(xx), nrow=length(xx), byrow=T),stringsAsFactors=FALSE)
  
  xx <- ifelse(
    test = xx == "A", yes = 2,
    no = ifelse(
      test = xx == "G", yes = 3,
      no = ifelse(
        test = xx == "C", yes = 3,
        no = ifelse(
          test = xx == "T", yes = 2,
          no = NA))))
  
  xx <- sapply(X = seq(from = 1, to = ncol(xx), by=3), FUN = function(x) {rowSums(xx[, x:(x+3-1)])} )
  colnames(xx) <- paste0("Pos_", 1: (NUMBER_OF_CODONS_TO_ANALYZE) )
  
  xx <- as_tibble(xx)
  
  operons_HBs_tibble <- bind_rows(operons_HBs_tibble,
                                  xx %>%
                                    mutate(Organism = ORGANISM_NAME) %>% 
                                    mutate(Operon_ID = paste0("Operon_", OPERON) ) %>% 
                                    mutate(ORF_position_within_operon = 1:n() ) %>% # name each ORF within an operon with consequtive numbers 
                                    select(Organism, Operon_ID, ORF_position_within_operon, Pos_1:paste0("Pos_", (NUMBER_OF_CODONS_TO_ANALYZE)) ) %>% 
                                    pivot_longer(-c(Organism, Operon_ID, ORF_position_within_operon), names_to = "Position", values_to = "HBs") %>% 
                                    mutate(Position2 = gsub(x = Position, pattern = "Pos_", replacement = "")) %>% 
                                    mutate(Position2 = as.numeric(Position2)),
                                  )
}

dir.create(path = "outputs/")
write_tsv(x = operons_HBs_tibble, path = paste0("outputs/", ORGANISM_NAME, "_operons_HBs_tibble.tsv") )
write_tsv(x = HBs_per_codon_tibble, path = paste0("outputs/", ORGANISM_NAME, "_HBs_per_codon_tibble.tsv") )

number_of_ORFs_per_ORF_group <- operons_HBs_tibble %>%
  group_by(Position2, ORF_position_within_operon) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  select(ORF_position_within_operon, n) %>% 
  unique() %>% 
  mutate(Organism = ORGANISM_NAME) %>% 
  mutate(ORF_position_within_operon_str = paste0("ORF_", ORF_position_within_operon))

write_tsv(x = number_of_ORFs_per_ORF_group, path = paste0("outputs/", ORGANISM_NAME, "_number_of_ORFs_per_ORF_group.tsv") )


#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
#----------------------------- figures ---------------------------------#
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#

Fig_2A <- number_of_ORFs_per_ORF_group %>% 
  ggplot(aes(x = ORF_position_within_operon, y = n)) +
  geom_col(fill=c(rep("black", 3), rep("grey", 16))) +
  geom_text(aes(label=n), size=1, nudge_y = 10) +
  geom_hline(yintercept = max(number_of_ORFs_per_ORF_group$n)/3, lty=2) +
  annotate(y = (max(number_of_ORFs_per_ORF_group$n)/3), x = 10, geom = "label", label=paste0(max(number_of_ORFs_per_ORF_group$n)," / 3"), size=1) +
  scale_x_continuous("", labels = 1:19, breaks = 1:19) +
  theme_new() +
  NULL

N_OF_ORFs = 6
set.seed(15)
my_comparisons <- list(c("ORF_1", "ORF_2"), c("ORF_1", "ORF_3"), c("ORF_2", "ORF_3"))

Fig_2B <- bind_rows(
  operons_HBs_tibble %>%
    filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
    filter(Position2 %in% 1:20) %>% 
    mutate(Codon_range = "Codon 1 to 20"),
  operons_HBs_tibble %>%
    filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
    filter(Position2 %in% 21:40) %>% 
    mutate(Codon_range = "Codon 21 to 40"),
  operons_HBs_tibble %>%
    filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
    filter(Position2 %in% 41:60) %>% 
    mutate(Codon_range = "Codon 41 to 60"),
  operons_HBs_tibble %>%
    filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
    filter(Position2 %in% 61:80) %>% 
    mutate(Codon_range = "Codon 61 to 80"),
  operons_HBs_tibble %>%
    filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
    filter(Position2 %in% 81:100) %>% 
    mutate(Codon_range = "Codon 81 to 100")
) %>% 
  mutate(ORF_position_within_operon_str = paste0("ORF_", ORF_position_within_operon)) %>%
  filter(ORF_position_within_operon <= N_OF_ORFs) %>% 
  ggplot(aes(x = ORF_position_within_operon_str, y = HBs)) +
  stat_summary(fun.data = "mean_cl_boot", fun.args = list(conf.int = 0.95, B = 100, reps = FALSE), geom = "errorbar", width=0.1) + # ?mean_cl_boot
  stat_summary(fun.y = "mean", geom="point", colour="black", size=2) +
  stat_summary(aes(colour=ORF_position_within_operon_str), fun.y = "mean", geom="point", size=1) +
  #stat_compare_means(label.y = c(7.64, 7.66, 7.64), comparisons = my_comparisons, method = "wilcox", label = "p.format", method.args = list(alternative = "less", paired = FALSE), size=2, tip.length = 0.001) + # Wilcoxon test: Compare two groups (non-parametric)
  stat_compare_means(label.y = 7.62, ref.group = "ORF_1", method = "wilcox", label = "p.format", method.args = list(alternative = "greater", paired = FALSE), size=2) + # Wilcoxon test: Compare two groups (non-parametric)
  stat_compare_means(label.y = 7.4, label.x = 2.0, method = "kruskal.test", colour="red", size=1.5) + # Kruskal-Wallis rank sum test: Compare multiple groups (non-parametric)
  geom_text(data = number_of_ORFs_per_ORF_group  %>%
              filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>%
              filter(ORF_position_within_operon <= N_OF_ORFs),
            aes(y = 7.38, label=paste0("n = ", n, " operons") ), nudge_x = 0, size=1.5) +
  ylab(label = "Hydrogen bonds (mean ± 95% C.I.)") +
  xlab(label = "ORF position within the operon") +
  facet_wrap(~Codon_range, ncol=5) +
  theme_new() +
  theme(legend.position = "none") +
  NULL

set.seed(15)
Fig_2C <- operons_HBs_tibble %>%
  filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
  filter(Position2 %in% 1:100) %>% 
  mutate(Codon_range = "Codon 1 to 100") %>% 
  mutate(ORF_position_within_operon_str = paste0("ORF_", ORF_position_within_operon)) %>%
  filter(ORF_position_within_operon <= N_OF_ORFs) %>% 
  ggplot(aes(x = ORF_position_within_operon_str, y = HBs)) +
  stat_summary(fun.data = "mean_cl_boot", fun.args = list(conf.int = 0.95, B = 100, reps = FALSE), geom = "errorbar", width=0.1) + # ?mean_cl_boot
  stat_summary(fun.y = "mean", geom="point", colour="black", size=2) +
  stat_summary(aes(colour=ORF_position_within_operon_str), fun.y = "mean", geom="point", size=1) +
  #stat_compare_means(label.y = c(7.58, 7.60, 7.58), comparisons = my_comparisons, method = "wilcox", label = "p.format", method.args = list(alternative = "less", paired = FALSE), size=2, tip.length = 0.001) + # Wilcoxon test: Compare two groups (non-parametric)
  stat_compare_means(label.y = 7.58, ref.group = "ORF_1", method = "wilcox", label = "p.format", method.args = list(alternative = "greater", paired = FALSE), size=2) + # Wilcoxon test: Compare two groups (non-parametric)
  stat_compare_means(label.y = 7.52, label.x = 2.0, method = "kruskal.test", colour="red", size=1.5) + # Kruskal-Wallis rank sum test: Compare multiple groups (non-parametric)
  geom_text(data = number_of_ORFs_per_ORF_group  %>%
              filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>%
              filter(ORF_position_within_operon <= N_OF_ORFs),
            aes(y = 7.51, label=paste0("n = ", n, " operons") ), nudge_x = 0, size=1.5) +
  ylab(label = "Hydrogen bonds (mean ± 95% C.I.)") +
  xlab(label = "ORF position within the operon") +
  facet_wrap(~Codon_range, ncol=4) +
  theme_new() +
  theme(legend.position = "none") +
  NULL

#------------- Hydrogen bonds per codon --------------------#
set.seed(15)
Fig_2D <- HBs_per_codon_tibble %>%
  filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
  filter(ORF_position_within_operon <= N_OF_ORFs) %>% 
  mutate(ORF_position_within_operon_str = paste0("ORF_", ORF_position_within_operon)) %>% 
  ggplot(aes(x = ORF_position_within_operon_str, y = HBs_per_codon)) +
  stat_summary(fun.data = "mean_cl_boot", fun.args = list(conf.int = 0.95, B = 100, reps = FALSE), geom = "errorbar", width=0.1) + # ?mean_cl_boot
  stat_summary(fun.y = "mean", geom="point", colour="black", size=2) +
  stat_summary(aes(colour=ORF_position_within_operon_str), fun.y = "mean", geom="point", size=1) +
  #stat_compare_means(label.y = c(7.595, 7.605, 7.595), comparisons = my_comparisons, method = "wilcox", label = "p.format", method.args = list(alternative = "less", paired = FALSE), size=2, tip.length = 0.001) + # Wilcoxon test: Compare two groups (non-parametric)
  stat_compare_means(label.y = 7.595, ref.group = "ORF_1", method = "wilcox", label = "p.format", method.args = list(alternative = "greater", paired = FALSE), size=2) + # Wilcoxon test: Compare two groups (non-parametric)
  stat_compare_means(label.y = 7.535, label.x = 2.0, method = "kruskal.test", colour="skyblue", size=1.5) + # Kruskal-Wallis rank sum test: Compare multiple groups (non-parametric)
  geom_text(data = number_of_ORFs_per_ORF_group  %>%
              filter(Organism == "Escherichia_coli_str_k_12_substr_mg1655") %>% 
              filter(ORF_position_within_operon <= N_OF_ORFs),
            aes(y = 7.525, label=paste0("n = ", n, " operons") ), nudge_x = 0, size=1.5) +
  ylab(label = "Hydrogen bonds per codon (mean ± 95% C.I.)") +
  xlab(label = "ORF position within the operon") +
  ggtitle(label = "Entire ORF") +
  #facet_wrap(~Organism, scales = "free") +
  theme_new() +
  theme(legend.position = "none") +
  NULL

Fig_2A_grid <- plot_grid(NULL, Fig_2A, NULL, ncol = 3, labels = c("","A", ""), label_size = 8)
Fig_2CD <- plot_grid(plotlist = list(NULL, Fig_2C, Fig_2D,  NULL), ncol = 4, rel_widths = c(0.1, 0.4, 0.4, 0.1), labels = c("", "C", "D", ""), label_size = 8)
Fig_2 <- plot_grid(Fig_2A_grid, Fig_2B, Fig_2CD, nrow = 3, rel_heights = c(0.3, 0.3, 0.4), labels = c("", "B", ""), label_size = 8)

ggsave(plot = Fig_2,
       filename = "outputs/Fig_2.pdf",
       device = "pdf",
       width = 180,
       height = 160,
       units = "mm",
       useDingbats = FALSE)







