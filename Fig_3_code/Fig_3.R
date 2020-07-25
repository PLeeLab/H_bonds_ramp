library("tidyverse")
library("Biostrings")
library("seqinr")
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
      #strip.background = element_blank(),
      #panel.grid.minor = element_blank(),
      #panel.grid.major = element_blank(),
      complete = TRUE
    )
}

expression_table <- read_tsv(file = "inputs/GSE45443_Transcripts_Samples_7_to_22.txt.gz") %>%
  mutate(Exp_1 = `Expression 1`) %>% 
  mutate(Exp_2 = `Expression 2`) %>% 
  mutate(Exp_3 = `Expression 3`) %>% 
  mutate(Exp_4 = `Expression 4`) %>% 
  mutate(Exp_5 = `Expression 5`) %>% 
  mutate(Exp_6 = `Expression 6`) %>%
  mutate(Exp_7 = `Expression 7`) %>% 
  mutate(Exp_8 = `Expression 8`) %>% 
  mutate(Exp_9 = `Expression 9`) %>%
  mutate(Exp_10 = `Expression 10`) %>% 
  mutate(Exp_11 = `Expression 11`) %>% 
  mutate(Exp_12 = `Expression 12`) %>%
  mutate(Exp_13 = `Expression 13`) %>% 
  mutate(Exp_14 = `Expression 14`) %>% 
  mutate(Exp_15 = `Expression 15`) %>%
  mutate(Exp_16 = `Expression 16`) %>%
  mutate(Gene = Name) %>% 
  select(Gene, starts_with("Exp_")) %>% 
  filter(Gene != "-")

#----- Fisrt group of expression
PERCENT_HIGH = 5
PERCENT_BOT = 13

#----- Second group of expression
# PERCENT_HIGH = 10
# PERCENT_BOT = 18

#----- Third group of expression
# PERCENT_HIGH = 15
# PERCENT_BOT = 23

#----- Fourth group of expression
# PERCENT_HIGH = 20
# PERCENT_BOT = 26

#----- Fifth group of expression
# PERCENT_HIGH = 25
# PERCENT_BOT = 30

#----- Sixth group of expression
# PERCENT_HIGH = 30
# PERCENT_BOT = 35


top_percent_genes_Exp_1 <- expression_table %>% 
  arrange(desc(Exp_1)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_2 <- expression_table %>% 
  arrange(desc(Exp_2)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_3 <- expression_table %>% 
  arrange(desc(Exp_3)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_4 <- expression_table %>% 
  arrange(desc(Exp_4)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_5 <- expression_table %>% 
  arrange(desc(Exp_5)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_6 <- expression_table %>% 
  arrange(desc(Exp_6)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_7 <- expression_table %>% 
  arrange(desc(Exp_7)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_8 <- expression_table %>% 
  arrange(desc(Exp_8)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_9 <- expression_table %>% 
  arrange(desc(Exp_9)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_10 <- expression_table %>% 
  arrange(desc(Exp_10)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_11 <- expression_table %>% 
  arrange(desc(Exp_11)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_12 <- expression_table %>% 
  arrange(desc(Exp_12)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_13 <- expression_table %>% 
  arrange(desc(Exp_13)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_14 <- expression_table %>% 
  arrange(desc(Exp_14)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_15 <- expression_table %>% 
  arrange(desc(Exp_15)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes_Exp_16 <- expression_table %>% 
  arrange(desc(Exp_16)) %>% 
  head(round(dim(expression_table)[1]*(PERCENT_HIGH/100))) %>%
  select(Gene) %>%
  as_vector()

top_percent_genes <- Reduce(intersect, list(top_percent_genes_Exp_1,
                                            top_percent_genes_Exp_2,
                                            top_percent_genes_Exp_3,
                                            top_percent_genes_Exp_4,
                                            top_percent_genes_Exp_5,
                                            top_percent_genes_Exp_6,
                                            top_percent_genes_Exp_7,
                                            top_percent_genes_Exp_8,
                                            top_percent_genes_Exp_9,
                                            top_percent_genes_Exp_10,
                                            top_percent_genes_Exp_11,
                                            top_percent_genes_Exp_12,
                                            top_percent_genes_Exp_13,
                                            top_percent_genes_Exp_14,
                                            top_percent_genes_Exp_15,
                                            top_percent_genes_Exp_16))



bot_percent_genes_Exp_1 <- expression_table %>% 
  arrange(Exp_1) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_2 <- expression_table %>% 
  arrange(Exp_2) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_3 <- expression_table %>% 
  arrange(Exp_3) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_4 <- expression_table %>% 
  arrange(Exp_4) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_5 <- expression_table %>% 
  arrange(Exp_5) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_6 <- expression_table %>% 
  arrange(Exp_6) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_7 <- expression_table %>% 
  arrange(Exp_7) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_8 <- expression_table %>% 
  arrange(Exp_8) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_9 <- expression_table %>% 
  arrange(Exp_9) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_10 <- expression_table %>% 
  arrange(Exp_10) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_11 <- expression_table %>% 
  arrange(Exp_11) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_12 <- expression_table %>% 
  arrange(Exp_12) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_13 <- expression_table %>% 
  arrange(Exp_13) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_14 <- expression_table %>% 
  arrange(Exp_14) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_15 <- expression_table %>% 
  arrange(Exp_15) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes_Exp_16 <- expression_table %>% 
  arrange(Exp_16) %>% 
  head(round(dim(expression_table)[1]*PERCENT_BOT/100)) %>% 
  select(Gene) %>%
  as_vector()

bot_percent_genes <- Reduce(intersect, list(bot_percent_genes_Exp_1,
                                            bot_percent_genes_Exp_2,
                                            bot_percent_genes_Exp_3,
                                            bot_percent_genes_Exp_4,
                                            bot_percent_genes_Exp_5,
                                            bot_percent_genes_Exp_6,
                                            bot_percent_genes_Exp_7,
                                            bot_percent_genes_Exp_8,
                                            bot_percent_genes_Exp_9,
                                            bot_percent_genes_Exp_10,
                                            bot_percent_genes_Exp_11,
                                            bot_percent_genes_Exp_12,
                                            bot_percent_genes_Exp_13,
                                            bot_percent_genes_Exp_14,
                                            bot_percent_genes_Exp_15,
                                            bot_percent_genes_Exp_16))

random_percent_genes <- sample(x = expression_table$Gene,
                               size = length(top_percent_genes),
                               replace = FALSE)

ORFeome <- readDNAStringSet(filepath = "inputs/Ecoli_ORFeome.fasta")
names(ORFeome) <- gsub(x = gsub(x = names(ORFeome), pattern = ".*.\\[gene=", replacement = ""), pattern = "] .*", replacement = "")

top_ORFs    <- ORFeome[names(ORFeome)%in%top_percent_genes]
#writeXStringSet(x = top_ORFs, filepath = "outputs/TOP.fasta", format = "fasta")
bot_ORFs <- ORFeome[names(ORFeome)%in%bot_percent_genes]
#writeXStringSet(x = bottom_ORFs, filepath = "outputs/BOTTOM.fasta", format = "fasta")
random_ORFs <- ORFeome[names(ORFeome)%in%random_percent_genes]
#writeXStringSet(x = bottom_ORFs, filepath = "outputs/RANDOM.fasta", format = "fasta")

dir.create(path = "outputs")

number_of_ORFs_per_cluster <- tibble(Percent = PERCENT_HIGH,
       A_Bottom = length(bot_ORFs),
       C_Top = length(top_ORFs),
       D_Random = length(random_ORFs))

number_of_ORFs_per_cluster

write_tsv(x = number_of_ORFs_per_cluster,
          path = paste0("outputs/number_of_ORFs_per_cluster_", PERCENT_HIGH,"_percent.tsv"))

NUMBER_OF_CODONS_TO_ANALYZE = 100

#--------------------- TOP ---------------------#
subset_CDS_file <- top_ORFs[width(top_ORFs) >= NUMBER_OF_CODONS_TO_ANALYZE * 3]
subset_CDS_file <- subseq(x = subset_CDS_file, start = 1, width = NUMBER_OF_CODONS_TO_ANALYZE * 3)
xx <- enframe(as.character(subset_CDS_file))$value
xx <- lapply(X=xx, FUN=seqinr::s2c)
xx <- data.frame(matrix(unlist(xx), nrow=length(xx), byrow=T),stringsAsFactors=FALSE)
Hydrogen_bonds_matrix <- ifelse(
  test = xx == "A", yes = 2,
  no = ifelse(
    test = xx == "G", yes = 3,
    no = ifelse(
      test = xx == "C", yes = 3,
      no = ifelse(
        test = xx == "T", yes = 2,
        no = NA))))
Hydrogen_bonds_matrix <- sapply(seq(1, ncol(Hydrogen_bonds_matrix), by=3), function(x) {rowSums(Hydrogen_bonds_matrix[, x:(x+3-1)])} )
Hydrogen_bonds_matrix <- as_tibble(Hydrogen_bonds_matrix, .name_repair = "minimal")
colnames(Hydrogen_bonds_matrix) <- c(1:NUMBER_OF_CODONS_TO_ANALYZE)
Hydrogen_bonds_matrix <- Hydrogen_bonds_matrix %>%
  gather(key = "Position", value = "Atoms")
Hydrogen_bonds_matrix$Position <- as.numeric(Hydrogen_bonds_matrix$Position)
Hydrogen_bonds_matrix$Element <- "Hydrogen bonds"

Elements_df_boot_top <- Hydrogen_bonds_matrix %>%
  group_by(Position) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(x = .$Atoms, conf.int = 0.95, B = 100, na.rm = TRUE, reps = FALSE)))) %>% 
  ungroup() %>% 
  mutate(Cluster = "C_Top")

#--------------------- BOTTOM ---------------------#
subset_CDS_file <- bot_ORFs[width(bot_ORFs) >= NUMBER_OF_CODONS_TO_ANALYZE * 3]
subset_CDS_file <- subseq(x = subset_CDS_file, start = 1, width = NUMBER_OF_CODONS_TO_ANALYZE * 3)
xx <- enframe(as.character(subset_CDS_file))$value
xx <- lapply(X=xx, FUN=seqinr::s2c)
xx <- data.frame(matrix(unlist(xx), nrow=length(xx), byrow=T),stringsAsFactors=FALSE)
Hydrogen_bonds_matrix <- ifelse(
  test = xx == "A", yes = 2,
  no = ifelse(
    test = xx == "G", yes = 3,
    no = ifelse(
      test = xx == "C", yes = 3,
      no = ifelse(
        test = xx == "T", yes = 2,
        no = NA))))
Hydrogen_bonds_matrix <- sapply(seq(1, ncol(Hydrogen_bonds_matrix), by=3), function(x) {rowSums(Hydrogen_bonds_matrix[, x:(x+3-1)])} )
Hydrogen_bonds_matrix <- as_tibble(Hydrogen_bonds_matrix, .name_repair = "minimal")
colnames(Hydrogen_bonds_matrix) <- c(1:NUMBER_OF_CODONS_TO_ANALYZE)
Hydrogen_bonds_matrix <- Hydrogen_bonds_matrix %>%
  gather(key = "Position", value = "Atoms")
Hydrogen_bonds_matrix$Position <- as.numeric(Hydrogen_bonds_matrix$Position)
Hydrogen_bonds_matrix$Element <- "Hydrogen bonds"

Elements_df_boot_bottom <- Hydrogen_bonds_matrix %>%
  group_by(Position) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(x = .$Atoms, conf.int = 0.95, B = 100, na.rm = TRUE, reps = FALSE)))) %>% 
  ungroup() %>% 
  mutate(Cluster = "A_Bottom")

#--------------------- RANDOM ---------------------#
subset_CDS_file <- random_ORFs[width(random_ORFs) >= NUMBER_OF_CODONS_TO_ANALYZE * 3]
subset_CDS_file <- subseq(x = subset_CDS_file, start = 1, width = NUMBER_OF_CODONS_TO_ANALYZE * 3)
xx <- enframe(as.character(subset_CDS_file))$value
xx <- lapply(X=xx, FUN=seqinr::s2c)
xx <- data.frame(matrix(unlist(xx), nrow=length(xx), byrow=T),stringsAsFactors=FALSE)
Hydrogen_bonds_matrix <- ifelse(
  test = xx == "A", yes = 2,
  no = ifelse(
    test = xx == "G", yes = 3,
    no = ifelse(
      test = xx == "C", yes = 3,
      no = ifelse(
        test = xx == "T", yes = 2,
        no = NA))))
Hydrogen_bonds_matrix <- sapply(seq(1, ncol(Hydrogen_bonds_matrix), by=3), function(x) {rowSums(Hydrogen_bonds_matrix[, x:(x+3-1)])} )
Hydrogen_bonds_matrix <- as_tibble(Hydrogen_bonds_matrix, .name_repair = "minimal")
colnames(Hydrogen_bonds_matrix) <- c(1:NUMBER_OF_CODONS_TO_ANALYZE)
Hydrogen_bonds_matrix <- Hydrogen_bonds_matrix %>%
  gather(key = "Position", value = "Atoms")
Hydrogen_bonds_matrix$Position <- as.numeric(Hydrogen_bonds_matrix$Position)
Hydrogen_bonds_matrix$Element <- "Hydrogen bonds"

Elements_df_boot_random <- Hydrogen_bonds_matrix %>%
  group_by(Position) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(x = .$Atoms, conf.int = 0.95, B = 100, na.rm = TRUE, reps = FALSE)))) %>% 
  ungroup() %>% 
  mutate(Cluster = "D_Random")

Elements_df_boot <- bind_rows(Elements_df_boot_top, Elements_df_boot_bottom, Elements_df_boot_random)
write_tsv(x = Elements_df_boot, path = paste0("outputs/HBs_per_expression_level_", PERCENT_HIGH,"_percent.tsv"))

#------------------------------ FIGURES ------------------------------#
expression_table <- read_tsv(file = "inputs/GSE45443_Transcripts_Samples_7_to_22.txt.gz") %>%
  mutate(Exp_1 = `Expression 1`) %>% 
  mutate(Exp_2 = `Expression 2`) %>% 
  mutate(Exp_3 = `Expression 3`) %>% 
  mutate(Exp_4 = `Expression 4`) %>% 
  mutate(Exp_5 = `Expression 5`) %>% 
  mutate(Exp_6 = `Expression 6`) %>%
  mutate(Exp_7 = `Expression 7`) %>% 
  mutate(Exp_8 = `Expression 8`) %>% 
  mutate(Exp_9 = `Expression 9`) %>%
  mutate(Exp_10 = `Expression 10`) %>% 
  mutate(Exp_11 = `Expression 11`) %>% 
  mutate(Exp_12 = `Expression 12`) %>%
  mutate(Exp_13 = `Expression 13`) %>% 
  mutate(Exp_14 = `Expression 14`) %>% 
  mutate(Exp_15 = `Expression 15`) %>%
  mutate(Exp_16 = `Expression 16`) %>%
  mutate(Gene = Name) %>% 
  select(Gene, starts_with("Exp_")) %>% 
  filter(Gene != "-")

Fig_3A <- plot_grid(
  expression_table %>% 
  pivot_longer(cols = -Gene, names_to = "Experiment", values_to = "Expression") %>% 
  filter(Experiment == "Exp_1") %>% 
  ggplot(aes(x = Experiment, y = reorder(Gene, Expression), fill=log2(Expression + 0.5))) +
  geom_tile() +
  ylab(label = "Genes (Ordered by expression)") +
  xlab(label = "Experiment") +
  theme_new() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none"),
  
  expression_table %>% 
    pivot_longer(cols = -Gene, names_to = "Experiment", values_to = "Expression") %>% 
    filter(Experiment == "Exp_2") %>% 
    ggplot(aes(x = Experiment, y = reorder(Gene, Expression), fill=log2(Expression + 0.5))) +
    geom_tile() +
    ylab(label = "") +
    xlab(label = "") +
    scale_fill_gradient(low = "#e0ecf4", high = "#4d004b", na.value = NA) +
    theme_new() +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"),
  
  expression_table %>% 
    pivot_longer(cols = -Gene, names_to = "Experiment", values_to = "Expression") %>% 
    filter(Experiment == "Exp_3") %>% 
    ggplot(aes(x = Experiment, y = reorder(Gene, Expression), fill=log2(Expression + 0.5))) +
    geom_tile() +
    ylab(label = "") +
    xlab(label = "") +
    scale_fill_gradient(low = "#fee0d2", high = "#67000d", na.value = NA) +
    theme_new() +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"),
  
  NULL,
  
  expression_table %>% 
    pivot_longer(cols = -Gene, names_to = "Experiment", values_to = "Expression") %>% 
    filter(Experiment == "Exp_15") %>% 
    ggplot(aes(x = Experiment, y = reorder(Gene, Expression), fill=log2(Expression + 0.5))) +
    geom_tile() +
    ylab(label = "") +
    xlab(label = "") +
    scale_fill_gradient(low = "#e5f5f9", high = "#00441b", na.value = NA) +
    theme_new() +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"),

  expression_table %>% 
    pivot_longer(cols = -Gene, names_to = "Experiment", values_to = "Expression") %>% 
    filter(Experiment == "Exp_16") %>% 
    ggplot(aes(x = Experiment, y = reorder(Gene, Expression), fill=log2(Expression + 0.5))) +
    geom_tile() +
    ylab(label = "") +
    xlab(label = "") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
    theme_new() +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none"),
  ncol = 6
  )


number_of_ORFs_per_cluster <- rbind(
  read_tsv(file = "outputs/number_of_ORFs_per_cluster_5_percent.tsv")  %>% mutate(Percent = "A_5%"),
  read_tsv(file = "outputs/number_of_ORFs_per_cluster_10_percent.tsv") %>% mutate(Percent = "B_10%"),
  read_tsv(file = "outputs/number_of_ORFs_per_cluster_15_percent.tsv") %>% mutate(Percent = "C_15%"),
  read_tsv(file = "outputs/number_of_ORFs_per_cluster_20_percent.tsv") %>% mutate(Percent = "D_20%"),
  read_tsv(file = "outputs/number_of_ORFs_per_cluster_25_percent.tsv") %>% mutate(Percent = "E_25%"),
  read_tsv(file = "outputs/number_of_ORFs_per_cluster_30_percent.tsv") %>% mutate(Percent = "F_30%")
) %>% pivot_longer(-Percent, names_to = "Cluster", values_to = "Number_of_ORFs")

all_HBs_per_gene_expression <- rbind(
  read_tsv(file = "outputs/HBs_per_expression_level_5_percent.tsv")  %>% mutate(Percent = "A_5%"),
  read_tsv(file = "outputs/HBs_per_expression_level_10_percent.tsv") %>% mutate(Percent = "B_10%"),
  read_tsv(file = "outputs/HBs_per_expression_level_15_percent.tsv") %>% mutate(Percent = "C_15%"),
  read_tsv(file = "outputs/HBs_per_expression_level_20_percent.tsv") %>% mutate(Percent = "D_20%"),
  read_tsv(file = "outputs/HBs_per_expression_level_25_percent.tsv") %>% mutate(Percent = "E_25%"),
  read_tsv(file = "outputs/HBs_per_expression_level_30_percent.tsv") %>% mutate(Percent = "F_30%")
  )

Fig_3B <- ggplot(data = all_HBs_per_gene_expression %>%
         #filter(Cluster != "D_Random") %>%
         filter(Cluster != "B_Mid") %>%
         #filter(Percent == "A_5%") %>%
         #filter(Cluster != "B_Mid" | Percent != "C_15%") %>%
         filter(Position >= 2 & Position <= 100),
       mapping = aes(x = Position, colour=Cluster)) +
  #geom_vline(xintercept = c(1), lty=2, colour="black") +
  #geom_pointrange(mapping = aes(y=Mean, ymin=Lower, ymax=Upper), alpha=1/1) +
  #geom_line(aes(y=Mean), alpha=1/1) +
  geom_point(aes(y=Mean), alpha=1/3) +
  geom_smooth(aes(y=Mean), se = TRUE, size=0.5) +
  geom_text(data = number_of_ORFs_per_cluster %>% filter(Cluster != "B_Mid"), aes(label=Number_of_ORFs), x = 50, y = 7.4, size=3) +
  xlab(label = "Codon position (5'-end)") +
  ylab(label = "HBs") +
  scale_color_brewer(palette = "Set1") +
  facet_grid(Percent~Cluster, scales="free") +
  theme_new() +
  theme(legend.position = "none") +
  NULL

Fig_3C <- all_HBs_per_gene_expression %>% 
  filter(Percent == "A_5%") %>% 
  filter(Cluster != c("B_Mid")) %>% 
  filter(Cluster != c("D_Random")) %>% 
  select(Position, Mean, Cluster) %>% 
  pivot_wider(id_cols = Position, names_from = Cluster, values_from = Mean) %>% 
  mutate(Delta = A_Bottom - C_Top) %>% 
  mutate(Delta_str = ifelse(test = Delta > 0, yes = "Positive", no = "Negative")) %>% 
  ggplot(aes(x = Position, y = Delta, colour=Delta_str)) +
  geom_hline(yintercept = 0, lty=1) +
  geom_point(size=0.5) +
  geom_linerange(aes(ymin = 0, ymax = Delta), size=0.25) +
  ylab(label = "Delta HBs (Bottom - Top)") +
  ylim(c(-0.6, 0.6)) +
  xlab(label = "Codon position (5'-end)") +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme_new() +
  theme(legend.position = "none") +
  NULL

Fig_3AC <- plot_grid(NULL, Fig_3A, NULL, Fig_3C, ncol = 1, labels = c("A", "", "","C"), label_size = 8, rel_heights = c(0.15, 0.3, 0.15, 0.4))
Fig_3 <- plot_grid(Fig_3AC, NULL, Fig_3B, ncol = 3, labels = c("", "", "B"), label_size = 8, rel_widths = c(0.4, 0.1, 0.5))

ggsave(plot = Fig_3,
       filename = "outputs/Fig_3.pdf",
       device = "pdf",
       width = 180,
       height = 120,
       units = "mm",
       useDingbats = FALSE)
