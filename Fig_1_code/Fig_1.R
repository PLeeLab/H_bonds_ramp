#-------- LIBRARIES -------#
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
#install.packages("seqinr")
#install.packages("reshape")
#install.packages("tidyverse")
#install.packages("cowplot")
#install.packages("corrr")
#install.packages("foreach")
#install.packages("doParallel")
#install.packages("doSNOW")
#install.packages("data.table")
library("Biostrings")
library("tidyverse")
library("seqinr")
library("foreach")
library("doParallel")
library("doSNOW")
library("data.table")
library("cowplot")
#library(package="benchmarkme")
#benchmarkme::get_cpu()
#-------------------------#

rm(list=ls()); # clear environment

theme_new <- function(){
  theme_bw() %+replace%
    theme(
      plot.title = element_text(   color = "black", size = 7, face = "bold", hjust = 0),
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

# Create a data frame (tibble) with the features of each codon
codon_features_df <- tibble()

# The For loop will iterate along the "words()" vector which contains all the codons and gather
# some information about each codon
for (CODON in words()) {
  codon_features_df <- rbind(codon_features_df,
                             tibble(amino_acid = aaa(translate(s2c(CODON))),
                                    AA = translate(s2c(CODON)),
                                    codon = toupper(CODON),
                                    aa_codon = paste0(aaa(translate(s2c(CODON))), "-", toupper(CODON)),
                                    ribonucl_1 = toupper(s2c(CODON)[1]),
                                    ribonucl_2 = toupper(s2c(CODON)[2]),
                                    ribonucl_3 = toupper(s2c(CODON)[3])))
}
# Let's have a look at the features dataset
glimpse(codon_features_df)
head(codon_features_df); tail(codon_features_df)

#----- Hydrogen bonds -----#
codon_features_df$Hydrogen_bonds_rb_1 <- ifelse(
  test = codon_features_df$ribonucl_1 == "A",
  yes = 2,
  no = ifelse(
    test = codon_features_df$ribonucl_1 == "G",
    yes = 3,
    no = ifelse(
      test = codon_features_df$ribonucl_1 == "C",
      yes = 3,
      no = ifelse(
        test = codon_features_df$ribonucl_1 == "T",
        yes = 2,
        no = NA
        #no = mean(c(3, 3, 2, 2))
      )
    )
  )
)

codon_features_df$Hydrogen_bonds_rb_2 <- ifelse(
  test = codon_features_df$ribonucl_2 == "A",
  yes = 2,
  no = ifelse(
    test = codon_features_df$ribonucl_2 == "G",
    yes = 3,
    no = ifelse(
      test = codon_features_df$ribonucl_2 == "C",
      yes = 3,
      no = ifelse(
        test = codon_features_df$ribonucl_2 == "T",
        yes = 2,
        no = NA
        #no = mean(c(3, 3, 2, 2))
      )
    )
  )
)

codon_features_df$Hydrogen_bonds_rb_3 <- ifelse(
  test = codon_features_df$ribonucl_3 == "A",
  yes = 2,
  no = ifelse(
    test = codon_features_df$ribonucl_3 == "G",
    yes = 3,
    no = ifelse(
      test = codon_features_df$ribonucl_3 == "C",
      yes = 3,
      no = ifelse(
        test = codon_features_df$ribonucl_3 == "T",
        yes = 2,
        no = NA
        #no = mean(c(3, 3, 2, 2))
      )
    )
  )
)

codon_features_df$Hydrogen_bonds_total <- codon_features_df$Hydrogen_bonds_rb_1 + codon_features_df$Hydrogen_bonds_rb_2 + codon_features_df$Hydrogen_bonds_rb_3

# Let's have a look at the final dataset
glimpse(codon_features_df)
head(codon_features_df)

Fig_1A <- codon_features_df %>%
  select(Hydrogen_bonds_total) %>%
  gather(key = "Element", value = "Number") %>% 
  mutate(Element = gsub("_total", "", Element)) %>% 
  ggplot(mapping = aes(x = Number, fill = Element)) +
  geom_bar(fill="grey") +
  xlab(label = "Number of hydrogen bonds") +
  ylab(label = "Number of codons") +
  theme_new() +
  theme(legend.position = "none") +
  NULL
  
Fig_1B <- codon_features_df %>% 
  mutate(Scaled_cost_Hydrogen_bonds_per_AA = ave(Hydrogen_bonds_total, amino_acid, FUN=scale)) %>% 
  select(codon, amino_acid, Scaled_cost_Hydrogen_bonds_per_AA) %>% 
  gather(key = "Element", value = "Relative_cost", -c(codon,amino_acid) ) %>% 
  mutate(Element = gsub("Scaled_cost_", "", Element)) %>% 
  mutate(Element = gsub("_per_AA", "", Element)) %>% 
  ggplot(mapping = aes(x = Element, y = codon)) +
  geom_tile(mapping = aes(fill=Relative_cost), colour="white") +
  xlab(label = "") +
  ylab(label = "Codon") +
  scale_fill_distiller(palette = "RdBu",
                       limits = c(-2,2),
                       name="Scaled number of hydrogen bonds") +
  facet_wrap(~amino_acid, scales = "free_y")+
  theme_new() +
  theme(legend.position = "right",
        axis.text.x = element_blank()) +
  NULL

Fig_1C <- codon_features_df %>% 
  mutate(Relative_cost_Hydrogen_bonds_per_AA = ave(Hydrogen_bonds_total, amino_acid, FUN=function(x){x/max(x)})) %>% 
  select(codon, amino_acid, Relative_cost_Hydrogen_bonds_per_AA) %>% 
  gather(key = "Element", value = "Relative_cost", -c(codon,amino_acid)) %>% 
  mutate(Element = gsub("Relative_cost_", "", Element)) %>% 
  mutate(Element = gsub("_per_AA", "", Element)) %>% 
  ggplot(mapping = aes(x = Element, y = Relative_cost)) +
  geom_jitter(width = 0.1, height = 0, alpha=1/2, shape=21, size=1, fill="black") +
  geom_boxplot(fill="transparent", outlier.alpha = 0, show.legend = FALSE, width=0.5) +
  xlab("") +
  ylab("Relative hydrogen bonds content") +
  theme_new() +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  NULL

Fig_1ABC_topRow  <- plot_grid(Fig_1A, Fig_1B, Fig_1C, labels = c('A', 'B', 'C'), ncol = 3, nrow = 1, rel_widths = c(1/5, 3/5, 1/5), label_size = 8)

#------ Position-dependent hydrogen bonding in E. coli ------------#

#-----------------------------------------------------------------#
#---------------------- [0] Initial setup ------------------------#
#-----------------------------------------------------------------#
PATH_TO_GENOME = "inputs/"
GENOMES = list.files(path = PATH_TO_GENOME, pattern = "_cds_from_genomic.fna.gz")
GENOME_NAMES = gsub(x = GENOMES, pattern = "_cds_from_genomic.fna.gz", replacement = "")
NUMBER_OF_CODONS_TO_ANALYZE = 250

Elements_df_boot_ALL   <- data.frame()
model_fitness_5end_ALL <- tibble()

#-----------------------------------------------------------------#
#------------------ [1] Analyzing REAL genome --------------------#
#-----------------------------------------------------------------#
CURRENT_GENOME = 1
CURRENT_GENOME_NAME = paste0("", GENOME_NAMES[CURRENT_GENOME])

subset_CDS_file <- readDNAStringSet(filepath = paste0(PATH_TO_GENOME, GENOMES[CURRENT_GENOME])) # PATH to the CDS FASTA file
  
# Remove CDSs that are not multiple of 3
subset_CDS_file <- subset_CDS_file[width(subset_CDS_file) %% 3 == 0]
#length(subset_CDS_file)
  
# Select only CDSs with length greater than NUMBER_OF_CODONS_TO_ANALYZE times 3 (e.g. 300 bases when 100 codons)
subset_CDS_file <- subset_CDS_file[width(subset_CDS_file) >= NUMBER_OF_CODONS_TO_ANALYZE * 3]
#length(subset_CDS_file)
  
# Select only the coding sequence region from the position # 1 to the position # NUMBER_OF_CODONS_TO_ANALYZE * 3
subset_CDS_file <- subseq(x = subset_CDS_file, start = 1, width = NUMBER_OF_CODONS_TO_ANALYZE * 3)
  
#----------------- PROCESS -----------------#
# Create a matrix for each element in which we will count the occurence of RNA bases (and hence codons) and their biosynthetic cost
# Note that the costs for T are counted as U as we are analyzing RNA biosynthesis cost
xx <- enframe(as.character(subset_CDS_file))$value
xx <- lapply(X=xx, FUN=seqinr::s2c)
xx <- data.frame(matrix(unlist(xx), nrow=length(xx), byrow=T),stringsAsFactors=FALSE)
#colnames(xx) <- paste0("Pos_", 1: (NUMBER_OF_CODONS_TO_ANALYZE*3) )
#head(xx)

# Transcripts - Using Uracil data instead of T
# Calculating the atoms of hydrogen invested in creating the hydrogen bonds between the base in the DNA and the newly transcribed mRNA
Hydrogen_bonds_matrix <- ifelse(
  test = xx == "A", yes = 2,
  no = ifelse(
    test = xx == "G", yes = 3,
    no = ifelse(
      test = xx == "C", yes = 3,
      no = ifelse(
        test = xx == "T", yes = 2,
        no = NA))))

# Just one check point
if (sum(is.na(c(Hydrogen_bonds_matrix))) != 0) {print("NA's found... CHECK YOUR SEQUENCE")} else {{print("Hydrogen_bonds_matrix ...calculation seems fine :D")}}

# Calculate the sum of elements in each codon given the nucleotides (n=3) in it
hydrogen_bonds_matrix_of_codons <- sapply(seq(1, ncol(Hydrogen_bonds_matrix), by=3), function(x) {rowSums(Hydrogen_bonds_matrix[, x:(x+3-1)])} )

Hydrogen_bonds_matrix <- as_tibble(hydrogen_bonds_matrix_of_codons, .name_repair = "minimal")
colnames(Hydrogen_bonds_matrix) <- c(1:NUMBER_OF_CODONS_TO_ANALYZE)
Hydrogen_bonds_matrix <- Hydrogen_bonds_matrix %>%
  gather(key = "Position", value = "Atoms")
Hydrogen_bonds_matrix$Position <- as.numeric(Hydrogen_bonds_matrix$Position)
Hydrogen_bonds_matrix$Element <- "Hydrogen bonds"

Elements_df_boot <- Hydrogen_bonds_matrix %>%
  group_by(Position) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(x = .$Atoms, conf.int = 0.95, B = 1000, na.rm = TRUE, reps = FALSE)))) %>% 
  ungroup()

Elements_df_boot$Genome <- CURRENT_GENOME_NAME

dir.create(path = "outputs", showWarnings = TRUE)

write.table(x = Elements_df_boot, file = "outputs/OUT_Elements_df_boot.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Fig_1E <- ggplot(data = Elements_df_boot %>%
                         filter(Position > 1 & Position <= 100), # 0 = no codon is excluded, 1 = only start codon is excluded,
                       mapping = aes(x = Position)) +
  geom_vline(xintercept = c(1), lty=2, colour="grey", lwd=0.25) +
  geom_pointrange(mapping = aes(y=Mean, ymin=Lower, ymax=Upper), alpha=1, size=0.1) +
  #geom_smooth(aes(y=Mean))+
  xlab(label = "Codon position (5'-end)") +
  ylab(label = "Hydrogen bonds per codon") +
  theme_new() +
  theme(#strip.background = element_blank(),
    legend.title = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.grid.major = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = "none") + #bottom
  guides(colour = guide_legend(nrow = 1))+
  NULL

Fig_1E_inset <- ggplot(data = Elements_df_boot %>%
                    filter(Position > 1 & Position <= 250), # 0 = no codon is excluded, 1 = only start codon is excluded,
                  mapping = aes(x = Position)) +
  geom_vline(xintercept = 1, lty=2, colour="grey", lwd=0.25) +
  geom_pointrange(mapping = aes(y=Mean, ymin=Lower, ymax=Upper), alpha=1, size=0.01) +
  #geom_smooth(aes(y=Mean))+
   xlab(label = "") +
   ylab(label = "") +
   theme_new() +
   theme(#strip.background = element_blank(),
     legend.title = element_blank(),
     #panel.grid.minor = element_blank(),
     #panel.grid.major = element_blank(),
     axis.text.x = element_text(colour = "black"),
     axis.text.y = element_text(colour = "black"),
     legend.position = "none") + #bottom
   guides(colour = guide_legend(nrow = 1))+
   NULL

counts_Hb6 <- colSums(hydrogen_bonds_matrix_of_codons == 6)
counts_Hb7 <- colSums(hydrogen_bonds_matrix_of_codons == 7)
counts_Hb8 <- colSums(hydrogen_bonds_matrix_of_codons == 8)
counts_Hb9 <- colSums(hydrogen_bonds_matrix_of_codons == 9)

if ( unique(counts_Hb6+counts_Hb7+counts_Hb8+counts_Hb9) != length(subset_CDS_file) ) {print("counts do no match! check the 'hydrogen_bonds_matrix_of_codons'")}

Hb_two_discrete_counts_of_codon_cost_per_position_df <- rbind(
  data.frame(Position = 1:ncol(hydrogen_bonds_matrix_of_codons),
             element = "Hydrogen bonds",
             cost = "Cheap (6Hb + 7Hb)",
             counts = counts_Hb6 + counts_Hb7),
  data.frame(Position = 1:ncol(hydrogen_bonds_matrix_of_codons),
             element = "Hydrogen bonds",
             cost = "Expensive (8Hb + 9Hb)",
             counts = counts_Hb8 + counts_Hb9)
)

Fig_1F <- ggplot(data = Hb_two_discrete_counts_of_codon_cost_per_position_df %>% 
                   filter(Position > 1 & Position < 100),
                 mapping = aes(x = Position, y = 100 * counts/length(subset_CDS_file), colour=cost)) +
  geom_vline(xintercept = 1, lty=2, colour="grey", lwd=0.25) +
  geom_line(aes(lty=cost), lwd=0.25) +
  geom_hline(yintercept = 50, lty=1, lwd=0.25) +
  #geom_smooth(method = "loess", formula = "y ~ x", show.legend = FALSE) +
  geom_point(shape=21, fill="white", size=0.75)+
  xlab(label = "Codon position (5'-end)")+
  ylab(label = "Percentage of all CDSs")+
  theme_new()+
  theme(#strip.background = element_blank(),
    legend.title = element_blank(),
    #panel.grid.minor = element_blank(),
    #panel.grid.major = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    legend.position = c(0.40, 0.85), legend.background = element_rect(fill = "transparent"))+
  guides(colour = guide_legend(nrow = 1))+
  NULL


Fig_1_DEF_bottomRow <- plot_grid(Fig_1E_inset, Fig_1E, Fig_1F, labels = c("D", "E", "F"), ncol = 3, nrow=1, label_size = 8)

ggsave(plot = plot_grid(Fig_1ABC_topRow, Fig_1_DEF_bottomRow, ncol = 1, nrow=2, rel_heights = c(4/7, 3/7)),
       filename = "outputs/Fig_1.pdf",
       device = "pdf",
       width = 180,
       height = 140,
       units = "mm",
       useDingbats = FALSE)

ggsave(plot = Fig_1E_inset,
       filename = "outputs/Fig_1E_inset.pdf",
       device = "pdf",
       width = 40,
       height = 40,
       units = "mm",
       useDingbats = FALSE)


