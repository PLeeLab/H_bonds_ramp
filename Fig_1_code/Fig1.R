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

GC_matrix <- ifelse(
  test = xx == "A", yes = 0,
  no = ifelse(
    test = xx == "G", yes = 1,
    no = ifelse(
      test = xx == "C", yes = 1,
      no = ifelse(
        test = xx == "T", yes = 0,
        no = NA))))

# Just one check point
if (sum(is.na(c(Hydrogen_bonds_matrix))) != 0) {print("NA's found... CHECK YOUR SEQUENCE")} else {{print("Hydrogen_bonds_matrix ...calculation seems fine :D")}}
if (sum(is.na(c(GC_matrix))) != 0) {print("NA's found... CHECK YOUR SEQUENCE")} else             {{print("GC_matrix ...............calculation seems fine :D")}}

# Calculate the sum of elements in each codon given the nucleotides (n=3) in it
hydrogen_bonds_matrix_of_codons <- sapply(seq(1, ncol(Hydrogen_bonds_matrix), by=3), function(x) {rowSums(Hydrogen_bonds_matrix[, x:(x+3-1)])} )
GC_matrix <-             sapply(seq(1, ncol(GC_matrix), by=3),             function(x) {rowSums(GC_matrix[, x:(x+3-1)])} )

#--------------- Counts for Expected matrix -----------------#
#--------------- Fig 2 -----------------#
Hydrogen_bonds_matrix <- as_tibble(hydrogen_bonds_matrix_of_codons, .name_repair = "minimal")
colnames(Hydrogen_bonds_matrix) <- c(1:NUMBER_OF_CODONS_TO_ANALYZE)
Hydrogen_bonds_matrix <- Hydrogen_bonds_matrix %>%
  gather(key = "Position", value = "Atoms")
Hydrogen_bonds_matrix$Position <- as.numeric(Hydrogen_bonds_matrix$Position)
Hydrogen_bonds_matrix$Element <- "Hydrogen bonds"

GC_matrix <- as_tibble(GC_matrix, .name_repair = "minimal")
colnames(GC_matrix) <- c(1:NUMBER_OF_CODONS_TO_ANALYZE)
GC_matrix <- GC_matrix %>%
  gather(key = "Position", value = "Atoms")
GC_matrix$Position <- as.numeric(GC_matrix$Position)
GC_matrix$Element <- "GC content"

Elements_df_boot <- rbind(Hydrogen_bonds_matrix, GC_matrix) %>%
  group_by(Position, Element) %>%
  do(data.frame(rbind(Hmisc::smean.cl.boot(x = .$Atoms, conf.int = 0.95, B = 1000, na.rm = TRUE, reps = FALSE)))) %>% 
  ungroup()

Elements_df_boot$Genome <- CURRENT_GENOME_NAME

dir.create(path = "outputs", showWarnings = TRUE)

write.table(x = Elements_df_boot, file = "outputs/OUT_Elements_df_boot.txt", sep = "\t", row.names = FALSE, quote = FALSE)

Fig_1B <- ggplot(data = Elements_df_boot %>%
                         filter(Element == "Hydrogen bonds") %>% # filter(Element == "GC content")
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

Fig_1B_inset <- ggplot(data = Elements_df_boot %>%
                    filter(Element == "Hydrogen bonds") %>% # filter(Element == "GC content")
                    filter(Position > 1 & Position <= 250), # 0 = no codon is excluded, 1 = only start codon is excluded,
                  mapping = aes(x = Position)) +
  geom_vline(xintercept = 1, lty=2, colour="grey", lwd=0.25) +
  geom_pointrange(mapping = aes(y=Mean, ymin=Lower, ymax=Upper), alpha=1, size=0.005) +
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

Fig_1C <- ggplot(data = Hb_two_discrete_counts_of_codon_cost_per_position_df %>% 
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

# Expensive codons (rich in HB) are depleted (selected against) at the 5'-end of coding sequences in order to reduce the cost of helicase ATP-dependent hydrolisis of the DNA-RNA hydrogen bonds in transcription elongation
# Selection against cost-reduced ribonucleotides/codons at the 5'-end of coding sequences

ggsave(plot = plot_grid(Fig_1B, Fig_1C, labels = c("B", "C"), ncol = 2, nrow=1, label_size = 8),
       filename = "outputs/Fig_1.pdf",
       device = "pdf",
       width = 120,
       height = 60,
       units = "mm",
       useDingbats = FALSE)

ggsave(plot = Fig_1B_inset,
       filename = "outputs/Fig_1_inset.pdf",
       device = "pdf",
       width = 40,
       height = 40,
       units = "mm",
       useDingbats = FALSE)


