library(tidyverse)
library(stringr)
library(ggsci)
library(dplyr)
library(scales)

# install.packages('~/iCloud/Desktop/script_test/svBreaks/', repos = NULL, type="source")

cleanTheme <- function(base_size = 12) {
  theme(
    # plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    # axis.text = element_text(size = 20),
    axis.title = element_text(size = 30)
  )
}

slideTheme <- function(base_size = 25) {
  theme(
    plot.title = element_text(hjust = 0.5, size = 50),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 50),
    strip.text = element_text(size = 20)
  )
}

getFreqs <- function(..., infile = '/Users/Nick_curie/Desktop/parserTest/filtered_231018/summary/merged/all_samples.txt'){
  af_data <-read.delim(infile, row.names=NULL)
  colnames(af_data) <- c("sample", "event", "source", "type", "chromosome1",	"bp1",	"chromosome2",	"bp2",	"split_reads",	"disc_reads",	"genotype",	"id",	"length",	"position",	"consensus",	"microhomology",	"configuration",	"allele_frequency",	"mechanism",	"log2",	"bp1_locus",	"bp2_locus",	"affected_genes",	"status",	"notes")

  af_data$'NA' <- NULL

  excluded_samples <- c( "A373R7", "A512R17", "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")

  af_data<-suppressWarnings(separate(af_data, bp1_locus, c("bp1_gene", "bp1_feature"), sep = ", ", extra = "merge"))
  af_data<-suppressWarnings(separate(af_data, bp2_locus, c("bp2_gene", "bp2_feature"), sep = ", ", extra = "merge"))
  af_data$hit <- paste(af_data$bp1_gene, af_data$bp2_gene, sep = '-')

  # af_data$status <- ifelse(is.na(af_data$status),TRUE, FALSE)

  af_data <- af_data %>%
    dplyr::filter(!sample %in% excluded_samples,
                  is.na(status)) %>%
    dplyr::select(sample, chromosome1, allele_frequency, bp1, bp2, hit)

    # dplyr::mutate(allele_frequency = as.numeric(as.character(allele_frequency))) %>%





  # af_data$allele_freq <- suppressWarnings(as.numeric(as.character(af_data$allele_freq)))


  # af_data$sample_shor <- ifelse(grepl('R', af_data$sample),
  #                               str_match(af_data$sample, ".*(R.*)")[2],
  #                               af_data$sample)
  af_data %>%
    mutate(sample_short = as.character(as.factor(ifelse(grepl('R', sample),
                                 str_match(sample, ".*(R.*)")[2],
                                 sample)))
           )

  sample_af_data <- af_data %>%
    dplyr::mutate(sample = factor(sample)) %>%
    dplyr::mutate(allele_freq = as.numeric(as.character(allele_frequency))) %>%
    dplyr::group_by(sample, hit) %>%
    dplyr::mutate(count = seq(n())) %>%
    dplyr::mutate(gene = paste(sample, hit, count, sep = "_")) %>%
    ungroup()

  sample_af_data <- sample_af_data %>%
    group_by(sample) %>%
    mutate(cell_fraction = ifelse(chromosome1 == "X", allele_frequency,
                                ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
    # mutate(sample_short = ifelse( grepl("R", sample), str_match(sample, ".*(R.*)")[2], sample)) %>%
    arrange(-allele_frequency) %>%
    ungroup()
  # af_data <- transform(af_data, gene = reorder(gene, -allele_freq))

  return(sample_af_data)
}

plotFreqs <- function(...){
  af_data <- getFreqs()
  
  geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }
  
  c <- unnest(d)
  
  af_data$colour <- ifelse(af_data$chromosome1 == 'X' | af_data$chromosome1 == 'Y', "sex", "autosome")
  af_data$colour <- ifelse(af_data$bp1 >= 2700000 & af_data$bp2 <= 3500000 & af_data$chromosome1 == "X", "notch", af_data$colour)

  cols <- c("#199E41FE", "#2366A1FE", "#BD2A2AFE")

  p <- ggplot(af_data)
  p <- p + geom_bar(aes(fct_reorder(gene, -cell_fraction), cell_fraction, fill=colour, group=gene), alpha = 0.75, stat = "identity", position='dodge')
  p <- p + facet_wrap(~sample, scales = "free_x", ncol=9)
  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      # axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      strip.text = element_text(size = 15)
      # legend.position = "top",
      # axis.text.y = element_text(size = 20)
    )
  p <- p + geom_hline(yintercept = 0.25)
  p <- p + scale_fill_jco()

  # mutPen <- paste("mutationPenetrance.pdf")
  # cat("Writing file", mutPen, "\n")
  # ggsave(mutPen, width = 15, height = 10)
  p

  # p <- p + scale_fill_manual(values = cols)
}


tumourEvolution <- function() {
  af_data <- getFreqs()

  af_data$colour <- ifelse(af_data$chromosome1 == 'X' | af_data$chromosome1 == 'Y', "sex", "autosome")
  af_data$colour <- ifelse(af_data$bp1 >= 2700000 & af_data$bp2 <= 3500000 & af_data$chromosome1 == "X", "notch", af_data$colour)


  cols <- c("#199E41FE", "#2366A1FE", "#BD2A2AFE")

  p <- ggplot(af_data)
  p <- p + geom_point(aes(fct_reorder(gene, -cell_fraction), cell_fraction, colour = colour), size = 2, stat = "identity")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p <- p + scale_fill_identity()

  p <- p + ylim(0, 1)

  p
}

mullerPlot <- function(){
  af_data <- getFreqs()
  Muller_df <- get_Muller_df(example_edges, example_pop_df)
}


plot_allele_freqs <- function(all_samples = '/Users/Nick_curie/Desktop/parserTest/filtered_231018/summary/merged/all_samples.txt'){
  all_data <- read.delim(all_samples, header = T)
  excluded_samples <- c("A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
  # excluded_samples <- c()
  geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }
  
  gene_hits <- all_data %>% 
    dplyr::filter(!sample %in% excluded_samples,
                  is.na(status)) %>% 
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>% 
    group_by(sample, event) %>% 
    dplyr::mutate(newCol = as.character(ifelse(any(geneIn('N', affected_genes)),'Notch', 'Other'))) %>% 
    mutate(cell_fraction = ifelse(chromosome1 == "X", allele_frequency,
                                  ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
    dplyr::group_by(sample) %>%
    # dplyr::mutate(psTime = rescale(-distinct(cell_fraction))) %>% 
    dplyr::ungroup() %>% 
    select(sample, newCol, allele_frequency, cell_fraction, affected_genes) %>% 
    droplevels()
  
  
  new_df <- gene_hits %>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(highest_n = ifelse(any(newCol == 'Notch'), max(cell_fraction[newCol=='Notch']), 0)) %>% 
    dplyr::distinct(sample, highest_n) %>% 
    ungroup()
  
  gene_hits$colour <- ifelse(gene_hits$newCol=='Notch', 'black', NA)
    
  p <- ggplot(new_df)
  p <- p + geom_errorbar(data=new_df, aes(fct_reorder(sample, -highest_n), ymin=-highest_n, ymax=-highest_n, colour=sample), alpha=.8,size=1.8, show.legend = FALSE)
  p <- p + geom_jitter(data=gene_hits, aes(sample, -cell_fraction, group=sample, colour = sample), size=2, width = 0.3, alpha = 0.9, show.legend = FALSE)
  p <- p + scale_y_continuous('Variant Allele Frequency', labels=seq(1,0,by=-.25))
  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      # axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
      axis.title.x=element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.y = element_text(size = 20)
      # legend.position = "top",
      # axis.text.y = element_text(size = 20)
    )
  # p <- p + facet_wrap(~sample)
  p
  
}



