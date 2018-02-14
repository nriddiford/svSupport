library(tidyverse)
library(stringr)

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

getFreqs <- function(infile = 'flat_vars.tsv'){
  af_data <- read.delim(infile, header = F)
  colnames(af_data) <- c("sample", "chrom", "bp1", "bp2", "allele_freq", "length", "bp1_feature", "bp2_feature", "genes")
  

  af_data<-suppressWarnings(separate(af_data, bp1_feature, c("bp1_gene", "bp1_feature"), sep = ", ", extra = "merge"))
  af_data<-suppressWarnings(separate(af_data, bp2_feature, c("bp2_gene", "bp2_feature"), sep = ", ", extra = "merge"))
  af_data$hit <- paste(af_data$bp1_gene, af_data$bp2_gene, sep = '-')
  
  af_data <- select(af_data, "sample", "chrom", "allele_freq", "chrom", "bp1", "bp2", "hit")
  
  af_data <- af_data %>%
    group_by(sample, hit) %>%
    mutate(count = seq(n())) %>%
    mutate(gene = ifelse(count == 1, as.character(hit), paste(hit, count, sep = "_"))) %>%
    arrange(hit, -allele_freq)
    
  af_data <- transform(af_data, gene = reorder(gene, -allele_freq))
  
  af_data$colour <- ifelse(af_data$chrom == 'X' | af_data$chrom == 'Y', "#199E41FE", "#2366A1FE")
  af_data$colour <- ifelse(af_data$bp1 >= 2700000 & af_data$bp2 <= 3500000 & af_data$chrom == "X", "#BD2A2AFE", af_data$colour)
  
  
  p <- ggplot(af_data)
  p <- p + geom_bar(aes(gene, allele_freq, fill = colour), alpha = 0.75, stat = "identity")
  p <- p + facet_wrap(~sample, scales = "free_x", ncol=9)
  p <- p + slideTheme() +
  # p <- p + theme(
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      # axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
      # legend.position = "top",
      # axis.text.y = element_text(size = 20)
    )  
  p <- p + scale_fill_identity()
  
  p
  
  
}

tumourEvolution <- function(sample=NA) {
  bp_data <- getData()
  
  bp_data <- bp_data %>%
    filter(!is.na(allele_freq)) %>%
    filter(genotype == "somatic_tumour") %>%
    droplevels()
  
  bp_data <- bp_data %>%
    group_by(gene, sample) %>%
    filter(bp_no != "bp2") %>%
    # filter(sample == 'HUM-7') %>%
    mutate(count = seq(n())) %>%
    mutate(gene2 = ifelse(count == 1, as.character(gene), paste(gene, count, sep = "_"))) %>%
    arrange(gene, -allele_freq) %>%
    transform(gene2 = reorder(gene2, -allele_freq))
  
  
  
  bp_data$colour <- ifelse(bp_data$bp >= 2700000 & bp_data$bp <= 3500000 & bp_data$chrom == "X", "#A52A2A", "gray37")
  
  p <- ggplot(bp_data)
  p <- p + geom_bar(aes(gene2, allele_freq, fill = colour), stat = "identity")
  # p <- p + scale_y_discrete("Allele frequency", breaks=seq(0,1,by=0.1), labels = seq(0,1,by=0.1), )
  p <- p + ylim(0, 1)
  p <- p + theme(
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1)
    # axis.text.x=element_blank(),
    # axis.ticks.x=element_blank()
  )
  p <- p + scale_fill_identity()
  p <- p + facet_wrap(~sample, scales = "free_x")
  p
}