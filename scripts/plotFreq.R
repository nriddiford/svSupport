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

getFreqs <- function(infile = 'variants_out.txt'){
  af_data <- read.delim(infile, header = F)
  colnames(af_data) <- c("sample", "chrom", "bp1", "bp2", "allele_freq", "type", "length", "bp1_feature", "bp2_feature", "genes")


  af_data<-suppressWarnings(separate(af_data, bp1_feature, c("bp1_gene", "bp1_feature"), sep = ", ", extra = "merge"))
  af_data<-suppressWarnings(separate(af_data, bp2_feature, c("bp2_gene", "bp2_feature"), sep = ", ", extra = "merge"))
  af_data$hit <- paste(af_data$bp1_gene, af_data$bp2_gene, sep = '-')

  af_data <- select(af_data, "sample", "chrom", "allele_freq", "chrom", "bp1", "bp2", "hit")

  af_data$allele_freq <- suppressWarnings(as.numeric(as.character(af_data$allele_freq)))
#
#   af_data %>%
#     mutate(sample = as.character(sample)) %>%
#     group_by(sample) %>%
#     mutate(sample_short = ifelse( grepl('R', sample), str_match(sample, ".*(R.*)")[2], sample)) %>%
#     ungroup()


  af_data <- af_data %>%
    group_by(sample, hit) %>%
    mutate(count = seq(n())) %>%
    mutate(gene = ifelse(count == 1, as.character(hit), paste(hit, count, sep = "_"))) %>%
    ungroup() %>%
    mutate(sample = as.character(sample)) %>%
    mutate(gene = factor(gene)) %>%
    group_by(sample) %>%
    mutate(allele_freq = ifelse(chrom == "X", allele_freq,
                                ifelse(allele_freq*2>1, 1, allele_freq*2))) %>%
    mutate(sample_short = ifelse( grepl("R", sample), str_match(sample, ".*(R.*)")[2], sample)) %>%
    arrange(-allele_freq)
  # af_data <- transform(af_data, gene = reorder(gene, -allele_freq))

  return(af_data)
}

plotFreqs <- function (){
  af_data <- getFreqs()

  af_data$colour <- ifelse(af_data$chrom == 'X' | af_data$chrom == 'Y', "sex", "autosome")
  af_data$colour <- ifelse(af_data$bp1 >= 2700000 & af_data$bp2 <= 3500000 & af_data$chrom == "X", "notch", af_data$colour)


  cols <- c("#199E41FE", "#2366A1FE", "#BD2A2AFE")

  p <- ggplot(af_data)
  p <- p + geom_bar(aes(fct_reorder(gene, allele_freq, .desc = TRUE), allele_freq, fill = colour), alpha = 0.75, stat = "identity")
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

  p <- p + scale_fill_manual(values = cols)

  p

}


tumourEvolution <- function() {
  af_data <- getFreqs()

  af_data$colour <- ifelse(af_data$chrom == 'X' | af_data$chrom == 'Y', "#199E41FE", "#2366A1FE")
  af_data$colour <- ifelse(af_data$bp1 >= 2700000 & af_data$bp2 <= 3500000 & af_data$chrom == "X", "#BD2A2AFE", af_data$colour)

  p <- ggplot(af_data)
  p <- p + geom_point(aes(fct_reorder(gene, allele_freq, .desc = TRUE), allele_freq, colour = colour), size = 2, stat = "identity")
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p <- p + scale_fill_identity()

  p <- p + ylim(0, 1)

  p
}

