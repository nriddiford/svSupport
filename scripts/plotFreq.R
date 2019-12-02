library(tidyverse)
library(stringr)
library(ggsci)
library(scales)
library(plotly)
library(forcats)
# install.packages('~/iCloud/Desktop/script_test/svBreaks/', repos = NULL, type="source")

excluded_samples <- c("A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7",
                      "A785-A788R9", "B241R41-2", "D050R01", "D050R03", "D050R05", "D050R07-1",
                      "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20",
                      "D050R22", "D050R24", "A373R7", "A512R17", "A373R7", "A512R17", "D106R11", "D106R1",
                      "D106R13", "D106R15", "D106R17", "D106R19", "D106R21", "D106R23", "D106R25", "D106R27",
                      "D106R29", "D106R3", "D106R31", "D106R33", "D106R5", "D106R7", "D106R9")

excluded_samples <- c("B241R41-2", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20",
                      "D050R22", "D050R24", "A373R7", "A512R17", "A373R7", "A512R17")

excluded_samples <- c("B241R41-2",  "A373R7", "A512R17", "D050R01", "D050R03", "D050R05", "D050R07-1", "D050R07-2", "D050R10", "D050R12", "D050R14", "D050R16", "D050R18", "D050R20", "D050R22", "D050R24")
excluded_samples <- c(excluded_samples, "A785-A788R1", "A785-A788R11", "A785-A788R3", "A785-A788R5", "A785-A788R7", "A785-A788R9")
excluding_females <- c(excluded_samples, "D106R11", "D106R1", "D106R13", "D106R15", "D106R17", "D106R19", "D106R21", "D106R3", "D106R5", "D106R7", "D106R9")




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


# readBed
#'
#' Simple function to read in bedfiles and return as a 3 col df
#' @keywords bed
#' @import dplyr
#' @export
readBed <- function(bed_in = NULL, annotated_file=FALSE){
  if(annotated_file){
    f <- read.delim(bed_in, header=T)
    df <- f %>%
      dplyr::mutate(chrom = chromosome,
                    start = pos,
                    end = start + 1) %>%
      dplyr::select(chrom, start, end)
  } else{
    df <- read.delim(bed_in, header=F)
    if(is.null(df$V3)){
      df$V3 <- df$V2 + 2
    }
    if(!(is.null(df$V4))){
      colnames(df[,c(1,2,3,4,5)]) <- c("chrom", "start", "end", "gene", "id")
      names(df)[1:5] <- c("chrom", "start", "end","gene", "id")
    } else {
      colnames(df[,c(1,2,3)]) <- c("chrom", "start", "end")
      names(df)[1:3] <- c("chrom", "start", "end")
    }
   
  }
  return(df)
}



plot_allele_freqs <- function(..., all_samples = '~/Desktop/script_test/alleleFreqs/all_samples_merged_14619.txt',
                              all_samples_snvs = '~/Desktop/script_test/mutationProfiles/data/annotated_snvs.txt',
                              all_samples_indels = '~/Desktop/script_test/mutationProfiles/data/annotated_indels.txt',
                              tes=F,
                              annotate_with = '~/Desktop/gene2bed/bed/dna_damage_repair.bed',
                              by='fraction', snvs=T, indels=T, drivers='N', show_subclonal=T, exclude=T){
  all_data <- read.delim(all_samples, header = T)
  
  geneIn <- function(gene, gene_list) {
    sapply(as.character(gene_list), function(x) gene %in% strsplit(x, ", ")[[1]], USE.NAMES=FALSE)
  }
  
  bed_file <- readBed(bed_in = annotate_with)
  colnames(bed_file) <- c("chrom", "start", "end", "gene", "id")
  
  # all_data <- all_data %>% dplyr::filter(type != 'DUP')
  
  gene_hits <- all_data %>% 
    dplyr::filter(...,
                  !status %in% c('F', 'aF')) %>% 
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>% 
    dplyr::mutate(type_decomposed = as.character(ifelse(str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>%
    group_by(sample, event) %>% 
    dplyr::mutate(gene_hit = as.character(ifelse(any(geneIn(drivers, affected_genes)),'Notch', 'Other'))) %>% 
    dplyr::mutate(special_hit = as.character(ifelse(any(geneIn(bed_file$gene, affected_genes)),'Hit', 'Other'))) %>% 
    
    dplyr::mutate(cell_fraction = ifelse(chromosome1 %in% c("X", "Y"), allele_frequency,
                                  ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
    dplyr::group_by(sample) %>%
    # dplyr::mutate(psTime = rescale(-distinct(cell_fraction))) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(sample, type_decomposed, gene_hit, allele_frequency, cell_fraction, affected_genes, special_hit) %>% 
    droplevels()
  
  
  new_df <- gene_hits %>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(highest_n = ifelse(any(gene_hit == 'Notch'), max(cell_fraction[gene_hit=='Notch']), 0)) %>% 
    dplyr::distinct(sample, highest_n) %>% 
    ungroup()
  
  if(snvs){
    snv_data <- read.delim(all_samples_snvs, header = T)
    colnames(snv_data) <- c('sample', 'chromosome',	'pos',	'ref',	'alt',	'trinuc',	'trans',	'decomp_trinuc',	'grouped_trans',	'allele_frequency',	'caller',	'variant_type',	'status',	'snpEff_anno',	'feature',	'gene',	'id')
    
    snv_hits <- snv_data %>% 
      dplyr::filter(...) %>% 
      dplyr::group_by(sample) %>% 
      dplyr::mutate(cell_fraction = ifelse(chromosome %in% c("X", "Y"), allele_frequency,
                                    ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
      dplyr::ungroup() %>% 
      dplyr::select(sample, allele_frequency, gene, cell_fraction) %>% 
      droplevels()
  }
  if(indels){
    indel_data <- read.delim(all_samples_indels, header = T)
    # colnames(snv_data) <- c('sample', 'chromosome',	'pos',	'ref',	'alt',	'trinuc',	'trans',	'decomp_trinuc',	'grouped_trans',	'allele_frequency',	'caller',	'variant_type',	'status',	'snpEff_anno',	'feature',	'gene',	'id')
    
    indel_hits <- indel_data %>% 
      dplyr::filter(...) %>% 
      dplyr::group_by(sample) %>% 
      dplyr::mutate(cell_fraction = ifelse(chromosome %in% c("X", "Y"), allele_frequency,
                                    ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
      dplyr::ungroup() %>% 
      dplyr::select(sample, allele_frequency, gene, cell_fraction) %>% 
      droplevels()
  }
  
  if(tes){
    merged_muts <- m %>%
      dplyr::mutate(shortName = ifelse(grepl('R', sample),
                                       sub(".*[R]", "R", sample),
                                       as.character(sample))
      )
    
    ## Get TE insertions per sample
    setwd("~/Desktop/mutCor/")
    te_data <- read.delim('te_insertions.txt', header = F)
    colnames(te_data) <- c('shortName', 'TE')
    te_data <- te_data %>% mutate(shortName = ifelse(shortName == "R41a", "R41-1", as.character(shortName)))
    te_data <- te_data %>% mutate(shortName = ifelse(shortName == "R41b", "R41-2", as.character(shortName)))
  }
  

  p <- ggplot(new_df)
  p <- p + geom_errorbar(data=new_df, aes(fct_reorder(sample, -highest_n), ymin=-highest_n, ymax=-highest_n, colour=sample), alpha=.8,size=1.8, show.legend = FALSE)
  
    if(by=='fraction'){
      if(show_subclonal){
        p <- p + geom_jitter(data=gene_hits, aes(sample, -cell_fraction, group=sample, colour = sample), size=2, width = 0.3, alpha = 0.9, show.legend = FALSE)
        if(snvs){
          p <- p + geom_jitter(data=snv_hits, aes(sample, -cell_fraction, group=sample, colour = sample), shape=4, size=.5, width = 0.3, alpha = 0.2, show.legend = FALSE)
        }
        if(indels){
          p <- p + geom_jitter(data=indel_hits, aes(sample, -cell_fraction, group=sample, colour = sample), shape=8, size=.5, width = 0.3, alpha = 0.2, show.legend = FALSE)
        }
        if(tes){
          p <- p + geom_jitter(data=indel_hits, aes(sample, -cell_fraction, group=sample, colour = sample), shape=8, size=.5, width = 0.3, alpha = 0.2, show.legend = FALSE)
        }
      }
      p <- p + scale_y_continuous('Fraction of cells with mutation', labels=seq(1,0,by=-.25))
    } else{
      if(show_subclonal){
        p <- p + geom_jitter(data=gene_hits, aes(sample, -allele_frequency, group=sample, colour = sample), size=2, width = 0.3, alpha = 0.9, show.legend = FALSE)
        if(snvs){
          p <- p + geom_jitter(data=snv_hits, aes(sample, -allele_frequency, group=sample, colour = sample), shape=4, size=.5, width = 0.3, alpha = 0.2, show.legend = FALSE)
        }
      }
      p <- p + scale_y_continuous('Variant Allele Frequency', labels=seq(1,0,by=-.25))
    }
  
  p <- p + slideTheme() +
    theme(
      panel.grid.major.y = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      # axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
      axis.title.x=element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.y = element_text(size = 20),
      legend.position = "top"
      # axis.text.y = element_text(size = 20)
    )
  p <- p + guides(sample = FALSE, group = FALSE)
  
  # ggplotly(p)
  print(p)
  return(new_df)
}


geneIn <- function(gene, gene_list) {
  # hits <- sapply(as.character(gene_list), function(x) tolower(gene) %in% tolower(strsplit(x, ", ")[[1]]), USE.NAMES=FALSE)
  affected_genes <- gene_list
  genes <- gene
  # if(!strsplit(affected_genes, ", ")[[1]]) return(FALSE)
  i <- which(tolower(genes) %in% tolower(strsplit(affected_genes, ", ")[[1]]))
  if(length(i)){
    l = list()
    for (h in 1:length(i)){
      l[h] = as.character(gene[[h]])
    }
    return(paste(l, collapse = ','))
  }
  return("Other")
  # i <- match(, genes)
  
  # return(sapply(as.character(gene_list), function(x) tolower(gene) %in% tolower(strsplit(x, ", ")[[1]]), USE.NAMES=FALSE))
}


get_SVs <- function(..., all_samples, bed_file, drivers){
  cat("Reading SV data...\n")
  
  all_data <- read.delim(all_samples, header = T)
  
  all_data$affected_genes <- as.character(all_data$affected_genes)

  gene_hits <- all_data %>%
    dplyr::filter(
      !status %in% c('F', 'aF')) %>%
    dplyr::rename(length = length.Kb.,
                  cn     = log2.cnv.) %>%
    dplyr::mutate(type_decomposed = as.character(ifelse(str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>%
    dplyr::group_by(sample, event) %>%
    dplyr::mutate(affected_genes = paste0(affected_genes, collapse = ", ")) %>% 
    # dplyr::mutate(gene_hit = ifelse(geneIn(drivers, affected_genes),'Notch', 'Other')) %>%
    # dplyr::mutate(special_hit = as.character(ifelse(type %in% c("COMPLEX", "DEL") && any(geneIn(bed_file$gene, affected_genes)), 'Hit', 'Other'))) %>%
    dplyr::mutate(gene_hit = as.character(geneIn(drivers, affected_genes)), "Other") %>% 
    dplyr::mutate(gene_hit = as.character(ifelse(gene_hit != "Other", "Notch", "Other"))) %>% 
    dplyr::mutate(special_hit = ifelse(type %in% c("COMPLEX", "DEL"), as.character(geneIn(bed_file$gene, affected_genes)), "Other")) %>% 
    
    dplyr::mutate(cell_fraction = ifelse(chromosome1 %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
    dplyr::mutate(class = 'SV') %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit) %>%
    droplevels()

  
  # all_data %>% 
  #   dplyr::filter(
  #                 !status %in% c('F', 'aF')) %>% 
  #   dplyr::rename(length = length.Kb.,
  #                 cn     = log2.cnv.) %>% 
  #   dplyr::mutate(type_decomposed = as.character(ifelse(str_detect(type, 'COMPLEX'), 'COMPLEX', as.character(type)))) %>%
  #   dplyr::group_by(sample, event) %>% 
  #   dplyr::mutate(gene_hit = as.character(geneIn(drivers, affected_genes))) %>% 
  #   dplyr::ungroup() %>% 
  #   dplyr::select(sample, allele_frequency, gene_hit, affected_genes) %>% 
  #   droplevels() %>% 
  #   View()


  return(as.data.frame(gene_hits))
}


get_SNVs <- function(..., all_samples_snvs, bed_file, drivers){
  cat("Reading SNV data...\n")
  
  snv_data <- read.delim(all_samples_snvs, header = T)
  colnames(snv_data) <- c('sample', 'chromosome',	'pos',	'ref',	'alt',	'trinuc',	'trans',	'decomp_trinuc',	'grouped_trans',	'allele_frequency',	'caller',	'variant_type',	'status',	'snpEff_anno',	'feature',	'gene',	'id')
  
  snv_hits <- snv_data %>% 
    dplyr::filter(...) %>% 
    dplyr::mutate(cell_fraction = ifelse(chromosome %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1, allele_frequency*2)),
                  class="SNV") %>%
    dplyr::mutate(gene_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & gene %in% drivers, 'Notch', 'Other'))) %>% 
    dplyr::mutate(special_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & tolower(gene) %in% bed_file$gene, as.character(gene), 'Other'))) %>% 
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit) %>% 
    droplevels()
  
  return(as.data.frame(snv_hits))
}


get_INDELs <- function(..., all_samples_indels, bed_file, drivers){
  cat("Reading INDEL data...\n")
  
  indel_data <- read.delim(all_samples_indels, header = T)
  # colnames(snv_data) <- c('sample', 'chromosome',	'pos',	'ref',	'alt',	'trinuc',	'trans',	'decomp_trinuc',	'grouped_trans',	'allele_frequency',	'caller',	'variant_type',	'status',	'snpEff_anno',	'feature',	'gene',	'id')
  
  indel_hits <- indel_data %>% 
    dplyr::filter(...) %>% 
    dplyr::mutate(cell_fraction = ifelse(chromosome %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1,allele_frequency*2)),
                  class="INDEL") %>%
    # dplyr::group_by(sample, chromosome, pos) %>% 
    dplyr::mutate(gene_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & gene %in% drivers, 'Notch', 'Other'))) %>% 
    dplyr::mutate(special_hit = as.character(ifelse(status %in% c('HIGH', 'MODERATE') & tolower(gene) %in% bed_file$gene, as.character(gene), 'Other'))) %>% 
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit) %>% 
    droplevels()
  
  return(as.data.frame(indel_hits))
}


get_TEs <- function(..., all_samples_tes, bed_file, drivers){
  cat("Reading TE data...\n")
  te_data <- read.delim(all_samples_tes, header = F)
  
  colnames(te_data) <- c("sample", "bp_no", "caller", "genotype", "chromosome1", "bp1", "gene1", "feature1", "chromosome2", "bp2", "gene2", "feature2", "type", "somefield", "allele_frequency", "confidence")
  
  te_data <- te_data %>% 
    dplyr::filter(...,
                  bp_no == 'bp1') %>% 
    tidyr::separate(sample, into = c("sample", "event"), sep = "_") %>% 
    dplyr::mutate(cell_fraction = ifelse(chromosome1 %in% c("X", "Y"), allele_frequency,
                                         ifelse(allele_frequency*2>1, 1,allele_frequency*2))) %>%
    dplyr::mutate(shortName = sample,
                  class = 'TE') %>% 
    dplyr::mutate(shortName = factor(shortName)) %>%
    dplyr::mutate(gene_hit = as.character(ifelse(!feature1 %in% c('intergenic', 'intron') & gene1 %in% drivers, 'Notch', 'Other'))) %>% 
    dplyr::mutate(special_hit = as.character(ifelse(!feature1 %in% c('intergenic', 'intron') & tolower(gene1) %in% bed_file$gene, as.character(gene1), 'Other'))) %>% 
    droplevels()
  
  te_data <- convert_names(x=te_data, short_to_long=T)
  te_data <- te_data %>% 
    dplyr::select(sample, allele_frequency, cell_fraction, class, gene_hit, special_hit)
  
  return(te_data)
}


highest_N <- function(gene_hits){
  new_df <- gene_hits %>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(highest_n = ifelse(any(gene_hit == 'Notch'), max(cell_fraction[gene_hit=='Notch']), 0)) %>% 
    dplyr::distinct(sample, highest_n) %>% 
    ungroup() %>% 
    as.data.frame()
  
  return(new_df)
}


convert_names <- function(x, short_to_long=F){
  name_conversion <- read.delim('/Users/Nick_curie/Desktop/samples_names_conversion.txt', header=F)
  colnames(name_conversion) <- c("sample", "name_short", "sex", "assay")
  
  if(short_to_long) {
    x$sample <- name_conversion$sample[ match(x$shortName, name_conversion$name_short) ]
  } else {
    x$shortName <- name_conversion$name_short[ match(x$sample, name_conversion$sample)]
  }
  return(x)
}
  

#   if(short_to_long){
#     for(s in x$shortName){
#         if(s %in% name_conversion$name_short){
#           cat("Found sample: ", as.character(name_conversion$sample[name_conversion$name_short==s]), s, "\n")
#           x$sample <- as.character(name_conversion$sample[name_conversion$name_short==s])
#         }
#       }
#   }
#   if(short_to_long){
#     x$sample <- name_conversion$sample[x$shortName]
#   } else {
#     x$shortName <- name_conversion$shortName[x$sample]
#   }
#   return(x)
# }


make_short_name <- function(x){
  df <- x %>%
    dplyr::mutate(shortName = ifelse(grepl('D', sample), sub(".*[R]", "D", sample),
                                     ifelse(grepl('R', sample), sub(".*[R]", "R", sample), as.character(sample))))
  return(df)
}


include <- c("A373R1", "A373R11", "A373R13", "A373R3", "A373R5", "A373R9", 
             "A512R19", "A512R21", "A512R23", "A573R25", "A573R27", "A573R29", 
             "A573R31", "A573R33", "B241R35", "B241R37", "B241R39", "B241R41-1", 
             "B241R43", "B241R45", "B241R49", "B241R51", "B241R53", "B241R55", 
             "B241R57", "B241R59", "B241R61", "B241R63", "D106R23", "D106R25", 
             "D106R27", "D106R29", "D106R31", "D106R33", "HUM-1", "HUM-4", 
             "HUM-7")

plot_tumour_evolution <- function(..., all_samples = '~/Desktop/script_test/alleleFreqs/all_samples_snvs_23719.txt',
                              all_samples_snvs = '~/Desktop/script_test/mutationProfiles/data/annotated_snvs.txt',
                              all_samples_indels = '~/Desktop/script_test/mutationProfiles/data/annotated_indels.txt',
                              all_samples_tes = '/Users/Nick_curie/Desktop/Analysis_pipelines/TEs/all_bps.txt',
                              annotate_with = '~/Desktop/gene2bed/bed/dna_damage_repair.bed', drivers=c('N', 'kuz', "Dl"),
                              show_sample = TRUE, snvs=TRUE, indels=TRUE, tes=TRUE, exclude=NULL){
  
  if(missing(exclude)) exclude <- c()
  
  bed_file <- readBed(bed_in = annotate_with)
  colnames(bed_file) <- c("chrom", "start", "end", "gene", "id")
  gene_hits <- get_SVs(all_samples=all_samples, bed_file=bed_file, drivers=drivers)
  Notch_svs <- highest_N(gene_hits)
  
  if(snvs) snv_hits <- get_SNVs( all_samples_snvs = all_samples_snvs, bed_file=bed_file, drivers=drivers )
  Notch_snvs <- highest_N(snv_hits)
  
  if(indels) indel_hits <- get_INDELs(all_samples_indels = all_samples_indels, bed_file=bed_file, drivers=drivers)
  Notch_indels <- highest_N(indel_hits)
  
  if(tes) te_hits <- get_TEs(all_samples_tes = all_samples_tes, bed_file=bed_file, drivers=drivers)
  Notch_tes <- highest_N(te_hits)
  
  df <- do.call("rbind", list(indel_hits, snv_hits, gene_hits, te_hits))
  
  # Get max Notch accross all mut types
  Notch_hits <- merge(Notch_svs,  Notch_snvs,   by = "sample", all = TRUE)
  Notch_hits <- merge(Notch_hits, Notch_indels, by = "sample", all = TRUE)
  colnames(Notch_hits) <- c('sample', 'sv', 'snv', 'indel')
  Notch_hits <- merge(Notch_hits, Notch_tes, by = "sample", all = TRUE)
  colnames(Notch_hits) <- c('sample', 'sv', 'snv', 'indel', 'te')
  
  Notch_hits[is.na(Notch_hits)] <- 0
  Notch_hits <- Notch_hits %>% 
    dplyr::mutate(highest_n=pmax(sv, snv, indel, te))
  
  # df <- make_short_name(x=df)
  # 
  # Notch_svs <- make_short_name(Notch_hits)
  
  df$sample <- factor(df$sample, levels(fct_reorder(Notch_hits$sample, -Notch_hits$highest_n)))
  
  df <- df %>% 
    dplyr::filter(..., !is.na(sample)) %>% 
    dplyr::mutate(class = factor(as.character(class))) %>% 
    droplevels()
  
  Notch_hits <- Notch_hits %>% 
    dplyr::filter(...) %>% 
    dplyr::mutate(sample=fct_reorder(sample, -highest_n)) %>%
    droplevels()
  
  # df <- plyr::join(df, Notch_svs, by='sample')
  
  special_hits <- df %>% 
    dplyr::filter(special_hit != "Other")
  
  special_hits$sample <- factor(special_hits$sample, levels(fct_reorder(Notch_hits$sample, -Notch_hits$highest_n)))
  
  special_hits %>% 
    dplyr::select(sample, class, special_hit, cell_fraction) %>% 
    print()
    
  p <- ggplot(Notch_hits)
  p <- p + geom_violin(data=df, aes(sample, -cell_fraction, fill = sample, colour=sample), alpha = 0.3, show.legend = FALSE,  adjust=0.3, scale='width')
  p <- p + geom_point(data=special_hits, aes(sample, -cell_fraction), size=2, shape=4, show.legend = FALSE)
  # p <- p + geom_violin(data=df[df$class=='INDEL',], aes(sample, -cell_fraction, colour = sample, fill = sample), position='dodge', size=0.1, alpha = 0.4, show.legend = FALSE)
  p <- p + geom_errorbar(data=Notch_hits, aes(fct_reorder(sample, -highest_n), ymin=-highest_n, ymax=-highest_n), alpha=0.6, size=0.7, show.legend = FALSE)
  
  p <- p + coord_flip()

  p <- p + labs(x = NULL)
  p <- p +  scale_y_continuous('Fraction of cells with mutation', labels=seq(1,0,by=-.25), expand = c(0, 0.01))

  p <- p + slideTheme() +
    theme(
      panel.grid.major.x = element_line(color = "grey80", size = 0.5, linetype = "dotted"),
      # axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1),
      axis.title.x = element_text(size = 15),
      axis.text = element_text(size = 15),
      panel.spacing = unit(2, "lines")
      )
 
  if(!show_sample) p <- p + theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
  p <- p + facet_wrap(~class, nrow=1)
  print(p)
  
  return(Notch_hits)
}


