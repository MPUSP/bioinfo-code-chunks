library(tidyverse)
library(circlize)
library(Biostrings)


# function to validate genomic coordinates against chromosome information
validate_genomic_input <- function(df, chrom_info) {
  chrom_bounds <- chrom_info %>%
    transmute(chr = name, chrom_start = start, chrom_end = end)
  invalid_chr <- df %>%
    anti_join(chrom_bounds, by = "chr")
  if (nrow(invalid_chr) > 0) {
    warning(
      "GFF file contains chromosomes that are missing from FASTA file:\n",
      paste0(invalid_chr$chr, collapse = "\n"),
      call. = FALSE
    )
    df <- df %>%
      semi_join(chrom_bounds, by = "chr")
  }
  out_of_bounds <- df %>%
    inner_join(chrom_bounds, by = "chr") %>%
    filter(start < chrom_start | end > chrom_end)
  if (nrow(out_of_bounds) > 0) {
    warning(
      "GFF file contains coordinates outside chromosome bounds in FASTA file:\n",
      paste0(out_of_bounds$chr, ":", out_of_bounds$start, "-", out_of_bounds$end, collapse = "\n"),
      call. = FALSE
    )
    df <- df %>%
      inner_join(chrom_bounds, by = "chr") %>%
      mutate(start = pmax(start, chrom_start), end = pmin(end, chrom_end)) %>%
      dplyr::select(-chrom_start, -chrom_end)
  }
  return(df)
}

# function to determine GC content and GC skew in sliding windows across the genome
calculate_gc_content <- function(genome_fasta, window = 1000) {
  df_gc_content <- data.frame()
  for (i in names(genome_fasta)) {
    chr_seq <- genome_fasta[[i]]
    chr_width <- length(chr_seq)
    if (chr_width < window) {
      message(paste0(
        "Chromosome ", i, " is shorter than window size (",
        chr_width, " < ", window, "). ",
        "Skipping GC content calculation for this chromosome."
      ))
      next
    }
    starts <- seq(1, chr_width - window + 1, by = window)
    views <- Views(chr_seq, start = starts, width = window)
    gc_freq <- letterFrequency(views, letters = c("G", "C")) / window
    g_content <- gc_freq[, "G"]
    c_content <- gc_freq[, "C"]
    gc_content <- g_content + c_content
    gc_skew <- ifelse(gc_content == 0, 0, (g_content - c_content) / gc_content)
    df_gc <- data.frame(
      chr = i,
      start = (starts - 1),
      end = (starts + window - 1),
      gc = gc_content,
      gc_skew = gc_skew
    )
    df_gc_content <- rbind(df_gc_content, df_gc)
  }
  return(df_gc_content)
}

# function to prepare input data
plot_circlize <- function(genome_fasta, genome_gff, extra = NULL, window = 1000) {
  # create summary df
  df_chrom <- data.frame(
    name = names(genome_fasta),
    start = rep(0, length(genome_fasta)),
    end = width(genome_fasta)
  )

  # create gene info df
  genes <- genome_gff[genome_gff$type == "gene"]
  df_genes <- tibble(
    chr = as.character(seqnames(genes)),
    start = start(genes),
    end = end(genes),
    strand = (as.numeric(strand(genes)) - 1.5) * -1 # -0.5 for "-" and +0.5 for "+"
  )

  # calculate GC content
  df_gc_content <- calculate_gc_content(genome_fasta, window = window)

  # final check to validate if genomic coordinates are OK
  df_genes <- validate_genomic_input(df_genes, df_chrom)
  df_gc_content <- validate_genomic_input(df_gc_content, df_chrom)

  # color vector for sectors (chromosomes)
  sector_colors <- RColorBrewer::brewer.pal(nrow(df_chrom), "Pastel1")

  # base track with cromosome information
  circos.clear()
  on.exit(circos.clear(), add = TRUE)
  circos.par("track.height" = 0.07, start.degree = 90)
  circos.genomicInitialize(df_chrom)
  circos.track(
    ylim = c(0, 1),
    bg.col = sector_colors,
    bg.border = NA, track.height = 0.07
  )

  # add gene track
  circos.genomicTrack(
    df_genes,
    ylim = c(0.5, -0.5),
    bg.border = sector_colors,
    panel.fun = function(region, value, ...) {
      circos.genomicRect(region, value,
        ytop.column = 1, ybottom = 0,
        border = NA,
        col = ifelse(value[[1]] > 0, "#c44040", "#3bb03b"), ...
      )
      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = grey(0.5))
    }
  )

  # add GC content track
  circos.genomicTrack(
    dplyr::select(df_gc_content, chr, start, end, gc),
    ylim = c(0.25, 0.75),
    bg.border = sector_colors,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, value, area = TRUE, border = NA)
    }
  )

  # add GC skew track
  circos.genomicTrack(
    dplyr::select(df_gc_content, chr, start, end, gc_skew),
    ylim = c(-0.35, 0.35),
    bg.border = sector_colors,
    panel.fun = function(region, value, ...) {
      circos.genomicLines(region, replace(value, value > 0, 0), area = TRUE, baseline = 0, border = NA, col = grey(0.8))
      circos.genomicLines(region, replace(value, value < 0, 0), area = TRUE, baseline = 0, border = NA, col = grey(0.6))
    }
  )

  # add extra tracks if provided
  if (!is.null(extra)) {
    for (track in extra) {
      if (!"ylim" %in% names(track)) {
        track$ylim <- c(0, 1)
      }
      if (!"height" %in% names(track)) {
        track$height <- 0.07
      }
      if (!"color" %in% names(track)) {
        track$color <- "#96389f"
      }
      circos.genomicTrack(
        track$data,
        ylim = track$ylim,
        bg.border = sector_colors,
        track.height = track$height,
        panel.fun = function(region, value, ...) {
          if (track$type == "points") {
            circos.genomicPoints(region, value, col = track$color, pch = 19, cex = 0.25, ...)
          } else if (track$type == "lines") {
            circos.genomicLines(region, value, col = track$color, ...)
          } else if (track$type == "rect") {
            circos.genomicRect(region, value,
              ytop.column = 1, ybottom = 0,
              border = NA, col = track$color, ...
            )
          } else {
            stop("Unsupported track type: ", track$type)
          }
        }
      )
    }
  }
}
