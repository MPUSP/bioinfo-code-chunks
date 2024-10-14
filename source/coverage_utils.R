# function to get coverage of arbitrary tracks;
# input can be any Granges or Galignment object that supports method "coverage"
get_coverage <- function(
    data, format = "bam", sample = "RNAseq",
    by_strand = TRUE, binning = FALSE, bin_width = 100, bin_fun = mean) {
  if (!by_strand) {
    message("Coverage by strand turned off, setting all ranges to '+'")
    strand(data) == "+"
  }
  if (all(strand(data) == "*")) {
    message("Input data is not strand-specific, setting all ranges to '+'")
    strand(data) == "+"
  }
  if (any(strand(data) == "*")) {
    message("Some ranges have undefined strand information, setting those to '+'")
    strand(data[strand(data) == "*"]) <- "+"
  }
  if (format == "bam") {
    df_coverage <- list(
      `+` = as.data.frame(coverage(data[strand(data) == "+"])),
      `-` = as.data.frame(coverage(data[strand(data) == "-"]) * -1)
    ) %>%
      dplyr::bind_rows(.id = "strand") %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(seqnames = factor(group_name), sample = sample) %>%
      dplyr::rename(score = value) %>%
      dplyr::select(-group, -group_name) %>%
      dplyr::group_by(strand) %>%
      dplyr::mutate(start = 0:(n() - 1), end = seq_len(n())) %>%
      dplyr::mutate(sample_strand = paste0(sample, " [", strand, "]")) %>%
      dplyr::select(seqnames, start, end, score, sample, strand, sample_strand) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  } else if (format == "BigWig") {
    df_coverage <- as.data.frame(data) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(
        start = start - 1,
        score = ifelse(strand == "-", sqrt(score^2) * -1, score),
        sample = sample,
        sample_strand = paste0(sample, " [", strand, "]")
      ) %>%
      dplyr::select(seqnames, start, end, score, sample, strand, sample_strand) %>%
      as.data.frame()
  } else {
    stop(paste0("the specified file format '", format, "' is unknown"))
  }
  # reduce complexity by calculating mean of 100 bp windows
  if (binning) {
    df_coverage <- df_coverage %>%
      mutate(
        bin = cut_width(start, width = bin_width, center = bin_width / 2, labels = FALSE) * bin_width
      ) %>%
      group_by(seqnames, sample, strand, sample_strand, bin) %>%
      summarize(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
      mutate(start = bin - (min(bin) - 1), end = bin) %>%
      as.data.frame()
  }
  return(df_coverage)
}

# function to plot a combined coverage track
track_coverage_combined <- function(
    df, start_coord = 0, end_coord = max(df$end),
    track_color = 1:6) {
  df %>%
    filter(start >= start_coord, end <= end_coord) %>%
    ggplot() +
    geom_area(position = "identity", aes(
      x = start,
      y = score,
      color = sample_strand,
      fill = sample_strand
    )) +
    lims(x = c(start_coord, end_coord)) +
    labs(x = "", y = "") +
    theme(legend.position.inside = c(0.9, 0.8), legend.key.size = unit(0.3, "cm")) +
    scale_color_manual(values = alpha(track_color, 0.6)) +
    scale_fill_manual(values = alpha(track_color, 0.4))
}

# function to plot separate coverage tracks using ggcoverage
track_coverage_separate <- function(
    df, start_coord = 0, end_coord = max(df$end),
    track_color = 1:6, facet.y.scale = "free", ...) {
  ggplot() +
    geom_coverage(
      data = filter(df, start >= start_coord, end <= end_coord),
      plot.type = "facet",
      group.key = "sample_strand",
      facet.key = "sample_strand",
      facet.color = track_color,
      color = track_color,
      facet.y.scale = facet.y.scale
    ) +
    lims(x = c(start_coord, end_coord)) +
    labs(x = "", y = "") +
    theme(
      panel.spacing.y = unit(-0.2, "pt"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# function to plot genes as arrows
track_genomic_features <- function(
    data, name_id, start_coord, end_coord, feature = NULL,
    arrow_width = 2.0, arrow_head = 1.0, arrow_type = "open",
    track_color = grey(0.4), track_name = NULL) {
  if (is.null(feature)) {
    data$feature <- "other"
  }
  data <- data %>%
    filter(end >= start_coord, start <= end_coord)
  if (!nrow(data)) {
    data <- add_row(data)
  }
  plt <- data %>%
    mutate(
      track = track_name,
      xstart = case_when(
        strand == "+" ~ start,
        strand == "-" ~ end
      ),
      xend = case_when(
        strand == "+" ~ end,
        strand == "-" ~ start,
      )
    ) %>%
    mutate(y = ifelse(strand == "+", 0.25, -0.25)) %>%
    ggplot() +
    geom_segment(
      aes(x = xstart, xend = xend, y = y, yend = y, color = feature),
      linewidth = arrow_width, lineend = "butt", linejoin = "mitre",
      arrow = arrow(angle = 30, length = unit(arrow_head, "points"), type = arrow_type)
    ) +
    ggrepel::geom_text_repel(aes(
      x = (xstart + xend) / 2, y = y,
      label = get(name_id), color = feature
    ), size = 2.5, nudge_y = 0.25) +
    coord_cartesian(
      ylim = c(-0.5, 0.75),
      xlim = c(start_coord, end_coord)
    ) +
    labs(x = "", y = "") +
    scale_color_manual(values = track_color) +
    theme(legend.position = "none")

  if (!is.null(track_name)) {
    plt <- plt +
      facet_wrap(~track, strip.position = "right") +
      theme(strip.background = element_rect(fill = grey(0.8), color = grey(0.3), ))
  }
  return(plt)
}

# function to find and replace offending values in GFF file
check_gff <- function(gff_file, pattern, replacement, suffix = "_cleaned") {
  gff <- read_lines(gff_file)
  dupl_names <- which(str_detect(gff, pattern))
  if (any(dupl_names)) {
    message(paste0("Found invalid attributes in ", length(dupl_names), " cases:"))
    gff[dupl_names] <- unname(sapply(dupl_names, function(x) {
      match <- str_match(gff[x], pattern)
      repl <- str_match(gff[x], replacement)
      message(str_glue("  Line: {x}, match: '{match}', replacement: '{repl}'"))
      str_replace(gff[x], pattern = pattern, replacement = repl)
    }))
    new_file <- str_replace(gff_file, "\\.gff$", paste0(suffix, ".gff"))
    message(paste0("exported cleaned GFF file to: ", new_file))
    write_lines(gff, file = new_file)
  } else {
    message("found no errors, no new GFF file exported")
  }
}
