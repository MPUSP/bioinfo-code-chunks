# function to get coverage of arbitrary tracks;
# input can be any Granges or Galignment object that supports method "coverage"
get_coverage <- function(data, Type = "RNAseq", by_strand = TRUE) {
  df_coverage <- list(
    `+` = as.data.frame(coverage(data[strand(data) == "+"])),
    `-` = as.data.frame(coverage(data[strand(data) == "-"]) * -1)
  ) %>%
    dplyr::bind_rows(.id = "Group") %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(seqnames = factor(group_name), Type = Type) %>%
    dplyr::rename(score = value) %>%
    dplyr::select(-group, -group_name) %>%
    dplyr::group_by(Group) %>%
    dplyr::mutate(start = 0:(n() - 1), end = seq_len(n())) %>%
    dplyr::select(seqnames, start, end, score, Type, Group) %>%
    dplyr::mutate(Type = paste0(Type, " [", Group, "]")) %>%
    dplyr::ungroup() %>%
    as.data.frame()
  return(df_coverage)
}

# function to plot a combined coverage track
track_coverage_combined <- function(
    df, start_coord = 0, end_coord = max(df$end),
    track_color = 1:6
) {
  df %>%
    filter(start >= start_coord, end <= end_coord) %>%
    ggplot() +
    geom_area(position = "identity", aes(
      x = start, y = score,
      color = Type, fill = Type
    )) +
    labs(x = "", y = "") +
    theme(legend.position = c(0.9, 0.8), legend.key.size = unit(0.3, "cm")) +
    scale_color_manual(values = alpha(track_color, 0.6)) +
    scale_fill_manual(values = alpha(track_color, 0.4))
}

# function to plot separate coverage tracks using ggcoverage
track_coverage_separate <- function(
    df, start_coord = 0, end_coord = max(df$end),
    track_color = 1:6
) {
  ggplot() +
    geom_coverage(
      data = filter(df, start >= start_coord, end <= end_coord),
      plot.type = "facet",
      group.key = "Type",
      facet.key = "Type",
      facet.color = track_color,
      color = track_color
    ) +
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
    track_color = grey(0.4), track_name = NULL
) {
  if (is.null(feature)) {
    data$feature = "other"
  }
  plt <- data %>%
    filter(end >= start_coord, start <= end_coord) %>%
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
    geom_text(aes(
      x = (xstart + xend) / 2, y = y,
      label = get(name_id), color = feature
    ), size = 2.5, nudge_y = 0.25) +
    coord_cartesian(
      ylim = c(-0.5, 0.75),
      xlim = c(start_coord, end_coord)
    ) +
    labs(x = "", y = "") +
    scale_color_manual(values = track_color) +
    theme(legend.position = 0)

  if (!is.null(track_name)) {
    plt <- plt +
      facet_wrap( ~ track, strip.position = "right") +
      theme(strip.background = element_rect(fill = grey(0.8), color = grey(0.3), ))
  }
  return(plt)
}
