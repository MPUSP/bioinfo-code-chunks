library(KEGGREST)
library(dplyr)
library(stringr)
library(tibble)

get_kegg_pathways <- function(id) {
  if (length(id) == 1 && all(str_detect(id, "^[a-z]{3}$"))) {
    org_id <- id
  } else if (length(id) >= 1 && all(str_detect(id, "^[a-z]{3}\\:.*"))) {
    org_id <- unique(str_extract(id, "^[a-z]{3}"))
  } else {
    warning("the supplied ID(s) do not match the pattern '<org-id>:<locus-tag>'")
    return(NULL)
  }
  df1 <- keggLink("pathway", id) %>%
    enframe(name = "locus_tag", value = "kegg_pathway_id") %>%
    mutate(kegg_pathway_id = str_remove(kegg_pathway_id, "path:| "))

  df2 <- lapply(org_id, function(org) {
    keggList("pathway", org) %>%
      enframe(name = "kegg_pathway_id", value = "kegg_pathway")
  }) %>%
    bind_rows() %>%
    filter(!duplicated(kegg_pathway))

  df <- left_join(by = "kegg_pathway_id", df1, df2) %>%
    mutate(
      locus_tag = str_remove(locus_tag, paste0(org_id, ":")),
      kegg_pathway_id = str_remove(kegg_pathway_id, "path:"),
      kegg_pathway = str_remove(kegg_pathway, " - .*$")
    )
  return(df)
}