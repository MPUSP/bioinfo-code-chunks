# bioinfo-code-chunks

Random selection of useful code chunks and examples for bioinformatics.

- Author: Michael Jahn, PhD
- Affiliation: Max Planck Unit for the Science of Pathogens (MPUSP), Berlin, Germany
- License: GPL-v3

This repository is a collection of reusable, self-contained code chunks and examples. The goal is to provide code examples for recurring tasks, that require mostly existing packages and therefore do not fulfill the requirements to become a new package.

### Structure

- `data`: raw data tables or files used for examples
- `docs`: rendered notebooks in `html` format, for example from R markdown
- `output`: output files saved for future reference or tests
- `pipeline`: R markdown or jupyter notebooks with examples
- `source`: script files

### Contents

- [Fast gene set enrichment (R)](https://MPUSP.github.io/bioinfo-code-chunks/fast_gene_set_enrichment.nb.html)
- [Hypergeotric test for gene enrichment (R)](https://MPUSP.github.io/bioinfo-code-chunks/hypergeometric-test.nb.html)
- [Quantify overlap between pathways (R)](https://MPUSP.github.io/bioinfo-code-chunks/quantify_overlap.nb.html)
- [Retrieve KEGG pathway information (R)](https://MPUSP.github.io/bioinfo-code-chunks/get_kegg_pathways.nb.html)
- [Plot sequence logos with `logomaker` (python)](https://MPUSP.github.io/bioinfo-code-chunks/plot_logos.html)
- [Plot coverage tracks (R)](https://MPUSP.github.io/bioinfo-code-chunks/plot_coverage.nb.html)

### Run examples

The pipelines collected in this repository are self-contained and executable. For R markdown and jupyter notebooks, the code _and_ the documentation are part of one and the same document. If you want to run the examples yourself, there are two possibilities:

1. copy and paste code chunks from the documents in `source/` or `pipeline/`
2. clone/download the repository, open and re-run notebooks in `pipeline/` (see below)

To download the repository on your local drive use `git clone` in a (linux) terminal:

```bash
cd /your-target-folder
git clone https://github.com/m-jahn/bioinfo-code-chunks.git
```

Open a pipeline with your favorite IDE ([Rstudio](https://posit.co/download/rstudio-desktop/), [jupyterlab](https://jupyter.org/), ...) and execute code (chunks) with the `Run` button. For R markdown notebooks, you can also open an interactive R session and render the document like this:

```bash
require(rmarkdown)
rmarkdown::render("document.Rmd")
```

### Import on the fly

It's possible to directly import functions from a github repo into running code.
For R code, simply `source` a function from github like this:

```
source("https://raw.githubusercontent.com/MPUSP/bioinfo-code-chunks/main/source/get_kegg_pathways.R")
```
