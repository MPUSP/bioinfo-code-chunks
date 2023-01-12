bioinfo-code-chunks
================================
Michael Jahn

Random selection of useful code chunks and examples for bioinformatics.

This repository is a collection of reusable, self-contained code chunks and examples. The goal is to provide code examples for recurring tasks, that require mostly existing packages and therefore do not fulfill the requirements to become a new package.

### Structure

- `data`: raw data tables or files used for examples
- `docs`: rendered notebooks in `html` format, for example from R markdown
- `source`: notebooks and raw script files

### Contents

- [Hypergeotric test for gene enrichment (R)](https://m-jahn.github.io/bioinfo-code-chunks/hypergeometric-test.nb.html)

### Run examples

The pipelines collected in this repository are self-contained and executable. For R markdown and jupyter notebooks, the code _and_ the documentation are part of one and the same document. If you want to run the examples yourself, there are two possibilities:

1. copy and paste code chunks from the documents in `source/`
2. clone/download the repository, open and re-run notebooks in `source/` (see below)

To download the repository on your local drive use `git clone` in a (linux) terminal:

``` bash
cd /your-target-folder
git clone https://github.com/m-jahn/bioinfo-code-chunks.git
```

Open a pipeline with your favorite IDE (Rstudio, jupyterlab, ...) and execute code (chunks) with the `Run` button. For R markdown notebooks, you can also open an interactive R session and render the document like this:

``` bash
require(rmarkdown)
rmarkdown::render("document.Rmd")
```