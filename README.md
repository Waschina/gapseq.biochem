# `gapseq.biochem`
Exploring gapseq's biochemistry database

## What is `gapseq.biochem`?

`gapseq.biochem` is a small R-package that aims to allow you to explore the
biochemistry database behind gapseq.

It is *not* the gapseq software itself. You can find all
information on gapseq here:

- Publication: [Genome Biology](https://doi.org/10.1186/s13059-021-02295-1)
- Documentation: [ReadTheDocs](https://gapseq.readthedocs.io/en/latest/)
- Source code: [GitHub](https://github.com/jotech/gapseq)

***Important note:*** The biochemistry database in gapseq is derived from 
ModelSEED (https://modelseed.org/). In order to give credits to the developers
and curators of ModelSEED, we strongly encourage you to cite also the
publication for the [ModelSEED database](https://doi.org/10.1093/nar/gkaa746)
when you are using gapseq.

## Installation

`gapseq.biochem` is in its development phase. The current development version
can be installed using:

```R
# install.packages("devtools")
devtools::install_github("Waschina/gapseq.biochem")
```
If you have not installed devtools yet, just un-comment the first line.
