# rSBM: Sequence-Based Mapping in R

This R-based SBM (rSBM) code was originally authored by Dr. Christopher Anderson (Univ. Rochester Medical Center). Please cite the following work:

> Anderson, C. S., Mccall, P. R., Stern, H. A., Yang, H., & Topham, D. J. (2018). Antigenic cartography of H1N1 influenza viruses using sequence-based antigenic distance calculation. BMC Bioinformatics, 19(1). doi:10.1186/s12859-018-2042-4

Sequence-based mapping (SBM) is a technique to *map* antigenic relationships between virus strains. Our work focuses on the antigenic relationships and relatedness between influenza A viruses (IAVs). This algorithm computes the Hamming distance between amino acid sequences to generate *n x n* square-distance matrix. This distance matrix can then be reduced with classic (metric) multidimensional scaling (MDS) to visualize the antigenic relatedness in a two-dimensional space (a plot with *x* and *y* axes).

## Prerequisites

1. FASTA file containing amino acid sequence, aligned with multiple sequence alignment (MSA) program of your choice (MUSCLE, Clustal Omega, etc).
2. `R` ([The R Project for Statistical Computing](https://www.r-project.org/), free and open source). The integrated development environment (IDE) RStudio is not necessary, but it would help a lot ([RStudio](https://rstudio.com/), free).
3. The following R packages:
   1. `Biostrings` from Bioconductor package repository.
   2. `seqinr` from CRAN.
   3. `stringdist` from CRAN.
   4. `maptools` from CRAN.
   5. `FactoMineR` from CRAN.
   6. `tidyverse` from CRAN (optional).

## Files Included In This Repository

1. `rSBM.R` main script.
2. `analysis.Rmd` RMarkdown file for an interactive notebook-style environment.
3. Folder `example` for publicly available IAV H3N2 and H1N1 sequences in FASTA format.
4. Tool to perform sequence alignment. User can use [MEGA X](https://www.megasoftware.net/) or [UGENE](http://ugene.net/) to achieve this.

## How To Run rSBM

The `rSBM.R` contains the functions to perform quality control (`HemagglutininQC()`), to generate distance matrix table (`matrix_generator()`), and to plot the distance (`plot_distance()`).

The RMarkdown `analysis.Rmd` has inline tutorial for you to follow. User is expected to know how to run R and RStudio.

## To-Do/Planned Features

- A separate function to return PC1 and PC2 values, then wrap to call with `ggplot2` for a much more pleasing plotting experience.
- A simple switch to enable user to input range of amino acid for epitope-specific distance calculation (user-defined).

## License

Please refer to each respective libraries for license. As for the `rSBM` algorithms and functions, we offer it under MIT License

```
Copyright (c) 2018-2020 Christopher Anderson, A. Karim Embong (Aizan), David J. Topham 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
