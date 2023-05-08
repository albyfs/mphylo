# mphylo
Multifurcated Phylogenetic Trees in R

## Description

[R](https://www.r-project.org) package [mphylo](https://github.com/albyfs/mphylo) offers a **MultiFurcating** version of the **Neighbor-Joining** method for reconstructing phylogenetic trees. Multifurcated phylogenetic trees can group more than two clusters when **tied distances** occur, and therefore they do not depend on the order of the input taxa.

This functionality is obtained with the function `mfnj`, which may be considered as a replacement for function `nj` (in package [ape](https://CRAN.R-project.org/package=ape)).


## Installation

There exist two main ways to install [mphylo](https://github.com/albyfs/mphylo):

- Installation from GitHub (you may need to install first [devtools](https://github.com/r-lib/devtools)):
    ```{r eval = FALSE}
    install.packages("devtools")
    library(devtools)
    install_github("albyfs/mphylo")
    ```
- Installation from source tarball:
    ```{r eval = FALSE}
    install.packages("mphylo_0.1.0.tar.gz", repos = NULL, type = "source")
    ```
Since [mphylo](https://github.com/albyfs/mphylo) includes C++ code, you may need to install first [Rtools](https://cran.r-project.org/bin/windows/Rtools/) in Windows, or [Xcode](https://developer.apple.com/xcode/) in MacOS.


## Author

- **Alberto Fernández**: Dept. Enginyeria Química, Universitat Rovira i Virgili, Tarragona (Spain). ([email](mailto:alberto.fernandez@urv.cat?subject=[mphylo])) ([ORCID](https://orcid.org/0000-0002-1241-1646)) ([Google Scholar](https://scholar.google.es/citations?user=AbH4r0IAAAAJ)) ([GitHub](https://github.com/albyfs))
