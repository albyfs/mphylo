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
    install.packages("mphylo-1.0.0.tar.gz", repos = NULL, type = "source")
    ```
Since [mphylo](https://github.com/albyfs/mphylo) includes C++ code, you may need to install first [Rtools](https://cran.r-project.org/bin/windows/Rtools/) in Windows, or [Xcode](https://developer.apple.com/xcode/) in MacOS.


## Tutorial

### Usage

```{r eval = FALSE}
mfnj(x, digits = NULL)
```

| Argument | Description |
| :--- | :--- |
| `x` | A structure of class `dist` containing non-negative distances. |
| `digits` | An integer value specifying the precision, i.e., the number of significant decimal digits to be used for the comparisons between distances. This is an important parameter, since equal distances at a certain precision may become different by increasing its value. Thus, it may be responsible of the existence of tied distances. If the value of this parameter is negative or `NULL` (default), then the precision is automatically set to that of the input distance with the largest number of significant decimal digits. |

### Result

An object of class `mfnj` that describes the multifurcated phylogenetic tree obtained. The object is a list with the following components:

| Component | Description |
| :--- | :--- |
| `call` | The call that produced the result. |
| `digits` | Number of significant decimal digits used as precision. |
| `size` | Number of taxa. |
| `labels` | Labels of the taxa. |
| `nwk` | A string describing the output phylogenetic tree in Newick format. |
| `polytomies` | Number of polytomies in the phylogenetic tree. |

### Example

```{r}
library(mphylo)
```

Let us first create a matrix of phylogenetic distances between taxa:

```{r}
m <- matrix(0, 9, 9)
m[lower.tri(m)] <- c(1.3, 4.3, 4.3, 2.7, 3.0, 1.7, 2.0,  8.7,
                          4.3, 4.3, 2.3, 3.0, 1.7, 2.0,  8.0,
                               0.7, 5.0, 1.3, 2.7, 3.0, 10.0,
                                    5.0, 1.3, 2.7, 3.0, 10.0,
                                         3.7, 2.3, 2.7, 10.0,
                                              2.0, 2.3,  8.7,
                                                   0.3,  9.0,
                                                         9.4)
x <- as.dist(m)
attr(x, "Labels") <- c("Abruzzo", "Pyrenees", "Kodiak", "Captive-3",
                       "Captive-4", "Captive-5", "Grizzly", "Polar-2", "Black")
```

Now we can reconstruct a phylogenetic tree from these distances:

```{r}
t <- mfnj(x, digits = 6)
```

The summary of the resulting phylogenetic tree is:

```{r}
summary(t)
```

And we can also plot this phylogenetic tree:

```{r}
plot(t)
```


## Reference

A. Fernández, N. Segura-Alabart, F. Serratosa. The MultiFurcating Neighbor-Joining Algorithm for Reconstructing Polytomic Phylogenetic Trees. _Journal of Molecular Evolution_ **91**, 773–779 (2023). DOI:[10.1007/s00239-023-10134-z](https://doi.org/10.1007/s00239-023-10134-z).


## Author

- **Alberto Fernández**: Dept. Enginyeria Química, Universitat Rovira i Virgili, Tarragona (Spain). ([email](mailto:alberto.fernandez@urv.cat?subject=[mphylo])) ([ORCID](https://orcid.org/0000-0002-1241-1646)) ([Google Scholar](https://scholar.google.es/citations?user=AbH4r0IAAAAJ)) ([GitHub](https://github.com/albyfs))
