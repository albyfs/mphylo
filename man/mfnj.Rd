\name{mfnj}
\alias{mfnj}
\alias{plot.mfnj}
\title{MultiFurcating Neighbor-Joining}
\description{
		A MultiFurcating version of the Neighbor-Joining method for reconstructing
		phylogenetic trees. Multifurcated phylogenetic trees can group more than two
		clusters when tied distances occur, and therefore they do not depend on the
		order of the input taxa.
}
\usage{
mfnj(x, digits = NULL)
}
\arguments{
    \item{x}{A structure of class \code{"dist"} containing non-negative
        distances.}
    \item{digits}{An integer value specifying the precision, i.e. the number of
        significant decimal digits to be used for the comparisons between
        distances. This is an important parameter, since equal distances at a
        certain precision may become different by increasing its value. Thus, it
        may be responsible of the existence of tied distances. If the value of
        this parameter is negative or \code{NULL} (default), then the precision
        is automatically set to that of the input distance with the largest
        number of significant decimal digits.}
}
\value{
    An object of class \code{"mfnj"} that describes the multifurcated
    phylogenetic tree obtained. The object is a list with the following
    components:
    \item{call}{The call that produced the result.}
    \item{digits}{Number of significant decimal digits used as precision.}
    \item{size}{Number of taxa.}
    \item{labels}{Labels of the taxa.}
    \item{nwk}{A string describing the output phylogenetic tree in Newick
        format.}
    \item{polytomies}{Number of polytomies in the phylogenetic tree.}

    Class \code{"mfnj"} has methods for the following generic functions:
    \code{\link{summary}} and \code{\link{plot}}.
}
\author{
    Alberto Fernandez \email{alberto.fernandez@urv.cat}.
}
\references{
    Randi, E.; Gentile, L.; Boscagli, G.; Huber, D.; Roth, H.U. (1994).
    Mitochondrial DNA sequence divergence among some west European brown bear
    (\emph{Ursus arctos} L.) populations. Lessons for conservation.
    \emph{Heredity}, \bold{73}(5): 480--489.
}
\seealso{
    \code{\link{dist}}.
}
\examples{
## Table 4 from Randi et al. (1994)
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

## Reconstruct phylogenetic tree from distances
t <- mfnj(x, digits = 6)
summary(t)
plot(t)
}
