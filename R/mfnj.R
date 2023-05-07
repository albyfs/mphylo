mfnj <- function(x, digits = NULL) {
	# Check parameters
	if (class(x) != "dist") {
		stop("'x' must be an object of class \"dist\"")
	}
	if (attr(x, "Size") < 3L) {
		stop("'x' must have at least 3 taxa")
	}
	if (anyNA(x)) {
		stop("NA values are not allowed in 'x'")
	}
	if (any(is.nan(x))) {
		stop("NaN values are not allowed in 'x'")
	}
	if (any(is.infinite(x))) {
		stop("Infinite values are not allowed in 'x'")
	}
	storage.mode(x) <- "double"
	if (min(x) < 0) {
		stop("Negative values are not allowed in 'x'")
	}
	size <- attr(x, "Size")
	labels <- attr(x, "Labels")
	if (is.null(labels)) {
		labels <- as.character(seq_len(size))
	}
	if (is.null(digits)) {
		digits <- -1L
	}
	# Reconstruct phylogenetic tree from distances
	lst <- rcppMfnj(labels=as.character(labels), x=as.numeric(x),
			digits=as.integer(digits))
	# Return object of class "mfnj"
	structure(list(
			call = match.call(),
			digits = lst$digits,
			size = size,
			labels = labels,
			nwk = lst$nwk,
			polytomies = lst$polytomies),
		class = "mfnj")
}

summary.mfnj <- function(object, ...) {
	# Print call
	cat("Call:\n", sep="")
	cl <- object$call
	cat(deparse(cl[[1L]]), "(x = ", deparse(cl$x), ",\n", sep="")
	cat("     digits = ", object$digits, ")\n\n", sep="")
	# Print size
	cat("Number of taxa: ", object$size, "\n\n", sep="")
	# Print labels
	cat("Labels:\n", sep="")
	print(object$labels, ...)
	# Print Newick
	cat("\nNewick tree:\n", sep="")
	cat(object$nwk, "\n\n", sep="")
	# Print polytomies
	cat("Number of polytomies: ", object$polytomies, "\n", sep="")
	invisible(object)
}

plot.mfnj <- function (x, ...) {
	phy <- ape::read.tree(text = x$nwk)
	ape::plot.phylo(phy, type = "unrooted", ...)
}
