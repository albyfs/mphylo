#include <algorithm>  // std::max, std::min
#include <cmath>  // std::floor, std::log10
#include <sstream>  // std::ostringstream
#include <string>  // std::string
#include <vector>  // std::vector

#include <Rcpp.h>

#include "Matrix.h"
#include "Phylogeny.h"

// [[Rcpp::export]]
Rcpp::List rcppMfnj(const Rcpp::StringVector& labels,
		const Rcpp::NumericVector& x, int digits = -1) {
	Matrix dist(Rcpp::as< std::vector<double> >(x));
	if (digits < 0) {
		digits = dist.precision();
	}
	// Check maximum precision
	double maxDist = std::max(dist.maxValue(), 1.0);
	int intDigits = 1 + (int)std::floor(std::log10(maxDist));
	int maxPrecision = MAX_DIGITS - intDigits - 1;
	digits = std::min(digits, maxPrecision);
	// Reconstruct phylogenetic tree from distances
	Phylogeny* phylo = new Phylogeny(dist, digits);
	phylo->reconstruct();
	// Save results
	Rcpp::List lst = Rcpp::List::create(
			Rcpp::Named("digits") = digits,
			Rcpp::Named("nwk") =
				phylo->getNewick(Rcpp::as< std::vector<std::string> >(labels)),
			Rcpp::Named("polytomies") = phylo->numPolytomies());
	delete phylo;
	return lst;
}
