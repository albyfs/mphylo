#include <algorithm>  // std::max, std::min
#include <cmath>  // std::abs, std::floor, std::log10, std::pow, std::round
#include <list>  // std::list
#include <queue>  // std::queue
#include <sstream>  // std::ostringstream
#include <string>  // std::string
#include <vector>  // std::vector

#include "Matrix.h"
#include "Merger.h"
#include "Phylogeny.h"

Phylogeny::Cluster::Cluster() {
	this->prevOTU = -1;
	this->nextOTU = -1;
	this->sumBranches = +INF;
}

Phylogeny::Phylogeny() {
	this->nTaxa = 0;
	this->nOTUs = 0;
	this->nPolytomies = 0;
	this->epsilon = std::pow(1.0, -(double)MAX_DIGITS);
	this->precision = 6;
	this->pow10precision = 1e6;
	this->firstOTU = -1;
	this->sMin = +INF;
}

Phylogeny::Phylogeny(const Matrix& dist, int precision) {
	this->nTaxa = dist.numRows();
	this->nOTUs = this->nTaxa;
	this->nPolytomies = 0;
	this->dist = Matrix(dist);
	this->sumBranches = Matrix(this->nTaxa);
	double maxDist = std::max(std::abs(dist.maxValue()), 1.0);
	int intDigits = 1 + (int)std::floor(std::log10(maxDist));
	int maxPrecision = MAX_DIGITS - intDigits - 1;
	this->epsilon = std::pow(10.0, -(double)(maxPrecision + 1));
	// 0 <= precision <= maxPrecision
	this->precision = std::max(precision, 0);
	this->precision = std::min(this->precision, maxPrecision);
	this->pow10precision = std::pow(10.0, (double)this->precision);
	// Initial partition of OTUs
	this->clusters = std::vector<Cluster>(this->nTaxa);
	for (int i = 0; i < this->nTaxa; i ++) {
	    this->clusters[i].prevOTU = i - 1;
	    this->clusters[i].nextOTU = i + 1;
		this->clusters[i].sumBranches = +INF;
	}
	this->firstOTU = 0;
	this->sMin = +INF;
	this->mergers.reserve(this->nTaxa - 1);
}

void Phylogeny::reconstruct() {
	// Repeat while there are OTUs to agglomerate
	while (this->nOTUs > 1) {
		sumBranchLengths();
		minimizeSumBranches();
		connectComponents();
		agglomerateOTUs();
		updateDistances();
		clearNearestNeighbors();
	}
	return;
}

int Phylogeny::numPolytomies() const {
	return this->nPolytomies;
}

std::vector<Merger> Phylogeny::getMergers() const {
	return this->mergers;
}

std::string Phylogeny::getNewick(const std::vector<std::string>& labels) const {
	std::vector<std::string> newick = labels;
	std::ostringstream oss;
	oss.setf(std::ios::fixed, std::ios::floatfield);  // Fixed precision
	oss.precision(std::max(this->precision, 0));  // Modify default precision
	for (int i = 0; i < (int)this->mergers.size(); i ++) {
		std::list< std::pair<int, double> > otus = this->mergers[i].getOTUs();
		std::list< std::pair<int, double> >::const_iterator it = otus.begin();
		std::pair<int, double> otu = *it;
		int j = otu.first;
		int jmin = j;
		double length = otu.second;
		oss.str("");  // clear oss
		oss << "(" << newick[j] << ":" << length;
		it ++;
		while (it != otus.end()) {
			otu = *it;
			j = otu.first;
			length = otu.second;
			oss << "," << newick[j] << ":" << length;
			it ++;
		}
		oss << ")";
		newick[jmin] = oss.str();
	}
	newick[0] += ";";
	return newick[0];
}

void Phylogeny::sumBranchLengths() {
	// R_i = sum_k D_ik
    std::vector<double> r = std::vector<double>(this->nTaxa, 0.0);
	int i = this->firstOTU;
	while (i < this->nTaxa) {
		int k = this->firstOTU;
		while (k < this->nTaxa) {
			if (k != i) {
				double dik = this->dist.value(i, k);
				r[i] += dik;
			}
			k = this->clusters[k].nextOTU;
		}
		i = this->clusters[i].nextOTU;
	}
	// S_ij = (N - 2) D_ij - R_i - R_j
	i = this->firstOTU;
	while (i < this->nTaxa) {
		int j = this->clusters[i].nextOTU;
		while (j < this->nTaxa) {
			double dij = this->dist.value(i, j);
			double sij = (this->nOTUs - 2) * dij - r[i] - r[j];
			this->sumBranches.setValue(i, j, sij);
			j = this->clusters[j].nextOTU;
		}
		i = this->clusters[i].nextOTU;
	}
	return;
}

void Phylogeny::minimizeSumBranches() {
	// Get the minimum sum of branch lengths and the corresponding list of OTUs
	this->sMin = +INF;
	int i = this->firstOTU;
	while (i < this->nTaxa) {
		// Get the minimum sum of branch lengths TO THE RIGHT of i
		double siMin = +INF;
		int jmin = -1;
		int j = this->clusters[i].nextOTU;
		while (j < this->nTaxa) {
		    double sij = precisionRound(this->sumBranches.value(i, j));
		    if (sij < siMin) {
		    	siMin = sij;
		    	jmin = j;
		    }
		    j = this->clusters[j].nextOTU;
		}
		this->clusters[i].sumBranches = siMin;
		if (jmin > -1) {
			// Add the first nearest neighbor TO THE RIGHT of i
			this->clusters[i].nearestNeighbors = {jmin};
			this->clusters[jmin].nearestNeighborOf.push_back(i);
		}
		if (siMin < this->sMin) {
			this->sMin = siMin;
			this->otusMin.clear();
			this->otusMin = {i};
		} else if (siMin == this->sMin) {
			this->otusMin.push_back(i);
    	}
		i = this->clusters[i].nextOTU;
	}
	return;
}

void Phylogeny::connectComponents() {
	// Complete nearest neighbors of minimum OTUs
	std::list<int>::iterator itmin = this->otusMin.begin();
	while (itmin != this->otusMin.end()) {
		int i = *itmin;
		// Include additional nearest neighbors TO THE RIGHT of i
		int j = this->clusters[i].nearestNeighbors.front();
		j = this->clusters[j].nextOTU;
		while (j < this->nTaxa) {
			double sij = precisionRound(this->sumBranches.value(i, j));
			if (sij == this->sMin) {
				this->clusters[i].nearestNeighbors.push_back(j);
				this->clusters[j].nearestNeighborOf.push_back(i);
			}
			j = this->clusters[j].nextOTU;
		}
		itmin ++;
	}
	// Connected components of minimum OTUs and their nearest neighbors
	this->connected = std::vector<bool>(this->nTaxa, false);
	itmin = this->otusMin.begin();
	while (itmin != this->otusMin.end()) {
		int i = *itmin;
		if (this->connected[i]) {
			// Remove from list of minimum OTUS if it is already connected
			itmin = this->otusMin.erase(itmin);
		} else {
			std::list<int> subsetI = connectedComponent(i);
			this->clusters[i].nearestNeighbors.clear();
			this->clusters[i].nearestNeighbors = subsetI;
			itmin ++;
		}
	}
	return;
}

std::list<int> Phylogeny::connectedComponent(int i) {
	std::list<int> subsetI;
	std::queue<int> q;
	q.push(i);
	while (!q.empty()) {
		int j = q.front();
		q.pop();
		if (!this->connected[j]) {
			this->connected[j] = true;
			if (!subsetI.empty()) {
				// If j is not the first, disconnect it from agglomerable OTUs
				disconnectOTU(j);
			}
			subsetI.push_back(j);
			double sj = precisionRound(this->clusters[j].sumBranches);
			if (sj == this->sMin) {
				std::list<int>::const_iterator itnn =
						this->clusters[j].nearestNeighbors.begin();
				while (itnn != this->clusters[j].nearestNeighbors.end()) {
					int k = *itnn;
					q.push(k);
					itnn ++;
				}
			}
			std::list<int>::iterator itnnof =
					this->clusters[j].nearestNeighborOf.begin();
			while (itnnof != this->clusters[j].nearestNeighborOf.end()) {
				int k = *itnnof;
				double sk = precisionRound(this->clusters[k].sumBranches);
				if (sk == this->sMin) {
					q.push(k);
				}
				itnnof ++;
			}
		}
	}
	subsetI.sort();
	return subsetI;
}

double Phylogeny::precisionRound(double value) const {
	// Add epsilon to avoid 0.49999999999999... being rounded to 0
	value += (value >= 0.0)? +this->epsilon : -this->epsilon;
	return std::round(value * this->pow10precision) / this->pow10precision;
}

void Phylogeny::disconnectOTU(int j) {
	int i = this->clusters[j].prevOTU;
	int k = this->clusters[j].nextOTU;
	if (i < 0) {
		this->firstOTU = k;
	} else {
		this->clusters[i].nextOTU = k;
	}
	if (k < this->nTaxa) {
		this->clusters[k].prevOTU = i;
	}
	this->clusters[j].prevOTU = -1;
	this->clusters[j].nextOTU = -1;
	return;
}

void Phylogeny::agglomerateOTUs() {
	std::list<int>::const_iterator itmin = this->otusMin.begin();
	while (itmin != this->otusMin.end()) {
		int i = *itmin;
		std::list<int> subsetI;
		std::list<int> subsetIc;
		splitOTUs(i, subsetI, subsetIc);
		std::vector<double> rI;
		std::vector<double> rIc;
		double sumRI;
		double sumRIc;
		sumDistances(subsetI, subsetIc, rI, rIc, sumRI, sumRIc);
		int nI = subsetI.size();
		int nIc = subsetIc.size();
		if ((nI > 2) && (this->nOTUs > 3)) {
			this->nPolytomies ++;
		}
		// Agglomerate OTUs into a new merger
		Merger merger;
		std::list<int>::const_iterator itI = subsetI.begin();
		while (itI != subsetI.end()) {
			int j = *itI;
			double length;
			if (nIc > 0) {  // there are still OTUs to agglomerate later
				length = sumRI / (double)(nI * (nI - 1)) + rIc[j] / (double)nIc
						- sumRIc / (double)(nI * nIc);
			} else {  // all remaining OTUs agglomerated together
				length = rI[j] / (double)(nI - 2)
						- sumRI / (double)((nI - 1) * (nI - 2));
			}
			merger.pushBackOTU(j, length);
			itI ++;
		}
		this->nOTUs -= nI - 1;
		if (this->nOTUs == 2) {  // there are only 2 remaining OTUs
			double length = newDistance(subsetI, subsetIc);
			std::list<int>::const_iterator itIc = subsetIc.begin();
			int j = *itIc;
			merger.pushFrontOTU(j, length);
			this->nOTUs -= 1;
		}
		this->mergers.push_back(merger);
		itmin ++;
	}
	return;
}

void Phylogeny::splitOTUs(int i, std::list<int>& subsetI,
		std::list<int>& subsetIc) const {
	subsetI = this->clusters[i].nearestNeighbors;
	std::list<int>::const_iterator it = this->otusMin.begin();
	while (it != this->otusMin.end()) {
		int j = *it;
		if (j != i) {
			std::list<int> subsetJ = this->clusters[j].nearestNeighbors;
			// Append subsetJ at the end of subsetIc
			subsetIc.splice(subsetIc.end(), subsetJ);
		}
		it ++;
	}
	int k = this->firstOTU;
	while (k < this->nTaxa) {
		if (!this->connected[k]) {
			subsetIc.push_back(k);
		}
		k = this->clusters[k].nextOTU;
	}
	return;
}

void Phylogeny::sumDistances(const std::list<int>& subsetI,
		const std::list<int>& subsetIc, std::vector<double>& rI,
		std::vector<double>& rIc, double& sumRI, double& sumRIc) const {
    rI = std::vector<double>(this->nTaxa, 0.0);
    rIc = std::vector<double>(this->nTaxa, 0.0);
    sumRI = 0.0;
    sumRIc = 0.0;
	std::list<int>::const_iterator iti = subsetI.begin();
	while (iti != subsetI.end()) {
		int i = *iti;
		std::list<int>::const_iterator itj = subsetI.begin();
		while (itj != subsetI.end()) {
			int j = *itj;
			if (j != i) {
				double dij = this->dist.value(i, j);
				rI[i] += dij;
				if (j > i) {
					sumRI += dij;
				}
			}
			itj ++;
		}
		std::list<int>::const_iterator itk = subsetIc.begin();
		while (itk != subsetIc.end()) {
			int k = *itk;
			double dik = this->dist.value(i, k);
			rIc[i] += dik;
			sumRIc += dik;
			itk ++;
		}
		iti ++;
	}
	return;
}

void Phylogeny::updateDistances() {
	std::list<int>::const_iterator iti = this->otusMin.begin();
	while (iti != this->otusMin.end()) {
		int i = *iti;
		std::list<int> subsetI = this->clusters[i].nearestNeighbors;
		std::list<int>::const_iterator itj = iti;
		itj ++;
		while (itj != this->otusMin.end()) {
			int j = *itj;
			std::list<int> subsetJ = this->clusters[j].nearestNeighbors;
			double dij = newDistance(subsetI, subsetJ);
			this->dist.setValue(i, j, dij);
			itj ++;
		}
		int k = this->firstOTU;
		while (k < this->nTaxa) {
			if (!this->connected[k]) {
				std::list<int> subsetK = {k};  // list with a single object
				double dik = newDistance(subsetI, subsetK);
				this->dist.setValue(i, k, dik);
			}
			k = this->clusters[k].nextOTU;
		}
		iti ++;
	}
	return;
}

double Phylogeny::newDistance(const std::list<int>& subsetI,
		const std::list<int>& subsetJ) const {
	int nI = subsetI.size();
	int nJ = subsetJ.size();
	double rIJ = sumDistancesBetween(subsetI, subsetJ);
	double dij = rIJ / (double)(nI * nJ);
	if (nI > 1) {
		double rII = sumDistancesWithin(subsetI);
		dij -= rII / (double)(nI * (nI - 1));
	}
	if (nJ > 1) {
		double rJJ = sumDistancesWithin(subsetJ);
		dij -= rJJ / (double)(nJ * (nJ - 1));
	}
	return dij;
}

double Phylogeny::sumDistancesBetween(const std::list<int>& subsetI,
		const std::list<int>& subsetJ) const {
    double rIJ = 0.0;
    std::list<int>::const_iterator iti = subsetI.begin();
    while (iti != subsetI.end()) {
    	int i = *iti;
    	std::list<int>::const_iterator itj = subsetJ.begin();
    	while (itj != subsetJ.end()) {
    		int j = *itj;
    		rIJ += this->dist.value(i, j);
    		itj ++;
    	}
    	iti ++;
    }
	return rIJ;
}

double Phylogeny::sumDistancesWithin(const std::list<int>& subsetI) const {
    double rII = 0.0;
    std::list<int>::const_iterator it1 = subsetI.begin();
    while (it1 != subsetI.end()) {
    	int i1 = *it1;
    	std::list<int>::const_iterator it2 = it1;
    	it2 ++;
    	while (it2 != subsetI.end()) {
    		int i2 = *it2;
    		rII += this->dist.value(i1, i2);
    		it2 ++;
    	}
    	it1 ++;
    }
	return rII;
}

void Phylogeny::clearNearestNeighbors() {
	int i = this->firstOTU;
	while (i < this->nTaxa) {
	  	this->clusters[i].nearestNeighbors.clear();
		this->clusters[i].nearestNeighborOf.clear();
		i = this->clusters[i].nextOTU;
	}
	return;
}
