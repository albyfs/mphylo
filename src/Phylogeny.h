#ifndef PHYLOGENY_H_
#define PHYLOGENY_H_

#include <list>  // std::list
#include <string>  // std::string
#include <vector>  // std::vector

#include "Matrix.h"
#include "Merger.h"

// Phylogenetic Tree
class Phylogeny {
public:
	Phylogeny();
	Phylogeny(const Matrix& dist, int precision);
    void reconstruct();
    int numPolytomies() const;
    std::vector<Merger> getMergers() const;
    std::string getNewick(const std::vector<std::string>& labels) const;
private:
    class Cluster {
    public:
    	Cluster();
    	int prevOTU;  // Previous agglomerable OTU
    	int nextOTU;  // Next agglomerable OTU
		double sumBranches;  // Sums of branch lengths NN TO THE RIGHT
		std::list<int> nearestNeighbors;  // Nearest neighbors TO THE RIGHT
		std::list<int> nearestNeighborOf;  // OTUs that current one is NN of
    };
    int nTaxa;  // Number of taxa
	int nOTUs;  // Number of OTUs still to agglomerate
	int nPolytomies;  // Number of polytomies
    Matrix dist;  // Distances between OTUs
	Matrix sumBranches;  // Sums of branch lengths
    double epsilon;  // Very small number
    int precision;  // Number of significant decimal digits
    double pow10precision;  // 10 to the power of significant decimal digits
    std::vector<Cluster> clusters;  // Clusters
    int firstOTU;  // First agglomerable OTU
	double sMin;  // Minimum sum of branch lengths
	std::list<int> otusMin;  // OTUS with the minimum sum of branch lengths
	std::vector<bool> connected;  // Connected components at the minimum sum
    std::vector<Merger> mergers;  // History of mergers
	void sumBranchLengths();
	void minimizeSumBranches();
    void connectComponents();
    std::list<int> connectedComponent(int i);
    double precisionRound(double value) const;
    void disconnectOTU(int j);
    void agglomerateOTUs();
    void splitOTUs(int i, std::list<int>& subsetI, std::list<int>& subsetIc)
    		const;
    void sumDistances(const std::list<int>& subsetI,
    		const std::list<int>& subsetIc, std::vector<double>& rI,
			std::vector<double>& rIc, double& sumRI, double& sumRIc) const;
    void updateDistances();
    double newDistance(const std::list<int>& subsetI,
    		const std::list<int>& subsetJ) const;
    double sumDistancesBetween(const std::list<int>& subsetI,
    		const std::list<int>& subsetJ) const;
    double sumDistancesWithin(const std::list<int>& subsetI) const;
    void clearNearestNeighbors();
};

#endif /* PHYLOGENY_H_ */
