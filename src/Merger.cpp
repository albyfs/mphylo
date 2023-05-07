#include "Merger.h"

#include <list>  // std::list
#include <utility>  // std::pair

Merger::Merger() {}

std::list< std::pair<int, double> > Merger::getOTUs() const {
	return this->otus;
}

void Merger::pushBackOTU(int i, double length) {
	std::pair<int, double> otu(i, length);
	this->otus.push_back(otu);
	return;
}

void Merger::pushFrontOTU(int i, double length) {
	std::pair<int, double> otu(i, length);
	this->otus.push_front(otu);
	return;
}
