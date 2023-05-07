#ifndef MERGER_H_
#define MERGER_H_

#include <list>  // std::list
#include <utility>  // std::pair

class Merger {
public:
    Merger();
    std::list< std::pair<int, double> > getOTUs() const;
    void pushBackOTU(int i, double length);
    void pushFrontOTU(int i, double length);
private:
    std::list< std::pair<int, double> > otus;  // OTUs merged and branch lengths
};

#endif /* MERGER_H_ */
