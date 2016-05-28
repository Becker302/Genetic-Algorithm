#ifndef MYCGA_H
#define MYCGA_H

#include "cgalib.cpp"

template <size_t N>
class myCGA : public CGAlib<N>{
public:
    explicit myCGA(unsigned int psize, unsigned int tsize, float pcross, float pmut);
    void switchProblem(unsigned short f);
protected:
    double assessFunction(const std::vector<double> fen);
    std::vector<double> decode(const std::bitset<N> gen);
    unsigned short func;
};

#endif // MYCGA_H
