#ifndef CGALIB_H
#define CGALIB_H

#include <vector>
#include <bitset>
#include <random>
#include <cmath>

enum SelectionType{Tournament = 0};
enum CrossType{OnePoint=6, TwoPoint=7, Uniform=8};
enum OptimType{Minimize = 10, Maximize = 11};

template <size_t N>
class CGAlib{
public:
    // chormosome size, tournament size, proba cruze, proba muta
    CGAlib(unsigned int psize, unsigned int tsize, float pcross, float pmut);

    void run(const unsigned long generations);
    void printPop();
    void setCrossType(CrossType);
    void setOptimType(OptimType Opt);
    std::vector<double> getActualBest();
    void abductionFlag(bool ab);
protected:
    CGAlib();
    void unsort(); //unsort the actual population
    void selection(); // select by certain policy the pairs to be mating and its puting in the pool
    void cross(); // crossover the singles in the pool
    void mutate(); // pass the mating pool to population and mutate for do an actual genereation
    void update_aptitudes();
    void abduct(); //ignorar, solo era curiosidad

    virtual double assessFunction(const std::vector<double> fen);
    virtual std::vector<double> decode(const std::bitset<N> gen);

    unsigned int pSize; //populate size
    unsigned int tSize; // tournament size
    float pCross; //P of cross
    float pMutation; //p of mutate
    bool abduction;

    SelectionType st;
    CrossType ct;
    OptimType ot;

    std::vector<double> aptitudes;
    std::vector<std::bitset<N> > population, pool;
};

#endif // CGALIB_H
