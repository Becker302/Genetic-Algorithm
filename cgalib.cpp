#include "cgalib.h"
#include <cassert>
#include <iostream>
#include <algorithm>
#include <ctime>

std::default_random_engine generator;

template <size_t N>
CGAlib<N>::CGAlib(unsigned int psize, unsigned int tsize, float pcross, float pmut):
    pSize(psize),
    tSize(tsize),
    pCross(pcross),
    pMutation(pmut),
    aptitudes(psize),
    population(psize, std::bitset<N>()),
    pool(0)
{
    assert(psize%2 == 0);
    std::cout<<"IN CONSTRUCTOR"<<std::endl;
    //std::default_random_engine generator;
    std::bernoulli_distribution bernoulli(0.5);

    for (unsigned int i = 0; i<pSize; i++){
        for(unsigned int j = 0; j<N; j++){
            if(bernoulli(generator)){
                population[i].flip(j);
            }
        }
        //std::cout<<population[i].to_string()<<std::endl;
    }
    st = SelectionType::Tournament;
    ct = CrossType::TwoPoint;
    ot = OptimType::Minimize;
    abduction = false;
}


//decode the genotype in bynary b1b2b3...bn to a fenotype vector X=(x1,x2,..,xn)e R
template <size_t N>
std::vector<double> CGAlib<N>::decode(std::bitset<N> gen){
    std::vector<double> vec;
    gen = gen;
    vec.push_back(1);
    return vec;
}

//get aptitude of vector X = (x1,x2,..,xn) with the function defined here (one dimension is X = (x1))
template <size_t N>
double CGAlib<N>::assessFunction(std::vector<double> fen){
    fen = fen;
    return 1;
}

template <size_t N>
void CGAlib<N>::unsort(){
    //std::cout<<"IN UNSORT"<<std::endl;
    std::random_shuffle(population.begin(), population.end());
}

template <size_t N>
void CGAlib<N>::update_aptitudes(){
    //std::cout<<"IN UPDATE APTITUDES"<<std::endl;
    std::vector<double> decoded;
    for(unsigned int i = 0; i < pSize; i++){
        decoded = decode(population[i]);
        //std::cout<<decoded[0]<<std::endl;
        aptitudes[i]=assessFunction(decoded);
        //std::cout<<aptitudes[i]<<std::endl;
    }
}

template <size_t N>
void CGAlib<N>::abductionFlag(bool ab){
    abduction = ab;
}

template <size_t N>
void CGAlib<N>::selection(){
    //std::cout<<"IN SELECTION"<<std::endl;
    if(st==SelectionType::Tournament){
        pool.clear();
        double max, min;
        unsigned int itmax, itmin;
        for(unsigned int i = 0; i<pSize; i++){
            max = aptitudes[i];
            itmax = i;
            min = aptitudes[i];
            itmin = i;
            for(unsigned int j=(i+1); j<(i+tSize); j++){
                if(j>=pSize){
                    if(aptitudes[j%pSize]>max){
                        itmax = j%pSize;
                        max = aptitudes[j%pSize];
                    }
                    if(aptitudes[j%pSize]<min) {
                        itmin = j%pSize;
                        min = aptitudes[j%pSize];
                    }
                }
                else{
                    if(aptitudes[j]>max){
                        itmax = j;
                        max = aptitudes[j];
                    }
                    if(aptitudes[j]<min) {
                        itmin = j;
                        min = aptitudes[j];
                    }
                }
            }
            if(ot == OptimType::Maximize){
                pool.push_back(population[itmax]);
            }
            if(ot == OptimType::Minimize){
                pool.push_back(population[itmin]);
            }
        }
        if(abduction){
            std::bernoulli_distribution alien(0.0161);
            if(alien(generator)) abduct();
        }
    }
}

template <size_t N>
void CGAlib<N>::cross(){
    //std::cout<<"IN CROSS"<<std::endl;
    //std::default_random_engine generator;
    if(ct==CrossType::OnePoint){
        std::uniform_int_distribution<int> distribution(1, N-1);
        std::bernoulli_distribution bernoulli(pCross);
        unsigned int crosspoint;
        for(unsigned int i=0;  i<pSize; i+=2){
            if(bernoulli(generator)){
                crosspoint = distribution(generator);
                std::bitset<N> tmp = population[i];
                for(unsigned int j = 0; j<crosspoint; j++){
                    pool[i][j]=pool[i+1][j];
                    pool[i+1][j]=tmp[j];
                }
            }
        }
    }
    if(ct==CrossType::TwoPoint){
        std::uniform_int_distribution<int> distribution(1, N-1);
        std::bernoulli_distribution bernoulli(pCross);
        unsigned int crosspoint1;
        unsigned int crosspoint2;
        for(unsigned int i=0; i<pSize; i+=2){
            if(bernoulli(generator)){
                crosspoint1 = distribution(generator);
                crosspoint2 = distribution(generator);
                while(crosspoint2 == crosspoint1){
                    crosspoint2 = distribution(generator);
                }
                if(crosspoint1>crosspoint2){
                    unsigned int tmp = crosspoint1;
                    crosspoint1 = crosspoint2;
                    crosspoint2 = tmp;
                }
                std::bitset<N> tmp=population[i];
                for(unsigned int j = crosspoint1; j<crosspoint2; j++){
                    pool[i][j]=pool[i+1][j];
                    pool[i+1][j]=tmp[j];
                }
            }
        }
    }
    if(ct == CrossType::Uniform){
        std::bernoulli_distribution distribution(0.5);
        for(unsigned int i = 0; i<pSize; i+=2){
            std::bitset<N> tmp = population[i];
            for(unsigned int j = 0; j<N; j++){
                if(distribution(generator)){
                    pool[i][j]=pool[i+1][j];
                    pool[i+1][j]=tmp[j];
                }
            }
        }
    }
}

template <size_t N>
void CGAlib<N>::mutate(){
    //std::cout<<"IN MUTATE"<<std::endl;
    //std::default_random_engine generator;
    std::uniform_int_distribution<int> aldis(0, N-1);
    std::bernoulli_distribution flipdis(pMutation);
    int allele;
    for(unsigned long i = 0; i<pSize; i++){
        population[i]=pool[i];
        if(flipdis(generator)){
            allele = aldis(generator);
            population[i][allele].flip();
        }
    }
}

template <size_t N>
void CGAlib<N>::run(const unsigned long generations){
    for(unsigned long i = 0; i<generations; i++){
        //std::cout<<"GENERATION "<<i<<std::endl;
        unsort();
        update_aptitudes();
        selection();
        cross();
        mutate();
        //std::cout<<std::endl;
    }
}

template <size_t N>
void CGAlib<N>::printPop(){
    for(unsigned long i = 0; i<pSize; i++){
        std::cout<<population[i]<<std::endl;
    }
}

template <size_t N>
void CGAlib<N>::setCrossType(CrossType type){
    ct = type;
}

template <size_t N>
void CGAlib<N>::setOptimType(OptimType Opt){
    ot == Opt;
}

template <size_t N>
std::vector<double> CGAlib<N>::getActualBest(){
    update_aptitudes();
    unsigned int itmax = 0 ,itmin = 0;
    double max = aptitudes[0];
    double min = aptitudes[0];
    std::vector<double> decoded;
    for(unsigned long i = 1; i<pSize; i++){
        if(aptitudes[i]>max){
            max = aptitudes[i];
            itmax = i;
        }
        if(aptitudes[i]<min){
            min = aptitudes[i];
            itmin = i;
        }
    }
    if(ot == Maximize){
        decoded = decode(population[itmax]);
    }
    else if(ot == Minimize){
        decoded = decode(population[itmin]);
    }

    return decoded;
}

template <size_t N>
void CGAlib<N>::abduct(){
    std::bitset<N> tmp;
    std::vector<double> aux = decode(population[0]);
    double actual = assessFunction(aux);
    double improve;
    for(unsigned int i = 0; i < N; i++){
        tmp = population[0];
        tmp.flip(i);
        aux = decode(tmp);
        improve = assessFunction(aux);
        if(ot == OptimType::Maximize){
            if(improve > actual){
                population[0].flip(i);
                actual = improve;
            }
        }
        else{
            if(improve < actual){
                population[0].flip(i);
                actual = improve;
            }
        }
    }

    aux = decode(population[pSize-1]);
    actual = assessFunction(aux);
    for(unsigned int i = 0; i < N; i++){
        tmp = population[pSize-1];
        tmp.flip(i);
        aux = decode(tmp);
        improve = assessFunction(aux);
        if(ot == OptimType::Maximize){
            if(improve > actual){
                population[pSize-1].flip(i);
                actual = improve;
            }
        }
        else{
            if(improve < actual){
                population[pSize-1].flip(i);
                actual = improve;
            }
        }
    }
}


