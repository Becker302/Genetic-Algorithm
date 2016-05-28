#include "mycga.h"
#include <cmath>
#include <cassert>
#define PI 3.14159265359

namespace {
    float GetFloat(std::bitset<32>& bits) {
        int *x,y;
        y = bits.to_ulong();
        x = &y;
        float *num = reinterpret_cast<float*>(x);
        return *num;
    }
    double GetDouble(std::bitset<64>& bits) {
        unsigned long int *x,y;
        y = bits.to_ulong();
        x = &y;
        double *num = reinterpret_cast<double*>(x);
        return *num;
    }
}

template <size_t N>
myCGA<N>::myCGA(unsigned int psize, unsigned int tsize, float pcross, float pmut):
    CGAlib<N>(psize, tsize, pcross, pmut),
    func(15)
{
    assert(N==128 || N==64);
}

//64 bits presition(float), 128 bits double
template <size_t N>
std::vector<double> myCGA<N>::decode(std::bitset<N> gen){
    std::vector<double> vec;

    if(N==128){
        std::bitset<64> bits1(0);
        std::bitset<64> bits2(0);
        double real1, real2;
        int counter = 31;
        for(unsigned int i = 0; i<64; i++){
            bits1[counter] = gen[i];
            counter--;
        }
        counter=31;
        for(unsigned int i = 64; i<N; i++){
            bits2[counter]=gen[i];
            counter--;
        }
        real1 = GetDouble(bits1);
        vec.push_back(real1);
        real2 = GetDouble(bits2);
        vec.push_back(real2);
    }
    if(N==64){
        std::bitset<32> bits1(0);
        std::bitset<32> bits2(0);
        double real1, real2;
        int counter = 31;
        for(unsigned int i = 0; i<32; i++){
            bits1[counter] = gen[i];
            counter--;
        }
        counter=31;
        for(unsigned int i = 32; i<N; i++){
            bits2[counter]=gen[i];
            counter--;
        }
        real1 = double(GetFloat(bits1));
        vec.push_back(real1);
        real2 = double(GetFloat(bits2));
        vec.push_back(real2);
    }
    return vec;
}

template <size_t N>
double myCGA<N>::assessFunction(std::vector<double> fen){
    double x = fen[0];
    double y = fen[1];
    double z;
    switch(func){
        case 0:
            //funcion de ackley (0,0)
            z = -20*exp(-0.2*sqrt(0.5*(x*x+y*y))) - exp(0.5*(cos(2*PI*x) + cos(2*PI*y))) + exp(1)+20;
        break;
        case 1:
            //funcion easom (pi,pi)
            z = -1*cos(x)*cos(y)*exp(-1*((x-PI)*(x-PI) + (y-PI)*(y-PI)));
        break;
        case 2:
            //funcion eggholder min (various in other domains)
            z = -(y + 45)*sin(sqrt(fabs(y + (x/2) +47))) - x*sin(sqrt(fabs(x-(y+47))));
        break;
        case 3:
            //levi (1,1)
            z = sin(3*PI*x)*sin(3*PI*x) + (x-1)*(x-1)*(1+sin(3*PI*y)*sin(3*PI*y)) + (y-1)*(y-1)*(1+sin(2*PI*y)*sin(2*PI*y));
        break;
        case 4:
            //schaffer N. 2 (0,0)
            z = 0.5 + ((sin(x*x -y*y)*sin(x*x -y*y) - 0.5)/((1+0.001*(x*x + y*y))*(1+0.001*(x*x + y*y))));
        break;
        case 5:
            //schaffer N. 4 (0,1.25313)
            z = 0.5 + ((cos(sin(fabs(x*x -y*y)))*cos(sin(fabs(x*x -y*y))) - 0.5)/((1+0.001*(x*x + y*y))*(1+0.001*(x*x + y*y))));
        break;
        case 6:
            //bukin (-10, 1)
            z = 100*sqrt(fabs(y-(0.01*x*x))) + 0.01*fabs(x+10);
        break;
        case 7:
            //McCormick (-0.54719, -1.54719)
            z = sin(x+y)+((x-y)*(x-y)) - 1.5*x + 2.5*y +1;
        break;
        default:
            //sfere (0,0)
            z = x*x + y*y;
        break;
    }
    return z;
}

template <size_t N>
void myCGA<N>::switchProblem(unsigned short f){
    func = f;
}
