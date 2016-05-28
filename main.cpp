/*
 *  Author: Ivan Martinez
 *  Assigment: Algoritmo Genetico Clasico
 *  Signature: Agentes Inteligentes
 *
 */

#include <iostream>
#include <vector>
#include "mycga.cpp"

/*
 * Se puede cambiar la funcion a optimizar
 * con la funcion
 * switchProblem(unsigned short funcion)
 * hay nueve funciones disponibles de wikipedia
 * para observar, en mycga.cpp estan sus nombres
 *
*/

int main(){
    //cromosoma (template) 2x64(presicion double) o 2x32(presicion float) bits,
    //parametros: poblacion(debe ser par) ,torneo, proba cruze, proba mutacion
    myCGA<64> solver(15000,8,0.9,0.1);
    std::vector<double> solution;

    std::cout<<"Funcion de Ackley"<<std::endl;
    solver.switchProblem(0);
    solver.run(150);
    solution = solver.getActualBest();

    std::cout<<"best x: "<<solution[0]<<std::endl;
    std::cout<<"best y: "<<solution[1]<<std::endl;

    return 0;
}
