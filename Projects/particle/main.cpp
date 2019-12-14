#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "SimulationDriver.h"

int main(int argc, char* argv[])
{
    std::cout << "main called" << std::endl;
    using T = float;
    constexpr int dim = 3;
    //  using TV = Eigen::Matrix<T,dim,1>;
    SimulationDriver<T,dim> driver;
    driver.run(240);

    return 0;
};
