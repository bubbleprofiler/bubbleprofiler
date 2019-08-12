#include <Eigen/Core>
#include <iostream>
#include "simple_2d_potential.hpp"

int main() {
    BubbleProfiler::Simple2DPotential potential(1.0);
    potential.init_effective_potential();
    Eigen::MatrixXd grid = potential.get_2d_potential_grid(10, -1, 1, -1, 1);
    std::cout << grid;
}
