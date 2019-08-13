#include <Eigen/Core>
#include <iostream>
#include "simple_2d_potential.hpp"

int main() {
    BubbleProfiler::Simple2DPotential potential(100.0);
    potential.init_effective_potential();
    std::shared_ptr<PhaseTracer::Effective_potential> effective_potential = potential.get_effective_potential();
    Eigen::VectorXd coord(2);
    coord << 100., 100.;
    std::cout << "V_tree(100,100,T=100) " << effective_potential->V_tree(coord, 100.) << "\n";
    std::cout << "V_tree(100,100,T=0) " << effective_potential->V_tree(coord, 0.) << "\n";
    std::cout << "V_daisy(100,100,T=100) " << effective_potential->V_daisy(coord, 100.) << "\n";
    std::cout << "V_one_loop(100,100,T=100) " << effective_potential->V_one_loop(coord, 100.) << "\n";
    std::cout << "V_one_loop(100,100,T=0) " << effective_potential->V_one_loop(coord, 0.) << "\n";
    std::cout << "V1T(100,100,T=100) " << effective_potential->V_finite_temp(coord, 100.) << "\n";
    std::cout << "V(100,100,T=100) " << effective_potential->V(coord, 100.) << "\n";
    
    // Eigen::MatrixXd grid = potential.get_2d_potential_grid(10, -101, 101, -101, 101);
    // std::cout << grid;
}
