#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <tuple>
#include "simple_2d_potential.hpp"
#include "gnuplot-iostream.h"

void print_grid(std::vector<std::tuple<double, double, double>> grid) {
    for (auto& t: grid) {
        std::cout << "x: " << std::get<0>(t) 
                  << ", y: " << std::get<1>(t) 
                  << ", z: " << std::get<2>(t) << '\n';
    }
}

void plot_grid(std::vector<std::tuple<double, double, double>> grid) {
    Gnuplot gp("tee /tmp/plot.gp | gnuplot -persist");

    // gp << "set dgrid3d\n";
    // gp << "unset surface\n";
    // gp << "set view map\n";
    // gp << "set contour\n";
    // gp << "set key outside\n";
    // gp << "set cntrparam cubicspline\n";
    // gp << "set cntrparam levels 50\n";
    // gp << "set pm3d interpolate 20,20\n";
    // gp << "set palette rgbformulae 33,13,10\n";
    // gp << "set xlabel 'x'\n";
    // gp << "set ylabel 'y'\n";
    
    // gp << "splot '-' using 1:2:3 notitle with lines lt 1\n";
    gp << "plot '-' using 1:2:3 with image\n";
    gp.send1d(grid);
}

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

    double v = 246.;
    double x_min = -v;
    double x_max = v;
    double y_min = -v;
    double y_max = v;
    int axis_size = 200;
    
    std::vector<std::tuple<double, double, double>> grid = potential.get_2d_potential_grid(axis_size, x_min, x_max, y_min, y_max);
    std::cout << "T = 0:\n";
    print_grid(grid);
    plot_grid(grid);

    potential.set_temperature(100.);

    grid = potential.get_2d_potential_grid(axis_size, x_min, x_max, y_min, y_max);
    std::cout << "T = 100:\n";
    print_grid(grid);
    plot_grid(grid);
}
