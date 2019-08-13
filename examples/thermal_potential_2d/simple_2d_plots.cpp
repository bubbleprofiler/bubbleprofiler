#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <tuple>
#include <limits>
#include <stdio.h>
#include "simple_2d_potential.hpp"
#include "gnuplot-iostream.h"

typedef std::tuple<double, double, double> data_row;

double find_minimum(std::vector<data_row> grid) {
    double min = std::numeric_limits<double>::infinity();

    for (auto& row : grid) {
        double z = std::get<2>(row);
        if (z < min) min = z;
    }

    return min;
}

void shift_to_zero(std::vector<data_row>& grid) {
    double min = find_minimum(grid);

    for (auto& row: grid) {
        row = std::make_tuple(
            std::get<0>(row),
            std::get<1>(row),
            std::get<2>(row) - min);
    }
}

void plot_grid(std::vector<data_row> grid, std::string title) {
    remove("/tmp/contour.txt");
    Gnuplot gp("tee /tmp/plot.gp | gnuplot -persist");
    gp << "reset\n";
    gp << "set dgrid3d 200,200,1\n";
    gp << "set contour base\n";
    gp << "unset surface\n";
    gp << "set cntrparam levels auto 100\n";
    gp << "set samples 250,2\n";
    gp << "set isosamples 2,250\n";
    gp << "set view map\n";
    gp << "set table '/tmp/contour.txt'\n";
    gp << "splot '-' using 1:2:(log10($3)) notitle\n";
    gp.send1d(grid);

    gp << "unset contour\n";
    gp << "unset table\n";
    gp << "set autoscale fix\n";

    gp << "set title '" << title << "'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    
    gp << "plot '-' u 1:2:3 w image not, '/tmp/contour.txt' u 1:2 lc black w l not\n";
    gp.send1d(grid);
}

int main() {
    BubbleProfiler::Simple2DPotential potential(100.0);
    potential.init_effective_potential();
    std::shared_ptr<PhaseTracer::Effective_potential> effective_potential = potential.get_effective_potential();

    double plot_scale = 2*246.;
    double x_min = -plot_scale;
    double x_max = plot_scale;
    double y_min = -plot_scale;
    double y_max = plot_scale;
    int axis_size = 200;


    double T = 0;
    double T_max = 100;
    double incr = 20;

    std::vector<data_row> grid;

    while (T <= T_max) {
        std::ostringstream title;
        title << "T = ";
        title << T;
        potential.set_temperature(T);
        grid = potential.get_2d_potential_grid(axis_size, x_min, x_max, y_min, y_max);
        shift_to_zero(grid);
        plot_grid(grid, title.str());
        T += incr;
    }
    
    // potential.set_temperature(100.);

    // std::vector<data_row> grid = potential.get_2d_potential_grid(axis_size, x_min, x_max, y_min, y_max);
    // std::cout << "T = 100:\n";
    // shift_to_zero(grid);
    // plot_grid(grid, "T = 100");

    // potential.set_temperature(0.);

    // grid.clear();
    // grid = potential.get_2d_potential_grid(axis_size, x_min, x_max, y_min, y_max);
    // std::cout << "T = 0:\n";
    // shift_to_zero(grid);
    // plot_grid(grid, "T = 0");


}
