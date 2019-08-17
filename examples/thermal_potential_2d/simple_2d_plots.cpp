#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <tuple>
#include <limits>
#include <stdio.h>
#include "simple_2d_potential.hpp"
#include "nlopt_optimizer.hpp"
#include "gnuplot-iostream.h"

typedef std::tuple<double, double, double> data_row;

struct Thermal_Vacua {
    Eigen::VectorXd true_vacuum;
    Eigen::VectorXd false_vacuum;
};

Thermal_Vacua locate_vacua(BubbleProfiler::Simple2DPotential potential,
    Eigen::VectorXd hint_false, Eigen::VectorXd hint_true, 
    double bbox_size) {
    
    Thermal_Vacua vacua;

    const auto v = [&potential](const Eigen::VectorXd& x) -> double {
        return potential(x);
    };

    Eigen::VectorXd lower_bounds(2);
    Eigen::VectorXd upper_bounds(2);

    BubbleProfiler::NLopt_optimizer optimizer(v, 2);
    optimizer.set_extremum_type(BubbleProfiler::NLopt_optimizer::Extremum_type::MIN);
    optimizer.set_max_time(1.0);

    // Locate false vacuum
    lower_bounds << hint_false(0) - bbox_size, hint_false(1) - bbox_size;
    upper_bounds << hint_false(0) + bbox_size, hint_false(1) + bbox_size;
    optimizer.set_lower_bounds(lower_bounds);
    optimizer.set_upper_bounds(upper_bounds);

    auto status = optimizer.optimize(hint_false);
    if (!BubbleProfiler::optimization_succeeded(status)) {
        std::cerr << "Error: unable to locate false vacuum." << std::endl;
        exit(EXIT_FAILURE);
    }

    vacua.false_vacuum = optimizer.get_extremum_location();

    // Locate true vacuum
    lower_bounds << hint_true(0) - bbox_size, hint_true(1) - bbox_size;
    upper_bounds << hint_true(0) + bbox_size, hint_true(1) + bbox_size;
    optimizer.set_lower_bounds(lower_bounds);
    optimizer.set_upper_bounds(upper_bounds);

    status = optimizer.optimize(hint_true);
    if (!BubbleProfiler::optimization_succeeded(status)) {
        std::cerr << "Error: unable to locate true vacuum." << std::endl;
        exit(EXIT_FAILURE);
    }

    vacua.true_vacuum = optimizer.get_extremum_location();

    return vacua;
}

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

void temperature_plots() {
    BubbleProfiler::Simple2DPotential potential(100.0);
    potential.init_effective_potential();

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
}

void bounce_action() {
    BubbleProfiler::Simple2DPotential potential(100.0);
    potential.init_effective_potential();
    potential.set_temperature(100.0);

    Eigen::VectorXd hint_false(2);
    hint_false(0) = 210.;
    hint_false(1) = -150.;

    Eigen::VectorXd hint_true(2);
    hint_true(0) = 300.;
    hint_true(1) = 350.;

    double bbox_size = 200;

    Thermal_Vacua vacua = locate_vacua(potential, hint_false, hint_true, bbox_size);

    std::cout << "False vacuum: " << vacua.false_vacuum << "\n";
    std::cout << "True vacuum: " << vacua.true_vacuum << "\n";
}

int main() {
    // temperature_plots();
    bounce_action();
}
