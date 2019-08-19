#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <tuple>
#include <limits>
#include <stdio.h>
#include "gnuplot-iostream.h"

#include "simple_2d_potential.hpp"
#include "nlopt_optimizer.hpp"
#include "profile_guesser.hpp"
#include "kink_profile_guesser.hpp"
#include "field_profiles.hpp"
#include "perturbative_profiler.hpp"

namespace BubbleProfiler {

typedef std::tuple<double, double, double> data_row;

struct Thermal_Vacua {
    Eigen::VectorXd true_vacuum;
    double val_true_vacuum;
    Eigen::VectorXd false_vacuum;
    double val_false_vacuum;
};

Thermal_Vacua locate_vacua(Simple2DPotential potential,
    Eigen::VectorXd hint_false, Eigen::VectorXd hint_true, 
    double bbox_size) {
    
    Thermal_Vacua vacua;

    const auto v = [&potential](const Eigen::VectorXd& x) -> double {
        return potential(x);
    };

    Eigen::VectorXd lower_bounds(2);
    Eigen::VectorXd upper_bounds(2);

    NLopt_optimizer optimizer(v, 2);
    optimizer.set_extremum_type(NLopt_optimizer::Extremum_type::MIN);
    optimizer.set_max_time(1.0);

    // Locate false vacuum
    lower_bounds << hint_false(0) - bbox_size, hint_false(1) - bbox_size;
    upper_bounds << hint_false(0) + bbox_size, hint_false(1) + bbox_size;
    optimizer.set_lower_bounds(lower_bounds);
    optimizer.set_upper_bounds(upper_bounds);

    auto status = optimizer.optimize(hint_false);
    if (!optimization_succeeded(status)) {
        std::cerr << "Error: unable to locate false vacuum." << std::endl;
        exit(EXIT_FAILURE);
    }

    vacua.false_vacuum = optimizer.get_extremum_location();
    vacua.val_false_vacuum = optimizer.get_extremum_value();

    // Locate true vacuum
    lower_bounds << hint_true(0) - bbox_size, hint_true(1) - bbox_size;
    upper_bounds << hint_true(0) + bbox_size, hint_true(1) + bbox_size;
    optimizer.set_lower_bounds(lower_bounds);
    optimizer.set_upper_bounds(upper_bounds);

    status = optimizer.optimize(hint_true);
    if (!optimization_succeeded(status)) {
        std::cerr << "Error: unable to locate true vacuum." << std::endl;
        exit(EXIT_FAILURE);
    }

    vacua.true_vacuum = optimizer.get_extremum_location();
    vacua.val_true_vacuum = optimizer.get_extremum_value();
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

void apply_cutoff(std::vector<data_row>& grid, double cutoff) {
    for (auto& row: grid) {
        row = std::make_tuple(
            std::get<0>(row),
            std::get<1>(row), 
            std::min(std::get<2>(row), cutoff));
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
    Simple2DPotential potential(90.0);
    potential.init_effective_potential();

    double plot_scale = 1.7*246.;
    // double x_min = -plot_scale;
    double x_min = 0;
    // double x_max = plot_scale;
    double x_max = plot_scale;
    double y_min = -plot_scale;
    double y_max = plot_scale;
    int axis_size = 200;


    double T = 115;
    double T_max = 125;
    double incr = 1;

    std::vector<data_row> grid;

    while (T <= T_max) {
        std::ostringstream title;
        title << "T = ";
        title << T;
        potential.set_temperature(T);
        grid = potential.get_2d_potential_grid(axis_size, x_min, x_max, y_min, y_max);
        shift_to_zero(grid);
        apply_cutoff(grid, 3e8);
        plot_grid(grid, title.str());
        T += incr;
    }
}

void bounce_action() {
    Simple2DPotential potential(90.0);
    potential.init_effective_potential();
    potential.set_temperature(105.0);

    Eigen::VectorXd hint_false(2);
    hint_false(0) = 208.;
    hint_false(1) = -128.;

    Eigen::VectorXd hint_true(2);
    hint_true(0) = 250.;
    hint_true(1) = 250.;

    double bbox_size = 100;

    Thermal_Vacua vacua = locate_vacua(potential, hint_false, hint_true, bbox_size);

    // Need to shift origin to false vacuum...
    potential.translate_origin(vacua.false_vacuum);
    Eigen::VectorXd shifted_true_vacuum = vacua.true_vacuum - vacua.false_vacuum;

    std::cout << "False vacuum: \n" << vacua.false_vacuum << "\n";
    std::cout << "V(false) = " << vacua.val_false_vacuum << "\n";
    std::cout << "True vacuum: \n" << vacua.true_vacuum << "\n";
    std::cout << "V(true) = " << vacua.val_true_vacuum << "\n";
    std::cout << "Shifted true vacuum: \n" << shifted_true_vacuum << "\n";

    std::shared_ptr<Kink_profile_guesser> kink_guesser
        = std::make_shared<Kink_profile_guesser>();
    std::shared_ptr<Profile_guesser> guesser(kink_guesser);

    Field_profiles ansatz = guesser->get_profile_guess(
        potential, shifted_true_vacuum, 2, -1.0, -1.0, 1.e-3, 100.);
    

}

}

int main() {
    // BubbleProfiler::temperature_plots();
    BubbleProfiler::bounce_action();
    return 0;
}