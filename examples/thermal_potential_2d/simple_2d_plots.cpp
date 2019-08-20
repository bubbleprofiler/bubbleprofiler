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
    double bbox_x, double bbox_y) {
    
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
    lower_bounds << hint_false(0) - bbox_x, hint_false(1) - bbox_y;
    upper_bounds << hint_false(0) + bbox_x, hint_false(1) + bbox_y;
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
    lower_bounds << hint_true(0) - bbox_x, hint_true(1) - bbox_y;
    upper_bounds << hint_true(0) + bbox_x, hint_true(1) + bbox_y;
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

void temperature_plot(Simple2DPotential potential) {
    potential.init_effective_potential();

    double plot_scale = 0.5*246.;
    double x_min = -plot_scale;
    double x_max = plot_scale;
    double y_min = -plot_scale;
    double y_max = plot_scale;
    // double x_min = -1;
    // double x_max = 1;
    // double y_min = -1;
    // double y_max = 1;
    int axis_size = 200;

    std::ostringstream title;
    title << "T = ";
    title << potential.get_temperature();
    std::vector<data_row> grid = potential.get_2d_potential_grid(axis_size, x_min, x_max, y_min, y_max);
    shift_to_zero(grid);
    // apply_cutoff(grid, 1e8);
    plot_grid(grid, title.str());
}

void bounce_action(Simple2DPotential potential, Eigen::VectorXd true_vacuum, Eigen::VectorXd false_vacuum) {

    // Eigen::VectorXd hint_false(2);
    // hint_false(0) = hint_false_x;
    // hint_false(1) = hint_false_y;

    // Eigen::VectorXd hint_true(2);
    // hint_true(0) = hint_true_x;
    // hint_true(1) = hint_true_y;

    // Thermal_Vacua vacua = locate_vacua(potential, hint_false, hint_true, 
    //     bbox_x, bbox_y);

    // Need to shift origin to false vacuum...
    potential.translate_origin(false_vacuum); 
    Eigen::VectorXd origin = Eigen::VectorXd::Zero(2);
    Eigen::VectorXd shifted_true_vacuum = true_vacuum - false_vacuum;

    std::cout << "False vacuum: \n" << false_vacuum << "\n";
    // std::cout << "V(false) = " << vacua.val_false_vacuum << "\n";
    std::cout << "True vacuum: \n" << true_vacuum << "\n";
    // std::cout << "V(true) = " << vacua.val_true_vacuum << "\n";
    std::cout << "Shifted true vacuum: \n" << shifted_true_vacuum << "\n";

    std::shared_ptr<Kink_profile_guesser> kink_guesser
        = std::make_shared<Kink_profile_guesser>();
    std::shared_ptr<Profile_guesser> guesser(kink_guesser);

    // Field_profiles ansatz = guesser->get_profile_guess(
    //     potential, shifted_true_vacuum, 2, -1.0, -1.0, 1.e-3, 1.);
   
    Field_profiles ansatz = guesser->get_profile_guess(
        potential, shifted_true_vacuum, 2, 1e-8, 1.3, 1.e-3, 1.);

    double alpha = kink_guesser->get_alpha();
    double domain_start = ansatz.get_domain_start();
    double domain_end = ansatz.get_domain_end();
    std::cout << "Alpha: " << alpha << '\n'
                << "domain start: " << domain_start << '\n'
                << "domain end " << domain_end << '\n';

    RK4_perturbative_profiler profiler;
    profiler.set_true_vacuum_loc(shifted_true_vacuum);
    profiler.set_false_vacuum_loc(origin);
    profiler.set_initial_guesser(guesser);

    auto convergence_tester = std::make_shared<Relative_convergence_tester>();
    convergence_tester->set_max_iterations(30);
    profiler.set_convergence_tester(convergence_tester);

    profiler.calculate_bubble_profile(potential);

    double action = profiler.get_euclidean_action();

    std::cout << "Action: " << action << '\n';
}

}

int main() {
    // defaults:
    // m1=120.,m2=50.,mu=25.,Y1=.1,Y2=.15,n=30
    double m1 = 100.;
    double m2 = 15.;
    double mu = 5.;
    double Y1 = .1;
    double Y2 = .05;
    
    double T = 143.25583145081288;

    // double hint_false_x = 0;
    // double hint_false_y = 0;
    // double hint_true_x = 80.3640942;
    // double hint_true_y = -20.10160559;

    // double bbox_x = 1;
    // double bbox_y = 1;

    Eigen::VectorXd true_vacuum(2);
    Eigen::VectorXd false_vacuum(2);

    true_vacuum << 80.3640942, -20.10160559;
    false_vacuum << 1.39969802e-05, -7.90122376e-06;

    BubbleProfiler::Simple2DPotential potential(m1, m2, mu, Y1, Y2);
    potential.set_temperature(T);
    potential.init_effective_potential();

    // BubbleProfiler::temperature_plot(potential);
    BubbleProfiler::bounce_action(potential, true_vacuum, false_vacuum);
    return 0;
}