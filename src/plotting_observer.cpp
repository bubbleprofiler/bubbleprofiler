#include "observers.hpp"
#include "field_profiles.hpp"
#include "gnuplot-iostream.h"

#include <stdio.h>
#include <iostream>
#include <Eigen/Core>

namespace BubbleProfiler {

typedef std::vector<std::vector<double> > gp_profile_data;

Plotting_observer::~Plotting_observer() {
    if (last_only) {
        plot_profile_data(profiles.back(), "final profiles");
    }
    else {
        for (unsigned int i = 0; i < profiles.size(); ++i) {
            std::ostringstream title;
            title << "profiles for iteration " << i;
            plot_profile_data(profiles[i], title.str());
        }
    }
}

void Plotting_observer::plot_profile_data(gp_profile_data profile_data, std::string plot_title) {
    std::ostringstream gp_command;
    gp_command << "plot '-' using 1:2 with lines title '" << fields[0] << "'";

    for (int i = 2; i <= n_fields; ++i) {
        gp_command << ", '' using 1:" <<  i + 1 << " with lines title '" << fields[i - 1] << "'";
    }
    gp_command << "\n";

    Gnuplot gp;
    gp << "set title '" << plot_title << "'\n";
    gp << gp_command.str();

    for (int i = 0; i < n_fields; ++i) {
        gp.send1d(profile_data);
    }
}

void Plotting_observer::operator()(
    const Field_profiles& profile, const Field_profiles& perturbation) {
    
    const Eigen::MatrixXd field_profiles = profile.get_field_profiles();
    const Eigen::VectorXd spatial_grid = profile.get_spatial_grid();

    const int n_grid_points = spatial_grid.size();
    const int n_fields = profile.get_number_of_fields();

    // Append radial coordinates to field profiles before plotting
    gp_profile_data profile_data;
    for (int i = 0; i < n_grid_points; ++i) {
        std::vector<double> row;
        row.push_back(spatial_grid(i));
        for (int j = 0; j < n_fields; ++j) {
            row.push_back(field_profiles(i, j));
        }
        profile_data.push_back(row);
    }

    profiles.push_back(profile_data);

    iteration_count++;
}

}