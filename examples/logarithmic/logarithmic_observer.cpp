/*
 * This file is part of BubbleProfiler.
 *
 * BubbleProfiler is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BubbleProfiler is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BubbleProfiler.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "logarithmic_observer.hpp"
#include "error.hpp"
#include "field_profiles.hpp"

#include <iomanip>

namespace BubbleProfiler {

Logarithmic_observer::Logarithmic_observer(
   const Logarithmic_potential& potential_,
   const std::string& output_file_)
   : potential(potential_)
   , output_file(output_file_)
{
   output_stream.open(output_file);
}

void Logarithmic_observer::write_file_header()
{
   if (!output_stream.is_open() || !output_stream.good()) {
      throw IO_error("unable to write to file '" + output_file + "'");
   }
   output_stream << "# "
                 << std::setw(12) << "iteration" << ' '
                 << std::setw(16) << "rho" << ' '
                 << std::setw(16) << "perturbation" << ' '
                 << std::setw(16) << "profile" << ' '
                 << std::setw(16) << "exact"
                 << std::endl;
}

void Logarithmic_observer::write_data(const Field_profiles& profile,
                                      const Field_profiles& perturbation,
                                      const Field_profiles& exact)
{
   if (!output_stream.is_open() || !output_stream.good()) {
      throw IO_error("unable to write to file '" + output_file + "'");
   }

   const auto coord_values = profile.get_spatial_grid();
   const auto numerical_values = profile.get_field_profiles();
   const auto perturbation_values = perturbation.get_field_profiles();
   const auto exact_values = exact.get_field_profiles();

   const int n_grid_points = coord_values.size();
   for (int i = 0; i < n_grid_points; ++i) {
      output_stream << "  "
                    << std::setw(12) << iteration_count << ' '
                    << std::setw(16) << std::setprecision(8)
                    << std::scientific << coord_values(i) << ' '
                    << std::setw(16) << std::setprecision(8)
                    << std::scientific << perturbation_values(i, 0) << ' '
                    << std::setw(16) << std::setprecision(8)
                    << std::scientific << numerical_values(i, 0) << ' '
                    << std::setw(16) << std::setprecision(8)
                    << std::scientific << exact_values(i, 0) << std::endl;
   }
}

void Logarithmic_observer::operator()(const Field_profiles& profile,
                                      const Field_profiles& perturbation)
{
   if (iteration_count == 0) {
      write_file_header();
   }

   const auto exact = potential.get_profile(profile.get_spatial_grid());
   write_data(profile, perturbation, exact);
   iteration_count++;
}

} // namespace BubbleProfiler
