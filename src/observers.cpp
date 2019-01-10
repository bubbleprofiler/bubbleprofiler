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

#include "observers.hpp"
#include "euclidean_action.hpp"
#include "field_profiles.hpp"

#include <Eigen/Core>

#include <iomanip>
#include <iostream>

namespace BubbleProfiler {

Profile_observer::Profile_observer(const std::vector<std::string>& fields_,
                                   const std::string& output_path_,
                                   const Potential& potential_)
   : fields(fields_)
   , output_path(output_path_)
   , potential(potential_)
{
   action_file.open(output_path + "/action.txt");
   profiles_file.open(output_path + "/field_profiles.txt");
   perturbations_file.open(output_path + "/perturbations.txt");
}


void Profile_observer::write_initial_action_to_file(
   const Field_profiles& profiles)
{
   action_file << "#" << ' '
               << std::setw(12) << "perturbation" << ' '
               << std::setw(16) << "action" << std::endl;

   write_action_to_file(profiles);
}

void Profile_observer::write_initial_profiles_to_file(
   const Field_profiles& profiles)
{
   logger.log_message(logging::Log_level::Trace, "Writing ansatz to file");

   profiles_file << "#" << ' '
                 << std::setw(12) << "perturbation" << ' '
                 << std::setw(16) << "rho";
   for (const std::string& name: fields) {
      profiles_file << ' ' << std::setw(16) << name;
   }
   profiles_file << std::endl;

   write_field_profiles_to_file(profiles_file, profiles);
}

void Profile_observer::write_initial_perturbations_to_file()
{
   perturbations_file << "#" << ' '
                 << std::setw(12) << "perturbation" << ' '
                 << std::setw(16) << "rho";
   for (const std::string& name: fields) {
      perturbations_file << ' ' << std::setw(16) << "eps_" + name;
   }
   perturbations_file << std::endl;
}

void Profile_observer::write_field_profiles_to_file(
   std::ofstream& stream, const Field_profiles& profiles)
{
   const auto coord_values = profiles.get_spatial_grid();
   const auto field_values = profiles.get_field_profiles();

   const int n_grid_points = coord_values.size();
   const int n_fields = fields.size();

   for (int i = 0; i < n_grid_points; ++i) {
      stream << "  "
                    << std::setw(12) << iteration_count << ' '
                    << std::setw(16) << std::setprecision(8)
                    << std::scientific << coord_values(i);
      for (int j = 0; j < n_fields; ++j) {
         stream << ' '
                       << std::setw(16) << std::setprecision(8)
                       << std::scientific << field_values(i, j);
      }
      stream << std::endl;
   }
}

void Profile_observer::write_action_to_file(const Field_profiles& profiles)
{
   const double action = calculate_action(potential, profiles);

   logger.log_message(logging::Log_level::Trace, "Action: "
                      + std::to_string(action));

   action_file << "  "
               << std::setw(12) << iteration_count << ' '
               << std::setw(16) << std::setprecision(8)
               << std::scientific << action << std::endl;
}

void Profile_observer::operator()(
   const Field_profiles& profiles, const Field_profiles& perturbations)
{
   if (iteration_count == 0) {
      write_initial_action_to_file(profiles);
      write_initial_profiles_to_file(profiles);
      write_initial_perturbations_to_file();
   } else {
      write_field_profiles_to_file(perturbations_file, perturbations);
      write_field_profiles_to_file(profiles_file, profiles);
      write_action_to_file(profiles);
   }
   iteration_count++;
}

} // namespace BubbleProfiler
