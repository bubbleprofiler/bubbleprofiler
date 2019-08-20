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

#ifndef BUBBLEPROFILER_OBSERVERS_HPP_INCLUDED
#define BUBBLEPROFILER_OBSERVERS_HPP_INCLUDED

#include "basic_logger.hpp"

#include <fstream>
#include <string>
#include <vector>

namespace BubbleProfiler {

class Field_profiles;
class Potential;

/*!
 * @class Dummy_observer
 * @brief A default observer that does not perform any actions when called
 *
 * This observer is used as a default if an observer is not explicitly
 * provided to routines taking such an option. Required to avoid code
 * duplication with  Perturbative_profiler::calculate_bubble_profile
 * overloads.
 */
struct Dummy_observer {
   void operator()(const Field_profiles&, const Field_profiles&) const {}
};

/*!
 * @class Profile_observer
 * @brief Observes the profile during iteration and writes out profile,
 *        perturbations and action for each step to file.
 */
class Profile_observer {
public:
   Profile_observer(const std::vector<std::string>& fields_,
                    const std::string& output_path_,
                    const Potential& potential_);
   Profile_observer(const Profile_observer&) = delete;
   Profile_observer& operator=(const Profile_observer&) = delete;

   void operator()(const Field_profiles& profile,
                   const Field_profiles& perturbation);
private:
   int iteration_count{0};
   std::vector<std::string> fields{};
   std::string output_path{};
   const Potential& potential;
   std::ofstream action_file;
   std::ofstream profiles_file;
   std::ofstream perturbations_file;
   logging::Basic_logger logger{};

   void write_initial_action_to_file(const Field_profiles& profiles);
   void write_initial_profiles_to_file(const Field_profiles& profiles);
   void write_initial_perturbations_to_file();

   void write_field_profiles_to_file(std::ofstream& stream,
                                     const Field_profiles& profiles);
   void write_action_to_file(const Field_profiles& profiles);
};

/*!
 * @class Plotting_observer
 * @brief Displays plots of field profiles at each step in the iteration.
 */
class Plotting_observer {
public:

   Plotting_observer(const std::vector<std::string>& fields_) : fields(fields_) {};

   Plotting_observer(const Plotting_observer&) = delete;
   Plotting_observer& operator=(const Plotting_observer&) = delete;

   void operator()(const Field_profiles& profile,
                   const Field_profiles& perturbation);

   private:
      int iteration_count{0};
      std::vector<std::string> fields{};
};

} // namespace BubbleProfiler

#endif
