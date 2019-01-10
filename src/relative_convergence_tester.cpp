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

/**
    @file
    @brief Judge whether perturbative solver converged
*/

#include "relative_convergence_tester.hpp"
#include "error.hpp"
#include "euclidean_action.hpp"
#include "field_profiles.hpp"
#include "math_wrappers.hpp"
#include "potential.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>

namespace BubbleProfiler {

Relative_convergence_tester::Relative_convergence_tester(
   double action_tol_, double fields_tol_)
   : action_relative_tol(action_tol_)
   , field_vals_relative_tol(fields_tol_)
{
   max_iterations = calculate_max_iterations();
}

Relative_convergence_tester::Relative_convergence_tester(double tol)
   : Relative_convergence_tester(tol, tol)
{
}

bool Relative_convergence_tester::is_converged(
   const Potential& potential, const Field_profiles& profiles)
{
   bool converged = false;
   const double current_action = calculate_action(potential, profiles);
   const double domain_start = profiles.get_domain_start();
   const Eigen::VectorXd current_field_vals = profiles.evaluate_at(domain_start);

   if (iteration_count > 0) {
      converged = check_action_converged(current_action) &&
         check_fields_converged(domain_start, current_field_vals);
   }

   old_action = current_action;
   old_field_vals = current_field_vals;
   ++iteration_count;

   if (converged) {
      logger.log_message(logging::Log_level::Trace, "Converged");
      logger.log_message(logging::Log_level::Trace, "Relative change in action < "
                         + std::to_string(action_relative_tol));
      logger.log_message(logging::Log_level::Trace, "Relative changes in fields at (r = "
                         + std::to_string(domain_start) + ") < "
                         + std::to_string(field_vals_relative_tol));
   }

   return converged;
}

int Relative_convergence_tester::calculate_max_iterations() const
{
   const auto min_tol = std::min(action_relative_tol, field_vals_relative_tol);
   return static_cast<int>(-std::log10(min_tol) * 10);
}

double Relative_convergence_tester::relative_difference(
   double x, double y) const
{
   const double ax = Abs(x);
   const double ay = Abs(y);
   const double largest = std::max(ax, ay);

   const double underflow = std::numeric_limits<double>::min();

   if (largest < underflow) {
      return 0.;
   }

   return Abs(x - y) / largest;
}

bool Relative_convergence_tester::check_action_converged(double action) const
{
   const double action_diff = relative_difference(action, old_action);

   logger.log_message(logging::Log_level::Trace, "Relative change in action = "
                      + std::to_string(action_diff));

   const bool action_converged = action_diff < action_relative_tol;

   std::stringstream log_str;
   log_str << "Action converged = " << (action_converged ? "true" : "false");
   logger.log_message(logging::Log_level::Trace, log_str.str());

   return action_converged;
}

bool Relative_convergence_tester::check_fields_converged(
   double domain_start, const Eigen::VectorXd& field_vals) const
{
   const int n_fields = field_vals.size();

   if (n_fields != old_field_vals.size()) {
      throw Setup_error(
         "Relative_convergence_tester::check_fields_converged: "
         "number of field values does not match previous number");
   }

   std::vector<double> field_diffs(n_fields);
   for (int i = 0; i < n_fields; ++i) {
      field_diffs[i] = relative_difference(field_vals(i), old_field_vals(i));
      logger.log_message(logging::Log_level::Trace, "Field_" + std::to_string(i)
                         + "(r = " + std::to_string(domain_start) + ") = "
                         + std::to_string(field_vals(i)));
   }
   const double max_diff = *std::max_element(field_diffs.cbegin(), field_diffs.cend());

   logger.log_message(logging::Log_level::Trace, "Maximum relative change in fields (at r = 0) = "
                      + std::to_string(max_diff));

   const bool fields_converged = max_diff < field_vals_relative_tol;
   std::stringstream log_str;
   log_str << "Fields converged = " << + (fields_converged ? "true" : "false");
   logger.log_message(logging::Log_level::Trace, log_str.str());
                      

   return fields_converged;
}

void Relative_convergence_tester::restart()
{
   iteration_count = 0.;
   old_action = 0.;
   old_field_vals = Eigen::VectorXd::Zero(old_field_vals.size());
}

} // namespace BubbleProfiler
