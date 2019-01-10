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

#include "shooting_profile_guesser.hpp"
#include "error.hpp"
#include "math_wrappers.hpp"
#include "potential.hpp"
#include "rotation.hpp"

#include <sstream>
#include <vector>

namespace BubbleProfiler {

namespace {

/*!
 * @brief Find a strict lower bound in a sorted grid
 *
 * This function finds the index of the element in the given
 * sorted grid that is strictly less than the requested value.
 * If no such element exists, -1 is returned.
 *
 * @param grid a sorted vector of values
 * @param value the value to locate a lower bound for
 * @return the index of the largest element strictly less than
 *         the given value, or -1 if no such element exists
 */
int find_strict_lower_bound(const Eigen::VectorXd& grid, double value)
{
   if (value < grid(0)) {
      return -1;
   } else if (value > grid(grid.size() - 1)) {
      return grid.size() - 1;
   }

   int closest = 0;
   grid.unaryExpr([value](double x) {
         return Abs(x - value); }).minCoeff(&closest);

   return grid(closest) < value ? closest : closest - 1;
}

/*!
 * @brief Find a strict upper bound in a sorted grid
 *
 * This function finds the index of the element in the given
 * sorted grid that is strictly greater than the requested value.
 * If no such element exists, -1 is returned.
 *
 * @param grid a sorted vector of values
 * @param value the value to locate an upper bound for
 * @return the index of the smallest element strictly greater than
 *         the given value, or -1 if no such element exists
 */
int find_strict_upper_bound(const Eigen::VectorXd& grid, double value)
{
   if (value < grid(0)) {
      return 0;
   } else if (value > grid(grid.size() - 1)) {
      return -1;
   }

   int closest = 0;
   grid.unaryExpr([value](double x) {
         return Abs(x - value); }).minCoeff(&closest);

   return grid(closest) > value ? closest : closest + 1;
}

Eigen::VectorXd extend_domain(const Eigen::VectorXd& grid, double new_start,
                              double new_end, int new_length)
{
   const int old_length = grid.size();

   if (old_length >= new_length) {
      throw Setup_error("old grid must be shorter than new grid");
   }

   const double old_start = grid(0);
   const double old_end = grid(old_length - 1);

   if (new_start > old_start || new_end < old_end) {
      throw Setup_error("new domain boundaries cannot be within old grid");
   }

   const double extra_length = old_start - new_start + new_end - old_end;
   const double ratio = (old_start - new_start) / extra_length;
   const int n_extra = new_length - old_length;
   const int n_front = static_cast<int>(std::floor(ratio * n_extra));
   const int n_back = n_extra - n_front;

   Eigen::VectorXd new_grid(new_length);
   new_grid.segment(n_front, old_length) = grid;

   if (n_front > 0) {
      const double front_step = (old_start - new_start) / n_front;
      for (int i = 0; i < n_front; ++i) {
         new_grid(i) = new_start + i * front_step;
      }
   }

   if (n_back > 0) {
      const double back_step = (new_end - old_end) / n_back;
      for (int i = 0; i < n_back; ++i) {
         new_grid(n_front + old_length + i) = old_end + i * back_step;
      }
   }

   return new_grid;
}

template <typename InputIt>
auto get_subset(int num, InputIt first, InputIt last)
   -> std::vector<typename std::iterator_traits<InputIt>::value_type>
{
   using value_type = typename std::iterator_traits<InputIt>::value_type;

   if (num == 1) {
      const auto halfway = std::distance(first, last) / 2;
      return std::vector<value_type>{ *(first + halfway) };
   } else if (num == 2) {
      return std::vector<value_type>{ *first, *(last - 1) };
   }

   const auto length = std::distance(first, last);
   const auto midpoint = length / 2;

   std::vector<value_type> subset;
   subset.reserve(num);

   if (num % 2 == 0) {
      const auto next_length = num / 2;
      const auto front
         = get_subset(next_length, first, first + midpoint);
      const auto back
         = get_subset(next_length, first + midpoint, last);
      subset.insert(std::end(subset), std::begin(front), std::end(front));
      subset.insert(std::end(subset), std::begin(back), std::end(back));
   } else {
      const auto next_length = (num - 1) / 2;
      const auto front
         = get_subset(next_length, first, first + midpoint);
      const auto back
         = get_subset(next_length, first + midpoint + 1, last);

      subset.insert(std::end(subset), std::begin(front), std::end(front));
      subset.push_back(*(first + midpoint));
      subset.insert(std::end(subset), std::begin(back), std::end(back));
   }

   return subset;
}

Eigen::VectorXd injective_coarsening(const Eigen::VectorXd& grid,
                                     int new_length)
{
   const int old_length = grid.size();

   Eigen::VectorXd new_grid(new_length);
   if (old_length % new_length == 0) {
      const int stride = old_length / new_length;
      for (int i = 0, j = 0; i < new_length; ++i, j += stride) {
         new_grid(i) = grid(j);
      }
   } else {
      std::vector<int> grid_indices(old_length);
      std::iota(std::begin(grid_indices), std::end(grid_indices), 0);
      const auto subset_indices = get_subset(new_length,
                                             std::begin(grid_indices),
                                             std::end(grid_indices));
      for (int i = 0; i < new_length; ++i) {
         new_grid(i) = grid(subset_indices[i]);
      }
   }

   return new_grid;
}

Eigen::VectorXd coarsen_grid(const Eigen::VectorXd& grid, int new_length)
{
   const int old_length = grid.size();

   if (old_length < new_length) {
      throw Setup_error("old grid must be finer than new grid");
   }

   return injective_coarsening(grid, new_length);
}

} // anonymous namespace

Field_profiles Shooting_profile_guesser::get_profile_guess(
   const Potential& potential, const Eigen::VectorXd& true_vacuum, int d,
   double domain_start, double domain_end, double initial_step_size,
   double interpolation_points_fraction)
{
   compute_vacuum_distance(potential, true_vacuum);
   calculate_potential_parameters(potential, true_vacuum);
   return calculate_field_profiles(d, domain_start, domain_end,
                                   initial_step_size,
                                   interpolation_points_fraction);
}

void Shooting_profile_guesser::compute_vacuum_distance(
   const Potential& potential, const Eigen::VectorXd& true_vacuum)
{
   dist_true_vacuum = true_vacuum.norm();

   if (dist_true_vacuum == 0) {
      throw Setup_error("Shooting_profile_guesser::get_profile_guess: "
                        "True and false vacua are coincident");
   }

   num_fields = potential.get_number_of_fields();
   Eigen::VectorXd origin(Eigen::VectorXd::Zero(num_fields));
   if (Abs(potential(true_vacuum) - potential(origin)) < 1.e-12) {
      throw Setup_error("Shooting_profile_guesser::get_profile_guess: "
                        "True and false vacua are degenerate");
   }
}

void Shooting_profile_guesser::calculate_potential_parameters(
   const Potential& potential, const Eigen::VectorXd& true_vacuum)
{
   // Make a copy of the potential for use in ansatz construction
   auto ansatz_potential = std::unique_ptr<Potential>(potential.clone());

   // Get change of basis matrix that orients first component
   // in direction of true vacuum
   cob_matrix = calculate_rotation_to_target(true_vacuum).transpose();

   std::stringstream log_str;
   log_str << "COB: " << cob_matrix;
   logger.log_message(logging::Log_level::Trace, log_str.str());

   // Apply the change of basis matrix. Note that we add
   // the distance to the true vacuum as a scale factor so that the
   // true vacuum is located at (1,0,0,...,0).
   ansatz_potential->apply_basis_change(dist_true_vacuum * cob_matrix);

   // Add a constant term to enforce V(false_vacuum) = 0
   ansatz_potential->add_constant_term(-1*ansatz_potential->operator()(
         Eigen::VectorXd::Zero(num_fields)));

   // New coordinates of the true vacuum
   Eigen::VectorXd _true_vacuum = Eigen::VectorXd::Zero(num_fields);
   _true_vacuum(0) = 1;

   // Calculate ansatz parameters
   logger.log_message(logging::Log_level::Trace, "d: "
                      + std::to_string(trial_dist));

   Eigen::VectorXd trial_vec = Eigen::VectorXd::Zero(num_fields);
   trial_vec(0) = trial_dist;

   double vtrial = ansatz_potential->operator()(trial_vec);
   double vmin = ansatz_potential->operator()(_true_vacuum);
   double r = vmin / vtrial;

   E = 2. * ((2. - trial_dist * trial_dist) * vmin
              - vtrial / (trial_dist * trial_dist))
      / ((trial_dist - 1.) * (trial_dist - 1.));
   alpha = (r * trial_dist * trial_dist * (trial_dist - 1.5) + 0.5)
      / (1. - r * trial_dist * trial_dist * (2. - trial_dist * trial_dist));

   logger.log_message(logging::Log_level::Trace, "Alpha: "
                      + std::to_string(alpha));
   logger.log_message(logging::Log_level::Trace, "|E|"
                      + std::to_string(Abs(E)));
}

void Shooting_profile_guesser::solve_one_dimensional(int d)
{
   const auto potential = [this](double phi) {
      return this->E * (
         -this->alpha * phi * phi * phi * phi
         + phi * phi * phi
         + 0.5 * (4. * this->alpha - 3.) * phi * phi);
   };

   const auto potential_first = [this](double phi) {
      return this->E * (-4. * this->alpha * phi * phi * phi
                        + 3. * phi * phi
                        + (4. * this->alpha - 3.) * phi);
   };

   const auto potential_second = [this](double phi) {
      return this->E * (-12. * this->alpha * phi * phi
                        + 6. * phi + (4. * this->alpha - 3.));
   };

   const double false_min = calculate_false_min_location();
   const double barrier = calculate_barrier_location();
   const double true_min = calculate_true_min_location();

   solver.solve(potential, potential_first, potential_second,
                false_min, true_min, barrier, d,
                Shooting::Solver_options::Compute_profile);
}

Field_profiles Shooting_profile_guesser::calculate_field_profiles(
   int d, double domain_start, double domain_end,
   double initial_step_size, double interpolation_points_fraction)
{
   solve_one_dimensional(d);

   Field_profiles initial_profile = solver.get_bubble_profile();
   const auto n_solution_points = initial_profile.get_number_of_grid_points();
   const double solution_start = initial_profile.get_domain_start();
   const double solution_end = initial_profile.get_domain_end();

   int num_grid_points = 0;
   if (domain_start < 0. && domain_end < 0.) {
      num_grid_points = 1 + static_cast<int>(
         std::ceil((solution_end - solution_start) / initial_step_size));
      if (num_grid_points == n_solution_points) {
         initial_profile.set_number_of_dimensions(d);
         initial_profile.set_interpolation_points_fraction(
            interpolation_points_fraction);
         return initial_profile;
      }
   } else {
      num_grid_points = 1 + static_cast<int>(
         std::ceil((domain_end - domain_start) / initial_step_size));
   }

   const auto initial_grid = initial_profile.get_spatial_grid();
   const Eigen::VectorXd rho
                      = update_spatial_grid(initial_grid, domain_start,
                                            domain_end, num_grid_points);

   Eigen::VectorXd ansatz(num_grid_points);
   for (int i = 0; i < num_grid_points; ++i) {
      const double r = rho(i);
      if (r < solution_start) {
         ansatz(i) = solver.quadratic_potential_solution(r);
      } else if (r >= solution_start && r <= solution_end) {
         ansatz(i) = initial_profile.evaluate_at(0, r);
      } else {
         ansatz(i) = get_large_distance_solution(r, initial_profile, d);
      }
   }

   // We still need to undo the change of basis
   Eigen::MatrixXd cob_reverse = cob_matrix.transpose();

   // Note - we're building all the profiles here, then
   // passing them over one at a time. This is inefficient,
   // but (I think) the path of least resistance right now.
   // Could be improved later.
   Eigen::MatrixXd m_profiles(num_grid_points, num_fields);

   Eigen::VectorXd temp_field_vec = Eigen::VectorXd::Zero(num_fields);

   for (int x = 0; x < num_grid_points; ++x) {
      // Ansatz basis field vec @ this radius
      temp_field_vec(0) = ansatz(x);

      // Transform it to the original field basis (and undo the
      // field coordinate scaling)
      Eigen::VectorXd orig_field_vec
         = cob_reverse * (dist_true_vacuum) * temp_field_vec;

      m_profiles.row(x) = orig_field_vec;
   }

   Field_profiles profiles(rho, m_profiles, interpolation_points_fraction);
   profiles.set_number_of_dimensions(d);

   return profiles;
}

Eigen::VectorXd Shooting_profile_guesser::update_spatial_grid(
   const Eigen::VectorXd& old_grid, double new_start, double new_end,
   int new_length) const
{
   const int old_length = old_grid.size();
   const double old_start = old_grid(0);
   const double old_end = old_grid(old_length - 1);

   double grid_start = new_start >= 0. ? new_start : old_start;
   double grid_end = new_end >= 0. ? new_end : old_end;

   if (grid_start > grid_end) {
      std::swap(grid_start, grid_end);
   }

   const int start_index = find_strict_upper_bound(old_grid, grid_start);
   const int end_index = find_strict_lower_bound(old_grid, grid_end);
   if (start_index < 0 || end_index < 0) {
      // either grid_start > largest value in old grid, or
      // grid_end < smallest value, and new grid lies entirely
      // outside old grid
      return Eigen::VectorXd::LinSpaced(new_length, grid_start, grid_end);
   }

   Eigen::VectorXd new_grid(new_length);
   new_grid(0) = grid_start;
   new_grid(new_length - 1) = grid_end;

   const int n_interior_points = end_index - start_index + 1;
   const Eigen::VectorXd interior_points =
                     old_grid.segment(start_index, n_interior_points);

   if (n_interior_points == new_length - 2) {
      new_grid.segment(1, new_length - 2) = interior_points;
   } else if (n_interior_points < new_length - 2) {
      new_grid.segment(1, new_length - 2) =
         extend_domain(interior_points, grid_start, grid_end, new_length - 2);
   } else {
         new_grid.segment(1, new_length - 2) =
            coarsen_grid(interior_points, new_length - 2);
   }

   return new_grid;
}

double Shooting_profile_guesser::calculate_false_min_location() const
{
   double local_min_loc;
   if (E == 0.) {
      local_min_loc = 0.;
   } else if (E > 0.) {
      if (alpha <= 0.375) {
         local_min_loc = 1.;
      } else if (alpha > 0.375 && alpha <= 0.75) {
         local_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else {
         local_min_loc = 0.;
      }
   } else {
      if (alpha <= 0.25) {
         local_min_loc = 0.;
      } else if (alpha > 0.25 && alpha <= 0.375) {
         local_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else if (alpha > 0.375 && alpha <= 0.5) {
         local_min_loc = 1.;
      } else if (alpha > 0.5 && alpha <= 0.75) {
         local_min_loc = 0.;
      } else {
         local_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      }
   }
   return local_min_loc;
}

double Shooting_profile_guesser::calculate_barrier_location() const
{
  double local_max_loc;
   if (E == 0.) {
      local_max_loc = 0.;
   } else if (E > 0.) {
      if (alpha <= 0.25) {
         local_max_loc = 0.;
      } else if (alpha > 0.25 && alpha <= 0.375) {
         local_max_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else if (alpha > 0.375 && alpha <= 0.5) {
         local_max_loc = 1.;
      } else if (alpha > 0.5 && alpha <= 0.75) {
         local_max_loc = 0.;
      } else {
         local_max_loc = 0.25 * (3. - 4. * alpha) / alpha;
      }
   } else {
      if (alpha <= 0.375) {
         // potential has one or two local maxima, or no barrier,
         // returns the smaller barrier if multiple local maxima
         local_max_loc = 1.;
      } else if (alpha > 0.375 && alpha <= 0.75) {
         local_max_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else {
         local_max_loc = 0.;
      }
   }
   return local_max_loc;
}

double Shooting_profile_guesser::calculate_true_min_location() const
{
   double global_min_loc;
   if (E == 0.) {
      global_min_loc = 0.;
   } else if (E > 0.) {
      if (alpha < 0.) {
         global_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else {
         // potential is unbounded from below
         throw Error("Shooting_profile_guesser::get_global_minimum_location: "
                     "global minimum does not exist");
      }
   } else {
      if (alpha <= 0.) {
         // potential is unbounded from below
         throw Error("Shooting_profile_guesser::get_global_minimum_location: "
                     "global minimum does not exist");
      } else if (alpha > 0. && alpha <= 0.25) {
         global_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else if (alpha > 0.25 && alpha <= 0.5) {
         global_min_loc = 0.;
      } else {
         global_min_loc = 1.;
      }
   }
   return global_min_loc;
}

double Shooting_profile_guesser::get_large_distance_solution(
   double rho, const Field_profiles& profile, int d) const
{
   const auto potential_second = [this](double phi) {
      return this->E * (-12. * this->alpha * phi * phi
                        + 6. * phi + (4. * this->alpha - 3.));
   };

   const double local_min = calculate_false_min_location();
   const double mass_squared = potential_second(local_min);

   const double rho_max = profile.get_domain_end();
   const double y = rho * Sqrt(Abs(mass_squared));
   const double y_max = rho_max * Sqrt(Abs(mass_squared));
   const double phip_max = profile.derivative_at(0, 1, rho_max);
   double correction = 0;
   if (mass_squared > 0) {
      if (d == 3) {
         correction = -phip_max * rho_max * rho_max * Exp(y_max - y) /
            (rho * (1. + y_max));
      } else {
         correction = -phip_max * rho_max * rho_max * BesselK(1, y) /
            (rho * (BesselK(1, y_max) + 0.5 * y_max * (
                       BesselK(0, y_max)
                       + BesselK(2, y_max))));
      }
   } else {
      throw Error("Shooting_profile_guesser::get_large_distance_solution: "
                  "cannot evaluate asymptotic solution for this potential");
   }

   return local_min + correction;
}

} // namespace BubbleProfiler
