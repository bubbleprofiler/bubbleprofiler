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

#ifndef BUBBLEPROFILER_GENERIC_PERTURBATIVE_PROFILER_HPP_INCLUDED
#define BUBBLEPROFILER_GENERIC_PERTURBATIVE_PROFILER_HPP_INCLUDED

/**
 * @file generic_perturbative_profiler.hpp
 * @brief contains the implementation of the Perturbative_profiler class
 */

#include "basic_logger.hpp"
#include "default_integration_policy.hpp"
#include "euclidean_action.hpp"
#include "error.hpp"
#include "field_profiles.hpp"
#include "gsl_interpolator.hpp"
#include "gsl_root_finder.hpp"
#include "kink_profile_guesser.hpp"
#include "observers.hpp"
#include "perturbations_ode_system.hpp"
#include "potential.hpp"
#include "relative_convergence_tester.hpp"

#include <Eigen/Core>

#include <algorithm>
#include <iterator>
#include <memory>
#include <numeric>
#include <sstream>
#include <tuple>

namespace BubbleProfiler {

/**
 * @class Perturbative_profiler
 * @brief Bounce solver using perturbative method
 *
 * This is the main interface class for the BubbleProfiler. The workflow
 * for using this class is as follows:
 * 1. Specify the domain size and discretization parameters via the set_
 *   methods.
 * 2. Supply a Profile_guesser class via set_initial_guesser to generate
 *   an initial ansatz solution.
 * 3. Supply a Profile_convergence_tester via set_convergence_tester to
 *   determine the stopping critera.
 * 4. Call calculate_bubble_profile to perform the calculation. Callers
 *   of this method supply:
 *   - A scalar field potential (Potential).
 *   - An algorithm for solving the ODE which determines the perturbative
 *     corrections (BVP_solver).
 *   - An Observer to monitor the output.
 * 5. Retrieve the resulting Field_profiles via get_bubble_profile.
 */
template <class Integration_policy>
class Perturbative_profiler : public Integration_policy {
public:

   //! @param start of the domain (radial coordinate)
   void set_domain_start(double start);

   //! @param end of the domain (radial coordinate)
   void set_domain_end(double end);

   //! @param n number of grid points per dimension.
   void set_initial_step_size(double h);

   //! @param f fraction of grid points to use in interpolating fields
   void set_interpolation_points_fraction(double f) {
      interpolation_points_fraction = f;
   }

   //! @param d number of dimensions (scalar fields)
   void set_number_of_dimensions(int d) { n_dims = d; }

   //! @param g Profile_guesser class to generate an initial ansatz solution
   void set_initial_guesser(std::shared_ptr<Profile_guesser> g) {
      guesser = g;
   }

   void set_root_finder(std::shared_ptr<Root_finder<Eigen::VectorXd> > rf) {
      root_finder = rf; }

   //! @param ct Profile_convergence_tester to determine stopping criterea
   void set_convergence_tester(std::shared_ptr<Profile_convergence_tester> ct) {
      convergence_tester = ct;
   }

   //! @param tv location of true vacuum in field space
   void set_true_vacuum_loc(const Eigen::VectorXd& tv) { true_vacuum = tv; }

   //! @param fv location of false vacuum in field space
   void set_false_vacuum_loc(const Eigen::VectorXd& fv) { false_vacuum = fv; }

   /*!
    * Calculate the bubble profile for the specified potential.
    * @tparam BVP_solver algorithm to solve BVP for perturbative corrections
    * @tparam Observer to report output
    */
   template <class Observer>
   void calculate_bubble_profile(Potential&, Observer&);

   /*!
    * @brief Calculate the bubble profile for the specified potential
    * @tparam BVP_solver algorithm to solve BVP for perturbative corrections
    */
   void calculate_bubble_profile(Potential&);

   //! Retrieve the field profiles after running calculate_bubble_profile
   const Field_profiles& get_bubble_profile() const { return profiles; }

   //! Retrieve the action after running calculate_bubble_profile
   double get_euclidean_action() const { return euclidean_action; }

private:
   double domain_start{-1.};
   double domain_end{-1.};
   double initial_step_size{1.e-2};
   int n_grid_points{1000};
   double interpolation_points_fraction{1.};
   int n_dims{3};
   int iteration{0};
   Field_profiles profiles{};
   double euclidean_action{0.};
   std::shared_ptr<Profile_guesser> guesser{
      std::make_shared<Kink_profile_guesser>()};
   std::shared_ptr<Profile_convergence_tester> convergence_tester{
      std::make_shared<Relative_convergence_tester>()};
   std::shared_ptr<Root_finder<Eigen::VectorXd> > root_finder{
      std::make_shared<GSL_root_finder<Eigen::Dynamic> >()};
   logging::Basic_logger logger{};

   Eigen::VectorXd true_vacuum{};
   Eigen::VectorXd false_vacuum{};

   void check_setup(Potential&) const;
   Field_profiles construct_initial_profiles(Potential&) const;
   int get_max_iterations() const;
   bool check_convergence(const Potential&) const;

   void integrate_ivp(const Perturbations_ODE_system&,
                      Eigen::VectorXd&, double) const;
   template <class Observer>
   void integrate_ivp(const Perturbations_ODE_system&,
                      Eigen::VectorXd&, double, Observer) const;
   std::tuple<Eigen::VectorXd,Eigen::MatrixXd> calculate_perturbation(
      Potential&);
   void integrate_system_to(const Perturbations_ODE_system&, double rho) const;

   void update_field_profiles(const Eigen::VectorXd&, const Eigen::MatrixXd&);
};

template <class Integration_policy>
void Perturbative_profiler<Integration_policy>::calculate_bubble_profile(
   Potential& potential)
{
   Dummy_observer observer;
   calculate_bubble_profile(potential, observer);
}

template <class Integration_policy>
template <class Observer>
void Perturbative_profiler<Integration_policy>::calculate_bubble_profile(
   Potential& potential, Observer& observer)
{
   check_setup(potential);

   // Set up the solution field profiles and initialize BVP
   profiles = construct_initial_profiles(potential);

   // Update domain boundaries
   domain_start = profiles.get_domain_start();
   domain_end = profiles.get_domain_end();

   // Compute required number of grid points to achieve the desired
   // initial step size
   n_grid_points = 1 + static_cast<int>(
      std::ceil((domain_end - domain_start) / initial_step_size));

   const auto n_fields = profiles.get_number_of_fields();
   Field_profiles perturbations(
      Eigen::MatrixXd::Zero(n_grid_points, n_fields), domain_start,
      domain_end, interpolation_points_fraction);

   // Perturbation loop - keep going until convergence criteria
   // are satisfied.
   observer(profiles, perturbations);
   const int max_iterations = get_max_iterations();

   iteration = 0;
   bool converged = false;
   while (iteration < max_iterations && !converged) {
      logger.log_message(logging::Log_level::Trace, "Beginning perturbation "
                         + std::to_string(iteration + 1));

      const auto result = calculate_perturbation(potential);
      perturbations = Field_profiles(std::get<0>(result), std::get<1>(result),
                                     interpolation_points_fraction);

      update_field_profiles(perturbations.get_spatial_grid(),
                            perturbations.get_field_profiles());

      observer(profiles, perturbations);

      converged = check_convergence(potential);

      iteration++;
   }

   if (!converged) {
      throw No_convergence_error(iteration);
   }

   Eigen::MatrixXd original_basis(profiles.get_field_profiles());
   original_basis.rowwise() += false_vacuum.transpose();

   profiles = Field_profiles(profiles.get_spatial_grid(), original_basis,
                             interpolation_points_fraction);
   profiles.set_number_of_dimensions(n_dims);

   potential.translate_origin(-false_vacuum);

   euclidean_action = calculate_action(potential, profiles);
}

template <class Integration_policy>
void Perturbative_profiler<Integration_policy>::set_domain_start(double start)
{
   if (start < 0.) {
      throw Domain_error("radial coordinates must be non-negative");
   }
   domain_start = start;
}

template <class Integration_policy>
void Perturbative_profiler<Integration_policy>::set_domain_end(double end)
{
   if (end < 0.) {
      throw Domain_error("radial coordinates must be non-negative");
   }
   domain_end = end;
}

template <class Integration_policy>
void Perturbative_profiler<Integration_policy>::set_initial_step_size(double h)
{
   if (h <= 0.) {
      throw Setup_error("integration step size must be positive");
   }
   initial_step_size = h;
}

template <class Integration_policy>
void Perturbative_profiler<Integration_policy>::check_setup(
   Potential& potential) const
{
   if (!convergence_tester) {
      throw Setup_error("Perturbative_profiler::check_setup: "
                        "no convergence tester provided");
   }

   if (!guesser) {
      throw Setup_error("Perturbative_profiler::check_setup: "
                        "no profile guesser provided");
   }

   if (!root_finder) {
      throw Setup_error("Perturbative_profiler::check_setup: "
                        "no root finder provided");
   }

   // Ensure false vacuum is located at the origin in field space
   potential.translate_origin(false_vacuum);
}

template <class Integration_policy>
Field_profiles Perturbative_profiler<Integration_policy>::construct_initial_profiles(
   Potential& potential) const
{
   // We will first make a shift of fields such that the false vacuum is at the origin;
   // this will remain the case for the rest of the calculation. The kink_profile_guesser
   // uses a different basis, but does not affect the main potential.
   Eigen::VectorXd true_vacuum_ = true_vacuum - false_vacuum;

   return guesser->get_profile_guess(potential, true_vacuum_, n_dims,
                                     domain_start, domain_end, initial_step_size,
                                     interpolation_points_fraction);
}

template <class Integration_policy>
int Perturbative_profiler<Integration_policy>::get_max_iterations() const
{
   return convergence_tester->get_max_iterations();
}

template <class Integration_policy>
bool Perturbative_profiler<Integration_policy>::check_convergence(
   const Potential& potential) const
{
   return convergence_tester->is_converged(potential, profiles);
}

template <class Integration_policy>
std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
Perturbative_profiler<Integration_policy>::calculate_perturbation(
   Potential& potential)
{
   const auto n_fields = profiles.get_number_of_fields();

   Perturbations_ODE_system system(potential, profiles, n_dims);

   const Eigen::VectorXd domain_end_values(profiles.evaluate_at(domain_end));
   Eigen::VectorXd domain_start_derivs(n_fields);
   for (auto i = 0; i < n_fields; ++i) {
      domain_start_derivs(i) = profiles.derivative_at(i, 1, domain_start);
   }

   const auto boundary_conditions =
      [n_fields, &domain_start_derivs, &domain_end_values]
      (const Eigen::VectorXd& start_state, const Eigen::VectorXd& end_state) {
      Eigen::VectorXd values(2 * n_fields);

      values.segment(0, n_fields)
      = end_state.segment(0, n_fields) + domain_end_values;
      values.segment(n_fields, n_fields)
      = start_state.segment(n_fields, n_fields) + domain_start_derivs;

      return values;
   };

   const double step_size = (domain_end - domain_start)
      / (n_grid_points - 1.);

   const auto calculate_errors
      = [this, &system, &boundary_conditions, step_size](
         const Eigen::VectorXd& start_state) {
      Eigen::VectorXd end_state(start_state);

      this->integrate_ivp(system, end_state, step_size);

      return boundary_conditions(start_state, end_state);
   };

   Eigen::VectorXd initial_guess(Eigen::VectorXd::Zero(2 * n_fields));
   initial_guess.segment(n_fields, n_fields) -= domain_start_derivs;

   const auto status = root_finder->find_root(calculate_errors, initial_guess);
   if (status != Root_finder_status::SUCCESS) {
      throw BVP_solver_error("Perturbative_profiler::calculate_perturbation: "
                             "unable to solve boundary value problem");
   }

   const Eigen::VectorXd domain_start_perturbations(
      root_finder->get_solution());
   Eigen::VectorXd domain_end_perturbations(domain_start_perturbations);

   std::vector<double> coord_values;
   std::vector<Eigen::VectorXd> perturbation_values;

   auto observer = [n_fields, &coord_values, &perturbation_values](
      const Eigen::VectorXd& state, double rho) {
      Eigen::VectorXd solution_values(
         Eigen::VectorXd::Map(state.data(), n_fields));
      perturbation_values.push_back(solution_values);
      coord_values.push_back(rho);
   };

   integrate_ivp(system, domain_end_perturbations, step_size, observer);

   // @todo reinstate final check on solution?

   const int n_perturbation_points = coord_values.size();
   const Eigen::VectorXd grid_points(
      Eigen::VectorXd::Map(coord_values.data(), n_perturbation_points));
   Eigen::MatrixXd perturbations(n_perturbation_points, n_fields);
   for (int i = 0; i < n_perturbation_points; ++i) {
      perturbations.row(i) = perturbation_values[i];
   }

   return std::make_tuple(grid_points, perturbations);
}

template <class Integration_policy>
void Perturbative_profiler<Integration_policy>::update_field_profiles(
   const Eigen::VectorXd& new_grid, const Eigen::MatrixXd& perturbations)
{
   const Eigen::VectorXd old_grid = profiles.get_spatial_grid();
   const Eigen::MatrixXd old_profiles = profiles.get_field_profiles();

   const int n_updated_points = new_grid.size();
   const int n_fields = perturbations.cols();

   // if the new grid to be used extends beyond the limits of the
   // previous grid, extend previous field profiles using constant
   // extrapolation
   const int n_old = old_grid.size();
   const double old_domain_start = old_grid(0);
   const double old_domain_end = old_grid(n_old - 1);

   const int n_before = (new_grid.array() < old_domain_start).count();
   const int n_after = (new_grid.array() > old_domain_end).count();
   const int extended_grid_size = n_before + n_old + n_after;
   Eigen::VectorXd extended_grid(extended_grid_size);
   Eigen::MatrixXd extended_profiles(extended_grid_size, n_fields);

   for (int i = 0; i < n_before; ++i) {
      extended_grid(i) = new_grid(i);
      extended_profiles.block(0, 0, n_before, n_fields).rowwise()
         = old_profiles.row(0);
   }
   extended_grid.segment(n_before, n_old) = old_grid;
   extended_profiles.block(n_before, 0, n_old, n_fields) = old_profiles;

   if (n_after > 0) {
      int offset = new_grid.size();
      while (new_grid(offset - 1) > old_domain_end) {
         --offset;
      }

      for (int i = 0; i < n_after; ++i) {
         extended_grid(n_before + n_old + i) = new_grid(offset + i);
         extended_profiles.block(n_before + n_old, 0, n_after, n_fields).rowwise()
            = old_profiles.row(n_old - 1);
      }
   }

   Eigen::MatrixXd interpolated_profiles(n_updated_points, n_fields);
   for (int i = 0; i < n_fields; ++i) {
      interpolated_profiles.col(i) = interpolate_f_at(
         extended_grid, extended_profiles.col(i), new_grid);
   }

   const Eigen::MatrixXd new_profiles = interpolated_profiles + perturbations;

   profiles = Field_profiles(new_grid, new_profiles,
                             interpolation_points_fraction);
   profiles.set_number_of_dimensions(n_dims);

   domain_start = new_grid(0);
   domain_end = new_grid(n_updated_points - 1);
}

template <class Integration_policy>
void Perturbative_profiler<Integration_policy>::integrate_ivp(
   const Perturbations_ODE_system& system,
   Eigen::VectorXd& perturbations, double step_size) const
{
   const auto dummy_observer = [](const Eigen::VectorXd&, double) {};
   integrate_ivp(system, perturbations, step_size, dummy_observer);
}

template <class Integration_policy>
template <class Observer>
void Perturbative_profiler<Integration_policy>::integrate_ivp(
   const Perturbations_ODE_system& system,
   Eigen::VectorXd& perturbations, double step_size, Observer observer) const
{
   Integration_policy::integrate_system(system, perturbations,
                                        domain_start, domain_end, step_size,
                                        observer);
}

} // namespace BubbleProfiler

#endif
