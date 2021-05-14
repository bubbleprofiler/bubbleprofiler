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
    @brief One-dimensional shooting method
    
    Solve a one-dimensional bounce problem with the shooting method.
*/

#include "shooting.hpp"
#include "error.hpp"
#include "math_wrappers.hpp"
#include "numeric.hpp"
#include "potential.hpp"

#include <boost/math/tools/minima.hpp>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <utility>

namespace BubbleProfiler {

double Shooting::get_euclidean_action() const
{
   if (!action_computed) {
      throw Error("Shooting::get_euclidean_action: action not computed");
   }
   return euclidean_action;
}

const Field_profiles& Shooting::get_bubble_profile() const
{
   if (!profile_computed) {
      throw Error("Shooting::get_bubble_profile: bubble profile not computed");
   }
   return profiles;
}

void Shooting::initialize_potential_parameters(
   const std::function<double(double)>& potential_,
   const std::function<double(double)>& potential_first_,
   const std::function<double(double)>& potential_second_,
   double false_min_, double true_min_, double barrier_, int dim_)
{
   potential = potential_;
   potential_first = potential_first_;
   potential_second = potential_second_;

   dim = dim_;
   barrier = barrier_;
   false_min = false_min_;
   true_min = true_min_;
   sign_min = Sign(false_min_ - true_min_);

   p_barrier = potential(barrier_);
   p_false_min = potential(false_min_);
   p_true_min = potential(true_min_);

   // Find settings for solver that depend on potential
   scale = bubble_scale();
   drho_guess = drho_frac * scale;
}

void Shooting::solve(
   const std::function<double(double)>& potential_,
   const std::function<double(double)>& potential_first_,
   const std::function<double(double)>& potential_second_,
   double false_min_, double true_min_, double barrier_, int dim_,
   unsigned int options)
{
   if (dim_ != 3 && dim_ != 4) {
      throw Setup_error("spacetime dimension must be 3 or 4");
   }

   action_computed = false;
   profile_computed = false;

   bool compute_action = (options & Solver_options::Compute_action) != 0;
   bool compute_profile = (options & Solver_options::Compute_profile) != 0;

   if (!compute_action && !compute_profile) {
      return;
   }

   initialize_potential_parameters(potential_, potential_first_,
                                   potential_second_, false_min_,
                                   true_min_, barrier_, dim_);

   // Build stepper for shooting
   shoot_stepper
      = make_controlled(shoot_ode_rel, shoot_ode_abs,
                        error_stepper_type());

   // Solve initial value problem
   lambda_sol = shooting();

   if (compute_action) {
      euclidean_action = action();
      action_computed = true;
   }

   if (compute_profile) {
      profiles = calculate_bubble_profile();
      profile_computed = true;
   }
}

void Shooting::solve(const Potential& potential_, double false_min_,
                     double true_min_, double barrier_, int dim_,
                     unsigned int options)
{
   if (potential_.get_number_of_fields() != 1) {
      throw Setup_error("potential must be one dimensional");
   }

   const auto potential_func = [&potential_](double phi) {
      Eigen::VectorXd fields(1);
      fields(0) = phi;
      return potential_(fields);
   };

   const auto potential_first_ = [&potential_](double phi) {
      Eigen::VectorXd fields(1);
      fields(0) = phi;
      return potential_.partial(fields, 0);
   };

   const auto potential_second_ = [&potential_](double phi) {
      Eigen::VectorXd fields(1);
      fields(0) = phi;
      return potential_.partial(fields, 0, 0);
   };

   solve(potential_func, potential_first_, potential_second_,
         false_min_, true_min_, barrier_, dim_, options);
}

void Shooting::ode(const state_type& x, state_type& dxdt, double rho)
{
   dxdt[0] = x[1];
   dxdt[1] = -(dim - 1.) / rho * x[1] + potential_first(x[0]);
}

double Shooting::calculate_roll_time(double lambda, double dVdphi,
                                        double d2Vdphi2) const
{
   const double change = evolve_change_rel * Abs(unmap(lambda)
                                                 - false_min);
   const double prefactor = Abs(dVdphi / d2Vdphi2);

   double y = 0.;
   if (d2Vdphi2 > 0.) {
      if (dim == 3) {
         if (change / prefactor > f_y_max) {
            // Asymptotic formula
            y = Log(2. * change / Abs(barrier - true_min)) + lambda;
         } else {
            y = asinch(1. + change / prefactor);
         }
      } else if (dim == 4) {
         if (change / prefactor > f_y_max) {
            // Asymptotic formula
            y = 0.5 * Log(27. * Pi / 16.)
                + Log(change / Abs(barrier - true_min)) + lambda;
         } else {
            y = approx_root_pos_4(0.5 * change / prefactor);
         }
      }
   } else {
      if (dim == 3) {
         y = asinc(1. - change / prefactor);
      } else if (dim == 4) {
         y = approx_root_neg_4(0.5 * change / prefactor);
      }
   }

   return y;
}

double Shooting::calculate_roll_velocity(double lambda, double y,
                                         double dVdphi, double d2Vdphi2) const
{
   const double change = evolve_change_rel * Abs(unmap(lambda)
                                                 - false_min);
   const double prefactor = Abs(dVdphi / d2Vdphi2);
   const double mass = Sqrt(Abs(d2Vdphi2));
   const double rho = y / mass;

   double v = 0.;
   if (d2Vdphi2 > 0.) {
      if (dim == 3) {
         const double sinch_y = 1. + change / prefactor;
         if (sinch_y > f_y_max) {
            v = sign_min * mass * change;  // Asymptotic formula
         } else if (sinch_y - 1. < f_y_min) {
            v = 2. * sign_min * change / rho;  // Asymptotic formula
         } else {
            v = sign_min * prefactor / rho * (Cosh(y) - sinch_y);
         }
      } else if (dim == 4) {
         const double f_y = 0.5 * change / prefactor;
         if (f_y > f_y_max) {
            v = sign_min * mass * change;  // Asymptotic formula
         } else if (f_y < f_y_min) {
            v = 2. * sign_min * change / rho;  // Asymptotic formula
         } else {
            v = sign_min * prefactor / rho * 2. * BesselI(2, y);
         }
      }
   } else {
      if (dim == 3) {
         const double sinc_y = 1. - change / prefactor;
         if (1. - sinc_y < f_y_min) {
            v = 2. * sign_min * change / rho;  // Asymptotic formula
         } else {
            v = - sign_min * prefactor / rho * (Cos(y) - sinc_y);
         }
      } else if (dim == 4) {
         const double f_y = 0.5 * change / prefactor;
         if (f_y < f_y_min) {
            v = 2. * sign_min * change / rho;  // Asymptotic formula
         } else {
            v = sign_min * prefactor / rho * 2. * BesselJ(2, y);
         }
      }
   }

   return v;
}

Shooting::evolve_type Shooting::evolve(double lambda) const
{
   const double phi_zero = unmap(lambda);
   const double mass_squared = potential_second(phi_zero);
   const double prime_start = potential_first(phi_zero);

   // Solve for \rho.
   const double mass = Sqrt(Abs(mass_squared));
   const double change = evolve_change_rel * Abs(phi_zero - false_min);

   const double y = calculate_roll_time(lambda, prime_start, mass_squared);
   const double rho = y / mass;

   // Find corresponding position and velocity.
   const double x = phi_zero + sign_min * change;
   const double v = calculate_roll_velocity(lambda, y, prime_start,
                                            mass_squared);

   if (std::isnan(rho) || std::isnan(x) || std::isnan(v)) {
      throw Numerical_error(
         "Could not evolve fields with approximate solution");
   }

   return {{rho, x, v}};
}

double Shooting::bubble_scale() const
{
   const double mass_squared = potential_second(barrier);
   if(mass_squared >= 0.0 || std::isnan(mass_squared) || std::isinf(mass_squared)  ) {
      throw Numerical_error(
         "Double derivative at the barrier (a local maximum) is not a negative finite number.");
   }
   return 2. * Pi / Sqrt(-mass_squared);
}

double Shooting::energy(double phi, double dot_phi) const
{
   return 0.5 * pow(dot_phi, 2) - (potential(phi) - p_false_min);
}

double Shooting::shoot(double lambda)
{
   constexpr double over = 1.;
   constexpr double under = -1.;

   const double phi_zero = unmap(lambda);

   // Reset action integral and profiles
   action_no_prefactor = 0.;
   rho_vals.clear();
   field_vals.clear();

   // Trivial cases - don't check them explicitly as may get stuck on
   // tops of maxima.
   if (std::isinf(lambda)) {
      return over;
   } else if (phi_zero == false_min) {
      return under;
   }

   // Check there is sufficient energy to make a shot before evolving
   // with approximate analytic solution etc.
   const double x = phi_zero;
   const double v = 0.;
   const bool no_energy = energy(x, v) < 0.;

   if (no_energy) {
      return under;
   }

   // Reset stepper
   shoot_stepper.reset();

   // Make ODE function with signature required by boost
   const auto ode_ = [this](const state_type& x, state_type& dxdt, double rho) {
      return this->ode(x, dxdt, rho);
   };

   // Evolve guess of \f$\phi_0\f$ forward using approximate solution,
   // resulting in initial conditions.
   evolve_type evolved = evolve(lambda);
   double rho = evolved[0];  // \f$\hat\rho\f$
   // \f$\phi(\hat\rho)\f$ and \f$\dot\phi(\hat\rho)\f$
   state_type y {{evolved[1], evolved[2]}};
   double drho = drho_guess;
   const double rho_max = rho + periods_max * scale;
   
   // Initialize integral
   double f_a = integrand(y[1], rho);

   // Initialize profile
   rho_vals = {rho};
   field_vals = {y[0]};

   for (int iter = 0; iter < iter_max; ++iter) {
      // Check for whether we reached an unacceptable value of \f$\rho\f$
      if (rho >= rho_max ||
          Abs(drho) <= std::numeric_limits<double>::epsilon()) {
         throw Numerical_error(
            "Could not determine whether shot was an undershot or "
            "overshot in shoot");
      }

      const double x = y[0];
      const double v = y[1];
      
      // Check whether we have gone past the true vacuum. This must come first,
      // as we could go past the true vacuum and then run out of energy etc.
      const bool overshot = sign_min * (x - false_min) > 0.;

      if (overshot) {
         return over;
      }

      // Check whether there is sufficient energy to make a shot
      const bool no_energy = energy(x, v) < 0.;

      if (no_energy) {
         return under;
      }

      // Check whether we are moving in the wrong direction
      const bool undershot = sign_min * v < 0.;

      if (undershot) {
         return under;
      }

      // Note \f$\rho\f$ - about to be updated in place
      const double rho_a = rho;

      // Step forward in \f$\rho\f$ and update \f$d\rho\f$ in place.
      // The while loop isn't strictly necessary but avoids repeating the
      // above checks.
      
      while (shoot_stepper.try_step(ode_, y, rho, drho)) {
      }
      
      // Contribution to action from trapezoid rule
      const double f_b = integrand(v, rho);

      // Don't use \f$d\rho\f$ as it is updated in place
      action_no_prefactor += 0.5 * (f_a + f_b) * (rho - rho_a);

      // Cache integrand
      f_a = f_b;
      
      // Cache profile
      rho_vals.push_back(rho);
      field_vals.push_back(y[0]);
   }

   throw Numerical_error("Could not determine whether shot was an undershot or "
                         "overshot in shoot");
}

double Shooting::unmap(double lambda) const
{
   return Exp(-lambda) * (barrier - true_min) + true_min;
}

double Shooting::shooting()
{
   const auto shoot_ = [this](double l) { return this->shoot(l); };
   std::pair<double, double> solution;
   boost::math::tools::eps_tolerance<double> stop(shoot_bisect_bits);
   double lower = 0.;
   double upper = bisect_lambda_max;
   double interval = bisect_lambda_max;

   // Try to find a root with bisection method between lower and upper.
   // If no solution is found, shift the range.

   for (int iter = 0; iter < iter_max; ++iter) {
      try {
         solution = boost::math::tools::bisect(shoot_, lower, upper, stop);
      } catch (boost::exception & e) {
         interval *= 2.;
         lower = upper;
         upper += interval;
         continue;
      }
      return 0.5 * (solution.first + solution.second);
   }

   throw Numerical_error("No solution found from bisect after reaching "
                         + std::to_string(iter_max)
                         + " attempts extending the range.");
}

double Shooting::integrand(double dot_phi, double rho) const
{
   return pow(dot_phi, 2) * pow(rho, dim - 1) / dim;
}

double Shooting::action()
{
   return interior(lambda_sol) + area_n_sphere(dim - 1) * action_no_prefactor;
}

Field_profiles Shooting::calculate_bubble_profile()
{
   // Use solution assuming quadratic potential inside bubble
   const std::size_t num_integration_points = rho_vals.size();
   const std::size_t num_interior_points = num_integration_points / 2;
   const double rho_step = rho_vals[0] / num_interior_points;

   Eigen::VectorXd rho_values(Eigen::VectorXd::Zero(num_interior_points
                                                    + num_integration_points));
   Eigen::MatrixXd field_values(
      Eigen::MatrixXd::Zero(num_interior_points + num_integration_points, 1));
   for (std::size_t i = 0; i < num_interior_points; ++i) {
      rho_values(i) = i * rho_step;
      try {
         field_values(i, 0) = quadratic_potential_solution(i * rho_step);
      } catch (const std::overflow_error & e) {
         field_values(i, 0) = 0;
         logger.log_message(logging::Log_level::Warning,
                            "Overflow error in bubble interior quadratic solution at "
                            "rho = " + std::to_string(rho_values(i)));
      }
   }

   rho_values.segment(num_interior_points, num_integration_points)
      = Eigen::VectorXd::Map(rho_vals.data(), num_integration_points);
   field_values.bottomRows(num_integration_points)
      = Eigen::VectorXd::Map(field_vals.data(), num_integration_points);

   return Field_profiles(rho_values, field_values);
}

double Shooting::quadratic_potential_solution(double rho) const
{
   return quadratic_potential_solution(rho, lambda_sol);
}

double Shooting::quadratic_potential_solution(double rho, double lambda) const
{
   const double phi_zero = unmap(lambda);

   // Trivial case in which we desire \f$\phi(\rho = 0)\f$. This would
   // otherwise result in nans in numerical evaluations of e.g.
   // \f$\lim_{x \to 0} \sinh(x) / x\f$.
   if (rho == 0.) {
      return phi_zero;
   }

   const double mass_squared = potential_second(phi_zero);

   double prime_start = potential_first(phi_zero);
   // @todo More appropriate threshold for Taylor expansion
   if (Abs(prime_start) <= std::numeric_limits<double>::epsilon()) {
      prime_start = potential_second(true_min) * (phi_zero - true_min);
   }
   const double prefactor = prime_start / mass_squared;

   const double y = rho * Sqrt(Abs(mass_squared));
   double delta_phi = 0.;
   if (mass_squared > 0) {
      if (dim == 3) {
         delta_phi = prefactor * (Sinh(y) / y - 1.);
      } else {
         delta_phi = prefactor * (2. * BesselI(1, y) / y - 1.);
      }
   } else {
      if (dim == 3) {
         delta_phi = prefactor * (Sin(y) / y - 1.);
      } else {
         delta_phi = prefactor * (2. * BesselJ(1, y) / y - 1.);
      }
   }

   return phi_zero + delta_phi;
}

double Shooting::interior(double lambda)
{
   const double phi_zero = unmap(lambda);
   const double mass_squared = potential_second(phi_zero);
   const double prime_start = potential_first(phi_zero);
   const double prefactor = Abs(prime_start / mass_squared);
   const double mass = Sqrt(Abs(mass_squared));
   const double change = evolve_change_rel * Abs(phi_zero - false_min);
   evolve_type evolved = evolve(lambda);
   const double rho = evolved[0];
   const double y = mass * rho;

   double value = 0.;
   if (mass_squared > 0.) {
      if (dim == 3) {
         if (y > y_max) {
            value = (change * change) / (6. * mass);  // Asymptotic formula
         } else {
            value = (prefactor * prefactor) / (6. * mass_squared * rho)
               * (1. + mass_squared * (rho * rho)
                  - Cosh(2. * y) + 0.5 * y * Sinh(2. * y));
         }
      } else if (dim == 4) {
         if (y > y_max) {
            value = 1.6875 * (change * change) / mass_squared;  // Asymptotic formula
         } else {
            const double bessel_1_y = BesselI(1, y);
            const double bessel_2_y = BesselI(2, y);
            value =  0.5 * (prefactor * prefactor) * rho / mass *
               (-y * (bessel_1_y * bessel_1_y)
                + 4. * bessel_1_y * bessel_2_y
                + y * (bessel_2_y * bessel_2_y));
         }
      }
   } else {
      if (dim == 3) {
         value = (prefactor * prefactor) / (6. * Abs(mass_squared) * rho)
            * (-1. + Abs(mass_squared) * (rho * rho)
               + Cos(2. * y) + 0.5 * y * Sin(2. * y));
      } else if (dim == 4) {
         const double bessel_1_y = BesselJ(1, y);
         const double bessel_2_y = BesselJ(2, y);
         value =  0.5 * (prefactor * prefactor) * rho / mass *
            (y * (bessel_1_y * bessel_1_y)
             - 4. * bessel_1_y * bessel_2_y
             + y * (bessel_2_y * bessel_2_y));
      }
   }

   return value;
}

}  // namespace BubbleProfiler
