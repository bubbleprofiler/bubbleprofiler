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

#ifndef BUBBLEPROFILER_SHOOTING_HPP_INCLUDED
#define BUBBLEPROFILER_SHOOTING_HPP_INCLUDED

/**
    @file
    @brief One-dimensional shooting method
    @example action.cpp
    @example logarithmic.cpp
    @example tabulate.cpp
    @example general_fubini.cpp
    @example scale.cpp
    @example thin_wall.cpp
*/

#include "field_profiles.hpp"
#include "basic_logger.hpp"

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

#include <functional>

namespace BubbleProfiler {

class Potential;

/** @brief Solve the one-dimensional problem using the shooting method */
class Shooting {
public:
   enum Solver_options : unsigned int {
      Compute_action = 0x01,
      Compute_profile = 0x02
   };

   /** Bits of precision in bisection of \f$\phi_0\f$ */
   void set_bisection_precision_bits(int b) { shoot_bisect_bits = b; }
   /** Relative tolerance for judging when\f$\phi \approx \phi_f\f$ */
   void set_action_arrived_rel(double tol) { action_arrived_rel = tol; }
   /** Absolute tolerance in ODE solver for a shot */
   void set_shooting_abs_tol(double tol) { shoot_ode_abs = tol; }
   /** Relative tolerance in ODE solver for a shot */
   void set_shooting_rel_tol(double tol) { shoot_ode_rel = tol; }
   /** Absolute tolerance in ODE solver for action integral */
   void set_action_abs_tol(double tol) { action_ode_abs = tol; }
    /** Relative tolerance in ODE solver for action integral */
   void set_action_rel_tol(double tol) { action_ode_rel = tol; }
   /** The initial ODE step size, \f$d\rho\f$, is this fraction of the bubble scale */
   void set_drho_frac(double frac) { drho_frac = frac; }
   /**
      The EOM is integrated analytically until
      \f[
      \frac{\left|\phi(\hat\rho) - \phi_0 \right|}{\left|\phi_f - \phi_0\right|},
      \f]
      changes by this amount.
    */
   void set_evolve_change_rel(double frac) { evolve_change_rel = frac; }
   /** Maximum value of \f$\lambda\f$ */
   void set_bisect_lambda_max(double l) { bisect_lambda_max = l; }
   /** Maximum number of iterations of shooting method */
   void set_max_iterations(int i) { iter_max = i; }
   /** Maximum number periods before the ODE integration aborts */
   void set_max_periods(double p) { periods_max = p; }
   /** Threshold for asymptotic approximation in analytic evolution */
   void set_f_y_max(double f) { f_y_max = f; }
  /** Threshold for asymptotic approximation in analytic evolution */
   void set_f_y_min(double f) { f_y_min = f; }
   /** Threshold for asymptotic approximation in analytic evolution */
   void set_y_max(double y) { y_max = y; }
   /**
      @overload
      Solver that accepts potential and first and second derivatives
   */
   void solve(const std::function<double(double)>& potential_,
              const std::function<double(double)>& potential_first_,
              const std::function<double(double)>& potential_second_,
              double false_min_, double true_min_, double barrier_,
              int dim_ = 3,
              unsigned int options = (Solver_options::Compute_action |
                                      Solver_options::Compute_profile));
    /**
       @overload
       Solver that accepts a Potential object
   */
   void solve(const Potential& potential, double false_min_,
              double true_min_, double barrier_, int dim_ = 3,
              unsigned int options = (Solver_options::Compute_action |
                                      Solver_options::Compute_profile));

   /** Bubble profile, i.e., \f$\phi(\rho)\f$ */
   const Field_profiles& get_bubble_profile() const;
   /** Euclidean action, \f$S\f$ */
   double get_euclidean_action() const;

   /** @overload */
   double quadratic_potential_solution(double) const;
   /** @returns Quadratic solution \f$\phi(\rho)\f$ */
   double quadratic_potential_solution(double, double) const;

private:
    logging::Basic_logger logger{};

   using evolve_type = boost::array<double, 3>;
   using state_type = boost::array<double, 2>;
   using error_stepper_type =
      boost::numeric::odeint::runge_kutta_dopri5<state_type>;
   using controlled_stepper_type =
      boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;

   /** Number of space-time dimensions, \f$d\f$ */
   int dim{3};

   /** Solution for \f$\lambda\f$ */
   double lambda_sol{0.};
   /** Solution for Euclidean action, \f$S\f$ */
   double euclidean_action{0.};
   Field_profiles profiles{};
   bool action_computed{false};
   bool profile_computed{false};

   // ODE steppers
   controlled_stepper_type shoot_stepper{};
   controlled_stepper_type action_stepper{};

   /** Guess of step size, \f$d\rho\f$ */
   double drho_guess{0.};
   /** Estimated scale of bubble */
   double scale{0.};

  /** Location of the barrier, \f$\phi_b\f$ */
   double barrier{0.};
   /** Location of the false minimum, \f$\phi_f\f$ */
   double false_min{0.};
   /** Location of the true minimum, \f$\phi_t\f$ */
   double true_min{1.};
   /** \f$\textrm{sign}\,\phi_f - \phi_t\f$ */
   double sign_min{1.};
   /** Potential at barrier, \f$V(\phi_b)\f$ */
   double p_barrier{0.};
   /** Potential at false minimum, \f$V(\phi_f)\f$ */
   double p_false_min{0.};
   /** Potential at true minimum, \f$V(\phi_t)\f$ */
   double p_true_min{0.};

   /** Potential, \f$V(\phi)\f$ */
   std::function<double(double)> potential{nullptr};
    /** Derivative of potential, \f$V^{\prime}(\phi)\f$ */
   std::function<double(double)> potential_first{nullptr};
   /** Second derivative of potential, \f$V^{\prime\prime}(\phi)\f$ */
   std::function<double(double)> potential_second{nullptr};

   // Tolerances and integration method settings
   int shoot_bisect_bits{5};
   double action_arrived_rel{1.e-3};
   double shoot_ode_abs{1.e-4};
   double shoot_ode_rel{1.e-4};
   double action_ode_abs{1.e-6};
   double action_ode_rel{1.e-6};
   double drho_frac{1.e-3};
   double evolve_change_rel{1.e-2};
   double bisect_lambda_max{5.};
   int iter_max{100000};
   double periods_max{1.e2};
   double f_y_max{1.e6};
   double f_y_min{1.e-3};
   double y_max{1.e1};

   void initialize_potential_parameters(
      const std::function<double(double)>& potential_,
      const std::function<double(double)>& potential_first_,
      const std::function<double(double)>& potential_second_,
      double false_min_, double true_min_, double barrier_, int dim_);

   /**
      @brief ODE for bounce solution

      Rewrite second-order as two coupled first order equations,
      \f[
      \dot\phi = \chi
      \f]
      and
      \f[
      \dot\chi = -\frac{2}{\rho} \phi + V^\prime(\phi)
      \f]

      @param x \f$(\phi, \chi)\f$
      @param dxdt \f$(\dot\phi, \dot\chi)\f$
      @param rho \f$\rho\f$
   */
   void ode(const state_type &x, state_type &dxdt, double rho);
   /**
      @brief Find the (dimensionless) time at which the field changes by
      a small amount, \f$m \hat\rho\f$

      @returns \f$m \hat\rho\f$
      @param lambda \f$\lambda\f$, related to \f$\phi_0\f$
      @param dVdphi \f$\left.\frac{dV}{d\phi}\right|_\lambda\f$
      @param d2Vdphi2 \f$\left.\frac{d^2V}{d\phi^2}\right|_\lambda\f$
   */
   double calculate_roll_time(double lambda, double dVdphi,
                                  double d2Vdphi2) const;
   /**
      @brief Find velocity \f$\hat\phi(\hat\rho)\f$

      @returns \f$\hat\phi(\hat\rho)\f$
      @param lambda \f$\lambda\f$, related to \f$\phi_0\f$
      @param y \f$y = m \hat\rho\f$
      @param dVdphi \f$\left.\frac{dV}{d\phi}\right|_\lambda\f$
      @param d2Vdphi2 \f$\left.\frac{d^2V}{d\phi^2}\right|_\lambda\f$
   */
   double calculate_roll_velocity(double lambda, double y,
                                  double dVdphi, double d2Vdphi2) const;
   /**
      @brief Uses an approximation solution to evolve the system forward
      until the field changes by a required amount

      The required amount is specified by
      `Shooting::rchange` using:
      \f[
      |\phi(\rho) - \phi_0| =
      Shooting::rchange |\phi_f - \phi_0|.
      \f]
      This method is used to avoid solving the ODE during the time
      period where friction dominates, which can be very slow, and
      avoids a singularity at \f$\rho =0\f$.

      @param lambda \f$\lambda\f$
      @returns \f$(\hat\rho, \phi(\hat\rho), \dot\phi(\hat\rho))\f$
    */
   evolve_type evolve(double lambda) const;
   /**
      @brief Characteristic scale of bubble

      Found from period of oscillation about barrier. We approximate the
      potential at the barrier by a quadratic such that
      \f[
      T = \frac{2\pi}{\sqrt{-V^{\prime\prime}(\phi_b)}}.
      \f]
   */
   double bubble_scale() const;

   /** @returns Energy at \f$(\phi, \dot\phi)\f$ */
   double energy(double phi, double dot_phi) const;

   /**
      @brief Performs a single shot of our shooting method

      Evolves field  from an initial guess for \f$\phi_0\f$, given by
      \f$\lambda\f$ and tests whether it under or over shoots the false vacuum.

      @param lambda \f$\lambda\f$
      @returns Whether overshot or undershot
   */
   double shoot(double lambda);

   /** @returns \f$\phi\f$ from \f$\lambda\f$  */
   double unmap(double lambda) const;

   /**
      @brief Solve for \f$\phi_0\f$ by bisection

      Bisect between over and under shots, beginning from the interval
      \f$(\phi_b, \phi_t)\f$.

      @returns \f$\lambda\f$ corresponding to \f$\phi_0\f$
   */
   double shooting();

   /**
      @brief Action integrand from kinetic term without pre-factor

      We use a trick that
      \f[
      I = \int d\rho \rho^{d-1} \left[ \frac12 \dot\phi^2 + V(\phi) - V(\phi_f)\right] \equiv  I_1 + I_2
      \f]
      but by the fact that action is stationary under \f$\rho \to a \rho\f$
      \f[
      d I_2 = (2 - d) I_1
      \f]
      such that
      \f[
      I = I_1 + I_2 = \frac{2}{d} I_1 = \frac{1}{d} \int d\rho \dot\phi^2 \rho^{d-1}
      \f]

      @param dot_phi \f$\dot\phi\f$
      @param rho \f$\rho\f$
      @returns Action without prefactor
   */
   double integrand(double dot_phi, double rho) const;

  /** @returns Action, \f$S\f$ */
   double action(double lambda);

   /** @returns Action without prefactor for interior of bubble */
   double interior(double lambda);

   /** @returns The bubble profile, \f$\phi(\rho)\f$ */
   Field_profiles calculate_bubble_profile(double lambda);
};

}  // namespace BubbleProfiler

#endif  // BUBBLEPROFILER_SHOOTING_HPP_INCLUDED
