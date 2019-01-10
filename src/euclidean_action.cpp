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
 * @file euclidean_action.cpp
 * @brief implementation of helper functions for calculating Euclidean action
 */

#include "euclidean_action.hpp"
#include "error.hpp"
#include "field_profiles.hpp"
#include "numeric.hpp"
#include "potential.hpp"

#include <cmath>

namespace BubbleProfiler {

/**
 * Given the bounce solution for \f$n\f$ scalar fields in \f$d\f$
 * dimensions defined on a finite domain, the contribution to the
 * action due to the field kinetic terms,
 * \f[
 *    S_T \equiv S_{d-1} \int_0^\infty d\rho \rho^{d-1} \sum_{i=1}^n
 *    \frac{1}{2} \dot{\phi}_i^2 ,
 * \f]
 * where \f$S_{d-1}\f$ is the surface area of the \f$d\f$-sphere,
 * is first calculated using numerical quadrature, with the
 * quadrature rule chosen according to the value of \c rule. The
 * adaptive quadrature routine attempts to compute an estimate
 * \f$\tilde{I}\f$ for the integral \f$I\f$ appearing above
 * that satisfies
 * \f[
 *    | \tilde{I} - I | \leq \max ( \epsilon_{\text{abs}},
 *    \epsilon_{\text{rel}} |I|) ,
 * \f]
 * with \f$\epsilon_{\text{abs}}\f$ and \f$\epsilon_{\text{rel}}\f$
 * given by \c abs_tol and \c rel_tol, respectively. The maximum
 * allowed number of subintervals used in this integration is
 * set by the value of \c max_intervals. The full value of the
 * action to be returned is then obtained from the value of
 * \f$S_T\f$ using the relation
 * \f[
 *    S_E = \frac{2}{d} S_T .
 * \f]
 *
 * @sa integrate_gsl_qag
 */
double calculate_kinetic_action(
   const Potential& /* potential */, const Field_profiles& profiles,
   std::size_t max_intervals, double rel_tol, double abs_tol,
   Integration_rule rule)
{
   const int n_fields = profiles.get_number_of_fields();
   const int d = profiles.get_number_of_dimensions();
   const double domain_start = profiles.get_domain_start();
   const double domain_end = profiles.get_domain_end();

   // finite domain correction
   double correction = 0;
   for (int i = 0; i < n_fields; ++i) {
      const double deriv = profiles.derivative_at(i, 1, domain_end);
      correction += 0.5 * deriv * deriv;
   }

   const auto integrand
      = [&profiles, n_fields, d, correction](double rho)
   {
      double result = 0.;
      for (int i = 0; i < n_fields; ++i) {
         const double deriv = profiles.derivative_at(i, 1, rho);
         result += 0.5 * deriv * deriv;
      }

      result -= correction;
      result *= std::pow(rho, d - 1);

      return result;
   };

   const auto result = integrate_gsl_qag(integrand, domain_start, domain_end,
                                         abs_tol, rel_tol, max_intervals, rule);

   const auto error = std::get<0>(result);
   if (error) {
      throw Numerical_error("calculation of action failed with error code: "
                            + std::to_string(error));
   }

   return 2. * area_n_sphere(d - 1) * std::get<1>(result) / d;
}

/**
 * Given the bounce solution for \f$n\f$ scalar fields in \f$d\f$
 * dimensions defined on a finite domain, the contribution to the
 * action due to the field potential terms,
 * \f[
 *    S_V \equiv S_{d-1} \int_0^\infty d\rho \rho^{d-1}
 *    V(\phi_i) ,
 * \f]
 * where \f$S_{d-1}\f$ is the surface area of the \f$d\f$-sphere
 * and the potential is assumed to vanish at the false vacuum,
 * is first calculated using numerical quadrature, with the
 * quadrature rule chosen according to the value of \c rule. The
 * adaptive quadrature routine attempts to compute an estimate
 * \f$\tilde{I}\f$ for the integral \f$I\f$ appearing above
 * that satisfies
 * \f[
 *    | \tilde{I} - I | \leq \max ( \epsilon_{\text{abs}},
 *    \epsilon_{\text{rel}} |I|) ,
 * \f]
 * with \f$\epsilon_{\text{abs}}\f$ and \f$\epsilon_{\text{rel}}\f$
 * given by \c abs_tol and \c rel_tol, respectively. The maximum
 * allowed number of subintervals used in this integration is
 * set by the value of \c max_intervals. The full value of the
 * action to be returned is then obtained from the value of
 * \f$S_V\f$ using the relation
 * \f[
 *    S_E = \frac{2}{2 - d} S_V .
 * \f]
 *
 * @sa integrate_gsl_qag
 */
double calculate_potential_action(
   const Potential& potential, const Field_profiles& profiles,
   std::size_t max_intervals, double rel_tol, double abs_tol,
   Integration_rule rule)
{
   const int d = profiles.get_number_of_dimensions();
   const double domain_start = profiles.get_domain_start();
   const double domain_end = profiles.get_domain_end();

   // finite domain correction
   double correction = potential(profiles.evaluate_at(domain_end));

   const auto integrand
      = [&potential, &profiles, d, correction](double rho)
   {
      double result = potential(profiles.evaluate_at(rho));
      result -= correction;
      result *= std::pow(rho, d - 1);

      return result;
   };

   const auto result = integrate_gsl_qag(integrand, domain_start, domain_end,
                                         abs_tol, rel_tol, max_intervals, rule);

   const auto error = std::get<0>(result);
   if (error) {
      throw Numerical_error("calculation of action failed with error code: "
                            + std::to_string(error));
   }

   return 2. * area_n_sphere(d - 1) * std::get<1>(result) / (2. - d);
}

/**
 * Given the bounce solution for \f$n\f$ scalar fields in \f$d\f$
 * dimensions defined on a finite domain, an approximate value for
 * the full Euclidean action
 * \f[
 *    S_E \equiv S_{d-1} \int_0^\infty d\rho \rho^{d-1} \left (
 *    \sum_{i=1}^n \frac{1}{2} \dot{\phi}_i^2 + V(\phi_i) \right ) ,
 * \f]
 * where \f$S_{d-1}\f$ is the surface area of the \f$d\f$-sphere
 * and the potential is assumed to vanish at the false vacuum,
 * is calculated using numerical quadrature, with the
 * quadrature rule chosen according to the value of \c rule. The
 * adaptive quadrature routine attempts to compute an estimate
 * \f$\tilde{I}\f$ for the integral \f$I\f$ appearing above
 * that satisfies
 * \f[
 *    | \tilde{I} - I | \leq \max ( \epsilon_{\text{abs}},
 *    \epsilon_{\text{rel}} |I|) ,
 * \f]
 * with \f$\epsilon_{\text{abs}}\f$ and \f$\epsilon_{\text{rel}}\f$
 * given by \c abs_tol and \c rel_tol, respectively. The maximum
 * allowed number of subintervals used in this integration is
 * set by the value of \c max_intervals.
 *
 * @sa integrate_gsl_qag
 */
double calculate_full_action(
   const Potential& potential, const Field_profiles& profiles,
   std::size_t max_intervals, double rel_tol,
   double abs_tol, Integration_rule rule)
{
   const int n_fields = profiles.get_number_of_fields();
   const int d = profiles.get_number_of_dimensions();
   const double domain_start = profiles.get_domain_start();
   const double domain_end = profiles.get_domain_end();

   // finite domain correction
   double correction = 0.;
   for (int i = 0; i < n_fields; ++i) {
      const double deriv = profiles.derivative_at(i, 1, domain_end);
      correction += 0.5 * deriv * deriv;
   }
   correction += potential(profiles.evaluate_at(domain_end));

   const auto integrand
      = [&potential, &profiles, n_fields, d, correction](double rho)
      {
         double result = 0.;
         for (int i = 0; i < n_fields; ++i) {
            const double deriv = profiles.derivative_at(i, 1, rho);
            result += 0.5 * deriv * deriv;
         }

         result += potential(profiles.evaluate_at(rho));
         result -= correction;
         result *= pow(rho, d - 1);

         return result;
      };

   const auto result = integrate_gsl_qag(integrand, domain_start, domain_end,
                                         abs_tol, rel_tol, max_intervals, rule);

   const auto error = std::get<0>(result);
   if (error) {
      throw Numerical_error("calculation of action failed with error code: "
                            + std::to_string(error));
   }

   return area_n_sphere(d - 1) * std::get<1>(result);
}

/**
 * The method used to evaluate the action is chosen based on the
 * value of \c use_kinetic, with a value of \c true indicating
 * that the action should be evaluated using only the contribution
 * from the field kinetic terms. If \c use_kinetic is \c false,
 * the action is evaluated using the contribution from the potential.
 *
 * Internally, this simply function calls
 * BubbleProfiler::calculate_kinetic_action or
 * BubbleProfiler::calculate_potential_action, as appropriate.
 *
 * @sa calculate_kinetic_action, calculate_potential_action
 */
double calculate_action(
   const Potential& potential, const Field_profiles& profiles,
   std::size_t max_intervals, double rel_tol,
   double abs_tol, Integration_rule rule,
   bool use_kinetic)
{
   return use_kinetic ?
      calculate_kinetic_action(potential, profiles, max_intervals,
                               rel_tol, abs_tol, rule) :
      calculate_potential_action(potential, profiles, max_intervals,
                                 rel_tol, abs_tol, rule);
}

} // namespace BubbleProfiler
