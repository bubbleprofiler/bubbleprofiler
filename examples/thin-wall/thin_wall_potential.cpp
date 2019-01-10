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

#include "thin_wall_potential.hpp"
#include "error.hpp"
#include "math_wrappers.hpp"

namespace BubbleProfiler {

Thin_wall_potential::Thin_wall_potential(double lambda_,
                                         double a_, double epsilon_)
   : lambda(lambda_)
   , a(a_)
   , epsilon(epsilon_)
{
   if (lambda <= 0.) {
      throw Setup_error("Thin_wall_potential: lambda must be positive");
   }

   const bool can_tunnel =
      Abs(epsilon) < 2. * lambda * a * a * a * a/ (3. * Sqrt(3));
   if (!can_tunnel) {
      throw Setup_error(
         "Thin_wall_potential: no barrier exists between minima");
   }
}

double Thin_wall_potential::operator()(
   const Eigen::VectorXd& coords) const
{
   if (coords.size() != 1) {
      throw Setup_error("Thin_wall_potential::partial: "
                        "number of coordinates must be one");
   }

   return operator()(coords(0));
}

double Thin_wall_potential::partial(
   const Eigen::VectorXd& coords, int i) const
{
   if (coords.size() != 1) {
      throw Setup_error("Thin_wall_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Thin_wall_potential::partial: invalid field index "
         + std::to_string(i));
   }

   return first_deriv(coords(0));
}

double Thin_wall_potential::partial(
   const Eigen::VectorXd& coords, int i, int j) const
{
   if (coords.size() != 1) {
      throw Setup_error("Thin_wall_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Thin_wall_potential::partial: invalid field index "
         + std::to_string(i));
   }

   if (j != 0) {
      throw Out_of_bounds_error(
         j, "Thin_wall_potential::partial: invalid field index "
         + std::to_string(j));
   }

   return second_deriv(coords(0));
}

double Thin_wall_potential::operator()(double coords) const
{
   const double phip = scale * coords + origin;

   return 0.125 * lambda * (phip * phip - a * a) * (phip * phip - a * a)
      + 0.5 * epsilon * (phip - a) / a;
}

double Thin_wall_potential::first_deriv(double coords) const
{
   const double phip = scale * coords + origin;

   return 0.5 * lambda * phip * phip * phip - 0.5 * lambda * a * a * phip
      + 0.5 * epsilon / a;
}

double Thin_wall_potential::second_deriv(double coords) const
{
   const double phip = scale * coords + origin;

   return 0.5 * lambda * (3. * phip * phip - a * a);
}

/**
 * The location of the local minimum is computed exactly
 * using the solution to the cubic equation that results
 * from imposing \f$ \partial V / \partial \phi = 0\f$.
 * For \f$\epsilon / (\lambda a^2) \ll 1\f$, one finds
 * that the local minimum is located approximately at
 * \f$\phi_f \approx a\f$.
 */
double Thin_wall_potential::get_local_minimum_location() const
{
   const double x = 1.5 * Sqrt(3.) * epsilon / (lambda * a * a * a * Abs(a));
   const double theta = ArcSin(x) / 3.;

   const double phip1 = 2. * Abs(a) * Sin(theta + 2. * Pi / 3.)
      / Sqrt(3);
   const double phip2 = 2. * Abs(a) * Sin(theta + 4. * Pi / 3.)
      / Sqrt(3);

   const double phi1 = (phip1 - origin) / scale;
   const double phi2 = (phip2 - origin) / scale;

   return (operator()(phi1) > operator()(phi2) ? phi1 : phi2);
}

/**
 * The location of the local maximum (i.e., the potential barrier)
 * is computed exactly using the solution to the cubic equation that results
 * from imposing \f$ \partial V / \partial \phi = 0\f$.
 * For \f$\epsilon / (\lambda a^2) \ll 1\f$, one finds
 * that the barrier is approximately located at
 * \f$\phi_b \approx \epsilon / (\lambda a^3)\f$.
 */
double Thin_wall_potential::get_local_maximum_location() const
{
   const double x = 1.5 * Sqrt(3.) * epsilon / (lambda * a * a * a * Abs(a));
   const double theta = ArcSin(x) / 3.;

   const double phip = 2. * Abs(a) * Sin(theta) / Sqrt(3);

   return (phip - origin) / scale;
}

/**
 * The location of the global minimum is computed exactly
 * using the solution to the cubic equation that results
 * from imposing \f$ \partial V / \partial \phi = 0\f$.
 * For \f$\epsilon / (\lambda a^2) \ll 1\f$, one finds
 * that the global minimum is located approximately at
 * \f$\phi_t \approx -a\f$.
 */
double Thin_wall_potential::get_global_minimum_location() const
{
   const double x = 1.5 * Sqrt(3.) * epsilon / (lambda * a * a * a * Abs(a));
   const double theta = ArcSin(x) / 3.;

   const double phip1 = 2. * Abs(a) * Sin(theta + 2. * Pi / 3.)
      / Sqrt(3);
   const double phip2 = 2. * Abs(a) * Sin(theta + 4. * Pi / 3.)
      / Sqrt(3);

   const double phi1 = (phip1 - origin) / scale;
   const double phi2 = (phip2 - origin) / scale;

   return (operator()(phi1) <= operator()(phi2) ? phi1 : phi2);
}

/**
 * An approximate value for the Euclidean action is calculated
 * using the analytic expression for the action derived in
 * the thin-wall approximation,
 * \f[
 *    S \approx \frac{8 \pi^2 \lambda^2 a^{12}}{3 \epsilon^3} .
 * \f]
 */
double Thin_wall_potential::get_thin_wall_action() const
{
   const double pisq = boost::math::double_constants::pi
      * boost::math::double_constants::pi;

   return 8. * pisq * lambda * lambda * std::pow(a, 12) /
      (3. * epsilon * epsilon * epsilon);
}

} // namespace BubbleProfiler
