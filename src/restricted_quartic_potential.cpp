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

#include "restricted_quartic_potential.hpp"
#include "error.hpp"

namespace BubbleProfiler {

Restricted_quartic_potential::Restricted_quartic_potential(double alpha_)
   : alpha(alpha_)
{
   if (alpha <= 0.5 || alpha >= 0.75) {
      throw Setup_error(
         "alpha must be between 0.5 and 0.75");
   }
}

Restricted_quartic_potential::Restricted_quartic_potential(
   double alpha_, double E_)
   : alpha(alpha_)
   , E(E_)
{

   if (alpha <= 0.5 || alpha >= 0.75) {
      throw Setup_error(
         "alpha must be between 0.5 and 0.75");
   }

   if (E <= 0.) {
      throw Setup_error(
         "E must be positive");
   }
}

double Restricted_quartic_potential::operator()(
   const Eigen::VectorXd& coords) const
{
   if (coords.size() != 1) {
      throw Setup_error("Restricted_quartic_potential::partial: "
                        "number of coordinates must be one");
   }

   return operator()(coords(0));
}

double Restricted_quartic_potential::partial(
   const Eigen::VectorXd& coords, int i) const
{
   if (coords.size() != 1) {
      throw Setup_error("Restricted_quartic_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Restricted_quartic_potential::partial: invalid field index "
         + std::to_string(i));
   }

   return first_deriv(coords(0));
}

double Restricted_quartic_potential::partial(const Eigen::VectorXd& coords,
                                    int i, int j) const
{
   if (coords.size() != 1) {
      throw Setup_error("Restricted_quartic_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Restricted_quartic_potential::partial: invalid field index "
         + std::to_string(i));
   }

   if (j != 0) {
      throw Out_of_bounds_error(
         j, "Restricted_quartic_potential::partial: invalid field index "
         + std::to_string(j));
   }

   return second_deriv(coords(0));
}


double Restricted_quartic_potential::operator()(double coords) const
{
   const double phip = scale * coords + origin;
   return 0.5 * (3. - 4. * alpha) * E * phip * phip
      - E * phip * phip * phip
      + alpha * E * phip * phip * phip * phip
      + offset;
}

double Restricted_quartic_potential::first_deriv(double coords) const
{
   const double phip = scale * coords + origin;
   return (3. - 4. * alpha) * E * phip - 3. * E * phip * phip
      + 4. * alpha * E * phip * phip * phip;
}

double Restricted_quartic_potential::second_deriv(double coords) const
{
   const double phip = scale * coords + origin;
   return (3. - 4. * alpha) * E - 6. * E * phip
      + 12. * alpha * E * phip * phip;
}

double Restricted_quartic_potential::get_global_minimum_location() const
{
   return 1 - origin;
}

double Restricted_quartic_potential::get_local_minimum_location() const
{
   return -origin;
}

double Restricted_quartic_potential::get_local_maximum_location() const
{
   return 0.25 * (3. - 4. * alpha) / alpha - origin;
}

} // namespace BubbleProfiler
