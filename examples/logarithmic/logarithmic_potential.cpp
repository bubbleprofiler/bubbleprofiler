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

#include "logarithmic_potential.hpp"
#include "error.hpp"
#include "math_wrappers.hpp"

namespace BubbleProfiler {

double Logarithmic_potential::operator()(
   const Eigen::VectorXd& coords) const
{
   if (coords.size() != 1) {
      throw Setup_error("Logarithmic_potential::partial: "
                        "number of coordinates must be one");
   }
   return this->operator()(coords(0));
}

double Logarithmic_potential::partial(
   const Eigen::VectorXd& coords, int i) const
{
   if (coords.size() != 1) {
      throw Setup_error("Logarithmic_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Logarithmic_potential::partial: invalid field index "
         + std::to_string(i));
   }

   return first_deriv(coords(0));
}

double Logarithmic_potential::partial(
   const Eigen::VectorXd& coords, int i, int j) const
{
   if (coords.size() != 1) {
      throw Setup_error("Logarithmic_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Logarithmic_potential::partial: invalid field index "
         + std::to_string(i));
   }

   if (j != 0) {
      throw Out_of_bounds_error(
         j, "Logarithmic_potential::partial: invalid field index "
         + std::to_string(j));
   }

   return second_deriv(coords(0));
}

double Logarithmic_potential::operator()(double coords) const
{
   const double phip = scale * coords + origin;

   if (Abs(phip) < std::numeric_limits<double>::min()) {
      return 0.;
   }

   return 0.5 * m * m * phip * phip * (1. - Log(phip * phip / (w * w)));
}

double Logarithmic_potential::first_deriv(double coords) const
{
   const double phip = scale * coords + origin;

   if (Abs(phip) < std::numeric_limits<double>::min()) {
      return 0.;
   }

   return -m * m * phip * Log(phip * phip / (w * w));
}

double Logarithmic_potential::second_deriv(double coords) const
{
   const double phip = scale * coords + origin;

   if (Abs(phip) < std::numeric_limits<double>::min()) {
      return 0.;
   }

   return -m * m * (2. + Log(phip * phip / (w * w)));
}

double Logarithmic_potential::get_local_minimum_location() const
{
   return -origin / scale;
}

double Logarithmic_potential::get_local_maximum_location() const
{
   const double phip = Abs(w);

   return (phip - origin) / scale;
}

double Logarithmic_potential::get_bounce_solution_at(double r) const
{
   const double phip = w * Exp(-0.5 * m * m * r *r + 2.);
   return (phip - origin) / scale;
}

Field_profiles Logarithmic_potential::get_profile(
   const Eigen::VectorXd& rho_values) const
{
   const auto n_grid_points = rho_values.size();
   Eigen::MatrixXd field_profile(Eigen::MatrixXd::Zero(n_grid_points, 1));
   for (int i = 0; i < n_grid_points; ++i) {
      field_profile(i, 0) = get_bounce_solution_at(rho_values(i));
   }

   return Field_profiles(rho_values, field_profile);
}

double Logarithmic_potential::get_action() const
{
   return 0.5 * Pi * Pi * Exp(4.) * w * w / (m * m);
}

} // namespace BubbleProfiler
