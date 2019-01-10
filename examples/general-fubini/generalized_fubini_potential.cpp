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

#include "generalized_fubini_potential.hpp"
#include "error.hpp"
#include "math_wrappers.hpp"

namespace BubbleProfiler {

Generalized_fubini_potential::Generalized_fubini_potential(double u_, double v_,
                                                           double m_)
   : u(u_), v(v_), m(m_)
{
   if (u <= 0.) {
      throw Setup_error("Generalized_fubini_potential: u must be positive");
   }
   if (v <= 0.) {
      throw Setup_error("Generalized_fubini_potential: v must be positive");
   }
   if (m <= 1.) {
      throw Setup_error(
         "Generalized_fubini_potential: m must be greater than 1");
   }
}

double Generalized_fubini_potential::operator()(
   const Eigen::VectorXd& coords) const
{
   if (coords.size() != 1) {
      throw Setup_error("Generalized_fubini_potential::partial: "
                        "number of coordinates must be one");
   }

   return operator()(coords(0));
}

double Generalized_fubini_potential::partial(
   const Eigen::VectorXd& coords, int i) const
{
   if (coords.size() != 1) {
      throw Setup_error("Generalized_fubini_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Generalized_fubini_potential::partial: invalid field index "
         + std::to_string(i));
   }

   return first_deriv(coords(0));
}

double Generalized_fubini_potential::partial(
   const Eigen::VectorXd& coords, int i, int j) const
{
   if (coords.size() != 1) {
      throw Setup_error("Generalized_fubini_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Generalized_fubini_potential::partial: invalid field index "
         + std::to_string(i));
   }

   if (j != 0) {
      throw Out_of_bounds_error(
         j, "Generalized_fubini_potential::partial: invalid field index "
         + std::to_string(j));
   }

   return second_deriv(coords(0));
}

double Generalized_fubini_potential::operator()(double coords) const
{
   const double phip = Abs(scale * coords + origin);

   const auto pow1 = 2. + 1. / m;
   const auto pow2 = 2. + 2. / m;

   return (4. * u * m * m * (m - 1.) * std::pow(phip, pow1)) /
      (2. * m + 1) - 2. * u * v * m * m * std::pow(phip, pow2);
}

double Generalized_fubini_potential::first_deriv(double coords) const
{
   const double phip = Abs(scale * coords + origin);

   const auto pow1 = 1. + 1. / m;
   const auto pow2 = 1. + 2. / m;

   return 4. * u * m * (m - 1.) * std::pow(phip, pow1)
      - 4. * u * v * m * (m + 1.) * std::pow(phip, pow2);
}

double Generalized_fubini_potential::second_deriv(double coords) const
{
   const double phip = Abs(scale * coords + origin);

   const auto pow1 = 1. / m;
   const auto pow2 = 2. / m;

   return 4. * u * (m * m - 1.) * std::pow(phip, pow1)
      - 4. * u * v * (m * m + 3. * m + 2.) * std::pow(phip, pow2);
}

double Generalized_fubini_potential::get_local_minimum_location() const
{
   return -origin / scale;
}

double Generalized_fubini_potential::get_local_maximum_location() const
{
   const double phip = std::pow(m - 1., m) / std::pow(v * (m + 1.), m);

   return (phip - origin) / scale;
}

double Generalized_fubini_potential::get_bounce_solution_at(double r) const
{
   const double phip = 1. / std::pow(u * r * r + v, m);
   return (phip - origin) / scale;
}

Field_profiles Generalized_fubini_potential::get_profile(
   const Eigen::VectorXd& rho_values) const
{
   const auto n_grid_points = rho_values.size();
   Eigen::MatrixXd field_profile(Eigen::MatrixXd::Zero(n_grid_points, 1));
   for (int i = 0; i < n_grid_points; ++i) {
      field_profile(i, 0) = get_bounce_solution_at(rho_values(i));
   }

   return Field_profiles(rho_values, field_profile);
}

double Generalized_fubini_potential::get_action() const
{
   return m * Pi * Pi / ((4. * m * m - 1.) * u * std::pow(v, 2. * m - 1.));
}

} // namespace BubbleProfiler
