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

#include "perturbations_ode_system.hpp"
#include "error.hpp"
#include "field_profiles.hpp"
#include "potential.hpp"

#include <cmath>
#include <sstream>

namespace BubbleProfiler {

Perturbations_ODE_system::Perturbations_ODE_system(
   Potential& potential_, Field_profiles& profiles_,
   int n_spacetime_dimensions_)
   : n_fields(potential_.get_number_of_fields())
   , n_spacetime_dimensions(n_spacetime_dimensions_)
   , potential(&potential_)
   , profiles(&profiles_)
{
}

void Perturbations_ODE_system::operator()(
   const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, double rho) const
{
   if (rho <= 0.) {
      throw Domain_error("Perturbations_ODE_system::operator(): "
                         "spatial coordinate must be positive");
   }

   if (!is_finite(x)) {
      std::ostringstream msg;
      msg << "Perturbations_ODE_system::operator(): "
          << "perturbations contained a nan at rho = " << rho;
      throw Numerical_error(msg.str());
   }

   if (dxdt.rows() != 2 * n_fields) {
      dxdt.resize(2 * n_fields, Eigen::NoChange);
   }

   dxdt.segment(0, n_fields) = x.segment(n_fields, n_fields);
   dxdt.segment(n_fields, n_fields) = calculate_mass_matrix(rho)
      * x.segment(0, n_fields) - ((n_spacetime_dimensions - 1) / rho)
      * x.segment(n_fields, n_fields) + calculate_inhomogeneities(rho);

   if (!is_finite(dxdt)) {
      std::ostringstream msg;
      msg << "Perturbations_ODE_system::operator(): "
          << "derivatives contained a nan at rho = " << rho;
      throw Numerical_error(msg.str());
   }
}

bool Perturbations_ODE_system::is_finite(const Eigen::VectorXd& x) const
{
#if EIGEN_VERSION_AT_LEAST(3,2,0)
   return x.allFinite();
#else
   const int n = x.size();
   for (int i = 0; i < n; ++i) {
      if (!std::isfinite(x(i))) {
         return false;
      }
   }
   return true;
#endif
}

Eigen::VectorXd Perturbations_ODE_system::calculate_inhomogeneities(
   double rho) const
{
   using Index = Eigen::VectorXd::Index;

   Eigen::VectorXd inh(Eigen::VectorXd::Zero(n_fields));
   const auto field_values = profiles->evaluate_at(rho);

   for (Index i = 0; i < n_fields; ++i) {
      inh(i) = potential->partial(field_values, i)
         - profiles->derivative_at(i, 2, rho)
         - ((n_spacetime_dimensions - 1.) / rho)
         * profiles->derivative_at(i, 1, rho);
   }

   return inh;
}

Eigen::MatrixXd Perturbations_ODE_system::calculate_mass_matrix(
   double rho) const
{
   using Index = Eigen::MatrixXd::Index;

   Eigen::MatrixXd mass_matrix(n_fields, n_fields);
   const auto field_values = profiles->evaluate_at(rho);

   for (Index j = 0; j < n_fields; ++j) {
      for (Index i = j; i < n_fields; ++i) {
         mass_matrix(i, j)
            = potential->partial(field_values, i, j);
      }
   }

   for (Index j = 0; j < n_fields; ++j) {
      for (Index i = 0; i < j; ++i) {
         mass_matrix(i, j) = mass_matrix(j, i);
      }
   }

   return mass_matrix;
}

} // namespace BubbleProfiler
