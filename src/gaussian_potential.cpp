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
 * @file gaussian_potential.cpp
 * @brief contains the implementation of the Gaussian_potential class
 */

#include "gaussian_potential.hpp"
#include "error.hpp"
#include "math_wrappers.hpp"

#include <cmath>

namespace BubbleProfiler {

Gaussian_potential::Gaussian_potential(
   double gamma_, double lambda_, std::size_t n_fields_)
   : gamma(gamma_)
   , lambda(lambda_)
   , n_fields(n_fields_)
{
   if (gamma <= 1) {
      throw Setup_error(
         "Gamma must be > 1 so that the origin is a false vacuum.");
   }

   origin = Eigen::VectorXd::Zero(n_fields);
   origin_translation = origin;
   mu = Eigen::VectorXd::Zero(n_fields);

   for (std::size_t i = 0; i < n_fields; ++i) {
      mu(i) = lambda / std::sqrt(n_fields);
   }

   basis_transform = Eigen::MatrixXd::Identity(n_fields, n_fields);
}

double Gaussian_potential::gaussian(const Eigen::VectorXd& coords,
                                    const Eigen::VectorXd& mu) const
{
   double exponent = -0.5*std::pow((coords - mu).norm(), 2);
   return std::pow(2 * Pi, -(n_fields / 2.0)) * std::exp(exponent);
}

double Gaussian_potential::gaussian_deriv(
   const Eigen::VectorXd& coords, const Eigen::VectorXd& mu, int i) const
{
   double exponent = -0.5*std::pow((coords - mu).norm(), 2);
   return -std::pow(2 * Pi, -(n_fields / 2.0))
      * (coords(i) - mu(i)) * std::exp(exponent);
}

double Gaussian_potential::gaussian_deriv2(
   const Eigen::VectorXd& coords, const Eigen::VectorXd& mu,
   int i, int j) const
{
   double exponent = -0.5*std::pow((coords - mu).norm(), 2);
   double delta_ij = i == j ? 1.0 : 0.0;
   return -std::pow(2 * Pi, -(n_fields / 2.0))
      *(delta_ij - (coords(i) - mu(i))*(coords(j) - mu(j)))*std::exp(exponent);
 }

double Gaussian_potential::operator()(const Eigen::VectorXd& coords) const
{
   Eigen::VectorXd transformed_coords =
      (basis_transform * coords) + origin_translation;
   return -(gaussian(transformed_coords, origin) +
            gamma*gaussian(transformed_coords, mu)) + constant_term;
}

double Gaussian_potential::partial(const Eigen::VectorXd &coords, int i) const
{
   Eigen::VectorXd transformed_coords =
      (basis_transform * coords) + origin_translation;

   // Vector of partials in original basis
   Eigen::VectorXd orig_partials(n_fields);
   for (std::size_t j = 0; j < n_fields; ++j) {
      orig_partials(j) = (-1.0)*(
         gaussian_deriv(transformed_coords, origin, j) +
         gamma*gaussian_deriv(transformed_coords, mu, j));
   }

   // Transform these to the current basis
   Eigen::VectorXd partials = (basis_transform.transpose()) * orig_partials;
   return partials(i);
}

double Gaussian_potential::partial(
   const Eigen::VectorXd& coords, int i, int j) const
{
   Eigen::VectorXd transformed_coords =
      (basis_transform * coords) + origin_translation;

   double partial2 = 0;

   for (std::size_t k = 0; k < n_fields; ++k) {
      for (std::size_t l = 0; l < n_fields; ++l) {
         double orig_partial = -(
            gaussian_deriv2(transformed_coords, origin, k, l) +
            gamma*gaussian_deriv2(transformed_coords, mu, k, l));
         partial2 += orig_partial *
            basis_transform.transpose()(i, k)*basis_transform.transpose()(j,l);
      }
   }

   return partial2;
}


void Gaussian_potential::translate_origin(const Eigen::VectorXd& translation)
{
   origin_translation = translation;
}

void Gaussian_potential::apply_basis_change(const Eigen::MatrixXd &new_basis)
{
   // Note we use the inverse transform on incoming coordinates!
   basis_transform = basis_transform * (new_basis.transpose());
}

void Gaussian_potential::add_constant_term(double constant)
{
   constant_term += constant;
}

std::size_t Gaussian_potential::get_number_of_fields() const
{
   return n_fields;
}

} // namespace BubbleProfiler
