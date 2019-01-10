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
 * @file field_profiles.cpp
 * @brief contains the implementation of the Field_profiles class
 */

#include "field_profiles.hpp"
#include "error.hpp"

#include <cmath>

namespace BubbleProfiler {

Field_profiles::Field_profiles(const Eigen::VectorXd& gridpoints_,
                               const Eigen::MatrixXd& profiles_,
                               double interpolation_points_fraction_)
   : gridpoints(gridpoints_)
   , profiles(profiles_)
   , splines(profiles_.cols())
{
   if (interpolation_points_fraction_ <= 0. ||
       interpolation_points_fraction_ > 1.) {
      throw Setup_error(
         "fraction of points used in interpolation must be between 0 and 1");
   }

   build_splines();
}

Field_profiles::Field_profiles(int num_fields_, double domain_start_,
                               double domain_end_, int n_steps_,
                               double interpolation_points_fraction_)
   : Field_profiles(Eigen::VectorXd::LinSpaced(
                       n_steps_, domain_start_, domain_end_),
                    Eigen::MatrixXd::Zero(n_steps_, num_fields_),
                    interpolation_points_fraction_)
{
}

Field_profiles::Field_profiles(const Eigen::MatrixXd& profiles_,
                               double domain_start_, double domain_end_,
                               double interpolation_points_fraction_)
   : Field_profiles(Eigen::VectorXd::LinSpaced(
                       profiles_.rows(), domain_start_, domain_end_),
                    profiles_, interpolation_points_fraction_)
{
}

void Field_profiles::set_interpolation_points_fraction(double f)
{
   if (f <= 0. || f > 1.) {
      throw Setup_error(
         "fraction of points used in interpolation must be between 0 and 1");
   }

   interpolation_points_fraction = f;
   must_update_coords = true;
}

void Field_profiles::set_field_profile(int i, const Eigen::VectorXd& p)
{
   if (p.rows() != profiles.rows()) {
      throw Setup_error("Field_profiles::set_field_profiles: "
                        "number of grid points does not match");
   }

   profiles.col(i) = p;
   build_spline_for_field(i);
}

void Field_profiles::initialize_interpolation_grid_values()
{
   const int n_steps = gridpoints.size();
   const int n_knots = static_cast<int>(
      std::ceil(interpolation_points_fraction * n_steps));

   if (n_knots > n_steps) {
      throw Setup_error("Field_profiles::initialize_interpolation_grid_values: "
                        "number of knots cannot exceed number of grid points");
   }

   const int stride = (n_steps - 1) / (n_knots - 1);

   if (interpolation_grid_values.size() != n_knots) {
      interpolation_grid_values.resize(n_knots);
   }

   interpolation_grid_values(0) = gridpoints(0);
   for (int i = 1; i < n_knots - 1; ++i) {
      interpolation_grid_values(i) = gridpoints(i * stride);
   }
   interpolation_grid_values(n_knots - 1) = gridpoints(n_steps - 1);
}

void Field_profiles::initialize_interpolation_field_values(int field)
{
   const int n_steps = gridpoints.size();
   const int n_knots = interpolation_grid_values.size();
   const int stride = (n_steps - 1) / (n_knots - 1);

   interpolation_field_values(0, field) = profiles(0, field);
   for (int i = 1; i < n_knots - 1; ++i) {
      interpolation_field_values(i, field) = profiles(i * stride, field);
   }
   interpolation_field_values(n_knots - 1, field) =
      profiles(n_steps - 1, field);
}

void Field_profiles::initialize_interpolation_field_values()
{
   const int n_steps = gridpoints.size();
   const int n_knots = interpolation_grid_values.size();
   const int stride = (n_steps - 1) / (n_knots - 1);

   if (interpolation_field_values.cols() != profiles.cols() ||
       interpolation_field_values.rows() != n_knots) {
      interpolation_field_values.resize(n_knots, profiles.cols());
   }

   interpolation_field_values.row(0) = profiles.row(0);
   for (int i = 1; i < n_knots - 1; ++i) {
      interpolation_field_values.row(i) = profiles.row(i * stride);
   }
   interpolation_field_values.row(n_knots - 1) = profiles.row(n_steps - 1);
}

void Field_profiles::build_spline_for_field(int field)
{
   if (must_update_coords) {
      build_splines();
   } else {
      initialize_interpolation_field_values(field);

      splines[field].construct_interpolant(
         interpolation_grid_values,
         interpolation_field_values.col(field),
         interpolation_grid_values.size());
   }
}

void Field_profiles::build_splines()
{
   if (must_update_coords) {
      initialize_interpolation_grid_values();
      must_update_coords = false;
   }

   initialize_interpolation_field_values();

   const int n_fields = profiles.cols();
   for (int i = 0; i < n_fields; ++i) {
      splines[i].construct_interpolant(interpolation_grid_values,
                                       interpolation_field_values.col(i),
                                       interpolation_grid_values.size());
   }
}

double Field_profiles::evaluate_at(int field, double rho) const
{
   if (field < 0 || field >= profiles.cols()) {
      throw Out_of_bounds_error(
         field,
         "Field_profiles::evaluate_at: invalid field index "
         + std::to_string(field));
   }

   return splines[field].evaluate_at(rho);
}

Eigen::VectorXd Field_profiles::evaluate_at(double rho) const
{
   const int n_fields = profiles.cols();
   Eigen::VectorXd result(n_fields);
   for (int i = 0; i < n_fields; ++i) {
      result(i) = splines[i].evaluate_at(rho);
   }
   return result;
}


double Field_profiles::derivative_at(int field, int order, double rho) const
{
   if (order != 1 && order != 2) {
      throw Not_implemented_error(
         "Field_profiles::derivative_at: "
         "only first or second derivatives may be computed");
   }

   if (field < 0 || field >= profiles.cols()) {
      throw Out_of_bounds_error(
         field,
         "Field_profiles::evaluate_at: invalid field index "
         + std::to_string(field));
   }

   return order == 1 ? splines[field].evaluate_deriv_at(rho) :
      splines[field].evaluate_second_deriv_at(rho);
}

} // namespace BubbleProfiler
