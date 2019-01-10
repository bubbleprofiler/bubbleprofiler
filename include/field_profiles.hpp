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

#ifndef BUBBLEPROFILER_FIELD_PROFILES_HPP_INCLUDED
#define BUBBLEPROFILER_FIELD_PROFILES_HPP_INCLUDED

/**
 * @file field_profiles.hpp
 * @brief contains the definition of the Field_profiles clas
 */

#include "gsl_interpolator.hpp"

#include <Eigen/Core>

#include <vector>

namespace BubbleProfiler {

/**
 * @class Field_profiles
 * @brief Discretized set of field profiles
 *
 * This class is a container for discretized field profiles, representing
 * a solution for or correction to the bubble profile. Methods are
 * provided to interrogate the dimensionality, size and discretization
 * parameters of the profile set. It is also possible to evaluate a point
 * in field space, either directly or via interpolation for off-grid points.
 * Finally, it is possible to add one field profile to another. This is
 * used to implement iterative corrections when calculating a profile solution.
 *
 * It should be noted that the profile data is assumed to be spherically
 * symmetric.
 */
class Field_profiles {
public:
   Field_profiles() = default;

   /**
    * @brief Create field profiles at given coordinates from data matrix
    * @param gridpoints_ Vector containing radial coordinate values
    * @param profiles_ Matrix containing field data
    * @param interpolation_points_fraction_ fraction of grid points to use
    *                                       for interpolation
    */
   Field_profiles(const Eigen::VectorXd& gridpoints_,
                  const Eigen::MatrixXd& profiles_,
                  double interpolation_points_fraction_ = 1.0);

   /**
    * @brief Create empty field profiles with specified parameters.
    * @param num_fields_ number of scalar fields
    * @param domain_start_ domain start (minimum radius)
    * @param domain_end_ domain end (maximum radius)
    * @param num_steps_ number of grid points per dimension
    * @param interpolation_points_fraction_ fraction of grid points to use
    *                                       for interpolation
    */
   Field_profiles(int num_fields_, double domain_start_,
                  double domain_end_, int num_steps_,
                  double interpolation_points_fraction_ = 1.0);

   /**
    * @brief Create field profiles from data matrix
    * @param profiles_ matrix containing field data
    * @param domain_start_ domain start (minimum radius)
    * @param domain_end_ domain end (maximum radius)
    * @param interpolation_points_fraction_ fraction of grid points to use
    *                                       for interpolation
    */
   Field_profiles(const Eigen::MatrixXd& profiles_,
                  double domain_start_, double domain_end_,
                  double interpolation_points_fraction_ = 1.0);

   ~Field_profiles() = default;
   Field_profiles(const Field_profiles&) = default;
   Field_profiles& operator=(const Field_profiles&) = default;
   Field_profiles(Field_profiles&&) = default;
   Field_profiles& operator=(Field_profiles&&) = default;

   void set_number_of_dimensions(int d) { n_spatial_dims = d; }
   void set_interpolation_points_fraction(double f);
   void set_field_profile(int i, const Eigen::VectorXd& p);

   int get_number_of_fields() const { return profiles.cols(); }
   int get_number_of_dimensions() const { return n_spatial_dims; }

   int get_number_of_grid_points() const { return gridpoints.size(); }
   double get_domain_start() const { return gridpoints(0); }
   double get_domain_end() const { return gridpoints(gridpoints.size() - 1); }

   /**
    * @brief Get a vector of the grid point coordinates.
    * @note these are the same for each dimension in the profile set.
    * @return Vector of grid point coordinates.
    */
   const Eigen::VectorXd& get_spatial_grid() const { return gridpoints; }

   /**
    * @brief Get the field profile data in matrix form.
    * @return field profile data matrix
    */
   const Eigen::MatrixXd& get_field_profiles() const { return profiles; }

   /**
    * @brief Get all field values at a given radius.
    * @param rho target radius
    * @return vector of field values
    */
   Eigen::VectorXd evaluate_at(double rho) const;

   /**
    * @brief Get a specific field value at a given radius.
    * @param field target field index
    * @param rho target radius
    * @return value of target field
    */
   double evaluate_at(int field, double rho) const;

   /**
    * @brief Take the radial derivative of a field at a given radius.
    * @param field target field index
    * @param order order of the derivative (1 or 2)
    * @param rho target radius
    * @return radial derivative of the field at target coordinate
    */
   double derivative_at(int field, int order, double rho) const;

private:
   Eigen::VectorXd gridpoints{};
   Eigen::MatrixXd profiles{};
   double n_spatial_dims{3};
   bool must_update_coords{true};
   double interpolation_points_fraction{1.0};
   Eigen::VectorXd interpolation_grid_values{};
   Eigen::MatrixXd interpolation_field_values{};
   std::vector<GSL_interpolator> splines{};

   void initialize_interpolation_grid_values();
   void initialize_interpolation_field_values(int);
   void initialize_interpolation_field_values();
   void build_spline_for_field(int);
   void build_splines();
};

} // namespace BubbleProfiler

#endif
