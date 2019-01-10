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

#ifndef BUBBLEPROFILER_PROFILE_GUESSER_HPP_INCLUDED
#define BUBBLEPROFILER_PROFILE_GUESSER_HPP_INCLUDED

#include "field_profiles.hpp"

#include <Eigen/Core>

namespace BubbleProfiler {

class Potential;

//! Abstract class to represent ansatz generators
class Profile_guesser {
public:
   virtual ~Profile_guesser() = default;

   /*!
    * @brief Calculate an initial guess for the bubble profile
    *
    * The calculated profiles should be evaluated at the requested number of
    * points, \c n_grid_points , between the minimum and maximum radial
    * coordinate values, \c domain_start and \c domain_end , respectively.
    * If a negative value is given for either of these, appropriate values
    * are to be guessed based on the given potential and vacuum location.
    *
    * @note The false vacuum is assumed to be at the origin in field space.
    *
    * @param potential the potential for which the profiles are to be computed
    * @param true_vacuum the location of the true vacuum in field space
    * @param n_dimensions the number of space-time dimensions
    * @param domain_start the minimum value of the radial coordinate at which
    *                     the profiles are evaluated.  If negative, the function
    *                     should guess an appropriate value.
    * @param domain_end the maximum value of the radial coordinate at which
    *                   the profiles are evaluated.  If negative, the function
    *                   should guess an appropriate value.
    * @param initial_step_size the initial step size to be used in constructing
    *                          the discretized solution
    * @param interpolation_points_fraction fraction of grid points to use
    *                                      for interpolation
    * @return an initial guess for the field profiles
    */
   virtual Field_profiles get_profile_guess(
      const Potential& potential, const Eigen::VectorXd& true_vacuum, int n_dimensions,
      double domain_start, double domain_end, double initial_step_size,
      double interpolation_points_fraction_) = 0;
};

} // namespace BubbleProfiler

#endif
