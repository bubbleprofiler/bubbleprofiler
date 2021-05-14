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

#ifndef BUBBLEPROFILER_ACTION_HPP_INCLUDED
#define BUBBLEPROFILER_ACTION_HPP_INCLUDED

/**
   @file
   @brief Functions for calculating action for one-dimensional potential
    parameterized by \f$E\f$ and \f$\alpha\f$

     The potential is
    \f[
    V(\phi) = E\left[\frac{-4 \alpha +3}{2} \phi^2 - \phi ^3 + \alpha \phi ^4\right].
    \f]

*/

namespace BubbleProfiler {

struct Shooting_settings {
   int shoot_bisect_bits{5};
   double shoot_ode_abs{1.e-4};
   double shoot_ode_rel{1.e-4};
   double drho_frac{1.e-3};
   double evolve_change_rel{1.e-2};
   double bisect_lambda_max{5};
   int iter_max{100000};
   double periods_max{1.e2};
   double f_y_max{1.e6};
   double f_y_min{1.e-3};
   double y_max{1.e1};
};

} // namespace BubbleProfiler

extern "C" double action(const double E,
                         const double alpha,
                         const int dim = 3,
                         BubbleProfiler::Shooting_settings settings
                         = BubbleProfiler::Shooting_settings());

#endif  // ACTION_HPP_
