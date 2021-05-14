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
    @file 
    @brief Functions for calculating action for one-dimensional potential
    parameterized by \f$E\f$ and \f$\alpha\f$. 
    
    The potential is 
    \f[
    V(\phi) = E\left[\frac{-4 \alpha +3}{2} \phi^2 - \phi ^3 + \alpha \phi ^4\right].
    \f]
    
    We must construct potentials etc for particular \f$E\f$ and \f$\alpha\f$ and
    pass them to the `Shooting` class.
*/

#include "action.hpp"
#include "shooting.hpp"

extern "C" double action(const double E,
                         const double alpha,
                         const int dim,
                         BubbleProfiler::Shooting_settings settings) {
    /**
        @returns Action, \f$S\f$
        @param E \f$E\f$
        @param alpha \f$\alpha\f$
        @param dim Number of space-time dimensions, \f$d\f$
        @param settings Settings for solver
    */
    using namespace BubbleProfiler;
    const double false_min = 0.;
    const double true_min = 1.;
    const double barrier = 0.75 / alpha - 1.;

    const auto potential = [E, alpha](double phi) {
      return -E * (-alpha * phi * phi * phi * phi
                   + phi * phi * phi
                   + 0.5 * (4. * alpha - 3.) * phi * phi);
    };

    const auto potential_first = [E, alpha](double phi) {
      return -E * (-4. * alpha * phi * phi * phi
                   + 3. * phi * phi
                   + (4. * alpha - 3.) * phi);
    };

    const auto potential_second = [E, alpha](double phi) {
      return -E * (-12. * alpha * phi * phi
                   + 6. * phi + (4. * alpha - 3.));
    };

    Shooting one_dim;
    one_dim.set_bisection_precision_bits(settings.shoot_bisect_bits);
    one_dim.set_shooting_abs_tol(settings.shoot_ode_abs);
    one_dim.set_shooting_rel_tol(settings.shoot_ode_rel);
    one_dim.set_drho_frac(settings.drho_frac);
    one_dim.set_bisect_lambda_max(settings.bisect_lambda_max);
    one_dim.set_max_iterations(settings.iter_max);
    one_dim.set_max_periods(settings.periods_max);
    one_dim.set_f_y_max(settings.f_y_max);
    one_dim.set_f_y_min(settings.f_y_min);
    one_dim.set_y_max(settings.y_max);

    one_dim.solve(potential, potential_first, potential_second,
                  false_min, true_min, barrier,
                  dim, Shooting::Solver_options::Compute_action);

    return one_dim.get_euclidean_action();
}
