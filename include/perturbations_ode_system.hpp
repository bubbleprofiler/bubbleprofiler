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

#ifndef BUBBLEPROFILER_PERTURBATIONS_ODE_SYSTEM_HPP_INCLUDED
#define BUBBLEPROFILER_PERTURBATIONS_ODE_SYSTEM_HPP_INCLUDED

#include <Eigen/Core>

namespace BubbleProfiler {

class Field_profiles;
class Potential;

/*!
 * @class Perturbations_ODE_system
 * @brief Class implementing the system of ODEs obeyed by profile perturbations
 *
 * This class provides for the calculation of the first-order system of
 * ODEs obeyed by the perturbations to the field profiles, and their first
 * derivatives, that arise when using the perturbative method.
 */
class Perturbations_ODE_system {
public:
   /*!
    * @brief Constructs the system of ODEs for the given potential
    * @param[in] potential_ the potential to calculate profile perturbations for
    * @param[in] profiles_ the current estimate for the bubble profile for which
    *                      perturbations are to be computed
    * @param[in] n_spacetime_dimensions_ the number of spacetime dimensions
    */
   Perturbations_ODE_system(Potential& potential_, Field_profiles& profiles_,
                            int n_spacetime_dimensions_);

   /*!
    * @brief Calculate the derivatives of the perturbations
    *
    * The state vector \c eps is assumed to contain the values of the
    * perturbations in the first \c n_fields entries, and the values
    * of their first derivatives in the remaining \c n_fields entries.
    *
    * @param[in] eps the values of the perturbations and their derivatives
    * @param[out] depsdr the values of the derivatives of the perturbations
    * @param[in] r the value of the radial coordinate
    */
   void operator()(const Eigen::VectorXd& eps, Eigen::VectorXd& depsdr,
                   double r) const;

private:
   /*! Number of fields in the potential */
   int n_fields{1};
   /*! Number of spacetime dimensions */
   int n_spacetime_dimensions{3};
   /*! Potential for which profile is calculated */
   Potential* potential{nullptr};
   /*! Current estimate for bubble profile to be perturbed */
   Field_profiles* profiles{nullptr};

   bool is_finite(const Eigen::VectorXd&) const;
   Eigen::VectorXd calculate_inhomogeneities(double) const;
   Eigen::MatrixXd calculate_mass_matrix(double) const;
};

} // namespace BubbleProfiler

#endif
