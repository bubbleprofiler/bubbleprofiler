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

#ifndef BUBBLEPROFILER_RELATIVE_CONVERGENCE_TESTER_HPP_INCLUDED
#define BUBBLEPROFILER_RELATIVE_CONVERGENCE_TESTER_HPP_INCLUDED

#include "basic_logger.hpp"
#include "profile_convergence_tester.hpp"

#include <Eigen/Core>

#include <vector>

namespace BubbleProfiler {

/*!
 * This class implements a convergence criterion based on two factors:
 *
 * 1. Relative change in Euclidean action
 * 2. Relative change in field values at the origin
 *
 * The user supplies relative tolerances for both factors. Note that in general
 * there is more than one scalar field; in which case the field with the
 * greatest relative change between iterations is used to test convergence.
 */
class Relative_convergence_tester : public Profile_convergence_tester {
public:
   /*!
    * @brief Construct relative convergence tester with default tolerances.
    */
   Relative_convergence_tester() = default;

   /*!
    * @brief Construct relative convergence tester with the same tolerance value for
    * change in action and change in field values at origin.
    *
    * The maximum number of iterations is set to
    * @f[
       -10 \log_{10}(\epsilon),
    * @f]
    * with \f$\epsilon\f$ equal to the supplied value of \c tol .
    *
    * @param tol relative change tolerance
    */
   explicit Relative_convergence_tester(double tol);

   /*!
    * @brief Construct relative convergence tester with different tolerance values for
    * change in action and change in field values at origin.
    *
    * The maximum number of iterations is set to
    * @f[
       \max(-10 \log_{10}(\epsilon)),
    * @f]
    * where \f$\epsilon\f$ is the smaller of \c action_tol and \c fields_tol .
    *
    * @param action_tol relative change tolerance for the action
    * @param fields_tol relative change tolerance for field values at origin
    */
   Relative_convergence_tester(double action_tol, double fields_tol);
   Relative_convergence_tester(const Relative_convergence_tester&) = default;
   Relative_convergence_tester(Relative_convergence_tester&&) = default;
   virtual ~Relative_convergence_tester() = default;
   Relative_convergence_tester& operator=(const Relative_convergence_tester&) = default;
   Relative_convergence_tester& operator=(Relative_convergence_tester&&) = default;

   /*!
    * @brief Check whether the candidate bubble solution has converged.
    *
    * This method will return \c true when the change in both metrics (i.e. the action
    * and field values at the origin) fall below the specified thresholds.
    *
    * @param potential the potential for which the bubble profile is being calculated
    * @param profiles the current estimate for the bubble field profiles
    * @return \c true if both the action and field profiles have converged
    */
   virtual bool is_converged(
         const Potential& potential, const Field_profiles& profiles) override;

   /*!
    * @brief Return the maximum number of iterations allowed.
    *
    * If this count is exceeded, before convergence is reached (as
    * indicated by \c is_converged returning \c true), the client may want to
    * issue a warning or implement some other error handling strategy.
    *
    * @return maximum number of iterations allowed
    */
   virtual int get_max_iterations() const override { return max_iterations; }

   /*!
    * @brief Restart the convergence tester to the start of the iteration.
    */
   virtual void restart() override;

   /*!
    * @brief Set the maximum number of iterations allowed.
    * @param it the maximum allowed number of iterations
    */
   void set_max_iterations(int it) { max_iterations = it; }

private:
   int iteration_count{0};
   int max_iterations{10};
   double action_relative_tol{1.e-4};
   double field_vals_relative_tol{1.e-4};
   double old_action{0.};
   Eigen::VectorXd old_field_vals{};
   logging::Basic_logger logger{};

   int calculate_max_iterations() const;
   double relative_difference(double, double) const;
   bool check_action_converged(double) const;
   bool check_fields_converged(double, const Eigen::VectorXd&) const;
};

} // namespace BubbleProfiler

#endif
