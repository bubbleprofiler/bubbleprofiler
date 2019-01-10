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

#ifndef BUBBLEPROFILER_DEFAULT_INTEGRATION_POLICY_HPP_INCLUDED
#define BUBBLEPROFILER_DEFAULT_INTEGRATION_POLICY_HPP_INCLUDED

#include "eigen_state_utils.hpp"
#include "integration_policy.hpp"

namespace BubbleProfiler {

using Constant_step_size_RK4 =
   Fixed_step_size_integrator<double, Eigen::VectorXd>;

using Constant_step_size_euler =
   Fixed_step_size_integrator<
   double, Eigen::VectorXd,
   boost::numeric::odeint::euler<
      Eigen::VectorXd, double, Eigen::VectorXd, double,
      State_algebra_dispatcher<Eigen::VectorXd>::algebra_type>
   >;

using Controlled_step_size_RK4 =
   Controlled_step_size_integrator<
   double, Eigen::VectorXd,
   boost::numeric::odeint::runge_kutta_cash_karp54<
      Eigen::VectorXd, double, Eigen::VectorXd, double,
      State_algebra_dispatcher<Eigen::VectorXd>::algebra_type>
   >;

using Controlled_step_size_RKD5 =
   Controlled_step_size_integrator<
   double, Eigen::VectorXd,
   boost::numeric::odeint::runge_kutta_dopri5<
      Eigen::VectorXd, double, Eigen::VectorXd, double,
      State_algebra_dispatcher<Eigen::VectorXd>::algebra_type>
   >;

} // namespace BubbleProfiler

#endif
