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

#ifndef BUBBLEPROFILER_PERTURBATIVE_PROFILER_HPP_INCLUDED
#define BUBBLEPROFILER_PERTURBATIVE_PROFILER_HPP_INCLUDED

#include "default_integration_policy.hpp"
#include "generic_perturbative_profiler.hpp"

namespace BubbleProfiler {

using Euler_perturbative_profiler = Perturbative_profiler<Constant_step_size_euler>;

using RK4_perturbative_profiler = Perturbative_profiler<Constant_step_size_RK4>;

using CRK4_perturbative_profiler = Perturbative_profiler<Controlled_step_size_RK4>;

using CRKD5_perturbative_profiler = Perturbative_profiler<Controlled_step_size_RKD5>;

} // namespace BubbleProfiler

#endif
