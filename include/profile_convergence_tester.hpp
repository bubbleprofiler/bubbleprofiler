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

#ifndef BUBBLEPROFILER_PROFILE_CONVERGENCE_TESTER_HPP_INCLUDED
#define BUBBLEPROFILER_PROFILE_CONVERGENCE_TESTER_HPP_INCLUDED

namespace BubbleProfiler {

class Potential;
class Field_profiles;

/*!
 * Abstract interface for a class used to determine whether a candidate
 * bubble solution has converged. Implementations are assumed to be stateful,
 * tracking relevant information about successive iterations until the
 * relevant convergence criterea are met.
 *
 * NOTE: Convergence testers may supply a maximum iteration count. However,
 * it is NOT guaranteed that implementations will return is_converged:TRUE when
 * this is exceeded. It is expected that instead, client classes will make use of
 * this value, and possibly provide a warning when the max count is exceeded
 * before the solution converges.
 */
class Profile_convergence_tester {
public:
   virtual ~Profile_convergence_tester() = default;
   virtual int get_max_iterations() const = 0;
   virtual bool is_converged(const Potential&, const Field_profiles&) = 0;

   //! Re-initialize the convergence tester.
   virtual void restart() = 0;
};

} // namespace BubbleProfiler

#endif
