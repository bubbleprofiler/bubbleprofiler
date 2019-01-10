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

#ifndef BUBBLEPROFILER_ROOT_FINDER_HPP_INCLUDED
#define BUBBLEPROFILER_ROOT_FINDER_HPP_INCLUDED

#include <functional>

namespace BubbleProfiler {

enum class Root_finder_status { SUCCESS, FAIL };

template <class Vector>
class Root_finder {
public:
   virtual ~Root_finder() = default;

   virtual const Vector& get_solution() const = 0;

   /*!
    * Run the root finding algorithm.
    * @tparam F Function data type
    * @param function Function to find root of
    * @param initial_guess Starting point for the algorithm
    * @return Status indicating whether the algorithm succeeded
    */
   template <class F>
   Root_finder_status find_root(F&& function, const Vector& guess);

protected:
   virtual Root_finder_status solve(
      const std::function<Vector(const Vector&)>&, const Vector&) = 0;
};

template <class Vector>
template <class F>
Root_finder_status Root_finder<Vector>::find_root(F&& f, const Vector& guess)
{
   std::function<Vector(const Vector&)> func{std::forward<F>(f)};
   return solve(func, guess);
}

} // namespace BubbleProfiler

#endif
