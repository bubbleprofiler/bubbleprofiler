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

#ifndef BUBBLEPROFILER_RAII_GUARD_HPP_INCLUDED
#define BUBBLEPROFILER_RAII_GUARD_HPP_INCLUDED

#include <utility>

namespace BubbleProfiler {

/*!
 * @class RAII_guard
 * @brief Carries out provided clean-up actions at destruction
 */
template <typename F>
class RAII_guard {
public:
   RAII_guard(F f_) : clean_up(std::move(f_)) {}
   RAII_guard(const RAII_guard&) = delete;
   RAII_guard(RAII_guard&&) noexcept = default;
   ~RAII_guard() { clean_up(); }
   RAII_guard& operator=(const RAII_guard&) = delete;
   RAII_guard& operator=(RAII_guard&&) noexcept = default;
private:
   F clean_up;
};

template <typename F>
constexpr RAII_guard<F> make_raii_guard(F f)
{
   return RAII_guard<F>(std::move(f));
}

} // namespace BubbleProfiler

#endif
