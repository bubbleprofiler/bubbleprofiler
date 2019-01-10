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

#ifndef BUBBLEPROFILER_LOGGING_MANAGER_HPP_INCLUDED
#define BUBBLEPROFILER_LOGGING_MANAGER_HPP_INCLUDED

#include "log_message.hpp"

namespace BubbleProfiler {

namespace logging {

class Logging_manager {
public:
   Logging_manager(const Logging_manager&) = delete;
   Logging_manager(Logging_manager&&) = delete;
   Logging_manager& operator=(const Logging_manager&) = delete;
   Logging_manager& operator=(Logging_manager&&) = delete;

   static Logging_manager& get_manager();

   void set_logging_enabled(bool flag) { logging_enabled = flag; }
   void set_minimum_log_level(Log_level level) { min_level = level; }

   void log(const Log_message& message);

protected:
   Logging_manager() = default;
   ~Logging_manager() = default;

private:
   bool logging_enabled{true};
   Log_level min_level{Log_level::Warning};

   bool above_min_level(Log_level) const;
};

} // namespace logging

} // namespace BubbleProfiler

#endif
