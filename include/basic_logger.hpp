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

#ifndef BUBBLEPROFILER_BASIC_LOGGER_HPP_INCLUDED
#define BUBBLEPROFILER_BASIC_LOGGER_HPP_INCLUDED

#include "log_message.hpp"
#include "logging_manager.hpp"

namespace BubbleProfiler {

namespace logging {

class Basic_logger {
public:
   void log_message(Log_level level, const std::string& msg) const;
};

} // namespace logging

} // namespace BubbleProfiler

#endif
