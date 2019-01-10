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

#include "log_message.hpp"
#include "error.hpp"

#include <ctime>
#include <sstream>

namespace BubbleProfiler {

namespace logging {

std::string to_log_level_string(Log_level level)
{
   switch(level) {
   case Log_level::Trace: return "trace";
   case Log_level::Debug: return "debug";
   case Log_level::Info: return "info";
   case Log_level::Warning: return "warning";
   case Log_level::Error: return "error";
   case Log_level::Fatal: return "fatal";
   default: throw Error("unrecognized log level");
   }

   throw Error("unrecognized log level");
}

std::string Default_log_message::get_formatted_time() const
{
   const auto time = Clock::to_time_t(log_time);

   char formatted[64];
   strftime(formatted, sizeof(formatted), "[%F %T %Z]", std::localtime(&time));

   return formatted;
}

std::string Default_log_message::get_log_entry() const
{
   std::stringstream entry;

   entry << get_formatted_time()
         << " [" << to_log_level_string(level) << "]\t"
         << message << '\n';

   return entry.str();
}

} // namespace logging

} // namespace BubbleProfiler
