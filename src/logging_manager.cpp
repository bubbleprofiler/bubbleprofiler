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

#include "logging_manager.hpp"
#include "error.hpp"

#include <iostream>

namespace BubbleProfiler {

namespace logging {

Logging_manager& Logging_manager::get_manager()
{
   static Logging_manager instance;

   return instance;
}

bool Logging_manager::above_min_level(Log_level level) const
{
   bool result = true;
   switch (min_level) {
   case Log_level::Trace: {
      break;
   }
   case Log_level::Debug: {
      if (level == Log_level::Trace) {
         result = false;
      }
      break;
   }
   case Log_level::Info: {
      if (level == Log_level::Trace ||
          level == Log_level::Debug) {
         result = false;
      }
      break;
   }
   case Log_level::Warning: {
      if (level == Log_level::Trace ||
          level == Log_level::Debug ||
          level == Log_level::Info) {
         result = false;
      }
      break;
   }
   case Log_level::Error: {
      if (level != Log_level::Error &&
          level != Log_level::Fatal) {
         result = false;
      }
      break;
   }
   case Log_level::Fatal: {
      if (level != Log_level::Fatal) {
         result = false;
      }
      break;
   }
   default: {
      throw Error("unrecognized log level");
   }
   }

   return result;
}

void Logging_manager::log(const Log_message& message)
{
   if (logging_enabled && above_min_level(message.get_log_level())) {
      std::cerr << message.get_log_entry();
   }
}

} // namespace logging

} // namespace BubbleProfiler
