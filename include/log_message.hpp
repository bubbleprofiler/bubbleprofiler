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

#ifndef BUBBLEPROFILER_LOG_MESSAGE_HPP_INCLUDED
#define BUBBLEPROFILER_LOG_MESSAGE_HPP_INCLUDED

#include <chrono>
#include <string>

namespace BubbleProfiler {

namespace logging {

enum class Log_level { Trace, Debug, Info, Warning, Error, Fatal };

std::string to_log_level_string(Log_level);

class Log_message {
public:
   virtual ~Log_message() = default;

   virtual Log_level get_log_level() const = 0;
   virtual std::string get_log_entry() const = 0;
};

class Default_log_message : public Log_message {
public:
   Default_log_message(Log_level level_, const std::string& message_)
      : level(level_), message(message_) {}
   virtual ~Default_log_message() = default;

   Log_level get_log_level() const override { return level; }
   std::string get_log_entry() const override;

private:
   using Clock = std::chrono::system_clock;
   using Time = std::chrono::time_point<Clock>;

   Log_level level{Log_level::Trace};
   std::string message{};
   Time log_time{Clock::now()};

   std::string get_formatted_time() const;
};

} // namespace logging

} // namespace BubbleProfiler

#endif
