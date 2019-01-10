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

/**
 * @file instream_profile_guesser.cpp
 * @brief contains the implementation of the Instream_profile_guesser class
 */

#include "instream_profile_guesser.hpp"
#include "error.hpp"
#include "field_profiles.hpp"
#include "math_wrappers.hpp"
#include "potential.hpp"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

namespace BubbleProfiler {

Instream_profile_guesser::Instream_profile_guesser(
   std::istream& input_stream_)
{
   read_profiles(input_stream_);
}

Instream_profile_guesser::Instream_profile_guesser(
   std::istream& input_stream_, const std::string& delimiters_)
   : delimiters(delimiters_)
{
   read_profiles(input_stream_);
}

Instream_profile_guesser::Instream_profile_guesser(
   std::istream& input_stream_, const std::string& delimiters_,
   const std::string& line_comment_start_)
   : delimiters(delimiters_)
   , line_comment_start(line_comment_start_)
{
   read_profiles(input_stream_);
}

void Instream_profile_guesser::reset()
{
   has_n_fields = false;
   spatial_grid.clear();
   field_values.clear();
}

void Instream_profile_guesser::read_profiles(std::istream& input_stream)
{
   reset();

   if (!input_stream.good()) {
      throw IO_error("cannot read from provided stream");
   }

   logger.log_message(logging::Log_level::Trace, "Reading from input stream");

   std::string line;
   has_n_fields = false;
   while (std::getline(input_stream, line)) {
      strip_comments(line);
      boost::trim(line);
      if (line.empty()) {
         continue;
      }
      process_line(line);
   }

   if (spatial_grid.empty()) {
      throw IO_error("no valid grid points found");
   }

   logger.log_message(logging::Log_level::Trace,
                      "Read " + std::to_string(spatial_grid.size())
                      + " points from input stream");
   logger.log_message(logging::Log_level::Trace,
                      "Found values for " + std::to_string(n_fields)
                      + " fields");
}

void Instream_profile_guesser::strip_comments(std::string& line) const
{
   const auto pos = line.find(line_comment_start);
   if (pos != std::string::npos) {
      line = line.substr(0, pos);
   }
}

void Instream_profile_guesser::process_line(const std::string& line)
{
   std::vector<std::string> fields;
   boost::split(fields, line, boost::is_any_of(delimiters),
                boost::token_compress_on);

   if (fields.empty()) {
      throw IO_error("no valid fields found in line: " + line);
   }

   const std::size_t current_n_fields = fields.size() - 1;
   if (current_n_fields == 0) {
      throw IO_error("no field values given in line: " + line);
   }

   if (has_n_fields) {
      if (current_n_fields != n_fields) {
         throw IO_error("incorrect number of fields in line: " + line);
      }
   } else {
      n_fields = current_n_fields;
      has_n_fields = true;
      field_values = std::vector<std::vector<double> >(n_fields);
   }

   try {
      const double grid_value = boost::lexical_cast<double>(fields[0]);
      spatial_grid.push_back(grid_value);

      for (std::size_t i = 0; i < n_fields; ++i) {
         const double field_value = boost::lexical_cast<double>(fields[i + 1]);
         field_values[i].push_back(field_value);
      }
   } catch (const boost::bad_lexical_cast& error) {
      throw IO_error("unable to read input data in line: " + line);
   }
}

Field_profiles Instream_profile_guesser::get_profile_guess(
   const Potential& potential,
   const Eigen::VectorXd& /* true_vacuum */,
   int n_dimensions, double domain_start, double domain_end,
   double initial_step_size, double interpolation_points_fraction)
{
   const std::size_t required_n_fields = potential.get_number_of_fields();
   if (required_n_fields != n_fields || !has_n_fields) {
      throw Setup_error(
         "number of fields in guess does not match number in potential");
   }

   const double stored_domain_start = spatial_grid.front();
   const double stored_domain_end = spatial_grid.back();

   if (domain_start >= 0.) {
      if (Abs(stored_domain_start - domain_start) > tolerance) {
         throw Setup_error(
            "stored domain start does not match requested domain start");
      }
   }

   if (domain_end >= 0.) {
      if (Abs(stored_domain_end - domain_end) > tolerance) {
         throw Setup_error(
            "stored domain end does not match requested domain end");
      }
   }

   // If definite values for the domain start and end, and the grid step
   // size are given, require that the calculated grid matches that read
   // from the stream.  This is imposed for simplicity, but could be
   // removed if the guessed profile is constructed by interpolating the
   // read in values at the desired grid points.
   if (domain_start >= 0. && domain_end >= 0.) {
      const std::size_t requested_n_grid_points = 1 + static_cast<std::size_t>(
         std::ceil((domain_end - domain_start) / initial_step_size));
      if (requested_n_grid_points != spatial_grid.size()) {
         throw Setup_error(
            "number of stored grid points does not match desired grid size");
      }
      const double rho_incr = (domain_end - domain_start) /
         (requested_n_grid_points - 1.);
      for (std::size_t i = 0; i < requested_n_grid_points; ++i) {
         const double current_rho = domain_start + i * rho_incr;
         if (Abs(spatial_grid[i] - current_rho) > tolerance) {
            throw Setup_error("spatial coordinate at grid point "
                              + std::to_string(i)
                              + " does not match required value");
         }
      }
   }

   Eigen::VectorXd grid(get_spatial_grid());
   Eigen::MatrixXd profiles(get_field_values());

   Field_profiles profile_guess(grid, profiles, interpolation_points_fraction);
   profile_guess.set_number_of_dimensions(n_dimensions);

   return profile_guess;
}

Eigen::VectorXd Instream_profile_guesser::get_spatial_grid() const
{
   Eigen::VectorXd grid(
      Eigen::VectorXd::Map(spatial_grid.data(), spatial_grid.size()));
   return grid;
}

Eigen::MatrixXd Instream_profile_guesser::get_field_values() const
{
   const auto n_grid_points = spatial_grid.size();
   Eigen::MatrixXd profiles(n_grid_points, n_fields);
   for (std::size_t i = 0; i < n_fields; ++i) {
      profiles.col(i) =
         Eigen::VectorXd::Map(field_values[i].data(), field_values[i].size());
   }
   return profiles;
}

} // namespace BubbleProfiler
