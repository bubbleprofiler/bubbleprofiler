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

#ifndef BUBBLEPROFILER_INSTREAM_PROFILE_GUESSER_HPP_INCLUDED
#define BUBBLEPROFILER_INSTREAM_PROFILE_GUESSER_HPP_INCLUDED

/**
 * @file instream_profile_guesser.hpp
 * @brief contains the definition of the Instream_profile_guesser class
 */

#include "basic_logger.hpp"
#include "profile_guesser.hpp"

#include <Eigen/Core>

#include <istream>
#include <string>
#include <vector>

namespace BubbleProfiler {

/**
 * @class Instream_profile_guesser
 * @brief Provides guess for field profiles read from an input stream
 */
class Instream_profile_guesser: public Profile_guesser {
public:
   /**
    * @brief reads a profile guess from the given input stream
    * @param input_stream_ the stream to read from
    */
   explicit Instream_profile_guesser(std::istream&);

   /**
    * @brief reads a profile guess from the given input stream
    * @param input_stream_ the stream to read from
    * @param delimiters_ a set of characters to use as field delimiters
    */
   Instream_profile_guesser(std::istream&, const std::string&);

   /**
    * @brief reads a profile guess from the given input stream
    * @param input_stream_ the stream to read from
    * @param delimiters_ a set of characters to use as field delimiters
    * @param line_comment_start_ characters indicating the start of a comment
    */
   Instream_profile_guesser(std::istream&, const std::string&,
                            const std::string&);
   virtual ~Instream_profile_guesser() = default;
   Instream_profile_guesser(const Instream_profile_guesser&) = default;
   Instream_profile_guesser& operator=(const Instream_profile_guesser&) = default;
   Instream_profile_guesser(Instream_profile_guesser&&) = default;
   Instream_profile_guesser& operator=(Instream_profile_guesser&&) = default;

   virtual Field_profiles get_profile_guess(
      const Potential&, const Eigen::VectorXd&, int,
      double, double, double, double) override;

   /**
    * @brief returns the number of grid points read from the stream
    * @return the number of grid points read
    */
   std::size_t get_number_of_grid_points() const { return spatial_grid.size(); }

   /**
    * @breif returns the number of fields for which values were read
    * @return the number of fields found in the input stream
    */
   std::size_t get_number_of_fields() const { return field_values.size(); }

   /**
    * @brief returns the values of the spacetime coordinate read from the stream
    * @return the grid of discretized coordinate values
    */
   Eigen::VectorXd get_spatial_grid() const;

   /**
    * @brief returns the values of the fields read from the stream
    * @return the field values at each grid point
    */
   Eigen::MatrixXd get_field_values() const;

   /**
    * @brief sets the spatial tolerance used when comparing grid points
    *
    * When using the read-in data to construct a set of guessed profiles,
    * the boundaries of the stored domain must agree to within the
    * specified tolerance with the values of \c domain_start
    * and \c domain_end provided to get_profile_guess() .
    *
    * @param t the tolerance to be used
    * @sa get_profile_guess()
    */
   void set_spatial_tolerance(double t) { tolerance = t; }

   /**
    * @brief set the characters to use as delimiters between fields
    * @param d a string containing the characters to be used as delimiters
    */
   void set_delimiters(const std::string& d) { delimiters = d; };

   /**
    * @brief sets the character sequence used to indicate a comment
    *
    * Lines in the provided input stream may contain a sequence of
    * one or more characters indicating a comment. All content in the
    * line following a comment sequence is ignored by the profile
    * guesser.
    *
    * @param c the string used to indicate the start of a comment
    */
   void set_line_comment(const std::string& c) { line_comment_start = c; };

   /**
    * @brief read a profile guess from the given stream
    * @param input_stream the stream to read from
    */
   void read_profiles(std::istream&);

private:
   std::string delimiters{" \t"};
   std::string line_comment_start{"#"};
   bool has_n_fields{false};
   std::size_t n_fields{0};
   double tolerance{1.e-6};
   std::vector<double> spatial_grid{};
   std::vector<std::vector<double> > field_values{};
   logging::Basic_logger logger{};

   void reset();
   void strip_comments(std::string&) const;
   void process_line(const std::string&);
};

} // namespace BubbleProfiler

#endif
