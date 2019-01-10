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
 * @file scale.cpp
 * @brief checks behaviour of action under changes in units
 *
 * This example program checks the behaviour of action under
 * changes of units in the case of \f$d = 4\f$ dimensions,
 * in which action is dimensionless.
 *
 * Usage:
 * @code
 *  ./scale <E> <alpha> <scale>
 * @endcode
 */

#include "action.hpp"
#include "error.hpp"
#include "logging_manager.hpp"
#include "shooting.hpp"

#include <iostream>

namespace BubbleProfiler {

struct Scale_options {
   double E{1.};
   double alpha{0.6};
   double scale{2.};
   bool verbose{false};
};

void print_usage()
{
   std::cout <<
      "Usage: scale [OPTION] <E> <alpha> <dim>\n\n"
      "Print action under a change of units.\n\n"
      "Example: scale 1.0 0.6 2.0\n\n"
      "Options:\n"
      "  -h, --help                      print this help message\n"
      "  -v, --verbose                   produce verbose output"
             << std::endl;
}

bool starts_with(const std::string& option, const std::string& prefix)
{
   return !option.compare(0, prefix.size(), prefix);
}

void parse_positional_args(const std::vector<std::string>& args,
                           Scale_options& options)
{
   const auto n_positional_args = args.size();

   if (n_positional_args == 0) {
      throw Setup_error(
         "value of E must be provided");
   } else if (n_positional_args == 1) {
      throw Setup_error(
         "value of alpha must be provided");
   } else if (n_positional_args == 2) {
      throw Setup_error(
         "scale factor must be provided");
   } else if (n_positional_args > 3) {
      throw Setup_error(
         "unrecognized command line argument '"
         + args[3] + "'");
   }

   try {
      options.E = std::stod(args[0]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value for E, '" + args[0] + "'");
   }

   try {
      options.alpha = std::stod(args[1]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value for alpha, '" + args[1] + "'");
   }

   try {
      options.scale = std::stod(args[2]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value for scale factor '" + args[2] + "'");
   }
}

Scale_options parse_cmd_line_args(int argc, const char* argv[])
{
   Scale_options options;

   std::vector<std::string> positional_args;
   int i = 1;
   bool finished_optionals{false};
   while (i < argc) {
      const std::string opt(argv[i++]);

      if (opt == "--") {
         finished_optionals = true;
         continue;
      }

      if (!starts_with(opt, "-") || finished_optionals) {
         positional_args.push_back(opt);
         continue;
      }

      if (opt == "-h" || opt == "--help") {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      if (opt == "-v" || opt == "--verbose") {
         options.verbose = true;
         continue;
      }

      throw Setup_error(
         "unrecognized command line argument '" + opt + "'");
   }

   parse_positional_args(positional_args, options);

   return options;
}

} // namespace BubbleProfiler

int main(int argc, const char* argv[])
{
   using namespace BubbleProfiler;

   if (argc == 1) {
      std::cerr << "Usage: scale [OPTION] <E> <alpha> <scale>\n"
                << "Try 'scale --help' for more information."
                << std::endl;
      return EXIT_FAILURE;
   }

   int exit_code = 0;
   try {
      Scale_options options = parse_cmd_line_args(argc, argv);

      if (options.verbose) {
         logging::Logging_manager::get_manager().set_minimum_log_level(
            logging::Log_level::Trace);
      }

      const double E = options.E;
      const double alpha = options.alpha;
      const double scale = options.scale;

      const double false_min = 0.;
      const double true_min = scale;
      const double barrier = (0.75 / alpha - 1.) * scale;

      const auto potential = [E, alpha, scale](double phi) {
         phi /= scale;
         return -pow(scale, 4) * E * (-alpha * phi * phi * phi * phi
                                      + phi * phi * phi
                                      + 0.5 * (4. * alpha - 3.) * phi * phi);
      };

      const auto potential_first = [E, alpha, scale](double phi) {
         phi /= scale;
         return -pow(scale, 3) * E * (-4. * alpha * phi * phi * phi
                                      + 3. * phi * phi
                                      + (4. * alpha - 3.) * phi);
      };

      const auto potential_second = [E, alpha, scale](double phi) {
         phi /= scale;
         return -pow(scale, 2) * E * (-12. * alpha * phi * phi
                                      + 6. * phi + (4. * alpha - 3.));
      };

      Shooting one_dim;
      one_dim.solve(potential, potential_first, potential_second,
                    false_min, true_min, barrier,
                    4, Shooting::Solver_options::Compute_action);

      std::cout << one_dim.get_euclidean_action() << std::endl;
   } catch (const Error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      exit_code = EXIT_FAILURE;
   } catch (...) {
      std::cerr << "Error: unrecognized error occurred."
                << std::endl;
      exit_code = EXIT_FAILURE;
   }

   return exit_code;
}
