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
    @file
    @brief Program that calculates action for a one-dimensional potential
    parameterized by \f$E\f$ and \f$\alpha\f$

    The potential is
    \f[
    V(\phi) = E\left[\frac{-4 \alpha +3}{2} \phi^2 - \phi ^3 + \alpha \phi ^4\right].
    \f]

    Usage:
    \code
    ./action <E> <alpha> <dim>
    \endcode
*/

#include "error.hpp"
#include "logging_manager.hpp"
#include "restricted_quartic_potential.hpp"
#include "shooting.hpp"

#include <iostream>
#include <string>
#include <vector>

namespace BubbleProfiler {

struct Action_options {
   double E{1.};
   double alpha{0.6};
   int dim{4};
   bool verbose{false};
};

void print_usage()
{
   std::cout <<
      "Usage: quartic [OPTION] <E> <alpha> <dim>\n\n"
      "Calculate action for one-dimensional quartic potential\n"
      "parameterized by E and alpha.\n\n"
      "Example: quartic 1.0 0.6 3\n\n"
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
                           Action_options& options)
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
         "number of dimensions must be provided");
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
      options.dim = std::stoi(args[2]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid number of dimensions '" + args[2] + "'");
   }
}

Action_options parse_cmd_line_args(int argc, const char* argv[])
{
   Action_options options;

   std::vector<std::string> positional_args;
   int i = 1;
   while (i < argc) {
      const std::string opt(argv[i++]);

      if (!starts_with(opt, "-")) {
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
      std::cerr << "Usage: quartic [OPTION] <E> <alpha> <dim>\n"
                << "Try 'quartic --help' for more information."
                << std::endl;
      return EXIT_FAILURE;
   }

   int exit_code = 0;
   try {
      Action_options options = parse_cmd_line_args(argc, argv);

      if (options.verbose) {
         logging::Logging_manager::get_manager().set_minimum_log_level(
            logging::Log_level::Trace);
      }

      Restricted_quartic_potential potential(options.alpha, options.E);
      const double false_min = potential.get_local_minimum_location();
      const double barrier = potential.get_local_maximum_location();
      const double true_min = potential.get_global_minimum_location();

      Shooting profiler;
      profiler.set_bisection_precision_bits(15);
      profiler.set_shooting_abs_tol(1.e-6);
      profiler.set_shooting_rel_tol(1.e-6);
      profiler.set_drho_frac(1.e-3);
      profiler.set_bisect_lambda_max(5);
      profiler.set_max_iterations(100000);
      profiler.set_max_periods(1.e2);
      profiler.set_f_y_max(1.e6);
      profiler.set_f_y_min(1.e-3);
      profiler.set_y_max(1.e1);

      profiler.solve(potential, false_min, true_min, barrier,
                    options.dim, Shooting::Solver_options::Compute_action);

      std::cout << "# E = " << options.E
                << ", alpha = " << options.alpha
                << ", dim = " << options.dim
                << ", action = " << profiler.get_euclidean_action()
                << std::endl;
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
