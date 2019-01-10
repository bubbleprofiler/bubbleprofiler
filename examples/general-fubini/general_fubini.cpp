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
    @brief Program that calculates action for a one-dimensional potential in case
    in which exact solution is known (the generalized Fubini potential).

     The potential is (see Eq. (6) of https://arxiv.org/abs/1412.3160)
    \f[
    V(\phi) = \frac{4 u m^2 (m-1)}{2 m+1}\phi^{(2m+1)/m} - 2 u v m^2 \phi^{(2m+2)/m}.
    \f]

    Usage:
    \code
    ./fubini <u> <v> <n>
    \endcode
*/

#include "error.hpp"
#include "field_profiles.hpp"
#include "generalized_fubini_observer.hpp"
#include "generalized_fubini_potential.hpp"
#include "gsl_root_finder.hpp"
#include "kink_profile_guesser.hpp"
#include "logging_manager.hpp"
#include "perturbative_profiler.hpp"
#include "shooting.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

namespace BubbleProfiler {

struct Fubini_options {
   double u{1.};
   double v{1.};
   double m{3.};
   bool calculate_profile{false};
   double global_min{-1.};
   std::string output_file{""};
   std::string perturbations_file{""};
   bool use_perturbative{false};
   bool verbose{false};
};

void print_usage()
{
   std::cout <<
      "Usage: fubini [OPTION] <u> <v> <m>\n\n"
      "Calculate action and field profile for one-dimensional binomial\n"
      "potential characterized by the parameters u, v, and m.\n\n"
      "Example: fubini -c -o 'output.txt' 1 1 3\n\n"
      "Options:\n"
      "  -c, --calculate-profile         calculate field profile as well\n"
      "  -g, --global-minimum=MIN        use MIN as location of true vacuum\n"
      "  -h, --help                      print this help message\n"
      "  -o, --output-file=FILE          write results to FILE\n"
      "  -P, --perturbations-file=FILE   write perturbations to FILE\n"
      "  -p, --perturbative              use the perturbative algorithm\n"
      "  -v, --verbose                   produce verbose output"
             << std::endl;
}

bool starts_with(const std::string& option, const std::string& prefix)
{
   return !option.compare(0, prefix.size(), prefix);
}

std::string get_option_value(const std::string& option,
                             const std::string& sep = "=")
{
   std::string value{""};
   const auto prefix_end = option.find(sep);

   if (prefix_end != std::string::npos) {
      value = option.substr(prefix_end + 1);
   }

   return value;
}

void parse_positional_args(const std::vector<std::string>& args,
                           Fubini_options& options)
{
   const auto n_positional_args = args.size();

   if (n_positional_args == 0) {
      throw Setup_error(
         "value of u must be provided");
   } else if (n_positional_args == 1) {
      throw Setup_error(
         "value of v must be provided");
   } else if (n_positional_args == 2) {
      throw Setup_error(
         "value of m must be provided");
   } else if (n_positional_args > 3) {
      throw Setup_error(
         "unrecognized command line argument '"
         + args[3] + "'");
   }

   try {
      options.u = std::stod(args[0]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value of u, '" + args[0] + "'");
   }

   try {
      options.v = std::stod(args[1]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value for v, '" + args[1] + "'");
   }

   try {
      options.m = std::stod(args[2]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value for m, '" + args[2] + "'");
   }
}

double parse_global_min(const std::string& value)
{
   if (value.empty()) {
      throw Setup_error(
         "'-g' or '--global-minimum' given but no value specified");
   }

   double global_min = 0.;
   try {
      global_min = std::stod(value);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid global minimum location '" + value + "'");
   }

   if (global_min <= 0.) {
      throw Setup_error(
         "explicitly specified true vacuum must be positive");
   }

   return global_min;
}

Fubini_options parse_cmd_line_args(int argc, const char* argv[])
{
   Fubini_options options;

   std::vector<std::string> positional_args;
   int i = 1;
   while (i < argc) {
      const std::string opt(argv[i++]);

      if (!starts_with(opt, "-")) {
         positional_args.push_back(opt);
         continue;
      }

      if (opt == "-c" || opt == "--calculate-profile") {
         options.calculate_profile = true;
         continue;
      }

      if (opt == "-g") {
         if (i == argc) {
            throw Setup_error(
               "'-g' given but no value for true vacuum provided");
         }
         const std::string value(argv[i++]);
         options.global_min = parse_global_min(value);
         continue;
      }

      if (starts_with(opt, "--global-minimum=")) {
         const std::string value = get_option_value(opt);
         options.global_min = parse_global_min(value);
         continue;
      }

      if (opt == "-h" || opt == "--help") {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      if (opt == "-o") {
         if (i == argc) {
            throw Setup_error(
               "'-o' given but no output file name provided");
         }
         const std::string filename(argv[i++]);
         if (starts_with(filename, "-") && filename != "-") {
            throw Setup_error(
               "'-o' given but no output file name provided");
         }
         options.output_file = filename;
         continue;
      }

      if (starts_with(opt, "--output-file=")) {
         const std::string filename = get_option_value(opt);
         if (filename.empty()) {
            throw Setup_error(
               "'--output-file=' given but no output file name provided");
         }
         options.output_file = filename;
         continue;
      }

      if (opt == "-P") {
         if (i == argc) {
            throw Setup_error(
               "'-P' given but no file name provided");
         }
         const std::string filename(argv[i++]);
         if (starts_with(filename, "-") && filename != "-") {
            throw Setup_error(
               "'-P' given but no file name provided");
         }
         options.perturbations_file = filename;
         continue;
      }

      if (starts_with(opt, "--perturbations-file=")) {
         const std::string filename = get_option_value(opt);
         if (filename.empty()) {
            throw Setup_error(
               "'--perturbations-file=' given but no file name provided");
         }
         options.perturbations_file = filename;
         continue;
      }

      if (opt == "-p" || opt == "--perturbative") {
         options.use_perturbative = true;
         continue;
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

std::tuple<double, Field_profiles> solve_shooting(
   const Generalized_fubini_potential& potential, double true_min,
   bool calculate_profile)
{
   const double false_min = potential.get_local_minimum_location();
   const double barrier = potential.get_local_maximum_location();

   const int dim = 4;

   unsigned int options
      = calculate_profile ?
      (Shooting::Solver_options::Compute_action |
       Shooting::Solver_options::Compute_profile)
      : Shooting::Solver_options::Compute_action;

   double action = 0.;
   Field_profiles profile;

   Shooting profiler;
   profiler.solve(potential, false_min, true_min, barrier, dim,
                  options);

   if ((options & Shooting::Solver_options::Compute_action) != 0) {
      action = profiler.get_euclidean_action();
   }

   if ((options & Shooting::Solver_options::Compute_profile) != 0) {
      profile = profiler.get_bubble_profile();
   }

   return std::make_tuple(action, profile);
}

std::tuple<double, Field_profiles> solve_perturbative(
   Generalized_fubini_potential& potential, double true_min,
   const std::string& observer_file)
{
   using Root_finder = GSL_root_finder<Eigen::Dynamic>;

   RK4_perturbative_profiler profiler;

   profiler.set_initial_step_size(1.e-3);

   profiler.set_initial_guesser(std::make_shared<Kink_profile_guesser>());

   const double action_tol = 1.e-3;
   const double fields_tol = 1.e-3;
   profiler.set_convergence_tester(
      std::make_shared<Relative_convergence_tester>(
         action_tol, fields_tol));

   auto root_finder = std::make_shared<Root_finder>();
   profiler.set_root_finder(root_finder);

   const int dim = 4;

   Eigen::VectorXd true_vacuum(1);
   true_vacuum << true_min;

   Eigen::VectorXd false_vacuum(1);
   false_vacuum << potential.get_local_minimum_location();

   profiler.set_number_of_dimensions(dim);
   profiler.set_true_vacuum_loc(true_vacuum);
   profiler.set_false_vacuum_loc(false_vacuum);

   double action = 0.;
   Field_profiles profile;

   if (observer_file.empty()) {
      profiler.calculate_bubble_profile(potential);
   } else {
      Generalized_fubini_observer observer(potential, observer_file);
      profiler.calculate_bubble_profile(potential, observer);
   }
   action = profiler.get_euclidean_action();
   profile = profiler.get_bubble_profile();

   return std::make_tuple(action, profile);
}

void write_profile(std::ostream& ostr, const Field_profiles& numerical_profile,
                   const Field_profiles& exact_profile)
{
   ostr << "# "
        << std::setw(16) << "rho" << ' '
        << std::setw(16) << "numeric" << ' '
        << std::setw(16) << "exact" << std::endl;

   const auto n_grid_points = numerical_profile.get_number_of_grid_points();
   if (n_grid_points != exact_profile.get_number_of_grid_points()) {
      throw Setup_error(
         "number of grid points do not match in numerical and exact solution");
   }

   const auto rho_values = numerical_profile.get_spatial_grid();
   const auto numerical_values = numerical_profile.get_field_profiles();
   const auto exact_values = exact_profile.get_field_profiles();

   for (int i = 0; i < n_grid_points; ++i) {
      ostr << "  "
           << std::setw(16) << std::setprecision(8)
           << std::scientific << rho_values(i) << ' '
           << std::setw(16) << std::setprecision(8)
           << std::scientific << numerical_values(i, 0) << ' '
           << std::setw(16) << std::setprecision(8)
           << std::scientific << exact_values(i, 0) << std::endl;
   }
}

void write_action(std::ostream& ostr, double numerical_action,
                  double exact_action)
{
   const double rel_diff = (numerical_action - exact_action)
      / (0.5 * (numerical_action + exact_action));

   ostr << "# Numerical action = "
        << std::setprecision(8) << std::scientific << numerical_action
        << std::endl;
   ostr << "# Exact action     = "
        << std::setprecision(8) << std::scientific << exact_action
        << std::endl;
   ostr << "# Relative error   = "
        << std::setprecision(8) << std::scientific << rel_diff
        << std::endl;
}

} // namespace BubbleProfiler

int main(int argc, const char* argv[])
{
   using namespace BubbleProfiler;

   if (argc == 1) {
      std::cerr << "Usage: fubini [OPTION] <u> <v> <m>\n"
                << "Try 'fubini --help' for more information."
                << std::endl;
      return EXIT_FAILURE;
   }

   int exit_code = 0;
   try {
      Fubini_options options = parse_cmd_line_args(argc, argv);

      if (options.global_min <= 0.) {
         options.global_min = 2.;
      }

      if (options.verbose) {
         logging::Logging_manager::get_manager().set_minimum_log_level(
            logging::Log_level::Trace);
      }

      Generalized_fubini_potential potential(options.u, options.v, options.m);

      // Sanity check on potential
      const double local_min = potential.get_local_minimum_location();
      if (potential(options.global_min) > potential(local_min)) {
         throw Setup_error(
            "true minimum of potential greater than false minimum");
      }

      std::tuple<double, Field_profiles> result;
      if (options.use_perturbative) {
         result = solve_perturbative(potential, options.global_min,
                                     options.perturbations_file);
      } else {
         result = solve_shooting(potential, options.global_min,
                                 options.calculate_profile);
      }

      std::ofstream output_file;
      if (!options.output_file.empty()) {
         output_file.open(options.output_file, std::ios::out);
      }

      std::ostream& output_stream
         = (options.output_file.empty() ? std::cout : output_file);

      if (options.calculate_profile) {
         const Field_profiles profiles = std::get<1>(result);
         const Field_profiles exact_profile
            = potential.get_profile(profiles.get_spatial_grid());
         write_profile(output_stream, std::get<1>(result), exact_profile);
      }

      const double exact_action = potential.get_action();
      write_action(output_stream, std::get<0>(result), exact_action);
   } catch (const Setup_error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      exit_code = EXIT_FAILURE;
   } catch (const BVP_solver_error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      exit_code = EXIT_FAILURE;
   } catch (const Numerical_error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      exit_code = EXIT_FAILURE;
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
