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
    @brief Program that calculates action for a one-dimensional potential in a
    thin-wall case,  in which exact solution is known.

     The potential is (see Eq. (4.5) of  https://journals.aps.org/prd/pdf/10.1103/PhysRevD.15.2929)
    \f[
    V(\phi) = \frac{\lambda}{8} (\phi^2 - a^2)^2 + \frac{\epsilon}{2 a} (\phi - a).
    \f]

    Usage:
    @code
    ./thin <lambda> <a> <delta>
     @endcode
*/

#include "error.hpp"
#include "field_profiles.hpp"
#include "gsl_root_finder.hpp"
#include "kink_profile_guesser.hpp"
#include "logging_manager.hpp"
#include "perturbative_profiler.hpp"
#include "shooting.hpp"
#include "thin_wall_observer.hpp"
#include "thin_wall_potential.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

namespace BubbleProfiler {

struct Thin_wall_options {
   double lambda{1.};
   double a{1.};
   double delta{0.01};
   int dim{4};
   bool calculate_profile{false};
   std::string output_file{""};
   std::string perturbations_file{""};
   bool use_perturbative{false};
   bool verbose{false};
};

void print_usage()
{
   std::cout <<
      "Usage: thin [OPTION] <lambda> <a> <delta>\n\n"
      "Calculate action and field profile for one-dimensional perturbed\n"
      "quartic potential characterized by the parameters lambda, a,\n"
      "and delta.\n\n"
      "Example: thin -c -o 'output.txt' 1 1 0.01\n\n"
      "Options:\n"
      "  -c, --calculate-profile         calculate field profile as well\n"
      "  -f, --finite-temperature        calculate using d = 3 dimensions\n"
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
                           Thin_wall_options& options)
{
   const auto n_positional_args = args.size();

   if (n_positional_args == 0) {
      throw Setup_error(
         "value of lambda must be provided");
   } else if (n_positional_args == 1) {
      throw Setup_error(
         "value of a must be provided");
   } else if (n_positional_args == 2) {
      throw Setup_error(
         "value of delta must be provided");
   } else if (n_positional_args > 3) {
      throw Setup_error(
         "unrecognized command line argument '"
         + args[3] + "'");
   }

   try {
      options.lambda = std::stod(args[0]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value of lambda, '" + args[0] + "'");
   }

   try {
      options.a = std::stod(args[1]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value for a, '" + args[1] + "'");
   }

   try {
      options.delta = std::stod(args[2]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid value for delta, '" + args[2] + "'");
   }
}

Thin_wall_options parse_cmd_line_args(int argc, const char* argv[])
{
   Thin_wall_options options;

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

      if (opt == "-f" || opt == "--finite-temperature") {
         options.dim = 3;
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
   const Thin_wall_potential& potential, int dim, bool calculate_profile)
{
   const double false_min = potential.get_local_minimum_location();
   const double barrier = potential.get_local_maximum_location();
   const double true_min = potential.get_global_minimum_location();

   unsigned int options
      = calculate_profile ?
      (Shooting::Solver_options::Compute_action |
       Shooting::Solver_options::Compute_profile)
      : Shooting::Solver_options::Compute_action;

   Shooting profiler;
   // Strict tolerances for thin-wall solutions
   profiler.set_bisection_precision_bits(30);
   profiler.set_shooting_abs_tol(1.e-15);
   profiler.set_shooting_rel_tol(1.e-15);
   profiler.set_action_abs_tol(1.e-15);
   profiler.set_action_rel_tol(1.e-15);
   profiler.set_evolve_change_rel(1.e-4);

   double action = 0.;
   Field_profiles profile;

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
   Thin_wall_potential& potential, int dim,
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

   Eigen::VectorXd true_vacuum(1);
   true_vacuum << potential.get_global_minimum_location();

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
      Thin_wall_observer observer(observer_file);
      profiler.calculate_bubble_profile(potential, observer);
   }
   action = profiler.get_euclidean_action();
   profile = profiler.get_bubble_profile();

   return std::make_tuple(action, profile);
}

void write_profile(std::ostream& ostr, const Field_profiles& numerical_profile)
{
   ostr << "# "
        << std::setw(16) << "rho" << ' '
        << std::setw(16) << "numeric" << std::endl;

   const auto n_grid_points = numerical_profile.get_number_of_grid_points();

   const auto rho_values = numerical_profile.get_spatial_grid();
   const auto numerical_values = numerical_profile.get_field_profiles();

   for (int i = 0; i < n_grid_points; ++i) {
      ostr << "  "
           << std::setw(16) << std::setprecision(8)
           << std::scientific << rho_values(i) << ' '
           << std::setw(16) << std::setprecision(8)
           << std::scientific << numerical_values(i, 0) << std::endl;
   }
}

void write_action(std::ostream& ostr, double numerical_action,
                  double thin_wall_action)
{
   const double rel_diff = (numerical_action - thin_wall_action)
      / (0.5 * (numerical_action + thin_wall_action));

   ostr << "# Numerical action = "
        << std::setprecision(8) << std::scientific << numerical_action
        << std::endl;
   ostr << "# Thin wall action = "
        << std::setprecision(8) << std::scientific << thin_wall_action
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
      std::cerr << "Usage: thin [OPTION] <lambda> <a> <delta>\n"
                << "Try 'thin --help' for more information."
                << std::endl;
      return EXIT_FAILURE;
   }

   int exit_code = 0;
   try {
      Thin_wall_options options = parse_cmd_line_args(argc, argv);

      if (options.verbose) {
         logging::Logging_manager::get_manager().set_minimum_log_level(
            logging::Log_level::Trace);
      }

      Thin_wall_potential potential(options.lambda, options.a, options.delta);

      std::tuple<double, Field_profiles> result;
      if (options.use_perturbative) {
         result = solve_perturbative(potential, options.dim,
                                     options.perturbations_file);
      } else {
         result = solve_shooting(potential, options.dim,
                                 options.calculate_profile);
      }

      std::ofstream output_file;
      if (!options.output_file.empty()) {
         output_file.open(options.output_file, std::ios::out);
      }

      std::ostream& output_stream
         = (options.output_file.empty() ? std::cout : output_file);

      if (options.calculate_profile) {
         write_profile(output_stream, std::get<1>(result));
      }

      const double thin_wall_action = potential.get_thin_wall_action();
      write_action(output_stream, std::get<0>(result), thin_wall_action);
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
