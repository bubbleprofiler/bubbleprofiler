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
    @brief Program that tabulates action for a one-dimensional potential.

    The potential is
    \f[
    V(\phi) = E\left[\frac{-4 \alpha +3}{2} \phi^2 - \phi ^3 + \alpha \phi ^4\right].
    \f]
    We pick \f$E = 1\f$ and vary \f$\alpha\f$.

    Usage:
    @code
    ./tabulate <dim> <min_alpha> <max_alpha> <step>
    @endcode
*/

#include "error.hpp"
#include "gsl_root_finder.hpp"
#include "kink_profile_guesser.hpp"
#include "logging_manager.hpp"
#include "perturbative_profiler.hpp"
#include "relative_convergence_tester.hpp"
#include "restricted_quartic_potential.hpp"
#include "shooting.hpp"

#include <chrono>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace BubbleProfiler {

namespace {

const std::string shooting_label = "shooting";
const std::string perturbative_label = "perturbative";

} // anonymous namespace

struct Tabulate_options {
   int dim{4};
   double E{1.};
   double min_alpha{0.51};
   double max_alpha{0.74};
   double alpha_step{0.01};
   bool use_shooting{true};
   bool use_perturbative{false};
   bool verbose{false};
};

struct Profiler_result {
   double action{0.};
   double time_in_ms{0.};
   int error{0};
   std::string info{""};
};

using Profiler_data = std::map<std::string, std::vector<Profiler_result> >;

void print_usage()
{
   std::cout <<
      "Usage: quartic_tabulate [OPTION] <dim> <min_alpha> <max_alpha> <alpha_step>\n\n"
      "Calculate action for one-dimensional quartic potentials with alpha\n"
      "varying between the given minimum and maximum values.\n\n"
      "Example: quartic_tabulate 4 0.51 0.74 0.1\n\n"
      "Options:\n"
      "  -a, --all-methods               print results for all algorithms\n"
      "  -h, --help                      print this help message\n"
      "  -p, --perturbative              use the perturbative algorithm\n"
      "  -v, --verbose                   produce verbose output"
             << std::endl;
}

bool starts_with(const std::string& option, const std::string& prefix)
{
   return !option.compare(0, prefix.size(), prefix);
}

void parse_positional_args(const std::vector<std::string>& args,
                           Tabulate_options& options)
{
   const auto n_positional_args = args.size();

   if (n_positional_args == 0) {
      throw Setup_error(
         "number of dimensions must be provided");
   } else if (n_positional_args == 1) {
      throw Setup_error(
         "minimum value of alpha must be provided");
   } else if (n_positional_args == 2) {
      throw Setup_error(
         "maximum value of alpha must be provided");
   } else if (n_positional_args == 3) {
      throw Setup_error(
         "alpha increment must be provided");
   } else if (n_positional_args > 4) {
      throw Setup_error(
         "unrecognized command line argument '"
         + args[4] + "'");
   }

   try {
      options.dim = std::stoi(args[0]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid number of dimensions '" + args[0] + "'");
   }

   try {
      options.min_alpha = std::stod(args[1]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid minimum value for alpha, '" + args[1] + "'");
   }

   try {
      options.max_alpha = std::stod(args[2]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid maximum value for alpha, '" + args[2] + "'");
   }

   try {
      options.alpha_step = std::stod(args[3]);
   } catch (const std::invalid_argument& e) {
      throw Setup_error(
         "invalid alpha increment, '" + args[3] + "'");
   }
}

Tabulate_options parse_cmd_line_args(int argc, const char* argv[])
{
   Tabulate_options options;

   std::vector<std::string> positional_args;
   int i = 1;
   while (i < argc) {
      const std::string opt(argv[i++]);

      if (!starts_with(opt, "-")) {
         positional_args.push_back(opt);
         continue;
      }

      if (opt == "-a" || opt == "--all-methods") {
         options.use_shooting = true;
         options.use_perturbative = true;
         continue;
      }

      if (opt == "-h" || opt == "--help") {
         print_usage();
         exit(EXIT_SUCCESS);
      }

      if (opt == "-p" || opt == "--perturbative") {
         options.use_shooting = false;
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

Profiler_result run_shooting_profiler(
   const Restricted_quartic_potential& potential, int dim)
{
   using Clock_t = std::chrono::high_resolution_clock;

   const double false_min = potential.get_local_minimum_location();
   const double barrier = potential.get_local_maximum_location();
   const double true_min = potential.get_global_minimum_location();

   Shooting profiler;
   profiler.set_bisection_precision_bits(5);
   profiler.set_action_arrived_rel(1.e-3);
   profiler.set_shooting_abs_tol(1.e-4);
   profiler.set_shooting_rel_tol(1.e-4);
   profiler.set_action_abs_tol(1.e-6);
   profiler.set_action_rel_tol(1.e-6);
   profiler.set_drho_frac(1.e-3);
   profiler.set_bisect_lambda_max(5);
   profiler.set_max_iterations(100000);
   profiler.set_max_periods(1.e2);
   profiler.set_f_y_max(1.e6);
   profiler.set_f_y_min(1.e-3);
   profiler.set_y_max(1.e1);

   Profiler_result result;
   const auto start_time = Clock_t::now();
   try {
      profiler.solve(potential, false_min, true_min, barrier, dim,
                     Shooting::Solver_options::Compute_action);
      result.action = profiler.get_euclidean_action();
   } catch (const Error& e) {
      result.error = 1;
      result.info = e.what();
   }
   const auto stop_time = Clock_t::now();
   std::chrono::duration<double, std::milli> t_in_ms = stop_time - start_time;
   result.time_in_ms = t_in_ms.count();

   return result;
}

Profiler_result run_perturbative_profiler(
   Restricted_quartic_potential potential, int dim)
{
   using Clock_t = std::chrono::high_resolution_clock;
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

   Profiler_result result;
   const auto start_time = Clock_t::now();
   try {
      profiler.calculate_bubble_profile(potential);
      result.action = profiler.get_euclidean_action();
   } catch (const Error& e) {
      result.error = 1;
      result.info = e.what();
   }
   const auto stop_time = Clock_t::now();
   std::chrono::duration<double, std::milli> t_in_ms = stop_time - start_time;
   result.time_in_ms = t_in_ms.count();

   return result;
}

Profiler_result run_profiler(
   const std::string& profiler,
   const Restricted_quartic_potential& potential, int dim)
{
   if (profiler == shooting_label) {
      return run_shooting_profiler(potential, dim);
   } else if (profiler == perturbative_label) {
      return run_perturbative_profiler(potential, dim);
   }

   throw Setup_error(
      "unrecognized profiler '" + profiler + "'");
}

void write_table_header(const std::vector<std::string>& profiler_labels)
{
   std::cout << "# "
             << std::setw(12) << "dim" << ' '
             << std::setw(16) << "E" << ' '
             << std::setw(16) << "alpha";
   for (const auto& profiler : profiler_labels) {
      std::cout << ' '
                << std::setw(16) << profiler + "_action" << ' '
                << std::setw(16) << profiler + "_time/ms" << ' '
                << std::setw(12) << profiler + "_error";
   }
   std::cout << std::endl;
}

void write_table(int dim, double E, const std::vector<double>& alpha_data,
                 const Profiler_data& profiler_data)
{
   std::vector<std::string> profilers;
   for (const auto& x : profiler_data) {
      profilers.push_back(x.first);
   }

   write_table_header(profilers);

   const std::size_t n_rows = alpha_data.size();
   for (std::size_t i = 0; i < n_rows; ++i) {
      std::cout << "  "
                << std::setw(12) << dim << ' '
                << std::setw(16) << std::setprecision(8)
                << std::scientific << E << ' '
                << std::setw(16) << std::setprecision(8)
                << std::scientific << alpha_data[i];

      std::string info = "";
      for (const auto& profiler : profilers) {
         const auto& data = profiler_data.at(profiler)[i];
         std::cout << ' '
                   << std::setw(16) << std::setprecision(8)
                   << std::scientific << data.action << ' '
                   << std::setw(16) << std::setprecision(8)
                   << std::scientific << data.time_in_ms << ' '
                   << std::setw(12) << data.error;

         if (!data.info.empty()) {
            if (info.empty()) {
               info = "\t# ";
            }
            info += profiler + ": " + data.info + ",";
         }
      }

      std::cout << info << std::endl;
   }
}

} // namespace BubbleProfiler

int main(int argc, const char* argv[])
{
   using namespace BubbleProfiler;

   if (argc == 1) {
      std::cerr << "Usage: quartic_tabulate [OPTION] <dim> <min_alpha> <max_alpha> <alpha_step>\n"
                << "Try 'quartic_tabulate --help' for more information."
                << std::endl;
      return EXIT_FAILURE;
   }

   int exit_code = 0;
   try {
      Tabulate_options options = parse_cmd_line_args(argc, argv);

      if (options.verbose) {
         logging::Logging_manager::get_manager().set_minimum_log_level(
            logging::Log_level::Trace);
      }

      if (std::abs(options.alpha_step)
          < std::numeric_limits<double>::epsilon()) {
         throw Setup_error(
            "alpha increment too small");
      }

      const bool forward_stepping = options.alpha_step > 0.;
      const bool start_gt_end = options.min_alpha > options.max_alpha;
      if ((forward_stepping && start_gt_end) ||
          (!forward_stepping && !start_gt_end)) {
         std::swap(options.min_alpha, options.max_alpha);
      }

      std::vector<double> alpha_data;
      Profiler_data profiler_data;
      if (options.use_shooting) {
         profiler_data[shooting_label] = std::vector<Profiler_result>();
      }

      if (options.use_perturbative) {
         profiler_data[perturbative_label] = std::vector<Profiler_result>();
      }

      double alpha = options.min_alpha;
      while (alpha <= options.max_alpha) {
         Restricted_quartic_potential potential(alpha, options.E);

         alpha_data.push_back(alpha);
         for (auto& x : profiler_data) {
            profiler_data[x.first].push_back(
               run_profiler(x.first, potential, options.dim));
         }

         alpha += options.alpha_step;
      }

      write_table(options.dim, options.E, alpha_data, profiler_data);

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
