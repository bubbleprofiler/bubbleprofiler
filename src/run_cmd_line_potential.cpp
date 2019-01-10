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
    @brief Command line interface to `BubbleProfiler`
    @todo
    Usage, e.g.,
    \code
    ./BubbleProfiler
    \endcode
*/

#include "error.hpp"
#include "gsl_root_finder.hpp"
#include "instream_profile_guesser.hpp"
#include "kink_profile_guesser.hpp"
#include "logging_manager.hpp"
#include "nlopt_optimizer.hpp"
#include "observers.hpp"
#include "perturbative_profiler.hpp"
#include "potential.hpp"
#include "relative_convergence_tester.hpp"
#include "shooting.hpp"
#include "shooting_profile_guesser.hpp"
#include "algebraic_potential.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

namespace BubbleProfiler {

enum class Integration_algorithm : int { Runge_kutta_4 = 0, Euler };

const std::string default_ansatz_file = "";
const double default_domain_start = -1.;
const double default_domain_end = -1.;
const double default_initial_step_size = 1.e-2;
const int default_interp_fraction = 1.0;
const Integration_algorithm default_algorithm =
   Integration_algorithm::Runge_kutta_4;
const double default_rtol_action = 1.e-3;
const double default_rtol_fields = 1.e-3;
const std::string default_output_file = "-";
const bool default_write_perturbations = true;
const bool default_force_output = false;
const double default_opt_timeout = 10.;
const int default_max_iterations = -1;
const bool default_assume_origin_false_vacuum = false;
const bool default_use_perturbative = false;
const bool default_write_profiles = false;
const bool default_shooting_ansatz = false;
const bool default_verbose = false;

// Options pertaining only to 1D shooting method
const int default_shoot_bisect_bits = 5;
const double default_action_arrived_rel = 1.e-3;
const double default_shoot_ode_abs = 1.e-4;
const double default_shoot_ode_rel = 1.e-4;
const double default_action_ode_abs = 1.e-6;
const double default_action_ode_rel = 1.e-6;
const double default_drho_frac = 1.e-3;
const double default_evolve_change_rel = 1.e-2;
const double default_bisect_lambda_max = 5.0;
const int default_iter_max = 100000;
const double default_periods_max = 1.e2;
const double default_f_y_max = 1.e6;
const double default_f_y_min = 1.e-3;
const double default_y_max = 1.e1;



struct Bubble_profiler_inputs {
   std::vector<std::string> fields;
   std::string potential;
   int n_dims;
   double domain_start;
   double domain_end;
   double initial_step_size;
   double interpolation_fraction;
   Integration_algorithm algorithm;
   double rtol_action;
   double rtol_fields;
   bool write_perturbations;
   bool force_output;
   std::string output_file;
   std::string output_path;
   std::vector<double> global_min;
   std::vector<double> local_min;
   std::vector<double> barrier;
   double opt_timeout;
   int max_iterations;
   std::string ansatz_file;
   bool assume_origin_false_vacuum;
   bool use_perturbative;
   bool write_profiles;
   bool shooting_ansatz;
   bool verbose;

   // Options pertaining only to 1D shooting method
   int shoot_bisect_bits;
   double action_arrived_rel;
   double shoot_ode_abs;
   double shoot_ode_rel;
   double action_ode_abs;
   double action_ode_rel;
   double drho_frac;
   double evolve_change_rel;
   double bisect_lambda_max;
   int iter_max;
   double periods_max;
   double f_y_max;
   double f_y_min;
   double y_max;
};

void validate(boost::any& v,
              const std::vector<std::string>& values,
              Integration_algorithm*, int)
{
   using namespace po;

   validators::check_first_occurrence(v);

   const std::string& val = validators::get_single_string(values);

   if (val == "runge-kutta-4") {
      v = boost::any(Integration_algorithm::Runge_kutta_4);
   } else if (val == "euler") {
      v = boost::any(Integration_algorithm::Euler);
   } else {
      throw validation_error(validation_error::invalid_option_value);
   }
}

Bubble_profiler_inputs parse_cmd_line_args(int argc, const char* argv[])
{
   Bubble_profiler_inputs input;

   po::options_description desc("run_cmd_line_potential options");
   desc.add_options()
      ("ansatz-file", po::value<std::string>(&input.ansatz_file)->default_value(default_ansatz_file),
       "Text file containing a predefined ansatz")
      ("false-vacuum-at-origin",
       po::bool_switch(&input.assume_origin_false_vacuum)->default_value(default_assume_origin_false_vacuum),
       "If given, assume the false location is at the origin of field space unless specified otherwise")
      ("barrier", po::value<std::vector<double> >()->multitoken(),
       "Location of barrier. Applies only when using the shooting method for 1-field problems."
       " If not supplied, barrier will be found numerically.")
      ("domain-start", po::value<double>(&input.domain_start)->default_value(default_domain_start),
       "Domain of integration start.  If negative, the profiler will attempt to guess a suitable value.")
      ("domain-end", po::value<double>(&input.domain_end)->default_value(default_domain_end),
       "Domain of integration end.  If negative, the profiler will attempt to guess a suitable value.")
      ("field", po::value<std::vector<std::string> >(&input.fields)->required(),
       "Specify a field name")
      ("force-output", po::bool_switch(&input.force_output)->default_value(default_force_output),
       "If given, overwrite existing output files")
      ("global-minimum", po::value<std::vector<double> >()->multitoken(),
       "Location of global minimum")
      ("initial-step-size", po::value<double>(&input.initial_step_size)->default_value(default_initial_step_size),
       "Approximate initial step size to use in solving ODE")
      ("help", "print help message")
      ("integration-method", po::value<Integration_algorithm>(&input.algorithm)->default_value(default_algorithm, "runge-kutta-4"),
       "Algorithm to use for integrating ODEs (valid options: runge-kutta-4, euler)")
      ("interpolation-fraction", po::value<double>(&input.interpolation_fraction)->default_value(default_interp_fraction),
       "Approximate fraction of grid points to use for interpolation")
      ("local-minimum", po::value<std::vector<double> >()->multitoken(),
       "Location of local minimum")
      ("max-iterations", po::value<int>(),
       "Perform at most this many iterations")
      ("opt-timeout", po::value<double>(&input.opt_timeout)->default_value(default_opt_timeout),
       "Timeout if using the optimizer to locate global minimum (0 for no timeout)")
      ("output-file", po::value<std::string>(&input.output_file)->default_value(default_output_file),
       "Name of output file.  If '-' is given, write to standard output")
      ("output-path", po::value<std::string>(&input.output_path)->default_value(""),
       "Output path for additional generated data files")
      ("perturbative", po::bool_switch(&input.use_perturbative)->default_value(default_use_perturbative),
       "If given, use perturbative method for one-dimensional potentials")
      ("potential", po::value<std::string>(&input.potential)->required(),
       "Specify a potential as GiNaCs string")
      ("n-dims", po::value<int>(&input.n_dims)->default_value(3),
       "number of spacetime dimensions")
      ("rtol-action", po::value<double>(&input.rtol_action)->default_value(default_rtol_action),
       "Relative tolerance in bounce action")
      ("rtol-fields", po::value<double>(&input.rtol_fields)->default_value(default_rtol_fields),
       "Relative tolerance in fields at r = 0")
      ("write-profiles",  po::bool_switch(&input.write_profiles)->default_value(default_write_profiles),
       "If given, write the field profiles to disk as well as the action")
      ("shooting-ansatz",  po::bool_switch(&input.shooting_ansatz)->default_value(default_shooting_ansatz),
       "If given, use the shooting method to calculate an ansatz")
      ("verbose", po::bool_switch(&input.verbose)->default_value(default_verbose),
       "If given, produce verbose logging output")
      ("shoot-bisect-bits", po::value<int>(&input.shoot_bisect_bits)->default_value(default_shoot_bisect_bits),
       "1D shooting: target relative precision in significant bits when bisecting")
      ("action-arrived-rel", po::value<double>(&input.action_arrived_rel)->default_value(default_action_arrived_rel),
       "1D shooting: relative tolerance for arriving at false vacuum when calculating action")
      ("shoot-ode-abs", po::value<double>(&input.shoot_ode_abs)->default_value(default_shoot_ode_abs),
       "1D shooting: absolute error tolerance for ODE integrator (bubble profile)")
      ("shoot-ode-rel", po::value<double>(&input.shoot_ode_rel)->default_value(default_shoot_ode_rel),
       "1D shooting: relative error tolerance for ODE integrator (bubble profile)")
      ("action-ode-abs", po::value<double>(&input.action_ode_abs)->default_value(default_action_ode_abs),
       "1D shooting: absolute error tolerance for ODE integrator (action integral)")
      ("action-ode-rel", po::value<double>(&input.action_ode_rel)->default_value(default_action_ode_rel),
       "1D shooting: relative error tolerance for ODE integrator (action integral)")
      ("drho-frac", po::value<double>(&input.drho_frac)->default_value(default_drho_frac),
       "1D shooting: initial step size relative to characteristic bubble scale")
      ("evolve-change-rel", po::value<double>(&input.evolve_change_rel)->default_value(default_evolve_change_rel),
       "1D shooting: evolve using approximate solution until field has changed by this much, relative to false vacuum")
      ("bisect-lambda-max", po::value<double>(&input.bisect_lambda_max)->default_value(default_bisect_lambda_max),
       "1D shooting: maximum value of bisection parameter lambda")
      ("iter-max", po::value<int>(&input.iter_max)->default_value(default_iter_max),
       "1D shooting: set maximum number of shooting iterations")
      ("periods-max", po::value<double>(&input.periods_max)->default_value(default_periods_max),
       "1D shooting: domain size in terms of characteristic oscillation periods")
      ("f-y-max", po::value<double>(&input.f_y_max)->default_value(default_f_y_max),
       "1D shooting: threshold for asymptotic approximation in analytic evolution")
      ("f-y-min", po::value<double>(&input.f_y_min)->default_value(default_f_y_min),
       "1D shooting: threshold for asymptotic approximation in analytic evolution")
      ("y-max", po::value<double>(&input.y_max)->default_value(default_y_max),
       "1D shooting: threshold for asymptotic approximation in analytic evolution")
      ;

   po::variables_map vm;
   po::store(
      po::parse_command_line(
         argc, argv, desc,
         po::command_line_style::unix_style ^ po::command_line_style::allow_short),
      vm);

   if (vm.count("help")) {
      std::cout << desc << std::endl;
      exit(EXIT_SUCCESS);
   }

   po::notify(vm);

   if (vm.count("max-iterations")) {
      const int max_iterations = vm["max-iterations"].as<int>();
      if (max_iterations < 0) {
         std::cerr << "Error: number of iterations must be non-negative"
                   << std::endl;
         exit(EXIT_FAILURE);
      }
      input.max_iterations = max_iterations;
   } else {
      input.max_iterations = default_max_iterations;
   }

   if (vm.count("global-minimum")) {
      input.global_min = vm["global-minimum"].as<std::vector<double> >();
   }
   if (vm.count("local-minimum")) {
      input.local_min = vm["local-minimum"].as<std::vector<double> >();
   }
   if (vm.count("barrier")) {
      input.barrier = vm["barrier"].as<std::vector<double> >();
   }

   return input;
}

/*!
 * @brief Write the value of the action to requested output stream
 *
 * @param action the value of the action
 * @param ostr the output stream to write to
 */
void write_action(double action, std::ostream& ostr)
{
   ostr << "# Action: " << action << std::endl;
}

/*!
 * @brief Write the given field profiles to requested output stream
 *
 * @param fields the names of the fields
 * @param profiles the values of the field profiles
 * @param ostr the output stream to write to
 */
void write_profiles(const std::vector<std::string>& fields,
                    const Field_profiles& profiles, std::ostream& ostr)
{
   ostr << "#" << ' '
        << std::left << std::setw(16) << "rho";

   for (const auto& name: fields) {
      ostr << ' ' << std::left << std::setw(16) << name;
   }
   ostr << std::endl;

   const auto coord_values = profiles.get_spatial_grid();
   const auto field_values = profiles.get_field_profiles();

   const int n_grid_points = coord_values.size();
   const int n_fields = fields.size();

   for (int i = 0; i < n_grid_points; ++i) {
      ostr << "  " << std::left
           << std::setw(16) << std::setprecision(8)
           << std::scientific << coord_values(i);
      for (int j = 0; j < n_fields; ++j) {
         ostr << ' ' << std::left
              << std::setw(16) << std::setprecision(8)
              << std::scientific << field_values(i, j);
      }
      ostr << std::endl;
   }
}

/*!
 * @brief Find the global minimum of the given potential
 * @param potential the potential to minimize
 * @param opt_timeout the allowed maximum amount of time to run for
 * @return the obtained location of the minimum
 */
Eigen::VectorXd find_global_min(const Potential& potential, double opt_timeout)
{
   const int n_fields = potential.get_number_of_fields();
   const auto v = [&potential](const Eigen::VectorXd& x) -> double {
      return potential(x);
   };
   NLopt_optimizer optimizer(v, n_fields);
   optimizer.set_extremum_type(NLopt_optimizer::Extremum_type::MIN);
   optimizer.set_lower_bounds(-1.e6);
   optimizer.set_upper_bounds(1.e6);
   optimizer.set_max_time(opt_timeout);

   const Eigen::VectorXd initial_guess(Eigen::VectorXd::Zero(n_fields));
   const auto status = optimizer.optimize(initial_guess);

   if (!optimization_succeeded(status)) {
      std::cerr << "Error: unable to locate global minimum." << std::endl;
      exit(EXIT_FAILURE);
   }

   return optimizer.get_extremum_location();
}

/*!
 * @brief Find barrier directly between two minima of a 1D potential
 *
 * In the one-dimensional case, this function attempts to find the
 * location of the barrier between two given minima, taken to be
 * the location of the maximum value of the potential between the two
 * minima.
 *
 * @param potential the potential to locate the barrier for
 * @param true_vacuum_loc the location of the global minimum
 * @param false_vacuum_loc the location of the local minimum
 * @param opt_timeout the allowed maximum amount of time to run for
 * @return the obtained location of the maximum
 */
Eigen::VectorXd find_one_dimensional_barrier(
   const Potential& potential,
   const Eigen::VectorXd& true_vacuum_loc,
   const Eigen::VectorXd& false_vacuum_loc,
   double opt_timeout)
{
   const int n_fields = potential.get_number_of_fields();
   if (n_fields != 1) {
      throw Setup_error("automatically locating potential barrier only "
                        "supported for single field case");
   }

   const auto v = [&potential](const Eigen::VectorXd& x) {
      return potential(x);
   };

   NLopt_optimizer optimizer(v, n_fields);

   // Don't need much precision for location of barrier
   optimizer.set_xtol_rel(1.e-2);
   optimizer.set_ftol_rel(1.e-2);

   optimizer.set_extremum_type(NLopt_optimizer::Extremum_type::MAX);
   optimizer.set_lower_bounds(std::min(true_vacuum_loc(0),
                                       false_vacuum_loc(0)));
   optimizer.set_upper_bounds(std::max(true_vacuum_loc(0),
                                       false_vacuum_loc(0)));
   optimizer.set_max_time(opt_timeout);

   Eigen::VectorXd initial_guess(0.5 * (true_vacuum_loc + false_vacuum_loc));
   const auto status = optimizer.optimize(initial_guess);

   if (!optimization_succeeded(status)) {
      std::cerr << "Error: unable to locate barrier. NLOPT status = "
                << status << std::endl;
      exit(EXIT_FAILURE);
   }

   return optimizer.get_extremum_location();
}

/*!
 * @brief Locate the requested critical points if not already provided
 *
 * If not empty, the values of \c input.global_min and \c input.local_min
 * are used to set the values of the global minimum and local minimum,
 * respectively.
 *
 * If \c input.global_min is empty, the function attempts to locate the
 * true vacuum by numerically optimizing the potential.  By default, the
 * function will not attempt to find the local minimum if
 * \c input.local_min is empty, and will exit with an error.  However,
 * if \c input.local_min is empty and \c input.assume_origin_false_vacuum
 * is set to \c true, then the location of the false vacuum will be set
 * to the origin in field space.
 *
 * @param potential the potential to find the critical points of
 * @param input settings to be used to configure the potential and solvers
 * @param true_vacuum_loc the vector in which to store the obtained
 *                        global minimum
 * @param false_vacuum_loc the vector in which to store the obtained
 *                         local minimum
 */
void initialize_extrema(const Potential& potential,
                        const Bubble_profiler_inputs& input,
                        Eigen::VectorXd& true_vacuum_loc,
                        Eigen::VectorXd& false_vacuum_loc)
{
   const std::size_t n_fields = potential.get_number_of_fields();

   if (input.local_min.empty()) {
      if (input.assume_origin_false_vacuum) {
         false_vacuum_loc = Eigen::VectorXd::Zero(n_fields);
      } else {
         throw Setup_error("location of false vacuum not specified");
      }
   } else {
      if (input.local_min.size() != n_fields) {
         throw Setup_error("number of fields does not match dimensions "
                           "of local minimum");
      }

      false_vacuum_loc = Eigen::VectorXd::Map(input.local_min.data(),
                                              input.local_min.size());
   }

   if (input.global_min.empty()) {
      true_vacuum_loc = find_global_min(potential, input.opt_timeout);
   } else {
      if (input.global_min.size() != n_fields) {
         throw Setup_error("number of fields does not match dimensions "
                           "of global minimum");
      }

      true_vacuum_loc = Eigen::VectorXd::Map(input.global_min.data(),
                                             input.global_min.size());
   }
}

/*!
 * @brief Locate the necessary critical points if not already provided
 *
 * If not empty, the values of \c input.global_min , \c input.local_min ,
 * and \c input.barrier are used to set the values of the global minimum,
 * local minimum, and barrier, respectively.
 *
 * If \c input.global_min is empty, the function attempts to locate the
 * true vacuum by numerically optimizing the potential.  By default, the
 * function will not attempt to find the local minimum if
 * \c input.local_min is empty, and will exit with an error.  However,
 * if \c input.local_min is empty and \c input.assume_origin_false_vacuum
 * is set to \c true, then the location of the false vacuum will be set
 * to the origin in field space.
 *
 * Locating the position of the barrier between two minima is only
 * supported for one-field potentials.  In these cases, if
 * \c input.barrier is empty, the local maximum of the potential is
 * obtained numerically; otherwise, the given value is used.
 *
 * @param potential the potential to find the critical points of
 * @param input settings to be used to configure the potential and solvers
 * @param true_vacuum_loc the vector in which to store the obtained
 *                        global minimum
 * @param false_vacuum_loc the vector in which to store the obtained
 *                         local minimum
 * @param barrier_loc the vector in which to store the obtained
 *                    barrier location
 */
void initialize_extrema(const Potential& potential,
                        const Bubble_profiler_inputs& input,
                        Eigen::VectorXd& true_vacuum_loc,
                        Eigen::VectorXd& false_vacuum_loc,
                        Eigen::VectorXd& barrier_loc)
{
   const std::size_t n_fields = potential.get_number_of_fields();

   if (n_fields != 1) {
      throw Setup_error("automatically locating potential barrier only "
                        "supported for single field case");
   }

   initialize_extrema(potential, input, true_vacuum_loc, false_vacuum_loc);

   if (input.barrier.empty()) {
      barrier_loc = find_one_dimensional_barrier(potential,
                                                 false_vacuum_loc,
                                                 true_vacuum_loc,
                                                 input.opt_timeout);
   } else {
      if (input.barrier.size() != n_fields) {
         throw Setup_error("number of fields does not match dimensions "
                           "of barrier location");
      }
      barrier_loc = Eigen::VectorXd::Map(input.barrier.data(),
                                         input.barrier.size());
   }
}

/*!
 * @brief Calculates the action and bubble profile using the shooting method
 *
 * @param input settings to be used to configure the potential and solver
 * @return a tuple containing the calculation action and, optionally, the
 *         field profiles
 */
std::tuple<double, Field_profiles> run_shooting_profiler(
   const Bubble_profiler_inputs& input)
{
   Algebraic_potential potential(input.fields, input.potential);
   const int n_fields = potential.get_number_of_fields();

   Eigen::VectorXd true_vacuum_loc(n_fields);
   Eigen::VectorXd false_vacuum_loc(n_fields);
   Eigen::VectorXd barrier(n_fields);

   initialize_extrema(potential, input, true_vacuum_loc, false_vacuum_loc,
                      barrier);

   unsigned int options = Shooting::Solver_options::Compute_action;
   if (input.write_profiles) {
      options = options | Shooting::Solver_options::Compute_profile;
   }
   Shooting solver;

   solver.set_bisection_precision_bits(input.shoot_bisect_bits);
   solver.set_action_arrived_rel(input.action_arrived_rel);
   solver.set_shooting_abs_tol(input.shoot_ode_abs);
   solver.set_shooting_rel_tol(input.shoot_ode_rel);
   solver.set_action_abs_tol(input.action_ode_abs);
   solver.set_action_rel_tol(input.action_ode_rel);
   solver.set_drho_frac(input.drho_frac);
   solver.set_evolve_change_rel(input.evolve_change_rel);
   solver.set_bisect_lambda_max(input.bisect_lambda_max);
   solver.set_max_iterations(input.iter_max);
   solver.set_max_periods(input.periods_max);
   solver.set_f_y_max(input.f_y_max);
   solver.set_f_y_min(input.f_y_min);
   solver.set_y_max(input.y_max);

   solver.solve(potential, false_vacuum_loc(0),
                true_vacuum_loc(0), barrier(0), input.n_dims, options);

   Field_profiles profile;
   if (input.write_profiles) {
      profile = solver.get_bubble_profile();
   }

   return std::make_tuple(solver.get_euclidean_action(), profile);
}

/*!
 * @brief Calculates the action and bubble profile using the perturbative method
 *
 * The provided BVP solver is used to solve the system of ODEs that
 * must be integrated to compute successive perturbations to the initial
 * solution ansatz.
 *
 * @param input settings to be used to configure the potential and solver
 * @return a tuple containing the calculation action and, optionally, the
 *         field profiles
 */
template <class Profiler>
std::tuple<double, Field_profiles> run_perturbative_profiler(
   const Bubble_profiler_inputs& input, Profiler& profiler)
{
   Algebraic_potential potential(input.fields, input.potential);

   const int n_fields = potential.get_number_of_fields();

   Eigen::VectorXd false_vacuum_loc(n_fields);
   Eigen::VectorXd true_vacuum_loc(n_fields);

   initialize_extrema(potential, input, true_vacuum_loc, false_vacuum_loc);

   if (input.domain_start >= 0.) {
      profiler.set_domain_start(input.domain_start);
   }
   if (input.domain_end >= 0.) {
      profiler.set_domain_end(input.domain_end);
   }
   profiler.set_initial_step_size(input.initial_step_size);
   profiler.set_interpolation_points_fraction(input.interpolation_fraction);
   profiler.set_false_vacuum_loc(false_vacuum_loc);
   profiler.set_true_vacuum_loc(true_vacuum_loc);
   profiler.set_number_of_dimensions(input.n_dims);

   auto root_finder = std::make_shared<GSL_root_finder<Eigen::Dynamic> >();
   profiler.set_root_finder(root_finder);

   std::shared_ptr<Profile_guesser> guesser;
   if (input.shooting_ansatz) {
      guesser = std::make_shared<Shooting_profile_guesser>();
   } else if (!input.ansatz_file.empty()) {
      std::ifstream ansatz_input(input.ansatz_file);
      guesser = std::make_shared<Instream_profile_guesser>(ansatz_input);
   } else {
      guesser = std::make_shared<Kink_profile_guesser>();
   }
   profiler.set_initial_guesser(guesser);

   auto convergence_tester = std::make_shared<Relative_convergence_tester>(
      input.rtol_action, input.rtol_fields);
   if (input.max_iterations >= 0) {
      convergence_tester->set_max_iterations(input.max_iterations);
   }
   profiler.set_convergence_tester(convergence_tester);

   if (!input.output_path.empty()) {
      Profile_observer observer(input.fields, input.output_path, potential);
      profiler.calculate_bubble_profile(potential, observer);
   } else {
      profiler.calculate_bubble_profile(potential);
   }

   Field_profiles profiles;
   if (input.write_profiles) {
      profiles = profiler.get_bubble_profile();
   }

   return std::make_tuple(profiler.get_euclidean_action(), profiles);
}

/*!
 * @brief Calculates the action and bubble profile using the perturbative method
 *
 * The BVP solver used in the calculation is determined from the provided
 * setting for the integration algorithm.
 *
 * @param input settings to be used to configure the potential and solver
 * @param bvp_solver the BVP solver to be used
 * @return a tuple containing the calculation action and, optionally, the
 *         field profiles
 */
std::tuple<double, Field_profiles> run_perturbative_profiler(
   const Bubble_profiler_inputs& input)
{
   switch (input.algorithm) {
   case Integration_algorithm::Runge_kutta_4: {
      RK4_perturbative_profiler profiler;
      return run_perturbative_profiler(input, profiler);
   }
   case Integration_algorithm::Euler: {
      Euler_perturbative_profiler profiler;
      return run_perturbative_profiler(input, profiler);
   }
   default: {
      throw Setup_error("unrecognized integration method");
   }
   }
}

/*!
 * @brief Calculates the action and bubble profile for the given potential
 *
 * @param input settings to be used to configure the potential and solver
 */
void run_profiler(const Bubble_profiler_inputs& input)
{
   const auto n_fields = input.fields.size();

   std::tuple<double, Field_profiles> result;
   if (n_fields == 1 && !input.use_perturbative) {
      result = run_shooting_profiler(input);
   } else {
      result = run_perturbative_profiler(input);
   }

   if (!input.output_file.empty()) {
      if (input.output_file == "-") {
         write_action(std::get<0>(result), std::cout);
         if (input.write_profiles) {
            write_profiles(input.fields, std::get<1>(result), std::cout);
         }
      } else {
         std::ofstream ostr(input.output_file);
         write_action(std::get<0>(result), ostr);
         if (input.write_profiles) {
            write_profiles(input.fields, std::get<1>(result), ostr);
         }
      }
   }
}

void configure_logging(const Bubble_profiler_inputs& input)
{
   auto& logging_manager = logging::Logging_manager::get_manager();

   if (input.verbose) {
      logging_manager.set_minimum_log_level(logging::Log_level::Trace);
   } else {
      logging_manager.set_minimum_log_level(logging::Log_level::Warning);
   }
}

} // namespace BubbleProfiler

int main(int argc, const char* argv[])
{
   using namespace BubbleProfiler;

   try {

      Bubble_profiler_inputs input = parse_cmd_line_args(argc, argv);

      configure_logging(input);

      std::cout << "Potential: " << input.potential << std::endl;
      for (const auto& f: input.fields) {
         std::cout << "Field: " << f << std::endl;
      }

      if (!input.output_file.empty()) {
         // check for existence of output file if not allowing overwriting
         if (!input.force_output && input.output_file != "-") {
            fs::path output_file(input.output_file);
            if (fs::exists(output_file)) {
               std::cerr << "Error: output file already exists.\n";
               std::cerr << "Please specify a different output file, or use\n";
               std::cerr << "the --force-output option.\n";
               return EXIT_FAILURE;
            }
         }
      }

      if (!input.output_path.empty()) {
         fs::path output_path(input.output_path);
         if (!fs::exists(output_path) || input.force_output) {
            fs::create_directory(output_path);
         } else if (fs::is_regular_file(output_path)) {
            throw Setup_error("Error: output path is a regular file");
         } else {
            std::cerr << "Error: output directory already exists.\n";
            std::cerr << "Please specify a different output path, or use\n";
            std::cerr << "the --force-output option.\n";
            return EXIT_FAILURE;
         }
      }

      run_profiler(input);

   } catch (const Setup_error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
   } catch (const Error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
   } catch (const std::runtime_error& e) {
      std::cerr << "runtime error:  " <<  e.what() << std::endl;
      return EXIT_FAILURE;
   } catch (const std::exception& e) {
      std::cerr << "Exception:  " <<  e.what() << std::endl;
      return EXIT_FAILURE;
   } catch (const std::string& a) {
      std::cerr << a;
      return EXIT_FAILURE;
   } catch (...) {
      std::cerr << "Unknown type of exception caught.\n";
      return EXIT_FAILURE;
   }

   return 0;
}
