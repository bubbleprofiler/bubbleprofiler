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

#include "gaussian_potential.hpp"
#include "profile_guesser.hpp"
#include "kink_profile_guesser.hpp"
#include "field_profiles.hpp"
#include "perturbative_profiler.hpp"

#include <boost/program_options.hpp>
#include <iostream>
#include <Eigen/Core>

namespace po = boost::program_options;

struct Gaussian_inputs {
   double gamma;
   double lambda;
   int n_fields;
   std::string output_path{""};
};

Gaussian_inputs parse_cmd_line_args(int argc, const char* argv[])
{
   Gaussian_inputs inputs;

   po::options_description desc("gaussian_potential_test options");

   desc.add_options()
           ("gamma", po::value<double>()->required(), "scale of secondary gaussian relative to primary")
           ("lambda", po::value<double>()->required(), "geometric distance between gaussians")
           ("n_fields", po::value<int>()->required(), "number of fields in potential")
           ("output", po::value<std::string>(), "output path for field profiles")
           ("help", "print help message");

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

   inputs.gamma = vm["gamma"].as<double>();
   inputs.lambda = vm["lambda"].as<double>();
   inputs.n_fields = vm["n_fields"].as<int>();

   if (vm.count("output")) {
      inputs.output_path = vm["output"].as<std::string>();
   }

   return inputs;
}

int main(int argc, const char* argv[])
{
   using namespace BubbleProfiler;

   Gaussian_inputs input = parse_cmd_line_args(argc, argv);

   Gaussian_potential potential(input.gamma, input.lambda, input.n_fields);

   std::shared_ptr<Kink_profile_guesser> kink_guesser
      = std::make_shared<Kink_profile_guesser>();
   std::shared_ptr<Profile_guesser> guesser(kink_guesser);

   Eigen::VectorXd false_vacuum(Eigen::VectorXd::Zero(input.n_fields));
   Eigen::VectorXd true_vacuum(Eigen::VectorXd::Zero(input.n_fields));

   for (int i = 0; i < input.n_fields; ++i) {
      true_vacuum(i) = input.lambda / std::sqrt(input.n_fields);
   }

   Field_profiles ansatz = guesser->get_profile_guess(
      potential, true_vacuum, input.n_fields, -1.0, -1.0, 1.e-3, 1);

   double alpha = kink_guesser->get_alpha();
   std::cout << "Alpha: " << alpha << '\n';

   RK4_perturbative_profiler profiler;

   profiler.set_true_vacuum_loc(true_vacuum);
   profiler.set_false_vacuum_loc(false_vacuum);
   profiler.set_initial_guesser(guesser);

   auto convergence_tester = std::make_shared<Relative_convergence_tester>();
   convergence_tester->set_max_iterations(30);
   profiler.set_convergence_tester(convergence_tester);

   // Build string of field names
   std::vector<std::string> field_names;
   for (int i = 0; i < input.n_fields; ++i) {
      field_names.push_back(std::to_string(i));
   }

   if (!input.output_path.empty()) {
      Profile_observer observer(field_names, input.output_path, potential);
      profiler.calculate_bubble_profile(potential, observer);
   }
   else {
      profiler.calculate_bubble_profile(potential);
   }

   std::cout << "Action: " << profiler.get_euclidean_action() << '\n';

   return 0;
}
