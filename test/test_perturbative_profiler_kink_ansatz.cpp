#include "catch/catch.hpp"

#include "basic_logger.hpp"
#include "field_profiles.hpp"
#include "gsl_root_finder.hpp"
#include "kink_profile_guesser.hpp"
#include "logging_manager.hpp"
#include "mock_potentials.hpp"
#include "relative_convergence_tester.hpp"
#include "perturbative_profiler.hpp"
#include "algebraic_potential.hpp"

#include <Eigen/Core>

#include <string>
#include <tuple>
#include <vector>

using namespace BubbleProfiler;

namespace {

logging::Basic_logger logger{};

std::tuple<bool, double, Field_profiles> run_profiler(
   Potential& potential,
   const Eigen::VectorXd& global_minimum,
   const Eigen::VectorXd& local_minimum,
   double domain_start, double domain_end, double step_size,
   double interp_fraction, double action_tol = 1.e-3, double fields_tol = 1.e-3,
   int d = 3)
{
   using Root_finder = GSL_root_finder<Eigen::Dynamic>;

   logging::Logging_manager::get_manager().set_minimum_log_level(
      logging::Log_level::Trace);
   logging::Basic_logger logger;

   auto root_finder = std::make_shared<Root_finder>();

   const auto observer = [](const Field_profiles&, const Field_profiles&) {};

   RK4_perturbative_profiler profiler;
   if (domain_start > 0) {
      profiler.set_domain_start(domain_start);
   }
   if (domain_end > 0) {
      profiler.set_domain_end(domain_end);
   }

   profiler.set_initial_step_size(step_size);
   profiler.set_interpolation_points_fraction(interp_fraction);
   profiler.set_initial_guesser(std::make_shared<Kink_profile_guesser>());
   profiler.set_convergence_tester(std::make_shared<Relative_convergence_tester>(
                                      action_tol, fields_tol));

   profiler.set_root_finder(root_finder);
   profiler.set_true_vacuum_loc(global_minimum);
   profiler.set_false_vacuum_loc(local_minimum);
   profiler.set_number_of_dimensions(d);

   bool found_solution = true;
   double action = 0.;
   Field_profiles profiles_solution;
   try {
      profiler.calculate_bubble_profile(potential, observer);
      action = profiler.get_euclidean_action();
      profiles_solution = profiler.get_bubble_profile();
   } catch (const Error& e) {
      found_solution = false;
      logger.log_message(logging::Log_level::Warning, e.what());
   } catch (...) {
      found_solution = false;
   }

   return std::make_tuple(found_solution, action, profiles_solution);
}

} // anonymous namespace

TEST_CASE("test kink ansatz for 1D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][1D potentials][polynomial potentials]")
{
   const int n_fields = 1;
   const std::vector<std::string> field_names{"x"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1.0;
   const double default_domain_end = -1.0;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 1D example 1")
   {
      const std::string potential_expr
         = "0.1*((-x + 2)^4 - 14*(-x + 2)^2 + 24*(-x + 2))";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 5.;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 1D example 2")
   {
      const std::string potential_expr
         = "x^4 - 14*x^2 + 24*x";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << -3.;
      local_minimum << 2.;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 1D example 3")
   {
      const std::string potential_expr
         = "-1. * ((4. * 0.55527638191 - 3.) / 2. * x^2 "
         "+ x^3 - 0.55527638191 * x^4)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on standard scaled potential")
   {
      const double alpha_min = 0.5141;
      const double alpha_max = 0.745;
      const int n_alpha_pts = 50;
      const double alpha_incr = (alpha_max - alpha_min) / (n_alpha_pts - 1.);

      const double E_min = -20.;
      const double E_max = -0.1;
      const int n_E_pts = 10;
      const double E_incr = (E_max - E_min) / (n_E_pts - 1.);

      int fails = 0;
      int total = n_alpha_pts * n_E_pts;

      for (int i = 0; i < n_alpha_pts; ++i) {
         for (int j = 0; j < n_E_pts; ++j) {
            const double alpha = alpha_min + i * alpha_incr;
            const double E = E_min + j * E_incr;
            Scaled_1D_potential potential(alpha, E);

            global_minimum << potential.get_global_minimum_location();
            local_minimum << potential.get_local_minimum_location();

            // Compute the polynomial coefficients for logging / reproducibility
            double c_quadratic = 0.5*E*(4*alpha - 3);
            double c_cubic = E;
            double c_quartic = (-1.0)*alpha;

            std::string log_potential =
                    std::to_string(c_quadratic) + "*x^2 " +
                    std::to_string(c_cubic) + "*x^3 " +
                    std::to_string(c_quartic) + "*x^4";

            logger.log_message(logging::Log_level::Trace, "Test potential: " + log_potential +
                    ", local_min: " + std::to_string(local_minimum(0)) + ", global_min: " +
                    std::to_string(global_minimum(0)));

            const auto result = run_profiler(
               potential, global_minimum, local_minimum,
               default_domain_start, default_domain_end, default_step_size,
               default_interp_fraction, default_action_tol, default_fields_tol);

            bool res = std::get<0>(result);
            if (!res) fails++;

            CHECK(res == true);
         }
      }
      logger.log_message(logging::Log_level::Trace, std::to_string(fails) +
              " out of " + std::to_string(total) + " tests failed.");

   }
}

TEST_CASE("test kink ansatz for 2D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][2D potentials][polynomial potentials]")
{
   const int n_fields = 2;
   const std::vector<std::string> field_names{"x", "y"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 2D example 1")
   {
      const std::string potential_expr
         = "(x^2 + y^2) * (1.8 * (x - 1)^2 + 0.2 * (y - 1)^2 - 0.4)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.0463723898865282536, 1.6634936652844771743;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 2D example 3")
   {
      const std::string potential_expr
         = "(x^2 + y^2) * (1.8 * (x - 1)^2 + 0.2 * (y - 1)^2 - 0.15)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.026807408218724, 1.3071345669072432;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, 7., default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 2D example 4")
   {
      const std::string potential_expr
         = "(x^2 + y^2) * (1.8 * (x - 1)^2 + 0.2 * (y - 1)^2 - 0.6)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.05502804185546, 1.884733882243574;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for 3D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][3D potentials][polynomial potentials]")
{
   const int n_fields = 3;
   const std::vector<std::string> field_names{"x", "y", "z"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 3D example 1")
   {
      const std::string potential_expr = "(-0.284821 + 0.684373 * (-1 + x)^2"
         " + 0.181928 * (-1 + y)^2 + 0.295089 * (-1 + z)^2)"
         " * (x^2 + y^2 + z^2)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.08187697579, 1.39800171207, 1.21288400399;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for 4D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][4D potentials][polynomial potentials]")
{
   const int n_fields = 4;
   const std::vector<std::string> field_names{"t", "x", "y", "z"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 4D example 1")
   {
      const std::string potential_expr = "(-0.258889 + 0.534808*(-1 + t)^2"
         " + 0.77023*(-1 + x)^2 + 0.838912*(-1 + y)^2"
         " + 0.00517238*(-1 + z)^2)*(t^2 + x^2 + y^2 + z^2)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.13435238329, 1.08960762926, 1.08167225882,
         -0.0889181410773;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for 5D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][5D potentials][polynomial potentials]")
{
   const int n_fields = 5;
   const std::vector<std::string> field_names{"s", "t", "x", "y", "z"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 5D example 1")
   {
      const std::string potential_expr = "(-0.658889 + 0.4747*(-1 + s)^2"
         " + 0.234808*(-1 + t)^2 + 0.57023*(-1 + x)^2 + 0.138912*(-1 + y)^2"
         " + 0.517238*(-1 + z)^2)*(s^2 + t^2 + x^2 + y^2 + z^2)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.14464229776, 1.34312068679, 1.11756180403,
         1.7600125205, 1.13118629431;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for 6D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][6D potentials][polynomial potentials]")
{
   const int n_fields = 6;
   const std::vector<std::string> field_names{"p", "s", "t", "x", "y", "z"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 6D example 1")
   {
      const std::string potential_expr = "(-0.658889 + 0.34234*(-1 + p)^2"
         " + 0.4747*(-1 + s)^2 + 0.234808*(-1 + t)^2 + 0.57023*(-1 + x)^2"
         " + 0.138912*(-1 + y)^2 + 0.517238*(-1 + z)^2)*(p^2 + s^2 + t^2"
         " + x^2 + y^2 + z^2)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.1939713966, 1.1327090551, 1.31037255877,
         1.10807369483, 1.66769539692, 1.12048004437;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for 7D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][7D potentials][polynomial potentials]")
{
   const int n_fields = 7;
   const std::vector<std::string> field_names{
      "r", "p", "s", "t", "x", "y", "z"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 7D example 1")
   {
      const std::string potential_expr = "(-0.658889 + 0.5233*(-1 + r)^2"
         " + 0.34234*(-1 + p)^2 + 0.4747*(-1 + s)^2 + 0.234808*(-1 + t)^2"
         " + 0.57023*(-1 + x)^2 + 0.138912*(-1 + y)^2"
         " + 0.517238*(-1 + z)^2)*(r^2 + p^2 + s^2 + t^2 + x^2 + y^2 + z^2)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.18014082812, 1.11093650195, 1.12369917273,
         1.2862525169, 1.10088532217, 1.60302762583, 1.11238278813;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for 8D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][8D potentials][polynomial potentials]")
{
   const int n_fields = 8;
   const std::vector<std::string> field_names{
      "q", "r", "p", "s", "t", "x", "y", "z"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 8D example 1")
   {
      const std::string potential_expr = "(-0.658889 + 0.2434*(-1 + q)^2"
         " + 0.5233*(-1 + r)^2 + 0.34234*(-1 + p)^2 + 0.4747*(-1 + s)^2"
         " + 0.234808*(-1 + t)^2 + 0.57023*(-1 + x)^2"
         " + 0.138912*(-1 + y)^2 + 0.517238*(-1 + z)^2)*(q^2"
         " + r^2 + p^2 + s^2 + t^2 + x^2 + y^2 + z^2)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.16287399482, 1.2453229329, 1.10087000448,
         1.11235744949, 1.25660315558, 1.09180625051, 1.52712119088,
         1.10217298316;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for 10D polynomial potentials",
          "[perturbative_profiler][kink_profile_guesser][10D potentials][polynomial potentials]")
{
   const int n_fields = 10;
   const std::vector<std::string> field_names{
      "b", "a", "q", "r", "p", "s", "t", "x", "y", "z"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = -1;
   const double default_domain_end = -1;
   const double default_step_size = 1.e-1;
   const int default_interp_fraction = 1.0;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 10D example 1")
   {
      const std::string potential_expr = "(-0.658889 + 0.28765*(-1 + b)^2"
         " + 0.7345*(-1 + a)^2 + 0.2434*(-1 + q)^2 + 0.5233*(-1 + r)^2"
         " + 0.34234*(-1 + p)^2 + 0.4747*(-1 + s)^2 + 0.234808*(-1 + t)^2"
         " + 0.57023*(-1 + x)^2 + 0.138912*(-1 + y)^2"
         " + 0.517238*(-1 + z)^2)*(a^2 + b^2 + q^2 + r^2 + p^2 + s^2"
         " + t^2 + x^2 + y^2 + z^2)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.06123721286, 1.17280494026, 1.14129805239,
         1.21084465797, 1.08813022635, 1.09803761364, 1.22025913169,
         1.08029470989, 1.43907456395, 1.08925529772;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol);

      REQUIRE(std::get<0>(result) == true);
   }
}

TEST_CASE("test kink ansatz for solvable logarithmic potential",
          "[perturbative_profiler][kink_profile_guesser][1D potentials][solvable potentials]")
{
   constexpr int n_fields = 1;
   const std::vector<std::string> field_names{"x"};

   const double default_action_tol = 1.e-3;
   const double default_fields_tol = 1.e-3;
   const double default_domain_start = 10.e-6;
   const double default_domain_end = -1.;
   const double default_step_size = 1.e-2;
   const int default_interp_fraction = 1.0;
   const int d = 4;

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on solvable logarithmic potential, d = 4")
   {
      const double action_rel_tol = 1.e-2;
      const double m = 0.5;
      const double w = 1.;

      logging::Logging_manager::get_manager().set_minimum_log_level(
      logging::Log_level::Trace);

      logger.log_message(logging::Log_level::Trace, "m =  "
                         + std::to_string(m));
      logger.log_message(logging::Log_level::Trace, "w =  "
                         + std::to_string(w));

      Solvable_logarithmic_potential potential(m, w);

      global_minimum << std::exp(2.) * w;
      local_minimum << potential.get_local_minimum_location();

      const auto result = run_profiler(
         potential, global_minimum, local_minimum,
         default_domain_start, default_domain_end, default_step_size,
         default_interp_fraction, default_action_tol, default_fields_tol, d);

      CHECK(std::get<0>(result));

      if (std::get<0>(result)) {
         const auto action = std::get<1>(result);
         const auto exact_action = potential.get_action();
         logger.log_message(logging::Log_level::Trace, "action: "
                            + std::to_string(action));
         logger.log_message(logging::Log_level::Trace, "exact_action: "
                            + std::to_string(exact_action));
         const double action_diff = std::abs(
            (action - exact_action)
            / (0.5 * (action + exact_action)));

         CHECK(action_diff < action_rel_tol);
      }
   }

   // The following scan currently contains points at which the
   // perturbative profiler fails to produce a solution:
   // SECTION("profiler succeeds on solvable logarithmic potential scan, d = 4")
   // {
   //    const double action_rel_tol = 1.e-2;

   //    const double m_min = 0.1;
   //    const double m_max = 10.;
   //    const int n_m_pts = 10;
   //    const double m_incr = (m_max - m_min) / (n_m_pts - 1.);

   //    const double w_min = 0.1;
   //    const double w_max = 10.;
   //    const int n_w_pts = 10;
   //    const double w_incr = (w_max - w_min) / (n_w_pts - 1.);

   //    int fails = 0;
   //    int total = n_m_pts * n_w_pts;
   //    bool check_1 = false;
   //    bool check_2 = false;

   //    for (int i = 0; i < n_m_pts; ++i) {
   //       for (int j = 0; j < n_w_pts; ++j) {
   //          const double m = m_min + i * m_incr;
   //          const double w = w_min + j * w_incr;

   //          logger.log_message(logging::Log_level::Trace, "m =  "
   //                             + std::to_string(m));
   //          logger.log_message(logging::Log_level::Trace, "w =  "
   //                             + std::to_string(w));

   //          Solvable_logarithmic_potential potential(m, w);

   //          // take as integration region 0 <= r <= r_max,
   //          // with r_max the smallest value of r such that
   //          // the exact solution is less than field_tol
   //          const double r_max
   //                  = std::sqrt(2. * (2. - std::log(field_tol / w)) / (m * m));

   //          global_minimum << std::exp(2.) * w;
   //          local_minimum << potential.get_local_minimum_location();

   //          const auto result = run_profiler(
   //             potential, global_minimum, local_minimum,
   //             default_domain_start, r_max, default_step_size,
   //             default_interp_fraction, default_action_tol, default_fields_tol, d);

   //          check_1 = std::get<0>(result);

   //          CHECK(check_1);

   //          if (std::get<0>(result)) {
   //             const auto action = std::get<1>(result);
   //             const auto exact_action = potential.get_action();
   //             logger.log_message(logging::Log_level::Trace, "action: "
   //                                + std::to_string(action));
   //             logger.log_message(logging::Log_level::Trace, "exact_action: "
   //                                + std::to_string(exact_action));
   //             const double action_diff = std::abs(
   //                (action - exact_action)
   //                / (0.5 * (action + exact_action)));

   //             check_2 = action_diff < action_rel_tol;

   //             CHECK(check_2);
   //          }

   //          if (!(check_1 && check_2)) fails++;
   //       }
   //    }
   //    logger.log_message(logging::Log_level::Trace, std::to_string(fails) +
   //                                                  " out of " + std::to_string(total) +
   //            " tests failed.");
   // }
}
