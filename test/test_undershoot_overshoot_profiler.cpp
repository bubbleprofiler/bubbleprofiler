#include "catch/catch.hpp"

#include "field_profiles.hpp"
#include "mock_potentials.hpp"
#include "shooting.hpp"
#include "algebraic_potential.hpp"

#include <string>
#include <vector>

using namespace BubbleProfiler;

namespace {

std::tuple<bool, double, Field_profiles> run_profiler(
   const Potential& potential,
   const Eigen::VectorXd& global_minimum,
   const Eigen::VectorXd& local_minimum,
   const Eigen::VectorXd& local_maximum,
   int n_spatial_dimensions,
   unsigned int options = (Shooting::Solver_options::Compute_action |
                           Shooting::Solver_options::Compute_profile))
{

   const double false_min = local_minimum(0);
   const double true_min = global_minimum(0);
   const double barrier = local_maximum(0);

   bool found_solution = true;
   double action = 0.;
   Field_profiles profiles_solution;
   try {
      Shooting solution;
      solution.set_bisection_precision_bits(10);
      solution.set_shooting_abs_tol(1.e-8);
      solution.set_shooting_rel_tol(1.e-8);
      solution.set_action_abs_tol(1.e-8);
      solution.set_action_rel_tol(1.e-8);
      solution.solve(potential, false_min, true_min, barrier,
                     n_spatial_dimensions, options);
      if ((options & Shooting::Solver_options::Compute_action) != 0) {
         action = solution.get_euclidean_action();
      }
      if ((options & Shooting::Solver_options::Compute_profile) != 0) {
         profiles_solution = solution.get_bubble_profile();
      }
   } catch (...) {
      found_solution = false;
   }

   return std::make_tuple(found_solution, action, profiles_solution);
}

} // anonymous namespace

TEST_CASE("test undershoot-overshoot solver for 1D polynomial potentials",
          "[undershoot_overshoot_profiler][1D potentials][polynomial potentials]")
{
   constexpr int n_fields = 1;
   const std::vector<std::string> field_names{"x"};

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_maximum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on 1D example 1, d = 3")
   {
      const std::string potential_expr
         = "0.1*((-x + 2)^4 - 14*(-x + 2)^2 + 24*(-x + 2))";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 5.;
      local_maximum << 1.;

      const int d = 3;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum, local_maximum, d);

      CHECK(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 1D example 1, d = 4")
   {
      const std::string potential_expr
         = "0.1*((-x + 2)^4 - 14*(-x + 2)^2 + 24*(-x + 2))";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 5.;
      local_maximum << 1.;

      const int d = 4;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum, local_maximum, d);

      CHECK(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 1D example 2, d = 3")
   {
      const std::string potential_expr
         = "x^4 - 14*x^2 + 24*x";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << -3;
      local_minimum << 2.;
      local_maximum << 1.;

      const int d = 3;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum, local_maximum, d);

      CHECK(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 1D example 2, d = 4")
   {
      const std::string potential_expr
         = "x^4 - 14*x^2 + 24*x";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << -3;
      local_minimum << 2.;
      local_maximum << 1.;

      const int d = 4;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum, local_maximum, d);

      CHECK(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 1D example 3, d = 3")
   {
      const std::string potential_expr
         = "-1. * ((4. * 0.55527638191 - 3.) / 2. * x^2 "
         "+ x^3 - 0.55527638191 * x^4)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.;
      local_minimum << 0.;
      local_maximum << 0.35067873303057423;

      const int d = 3;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum, local_maximum, d);

      CHECK(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on 1D example 3, d = 4")
   {
      const std::string potential_expr
         = "-1. * ((4. * 0.55527638191 - 3.) / 2. * x^2 "
         "+ x^3 - 0.55527638191 * x^4)";

      Algebraic_potential potential(field_names, potential_expr);

      global_minimum << 1.;
      local_minimum << 0.;
      local_maximum << 0.35067873303057423;

      const int d = 4;

      const auto result = run_profiler(
         potential, global_minimum, local_minimum, local_maximum, d);

      CHECK(std::get<0>(result) == true);
   }

   SECTION("profiler succeeds on standard scaled potential, d = 3")
   {
      const double alpha_min = 0.51;
      const double alpha_max = 0.745;
      const int n_alpha_pts = 50;
      const double alpha_incr = (alpha_max - alpha_min) / (n_alpha_pts - 1.);

      const double E_min = -20.;
      const double E_max = -0.1;
      const int n_E_pts = 10;
      const double E_incr = (E_max - E_min) / (n_E_pts - 1.);

      const int d = 3;

      for (int i = 0; i < n_alpha_pts; ++i) {
         for (int j = 0; j < n_E_pts; ++j) {
            const double alpha = alpha_min + i * alpha_incr;
            const double E = E_min + j * E_incr;
            Scaled_1D_potential potential(alpha, E);

            global_minimum << potential.get_global_minimum_location();
            local_minimum << potential.get_local_minimum_location();
            local_maximum << potential.get_local_maximum_location();

            const auto result = run_profiler(
               potential, global_minimum, local_minimum, local_maximum, d);

            CHECK(std::get<0>(result) == true);
         }
      }
   }

   SECTION("profiler succeeds on standard scaled potential, d = 4")
   {
      const double alpha_min = 0.51;
      const double alpha_max = 0.745;
      const int n_alpha_pts = 50;
      const double alpha_incr = (alpha_max - alpha_min) / (n_alpha_pts - 1.);

      const double E_min = -20.;
      const double E_max = -0.1;
      const int n_E_pts = 10;
      const double E_incr = (E_max - E_min) / (n_E_pts - 1.);

      const int d = 4;

      for (int i = 0; i < n_alpha_pts; ++i) {
         for (int j = 0; j < n_E_pts; ++j) {
            const double alpha = alpha_min + i * alpha_incr;
            const double E = E_min + j * E_incr;
            Scaled_1D_potential potential(alpha, E);

            global_minimum << potential.get_global_minimum_location();
            local_minimum << potential.get_local_minimum_location();
            local_maximum << potential.get_local_maximum_location();

            const auto result = run_profiler(
               potential, global_minimum, local_minimum, local_maximum, d);

            CHECK(std::get<0>(result) == true);
         }
      }
   }

}

TEST_CASE("test undershoot-overshoot solver for generalized Fubini potential",
          "[undershoot_overshoot_profiler][1D potentials][solvable potentials]")
{
   constexpr int n_fields = 1;
   const std::vector<std::string> field_names{"x"};

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_maximum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on generalized Fubini potential, d = 4")
   {
      const double u_min = 1.;
      const double u_max = 3.;
      const int n_u_pts = 10;
      const double u_incr = (u_max - u_min) / (n_u_pts - 1.);

      const double v_min = 1.;
      const double v_max = 3.;
      const int n_v_pts = 10;
      const double v_incr = (v_max - v_min) / (n_v_pts - 1.);

      const double n_min = 1.5;
      const double n_max = 4.;
      const int n_n_pts = 10;
      const double n_incr = (n_max - n_min) / (n_n_pts - 1.);

      const double action_rel_tol = 1.e-2;

      const int d = 4;

      for (int i = 0; i < n_u_pts; ++i) {
         for (int j = 0; j < n_v_pts; ++j) {
            for (int k = 0; k < n_n_pts; ++k) {
               const double u = u_min + i * u_incr;
               const double v = v_min + j * v_incr;
               const double n = n_min + k * n_incr;

               Generalized_fubini_potential potential(u, v, n);

               global_minimum << 2.;
               local_minimum << potential.get_local_minimum_location();
               local_maximum << potential.get_local_maximum_location();

               const auto result = run_profiler(
                  potential, global_minimum, local_minimum, local_maximum, d,
                  Shooting::Solver_options::Compute_action);

               CHECK(std::get<0>(result) == true);

               if (std::get<0>(result)) {
                  const auto action = std::get<1>(result);
                  const auto exact_action = potential.get_action();
                  const double action_diff = std::abs(
                     (action - exact_action)
                     / (0.5 * (action + exact_action)));

                  CHECK(action_diff < action_rel_tol);
               }
            }
         }
      }
   }
}

TEST_CASE("test undershoot-overshoot solver for solvable logarithmic potential",
          "[undershoot_overshoot_profiler][1D potentials][solvable potentials]")
{
   constexpr int n_fields = 1;
   const std::vector<std::string> field_names{"x"};

   Eigen::VectorXd global_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_minimum(Eigen::VectorXd::Zero(n_fields));
   Eigen::VectorXd local_maximum(Eigen::VectorXd::Zero(n_fields));

   SECTION("profiler succeeds on solvable logarithmic potential, d = 4")
   {
      const double m_min = 0.1;
      const double m_max = 10.;
      const int n_m_pts = 10;
      const double m_incr = (m_max - m_min) / (n_m_pts - 1.);

      const double w_min = 0.1;
      const double w_max = 10.;
      const int n_w_pts = 10;
      const double w_incr = (w_max - w_min) / (n_w_pts - 1.);

      const double action_rel_tol = 1.e-2;

      const int d = 4;

      for (int i = 0; i < n_m_pts; ++i) {
         for (int j = 0; j < n_w_pts; ++j) {
            const double m = m_min + i * m_incr;
            const double w = w_min + j * w_incr;

            Solvable_logarithmic_potential potential(m, w);

            global_minimum << 10. * w;
            local_minimum << potential.get_local_minimum_location();
            local_maximum << potential.get_local_maximum_location();

            const auto result = run_profiler(
               potential, global_minimum, local_minimum, local_maximum, d,
               Shooting::Solver_options::Compute_action);

            CHECK(std::get<0>(result) == true);

            if (std::get<0>(result)) {
               const auto action = std::get<1>(result);
               const auto exact_action = potential.get_action();
               const double action_diff = std::abs(
                  (action - exact_action)
                  / (0.5 * (action + exact_action)));

               CHECK(action_diff < action_rel_tol);
            }
         }
      }
   }
}
