#include "catch/catch.hpp"

#include "algebraic_potential.hpp"
#include "error.hpp"
#include "euclidean_action.hpp"
#include "field_profiles.hpp"
#include "kink_profile_guesser.hpp"
#include "mock_potentials.hpp"

#include <cassert>
#include <cmath>
#include <string>
#include <vector>

using namespace BubbleProfiler;

TEST_CASE("test action calculated matches analytic value",
          "[euclidean_action][1D potentials]")
{
   SECTION("matches analytic result for Fubini potential in d = 4 dimensions")
   {
      const double abs_tol = 1.e-4;

      const int d = 4;

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

      const double field_tol = 1.e-7;
      const double step_size = 1.e-2;
      const double r_min = 0.;

      const std::size_t max_intervals = 1000;
      const double integration_abs_tol = 1.e-2 * abs_tol;
      const double integration_rel_tol = 1.e-2 * abs_tol;

      for (int i = 0; i < n_u_pts; ++i) {
         for (int j = 0; j < n_v_pts; ++j) {
            for (int k = 0; k < n_n_pts; ++k) {
               const double u = u_min + i * u_incr;
               const double v = v_min + j * v_incr;
               const double n = n_min + k * n_incr;

               Generalized_fubini_potential potential(u, v, n);

               // take as integration region 0 <= r <= r_max,
               // with r_max the smallest value of r such that
               // the exact solution is less than field_tol
               const double r_max
                  = std::sqrt((std::pow(field_tol, -1. / n) - v) / u);
               const int n_grid_points = 1 + static_cast<int>(
                  std::ceil(r_max / step_size));
               const Eigen::VectorXd grid(
                  Eigen::VectorXd::LinSpaced(n_grid_points, r_min, r_max));

               const Eigen::VectorXd solution(
                  potential.get_bounce_solution_at(grid));
               Field_profiles profile(grid, solution);
               profile.set_number_of_dimensions(d);

               const double exact_action = potential.get_action();

               const double potential_action
                  = calculate_potential_action(
                     potential, profile,
                     max_intervals, integration_rel_tol,
                     integration_abs_tol);
               const double kinetic_action
                  = calculate_kinetic_action(
                     potential, profile,
                     max_intervals, integration_rel_tol,
                     integration_abs_tol);
               const double action
                  = calculate_action(
                     potential, profile,
                     max_intervals, integration_rel_tol,
                     integration_abs_tol);

               // area_n_sphere(3) = 19.7392
               const double max_abs_err = 20. * integration_abs_tol;
               const double max_rel_err = 20. * integration_rel_tol * std::abs(exact_action);

               const bool potential_action_close =
                  (std::abs(potential_action - exact_action) < max_abs_err ||
                   std::abs(potential_action - exact_action) < max_rel_err);
               const bool kinetic_action_close =
                  (std::abs(kinetic_action - exact_action) < max_abs_err ||
                   std::abs(kinetic_action - exact_action) < max_rel_err);
               const bool full_action_close =
                  (std::abs(action - exact_action) < max_abs_err ||
                   std::abs(action - exact_action) < max_rel_err);

               CHECK(potential_action_close);
               CHECK(kinetic_action_close);
               CHECK(full_action_close);
            }
         }
      }
   }

   SECTION("matches analytic result for logarithmic potential in d = 4 dimensions")
   {
      const double abs_tol = 1.e-3;

      const int d = 4;

      const double m_min = 0.1;
      const double m_max = 10.;
      const int n_m_pts = 10;
      const double m_incr = (m_max - m_min) / (n_m_pts - 1.);

      const double w_min = 0.1;
      const double w_max = 10.;
      const int n_w_pts = 10;
      const double w_incr = (w_max - w_min) / (n_w_pts - 1.);

      const double field_tol = 1.e-10;
      const double step_size = 1.e-2;
      const double r_min = 0.;

      const std::size_t max_intervals = 1000;
      const double integration_abs_tol = 1.e-3 * abs_tol;
      const double integration_rel_tol = 1.e-3 * abs_tol;

      for (int i = 0; i < n_m_pts; ++i) {
         for (int j = 0; j < n_w_pts; ++j) {
            const double m = m_min + i * m_incr;
            const double w = w_min + j * w_incr;

            Solvable_logarithmic_potential potential(m, w);

            // take as integration region 0 <= r <= r_max,
            // with r_max the smallest value of r such that
            // the exact solution is less than field_tol
            const double r_max
               = std::sqrt(2. * (2. - std::log(field_tol / w)) / (m * m));
            const int n_grid_points = 1 + static_cast<int>(
               std::ceil(r_max / step_size));
            const Eigen::VectorXd grid(
               Eigen::VectorXd::LinSpaced(n_grid_points, r_min, r_max));

            const Eigen::VectorXd solution(
               potential.get_bounce_solution_at(grid));
            Field_profiles profile(grid, solution);
            profile.set_number_of_dimensions(d);

            const double exact_action = potential.get_action();

            const double potential_action
               = calculate_potential_action(
                  potential, profile,
                  max_intervals, integration_rel_tol,
                  integration_abs_tol);
            const double kinetic_action
               = calculate_kinetic_action(
                  potential, profile,
                  max_intervals, integration_rel_tol,
                  integration_abs_tol);
            const double action
               = calculate_action(
                  potential, profile,
                  max_intervals, integration_rel_tol,
                  integration_abs_tol);

            // area_n_sphere(3) = 19.7392
            const double max_abs_err = 20. * integration_abs_tol;
            const double max_rel_err = 20. * integration_rel_tol * std::abs(exact_action);

            const bool potential_action_close =
               (std::abs(potential_action - exact_action) < max_abs_err ||
                std::abs(potential_action - exact_action) < max_rel_err);
            const bool kinetic_action_close =
               (std::abs(kinetic_action - exact_action) < max_abs_err ||
                std::abs(kinetic_action - exact_action) < max_rel_err);
            const bool full_action_close =
               (std::abs(action - exact_action) < max_abs_err ||
                std::abs(action - exact_action) < max_rel_err);

            CHECK(potential_action_close);
            CHECK(kinetic_action_close);
            CHECK(full_action_close);
         }
      }
   }
}

TEST_CASE("test action calculation for two field case",
          "[euclidean_action][2D potentials]")
{
   constexpr int num_fields = 2;
   const std::vector<std::string> field_names{"x", "y"};
   const std::string potential_expr = "(x^2 + y^2) * (1.8 * (x - 1)^2 + 0.2 * (y - 1)^2 - 0.4)";

   SECTION("correct action for reference potential")
   {
      Algebraic_potential potential(field_names, potential_expr);

      Eigen::VectorXd global_min(num_fields);
      Eigen::VectorXd local_min(num_fields);

      global_min << 1.0463723898865282536, 1.6634936652844771743;
      local_min << 0., 0.;

      Kink_profile_guesser ansatz;

      const int n_dimensions = 3;
      const double domain_start = 1.e-8;
      const double domain_end = 9.;
      const double step_size = 1.e-3;
      const double interp_fraction = 1.0;
      Field_profiles profiles =
         ansatz.get_profile_guess(potential, global_min - local_min,
                                  n_dimensions, domain_start, domain_end,
                                  step_size, interp_fraction);

      const double action = calculate_full_action(potential, profiles);

      CHECK(std::abs(action - 30.924575383899476577) < 1.e-3);
   }
}
