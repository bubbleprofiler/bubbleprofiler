#include "catch/catch.hpp"

#include "euclidean_action.hpp"
#include "field_profiles.hpp"
#include "mock_potentials.hpp"
#include "algebraic_potential.hpp"
#include "error.hpp"

#define private public
#include "kink_profile_guesser.hpp"
#undef private

#include <Eigen/Core>

#include <string>
#include <vector>
#include <iostream>

using namespace BubbleProfiler;

TEST_CASE("test simple ansatz construction for 1D polynomial potentials",
          "[kink_profile_guesser][1D potentials][polynomial potentials]")
{
   constexpr int num_fields = 1;
   constexpr int n_dimensions = 3;
   constexpr double default_tol = 1.e-7;
   constexpr double default_domain_start = 1.e-8;
   constexpr double default_domain_end = 9.;
   constexpr double default_step_size = 1.e-3;
   constexpr double interp_fraction = 1.0;

   Eigen::VectorXd true_vacuum_loc(Eigen::VectorXd::Ones(num_fields));
   Eigen::VectorXd false_vacuum_loc(Eigen::VectorXd::Zero(num_fields));

   SECTION("correctly handles standard scaled potential")
   {
      const double alpha_min = 0.515;
      const double alpha_max = 0.745;
      const int n_alpha_pts = 10;
      const double alpha_incr = (alpha_max - alpha_min) / (n_alpha_pts - 1.);

      const double E_min = -20.;
      const double E_max = -0.1;
      const int n_E_pts = 10;
      const double E_incr = (E_max - E_min) / (n_E_pts - 1.);

      for (int i = 0; i < n_alpha_pts; ++i) {
         for (int j = 0; j < n_E_pts; ++j) {
            const double alpha = alpha_min + i * alpha_incr;
            const double E = E_min + j * E_incr;
            Scaled_1D_potential potential(alpha, E);

            true_vacuum_loc(0) = potential.get_global_minimum_location();
            false_vacuum_loc(0) = potential.get_local_minimum_location();
            potential.translate_origin(false_vacuum_loc);

            Kink_profile_guesser ansatz;

            Field_profiles profiles =
               ansatz.get_profile_guess(
                  potential, true_vacuum_loc - false_vacuum_loc, n_dimensions,
                  default_domain_start, default_domain_end, default_step_size,
                  interp_fraction);

            const double alpha_expected = potential.get_fitted_alpha();
            const double E_expected = potential.get_fitted_E();

            REQUIRE(std::abs(alpha_expected - alpha) < 1.e-12);
            REQUIRE(std::abs(E_expected - E) < 1.e-12);

            CHECK(std::abs(ansatz.alpha - alpha_expected) < default_tol);
            CHECK(std::abs(ansatz.aE - std::abs(E_expected)) < default_tol);
         }
      }
   }

   SECTION("correctly handles global minimum at the origin")
   {
      const double alpha_min = 0.27;
      const double alpha_max = 0.47;
      const int n_alpha_pts = 5;
      const double alpha_incr = (alpha_max - alpha_min) / (n_alpha_pts - 1.);

      const double E_min = -20.;
      const double E_max = -0.1;
      const int n_E_pts = 10;
      const double E_incr = (E_max - E_min) / (n_E_pts - 1.);

      for (int i = 0; i < n_alpha_pts; ++i) {
         for (int j = 0; j < n_E_pts; ++j) {
            const double alpha = alpha_min + i * alpha_incr;
            const double E = E_min + j * E_incr;
            Scaled_1D_potential potential(alpha, E);

            true_vacuum_loc(0) = potential.get_global_minimum_location();
            false_vacuum_loc(0) = potential.get_local_minimum_location();

            potential.translate_origin(false_vacuum_loc);

            REQUIRE(std::abs(true_vacuum_loc(0)) < 1.e-10);

            Kink_profile_guesser ansatz;

            Field_profiles profiles =
               ansatz.get_profile_guess(
                  potential, true_vacuum_loc - false_vacuum_loc, n_dimensions,
                  default_domain_start, default_domain_end, default_step_size,
                  interp_fraction);

            const double alpha_expected = potential.get_fitted_alpha();
            const double E_expected = potential.get_fitted_E();

            CHECK(std::abs(ansatz.alpha - alpha_expected) < default_tol);
            CHECK(std::abs(ansatz.aE - std::abs(E_expected)) < default_tol);
         }
      }
   }

   SECTION("correctly handles global minimum at negative field values")
   {
      const double alpha_min = -0.5;
      const double alpha_max = -0.01;
      const int n_alpha_pts = 10;
      const double alpha_incr = (alpha_max - alpha_min) / (n_alpha_pts - 1.);

      const double E_min = 0.1;
      const double E_max = 20.;
      const int n_E_pts = 10;
      const double E_incr = (E_max - E_min) / (n_E_pts - 1.);

      for (int i = 0; i < n_alpha_pts; ++i) {
         for (int j = 0; j < n_E_pts; ++j) {
            const double alpha = alpha_min + i * alpha_incr;
            const double E = E_min + j * E_incr;
            Scaled_1D_potential potential(alpha, E);

            true_vacuum_loc(0) = potential.get_global_minimum_location();
            false_vacuum_loc(0) = potential.get_local_minimum_location();

            potential.translate_origin(false_vacuum_loc);

            REQUIRE(true_vacuum_loc(0) < 0.);

            Kink_profile_guesser ansatz;

            Field_profiles profiles =
               ansatz.get_profile_guess(
                  potential, true_vacuum_loc - false_vacuum_loc, n_dimensions,
                  default_domain_start, default_domain_end, default_step_size,
                  interp_fraction);

            const double alpha_expected = potential.get_fitted_alpha();
            const double E_expected = potential.get_fitted_E();

            CHECK(std::abs(ansatz.alpha - alpha_expected) < default_tol);
            CHECK(std::abs(ansatz.aE - std::abs(E_expected)) < default_tol);
         }
      }
   }

   SECTION("correctly handles local minimum at negative field values")
   {
      const double alpha_min = 0.76;
      const double alpha_max = 1.;
      const int n_alpha_pts = 5;
      const double alpha_incr = (alpha_max - alpha_min) / (n_alpha_pts - 1.);

      const double E_min = -20.;
      const double E_max = -0.1;
      const int n_E_pts = 5;
      const double E_incr = (E_max - E_min) / (n_E_pts - 1.);

      for (int i = 0; i < n_alpha_pts; ++i) {
         for (int j = 0; j < n_E_pts; ++j) {
            const double alpha = alpha_min + i * alpha_incr;
            const double E = E_min + j * E_incr;
            Scaled_1D_potential potential(alpha, E);

            true_vacuum_loc(0) = potential.get_global_minimum_location();
            false_vacuum_loc(0) = potential.get_local_minimum_location();
            potential.translate_origin(false_vacuum_loc);

            REQUIRE(false_vacuum_loc(0) < 0.);

            Kink_profile_guesser ansatz;

            Field_profiles profiles =
               ansatz.get_profile_guess(
                  potential, true_vacuum_loc - false_vacuum_loc, n_dimensions,
                  default_domain_start, default_domain_end, default_step_size,
                  interp_fraction);

            const double alpha_expected = potential.get_fitted_alpha();
            const double E_expected = potential.get_fitted_E();

            CHECK(std::abs(ansatz.alpha - alpha_expected) < default_tol);
            CHECK(std::abs(ansatz.aE - std::abs(E_expected)) < default_tol);
         }
      }
   }

   SECTION("produces error when given degenerate minima")
   {
      const double alpha = 0.5;
      const double E = -2.;

      Scaled_1D_potential potential(alpha, E);

      // NB the Scaled_1D_potential in this case identifies
      // the minimum at the origin as the global minimum, thus
      // for the purpose of the test we swap them
      true_vacuum_loc(0) = potential.get_local_minimum_location();
      false_vacuum_loc(0) = potential.get_global_minimum_location();

      REQUIRE(std::abs(potential(true_vacuum_loc)
                       - potential(Eigen::VectorXd::Zero(num_fields))) < 1.e-12);

      Kink_profile_guesser ansatz;
      CHECK_THROWS(
         ansatz.get_profile_guess(
            potential, true_vacuum_loc - false_vacuum_loc, n_dimensions,
            default_domain_start, default_domain_end, default_step_size,
            interp_fraction));
   }

   SECTION("produces error when given coincident minima")
   {
      const double alpha = 0.6;
      const double E = -2.;

      Scaled_1D_potential potential(alpha, E);

      true_vacuum_loc(0) = 0.;
      false_vacuum_loc(0) = 0.;

      Kink_profile_guesser ansatz;

      CHECK_THROWS_AS(
         ansatz.get_profile_guess(
            potential, true_vacuum_loc - false_vacuum_loc, n_dimensions,
            default_domain_start, default_domain_end, default_step_size,
            interp_fraction), Setup_error);
   }
}

// TEST_CASE("test simple ansatz construction for 2D polynomial potentials",
//           "[kink_profile_guesser][2D potentials][polynomial potentials]")
// {
//    SECTION("correct ansatz for reference potential")
//    {
//       Algebraic_potential potential(
//          std::vector<std::string>({"x", "y"}),
//          "(x^2 + y^2)*(1.8*(x - 1)^2 + 0.2*(y-1)^2 - 0.4)");

//       Eigen::VectorXd true_vacuum_loc(2);
//       Eigen::VectorXd false_vacuum_loc(2);

//       true_vacuum_loc << 1.0463723898865282536, 1.6634936652844771743;
//       false_vacuum_loc << 0., 0.;

//       potential.translate_origin(false_vacuum_loc);

//       Kink_profile_guesser ansatz;

//       const int n_dimensions = 3;
//       const double domain_start = 1.e-8;
//       const double domain_end = 9.;
//       const double step_size = 1.e-3;
//       const int n_knots = 100;
//       Field_profiles profiles =
//          ansatz.get_profile_guess(potential, true_vacuum_loc - false_vacuum_loc,
//                                   n_dimensions, domain_start, domain_end,
//                                   step_size, n_knots);
//    }
// }
