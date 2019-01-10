#include "catch/catch.hpp"

#include "gsl_root_finder.hpp"

#include <Eigen/Core>

using namespace BubbleProfiler;

TEST_CASE("test GSL wrappers for linear systems of one variable",
          "[gsl_root_finder]")
{
   constexpr std::size_t n_dims = 1;

   SECTION("correct solution for simple 1D system with unique solution")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument& x) -> Function_value {
         Function_value result(n_dims);

         result(0) = 2. * x(0) + 3. - 5.;

         return result;
      };

      const Function_argument initial_guess(Eigen::VectorXd::Zero(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      REQUIRE(status == Root_finder_status::SUCCESS);

      const Function_argument solution = root_finder.get_solution();

      Function_argument exact_solution(n_dims);
      exact_solution << 1.;

      CHECK(std::abs(solution(0) - exact_solution(0)) < 1.e-10);
   }

   SECTION("reports failure to find solution in 1D system with no solution")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument&) -> Function_value {
         Function_value result(n_dims);

         result(0) = 1.;

         return result;
      };

      const Function_argument initial_guess(Eigen::VectorXd::Zero(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      CHECK(status == Root_finder_status::FAIL);
   }

   SECTION("reports failure to find solution in 1D system"
           " with infinite solutions")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument&) -> Function_value {
         Function_value result(n_dims);

         result(0) = 0.;

         return result;
      };

      const Function_argument initial_guess(Eigen::VectorXd::Zero(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      CHECK(status == Root_finder_status::FAIL);
   }

   SECTION("successfully finds solution in 1D non-linear system")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument& x) -> Function_value {
         Function_value result(n_dims);

         result(0) = std::sqrt(x(0)) - 2.;

         return result;
      };

      const Function_argument initial_guess(Eigen::VectorXd::Zero(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      REQUIRE(status == Root_finder_status::SUCCESS);

      const Function_argument solution = root_finder.get_solution();

      Function_argument exact_solution(n_dims);
      exact_solution << 4.;

      CHECK(std::abs(solution(0) - exact_solution(0)) < 1.e-5);
   }
}

TEST_CASE("test GSL wrappers for linear systems of two variables",
          "[gsl_root_finder]")
{
   constexpr std::size_t n_dims = 2;

   SECTION("correct solution for simple 2D system with unique solution")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument& x) -> Function_value {
         Function_value result(n_dims);

         result(0) = x(0) + x(1) - 2.;
         result(1) = x(0) + 2. * x(1) - 1.;

         return result;
      };

      const Function_argument initial_guess(Eigen::VectorXd::Zero(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      REQUIRE(status == Root_finder_status::SUCCESS);

      const Function_argument solution = root_finder.get_solution();

      Function_argument exact_solution(n_dims);
      exact_solution << 3., -1.;

      CHECK(std::abs(solution(0) - exact_solution(0)) < 1.e-10);
      CHECK(std::abs(solution(1) - exact_solution(1)) < 1.e-10);
   }

   SECTION("reports failure to find solution in 2D system with no solution")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument& x) -> Function_value {
         Function_value result(n_dims);

         result(0) = x(0) + x(1) - 2.;
         result(1) = 2. * x(0) + 2. * x(1) - 1.;

         return result;
      };

      const Function_argument initial_guess(Eigen::VectorXd::Zero(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      CHECK(status == Root_finder_status::FAIL);
   }

   SECTION("finds approximate solution to a 2D system with infinite solutions")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument& x) -> Function_value {
         Function_value result(n_dims);

         result(0) = x(0) + x(1) - 2.;
         result(1) = 2. * x(0) + 2. * x(1) - 4.;

         return result;
      };

      const Function_argument initial_guess(Eigen::VectorXd::Zero(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      CHECK(status == Root_finder_status::SUCCESS);
   }

   SECTION("succesfully finds a solution to a 2D non-linear system")
   {
      using Root_finder = GSL_root_finder<n_dims>;
      using Function_argument = Root_finder::Function_argument;
      using Function_value = Root_finder::Function_value;

      Root_finder root_finder;

      const auto f = [n_dims](const Function_argument& x) -> Function_value {
         Function_value result(n_dims);

         result(0) = x(0) * x(0) + x(1) * x(1) - 4.;
         result(1) = 0.25 * (x(0) - 1.) * (x(0) - 1.)
         + x(1) * x(1) / 9. - 1.;

         return result;
      };

      const Function_argument initial_guess(-Eigen::VectorXd::Ones(n_dims));

      const auto status = root_finder.find_root(f, initial_guess);

      REQUIRE(status == Root_finder_status::SUCCESS);

      const Function_argument solution = root_finder.get_solution();

      Function_value solution_value(f(solution));

      CHECK(std::abs(solution_value(0)) < 1.e-4);
      CHECK(std::abs(solution_value(1)) < 1.e-4);
   }
}
