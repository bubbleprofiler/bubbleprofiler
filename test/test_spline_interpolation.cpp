#include "catch/catch.hpp"

#include "gsl_interpolator.hpp"
#include "error.hpp"

#include <Eigen/Core>

using namespace BubbleProfiler;

TEST_CASE("test cubic spline interpolation")
{
   SECTION("throws error when given fewer y values than x values")
   {
      Eigen::VectorXd x(5);
      Eigen::VectorXd y(4);

      x << 1., 2., 3., 4., 5.;
      y << 1., 2., 3., 4.;

      GSL_interpolator interpolator;
      CHECK_THROWS_AS(interpolator.construct_interpolant(x, y), Setup_error);
   }

   SECTION("throws error when given fewer x values than y values")
   {
      Eigen::VectorXd x(4);
      Eigen::VectorXd y(5);

      x << 1., 2., 3., 4.;
      y << 1., 2., 3., 4., 5.;

      GSL_interpolator interpolator;
      CHECK_THROWS_AS(interpolator.construct_interpolant(x, y), Setup_error);
   }

   SECTION("throws error when only one point is given")
   {
      Eigen::VectorXd x(1);
      Eigen::VectorXd y(1);

      x(0) = 1.;
      y(0) = 1.;

      GSL_interpolator interpolator;
      CHECK_THROWS_AS(interpolator.construct_interpolant(x, y), Setup_error);
   }

   SECTION("spline interpolates through given points")
   {
      Eigen::VectorXd x(5);
      Eigen::VectorXd y(5);

      x << 1., 2., 3., 4., 5.;
      y << 1.5, 0.5, 2.3, 3., 1.8;

      GSL_interpolator interpolator;
      interpolator.construct_interpolant(x, y);

      const int n_points = x.size();
      for (int i = 0; i < n_points; ++i) {
         const double value = interpolator.evaluate_at(x(i));
         CHECK(std::abs(value - y(i)) < 1.e-5);
      }
   }
}
