#include "catch/catch.hpp"

#include "rotation.hpp"
#include "error.hpp"

#include <Eigen/Dense>

using namespace BubbleProfiler;

TEST_CASE("test two-dimensional rotations", "[rotation]")
{
   const int num_dims = 2;
   const double tolerance = 1.e-10;

   SECTION("basic rotation test")
   {
      const Eigen::VectorXd target{Eigen::VectorXd::Ones(num_dims)};
      const Eigen::MatrixXd rotation_matrix(
         calculate_rotation_to_target(target));

      const double oneOverSqrt2 = 1. / std::sqrt(2.);
      const Eigen::MatrixXd rotation_magnitudes(rotation_matrix.cwiseAbs());

      REQUIRE(fabs(rotation_magnitudes(0,0) - oneOverSqrt2) < tolerance);
      REQUIRE(fabs(rotation_magnitudes(0,1) - oneOverSqrt2) < tolerance);
      REQUIRE(fabs(rotation_magnitudes(1,0) - oneOverSqrt2) < tolerance);
      REQUIRE(fabs(rotation_magnitudes(1,1) - oneOverSqrt2) < tolerance);

      REQUIRE(fabs(rotation_matrix.determinant() - 1.) < tolerance);
   }

   SECTION("basis already aligned gives trivial rotation")
   {
      Eigen::VectorXd target{Eigen::VectorXd::Zero(num_dims)};
      target(0) = 1.;

      Eigen::MatrixXd rotation_matrix(calculate_rotation_to_target(target));

      REQUIRE(rotation_matrix.isIdentity());
   }

   SECTION("throws exception when target is origin")
   {
      Eigen::VectorXd target{Eigen::VectorXd::Zero(num_dims)};

      REQUIRE_THROWS_AS(calculate_rotation_to_target(target), Setup_error);
   }
}
