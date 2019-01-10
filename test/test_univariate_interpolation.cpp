#include "catch/catch.hpp"

#include "univariate_interpolation.hpp"

#include <array>
#include <cmath>

using namespace BubbleProfiler;

TEST_CASE("test univariate linear Lagrange interpolation with std::array",
          "[univariate_interpolation][lagrange_interpolation][linear_interpolation][interpolation]")
{
   constexpr double tol = 1.e-10;

   SECTION("returns function values when evaluated at data points")
   {
      constexpr std::array<double,6> xdata
         = {{-1.2, 0., 3.2, 4.3, 6.1, 6.2}};
      constexpr std::array<double,6> ydata
         = {{0.991458, 0., -0.727878, -0.351859, -0.469842, 0.674944}};

      for (std::size_t i = 0, num_test_points = xdata.size();
           i < num_test_points; ++i) {
         const double interpolated_value
            = linear_lagrange_interpolation_at_point(xdata[i], xdata, ydata);
         CHECK(std::abs(interpolated_value - ydata[i]) < tol);
      }
   }

   SECTION("correctly interpolates inside the data range")
   {
      constexpr std::array<double,5> xdata = {{1., 2., 3., 4., 5.}};
      constexpr std::array<double,5> ydata = {{1., 4., 9., 16., 25.}};

      constexpr std::array<double,4> test_x_values = {{1.2, 2.5, 3.1, 4.9}};
      constexpr std::array<double,4> expected_y_values = {{1.6, 6.5, 9.7, 24.1}};

      for (std::size_t i = 0, num_test_points = test_x_values.size();
           i < num_test_points; ++i) {
         const double obtained_y_value = linear_lagrange_interpolation_at_point(
            test_x_values[i], xdata, ydata);
         CHECK(std::abs(obtained_y_value - expected_y_values[i]) < tol);
      }
   }

   SECTION("performs linear extrapolation for value below data range")
   {
      constexpr std::array<double,3> xdata = {{-2.1, 1.1, 3.2}};
      constexpr std::array<double,3> ydata = {{10.6281, -0.9559, 84.3776}};

      constexpr double test_x_value = -4.2;
      constexpr double expected_y_value = 18.2301;
      const double obtained_y_value = linear_lagrange_interpolation_at_point(
         test_x_value, xdata, ydata);

      CHECK(std::abs(obtained_y_value - expected_y_value) < tol);
   }

   SECTION("performs linear extrapolation for value above data range")
   {
      constexpr std::array<double,3> xdata = {{2.3, 4.7, 5.}};
      constexpr std::array<double,3> ydata = {{0.832909, 1.54756, 1.60944}};

      constexpr double test_x_value = 6.2;
      constexpr double expected_y_value = 1.85696;
      const double obtained_y_value = linear_lagrange_interpolation_at_point(
         test_x_value, xdata, ydata);

      CHECK(std::abs(obtained_y_value - expected_y_value) < tol);
   }
}

TEST_CASE("test univariate quadratic Lagrange interpolation with std::array",
          "[univariate_interpolation][lagrange_interpolation][quadratic_interpolation][interpolation]")
{
   constexpr double tol = 1.e-10;

   SECTION("returns function values when evaluated at data points")
   {
      constexpr std::array<double,6> xdata
         = {{-5.4, -4., -3.5, -2.1, -1.3, -1.1}};
      constexpr std::array<double,6> ydata
         = {{215.784, 96., 67.375, 18.081, 5.577, 3.751}};

      for (std::size_t i = 0, num_test_points = xdata.size();
           i < num_test_points; ++i) {
         const double interpolated_value
            = quadratic_lagrange_interpolation_at_point(xdata[i], xdata, ydata);
         CHECK(std::abs(interpolated_value - ydata[i]) < tol);
      }
   }

   SECTION("correctly interpolates inside the data range")
   {
      constexpr std::array<double,5> xdata = {{1.5, 3.1, 5.8, 6.7, 7.}};
      constexpr std::array<double,5> ydata
         = {{1.225, 24.281, 183.932, 287.693, 329.3}};

      constexpr std::array<double,4> test_x_values = {{1.7, 3.5, 5.85, 6.8}};
      constexpr std::array<double,4> expected_y_values
         = {{1.195, 38.365, 189.0335, 301.172}};

      for (std::size_t i = 0, num_test_points = test_x_values.size();
           i < num_test_points; ++i) {
         const double obtained_y_value
            = quadratic_lagrange_interpolation_at_point(
               test_x_values[i], xdata, ydata);

         CHECK(std::abs(obtained_y_value - expected_y_values[i]) < tol);
      }
   }

   SECTION("performs quadratic extrapolation for value below data range")
   {
      constexpr std::array<double,3> xdata = {{-2.1, 1.1, 3.2}};
      constexpr std::array<double,3> ydata = {{10.6281, -0.9559, 84.3776}};

      constexpr double test_x_value = -4.2;
      constexpr double expected_y_value = 111.1656;
      const double obtained_y_value = quadratic_lagrange_interpolation_at_point(
         test_x_value, xdata, ydata);

      CHECK(std::abs(obtained_y_value - expected_y_value) < tol);
   }

   SECTION("performs quadratic extrapolation for value above data range")
   {
      constexpr std::array<double,3> xdata = {{2.3, 4.7, 5.}};
      constexpr std::array<double,3> ydata = {{0.832909, 1.54756, 1.60944}};

      constexpr double test_x_value = 6.2;
      constexpr double expected_y_value = 1.7959569444444443;
      const double obtained_y_value = quadratic_lagrange_interpolation_at_point(
         test_x_value, xdata, ydata);

      CHECK(std::abs(obtained_y_value - expected_y_value) < tol);
   }
}
