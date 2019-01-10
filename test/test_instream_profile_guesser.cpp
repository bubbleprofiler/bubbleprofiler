#include "catch/catch.hpp"

#include "instream_profile_guesser.hpp"

#include <cmath>
#include <sstream>

using namespace BubbleProfiler;

TEST_CASE("test reads valid input correctly",
          "[instream_profile_guesser]")
{
   const std::size_t n_lines = 4;
   const std::size_t n_fields = 2;
   const std::string test_input = "\
1.0 1.2 3.4\n     \
2.0\t9.0e-1 0.5\n \
3.0 0.1 -0.1\n    \
4.0   0.001 0.0\n \
";

   SECTION("reads the correct number of lines")
   {
      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream);

      CHECK(guesser.get_number_of_grid_points() == n_lines);
   }

   SECTION("reads the correct number of fields")
   {
      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream);

      CHECK(guesser.get_number_of_fields() == n_fields);
   }

   SECTION("reads the correct grid values")
   {
      const double tolerance = 1.e-12;

      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream);

      REQUIRE(guesser.get_number_of_grid_points() == n_lines);

      const Eigen::VectorXd read_grid(guesser.get_spatial_grid());
      CHECK(std::abs(read_grid[0] - 1.) < tolerance);
      CHECK(std::abs(read_grid[1] - 2.) < tolerance);
      CHECK(std::abs(read_grid[2] - 3.) < tolerance);
      CHECK(std::abs(read_grid[3] - 4.) < tolerance);
   }

   SECTION("reads the correct field values")
   {
      const double tolerance = 1.e-12;

      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream);

      REQUIRE(guesser.get_number_of_grid_points() == n_lines);
      REQUIRE(guesser.get_number_of_fields() == n_fields);

      const Eigen::MatrixXd read_fields(guesser.get_field_values());
      CHECK(std::abs(read_fields(0,0) - 1.2) < tolerance);
      CHECK(std::abs(read_fields(1,0) - 0.9) < tolerance);
      CHECK(std::abs(read_fields(2,0) - 0.1) < tolerance);
      CHECK(std::abs(read_fields(3,0) - 0.001) < tolerance);

      CHECK(std::abs(read_fields(0,1) - 3.4) < tolerance);
      CHECK(std::abs(read_fields(1,1) - 0.5) < tolerance);
      CHECK(std::abs(read_fields(2,1) + 0.1) < tolerance);
      CHECK(std::abs(read_fields(3,1)) < tolerance);
   }
}

TEST_CASE("test handles comments correctly",
          "[instream_profile_guesser]")
{
   SECTION("handles default comment character at start of line correctly")
   {
      const std::string test_input = "\
# r x y z\n \
0.0 1.0 2.0\t-3.0\n \
";
      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream);

      CHECK(guesser.get_number_of_grid_points() == 1);
      CHECK(guesser.get_number_of_fields() == 3);
   }

   SECTION("handles default comment character at end of line correctly")
   {
      const std::string test_input = "\
0.0 1.0 -1.0 3.0e-02 # comment one\n  \
1.1 1.3 0.6   0.3\t#comment two\n     \
";
      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream);

      CHECK(guesser.get_number_of_grid_points() == 2);
      CHECK(guesser.get_number_of_fields() == 3);
   }

   SECTION("handles alternative comment string at start of line correctly")
   {
      const std::string test_input = "\
// r x y\n \
0.0 -0.3 0.5\n \
";
      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream, " \t", "//");

      CHECK(guesser.get_number_of_grid_points() == 1);
      CHECK(guesser.get_number_of_fields() == 2);
   }

   SECTION("handles alternative comment string at end of line correctly")
   {
      const std::string test_input = "\
0.0 -0.3 0.5 // comment one\n \
1.5 0.4 -0.4\n \
2.0 0.1 -0.4 // comment two\n \
";
      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream, " \t", "//");

      CHECK(guesser.get_number_of_grid_points() == 3);
      CHECK(guesser.get_number_of_fields() == 2);
   }
}

TEST_CASE("test handles alternative delimiters correctly",
          "[instream_profile_guesser]")
{
   SECTION("test handles single alternative delimiter")
   {
      const double tolerance = 1.e-12;
      const std::string test_input = "\
# r x y z\n \
0.0,1.3,-0.3,0.2\n \
1.1,1.4,-0.8,3e-2\n \
";
      std::stringstream input_stream(test_input);
      Instream_profile_guesser guesser(input_stream, ",");

      CHECK(guesser.get_number_of_grid_points() == 2);
      CHECK(guesser.get_number_of_fields() == 3);

      const Eigen::VectorXd grid(guesser.get_spatial_grid());
      CHECK(std::abs(grid(0)) < tolerance);
      CHECK(std::abs(grid(1) - 1.1) < tolerance);

      const Eigen::MatrixXd field_values(guesser.get_field_values());
      CHECK(std::abs(field_values(0,0) - 1.3) < tolerance);
      CHECK(std::abs(field_values(1,0) - 1.4) < tolerance);
      CHECK(std::abs(field_values(0,1) + 0.3) < tolerance);
      CHECK(std::abs(field_values(1,1) + 0.8) < tolerance);
      CHECK(std::abs(field_values(0,2) - 0.2) < tolerance);
      CHECK(std::abs(field_values(1,2) - 0.03) < tolerance);
   }
}
