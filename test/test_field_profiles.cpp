#include "catch/catch.hpp"

#include "field_profiles.hpp"

#include <Eigen/Core>

using namespace BubbleProfiler;

TEST_CASE("test field profile getters and setters", "[basic_field_profiles]")
{
   SECTION("domain boundaries constructed correctly")
   {
      const int n_fields = 2;
      const double domain_start = 1.1;
      const double domain_end = 9.8;
      const int n_steps = 100;

      Field_profiles profiles_one(n_fields, domain_start, domain_end,
                                  n_steps);

      CHECK(profiles_one.get_domain_start() == domain_start);
      CHECK(profiles_one.get_domain_end() == domain_end);

      Eigen::MatrixXd numerical_profiles(Eigen::MatrixXd::Zero(50,2));
      Field_profiles profiles_two(numerical_profiles, domain_start,
                                  domain_end);

      CHECK(profiles_two.get_domain_start() == domain_start);
      CHECK(profiles_two.get_domain_end() == domain_end);
   }

   SECTION("spacetime dimensions set correctly")
   {
      const int n_dims = 3;

      Field_profiles profiles;
      profiles.set_number_of_dimensions(n_dims);

      CHECK(profiles.get_number_of_dimensions() == n_dims);
   }

   SECTION("number of fields constructed correctly")
   {
      const int n_fields = 2;
      const double domain_start = 1.1;
      const double domain_end = 9.8;
      const int n_steps = 100;

      Field_profiles profiles_one(n_fields, domain_start, domain_end,
                                  n_steps);

      CHECK(profiles_one.get_number_of_fields() == n_fields);

      const int n_fields_numerical = 4;
      Eigen::MatrixXd numerical_profiles(
         Eigen::MatrixXd::Zero(50,n_fields_numerical));
      Field_profiles profiles_two(numerical_profiles, domain_start,
                                  domain_end);

      CHECK(profiles_two.get_number_of_fields() == n_fields_numerical);
   }

   SECTION("number of grid points constructed correctly")
   {
      const int n_fields = 2;
      const double domain_start = 1.1;
      const double domain_end = 9.8;
      const int n_steps = 100;

      Field_profiles profiles_one(n_fields, domain_start, domain_end,
                                  n_steps);

      CHECK(profiles_one.get_number_of_grid_points() == n_steps);

      const int n_fields_numerical = 4;
      const int n_steps_numerical = 200;
      Eigen::MatrixXd numerical_profiles(
         Eigen::MatrixXd::Zero(n_steps_numerical,n_fields_numerical));
      Field_profiles profiles_two(numerical_profiles, domain_start,
                                  domain_end);

      CHECK(profiles_two.get_number_of_grid_points() == n_steps_numerical);
   }
}
