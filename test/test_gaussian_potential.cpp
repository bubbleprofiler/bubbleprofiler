#include "catch/catch.hpp"
#include "gaussian_potential.hpp"

#include <vector>
#include <Eigen/Core>
#include <iostream>

using namespace BubbleProfiler;

TEST_CASE("Check values and derivatives for test potential") {
   double tolerance = 1.e-7;
   double gamma = 2.0;
   double lambda = 3;
   int n_fields = 2;

   Eigen::VectorXd points_x1(4);
   Eigen::VectorXd points_x2(4);
   Eigen::VectorXd data_v(4);
   Eigen::VectorXd data_dvx1(4);
   Eigen::VectorXd data_dvx2(4);
   Eigen::VectorXd data_dvx1x2(4);

   points_x1 << 0.5,0.0,1.0,1.0;
   points_x2 << 0.5,1.0,0.5,1.0;
   data_v << -0.1469235944499544, -0.11442421040929, -0.1307939799680684, -0.149078454823773;
   data_dvx1 << 0.02472743188255109, -0.03795436189116938, 0.03405227351112536, -0.04296175545580584;
   data_dvx2 << 0.02472743188255109, 0.07646984851812064, -0.03134471647290885, -0.04296175545580584;
   data_dvx1x2 << -0.09137773390141901, -0.04255899811539305, -0.1255045802764081, -0.1723768391121575;

   Gaussian_potential potential = Gaussian_potential(gamma, lambda, n_fields);

   Eigen::VectorXd coord(2);

   for (int i = 0; i < 4; ++i) {
      coord << points_x1(i), points_x2(i);

      CHECK(std::abs(potential(coord) - data_v(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 0) - data_dvx1(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 1) - data_dvx2(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 0, 1) - data_dvx1x2(i)) < tolerance);
   }
}

TEST_CASE("Check values and derivatives for test potential after translation") {
   double tolerance = 1.e-7;
   double gamma = 2.0;
   double lambda = 3;
   int n_fields = 2;

   Eigen::VectorXd points_x1(4);
   Eigen::VectorXd points_x2(4);
   Eigen::VectorXd data_v(4);
   Eigen::VectorXd data_dvx1(4);
   Eigen::VectorXd data_dvx2(4);
   Eigen::VectorXd data_dvx1x2(4);

   points_x1 << -3.5, -3, -2, -2;
   points_x2 << -3.5, -2, -3, -2;
   data_v << -0.1242801156644124, -0.11442421040929, -0.11442421040929, -0.149078454823773;
   data_dvx1 << -0.06284035097790074, -0.03795436189116938, 0.07646984851812064, -0.04296175545580584;
   data_dvx2 << -0.06284035097790074, 0.07646984851812064, -0.03795436189116938, -0.04296175545580584;
   data_dvx1x2 << -0.03325586815821489, -0.04255899811539305, -0.04255899811539305, -0.1723768391121575;

   Gaussian_potential potential = Gaussian_potential(gamma, lambda, n_fields);
   Eigen::VectorXd translation(2);
   translation << 3.0, 3.0;
   potential.translate_origin(translation);

   Eigen::VectorXd coord(2);

   for (int i = 0; i < 4; ++i) {
      coord << points_x1(i), points_x2(i);

      CHECK(std::abs(potential(coord) - data_v(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 0) - data_dvx1(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 1) - data_dvx2(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 0, 1) - data_dvx1x2(i)) < tolerance);
   }
}

TEST_CASE("Check values and derivatives for test potential after translation and basis change") {
   double tolerance = 1.e-7;
   double gamma = 2.0;
   double lambda = 3;
   int n_fields = 2;

   Eigen::VectorXd points_x1(4);
   Eigen::VectorXd points_x2(4);
   Eigen::VectorXd data_v(4);
   Eigen::VectorXd data_dvx1(4);
   Eigen::VectorXd data_dvx2(4);
   Eigen::VectorXd data_dvx1x2(4);

   points_x1 << 1., 1., 1., 1.;
   points_x2 << -4.0, -2.0, 0.0, 1.0;
   data_v << -0.09804459053733929, -0.1527348898350441, -0.08921783744495128, -0.01561627151011979;
   data_dvx1 << 0.09804459053733929, 0.1527348898350441, 0.08921783744495129, 0.01561627151011979;
   data_dvx2 << 0.01085220487429217, -0.0922509396971322, 0.1109014539420692, 0.03502199735755864;
   data_dvx1x2 << -0.01085220487429217, 0.09225093969713223, -0.1109014539420693, -0.03502199735755864;

   Gaussian_potential potential = Gaussian_potential(gamma, lambda, n_fields);
   Eigen::VectorXd translation(2);
   translation << 3.0, 3.0;
   potential.translate_origin(translation);

   // NB this is the inverse of the actual basis transformation, to
   // fit with e.g. the convention used in kink_profile_guesser
   Eigen::MatrixXd cob_matrix(2,2);
   cob_matrix << 1., -1., 1., 1.;
   cob_matrix = (1./std::sqrt(2)) * cob_matrix;

   potential.apply_basis_change(cob_matrix);

   Eigen::VectorXd coord(2);

   for (int i = 0; i < 4; ++i) {
      coord << points_x1(i), points_x2(i);

      CHECK(std::abs(potential(coord) - data_v(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 0) - data_dvx1(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 1) - data_dvx2(i)) < tolerance);
      CHECK(std::abs(potential.partial(coord, 0, 1) - data_dvx1x2(i)) < tolerance);
   }
}
