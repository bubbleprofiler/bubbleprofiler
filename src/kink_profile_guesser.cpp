/*
 * This file is part of BubbleProfiler.
 *
 * BubbleProfiler is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BubbleProfiler is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with BubbleProfiler.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file kink_profile_guesser.cpp
 * @brief contains the implementation of the Kink_profile_guesser class
 */

#include "kink_profile_guesser.hpp"
#include "error.hpp"
#include "field_profiles.hpp"
#include "math_wrappers.hpp"
#include "potential.hpp"
#include "rotation.hpp"
#include "univariate_interpolation.hpp"

#include <boost/math/tools/roots.hpp>

namespace BubbleProfiler {

const std::array<double,50> Kink_profile_guesser::alpha_grid_thin_3d = {{
   0.5001, 0.500527, 0.500953, 0.50138, 0.501806, 0.502233, 0.502659,
   0.503086, 0.503512, 0.503939, 0.504365, 0.504792, 0.505218, 0.505645,
   0.506071, 0.506498, 0.506924, 0.507351, 0.507778, 0.508204, 0.508631,
   0.509057, 0.509484, 0.50991, 0.510337, 0.510763, 0.51119, 0.511616,
   0.512043, 0.512469, 0.512896, 0.513322, 0.513749, 0.514176, 0.514602,
   0.515029, 0.515455, 0.515882, 0.516308, 0.516735, 0.517161, 0.517588,
   0.518014, 0.518441, 0.518867, 0.519294, 0.51972, 0.520147, 0.520573,
   0.521}};

const std::array<double,46> Kink_profile_guesser::alpha_grid_3d = {{
   0.52, 0.525, 0.53, 0.535, 0.54, 0.545, 0.55, 0.555, 0.56, 0.565,
   0.57, 0.575, 0.58, 0.585, 0.59, 0.595, 0.600, 0.605, 0.61, 0.615,
   0.62, 0.625, 0.63, 0.635, 0.64, 0.645, 0.65, 0.655, 0.66, 0.665,
   0.67, 0.675, 0.68, 0.685, 0.69, 0.695, 0.7, 0.705, 0.71, 0.715,
   0.72, 0.725, 0.73, 0.735, 0.74, 0.745
}};

const std::array<double,46> Kink_profile_guesser::phi0_grid_3d = {{
   0.999937813492861, 0.9998455972594864, 0.999652910051955,
   0.9992783058948789, 0.9986121924679547, 0.9975330072377783,
   0.9959228817337971, 0.9936814248400437, 0.9907300277928794,
   0.9870109949684415, 0.982483019871098, 0.9771146067103315,
   0.970878370010354, 0.9637457498588253, 0.9556829606935944,
   0.9466496322134932, 0.9365966418954739, 0.9254664043974916,
   0.9131939076638974, 0.8997081810874542, 0.8849327027971878,
   0.8687916528072536, 0.8512041830047752, 0.8320951904597592,
   0.8113942473212112, 0.7890322303569218, 0.7649566854809807,
   0.739107381734247, 0.7114770571965077, 0.6820156405313316,
   0.6507569028925989, 0.6176739734302936, 0.582774829968536,
   0.5461906624510391, 0.5080001141498214, 0.46818937520654036,
   0.4271774192183204, 0.3846163877046799, 0.3412559463076824,
   0.296820457223674, 0.2527890892404611, 0.20800300649099204,
   0.16344942922240233, 0.1202112946983806, 0.07758879848049229,
   0.037771036089571164
}};

const std::array<double,46> Kink_profile_guesser::delta_grid_3d = {{
   16.922965327830823, 13.571272753610996, 11.331058461502998,
   9.72627170004579, 8.519044051949859, 7.5773396148807715,
   6.8220259834118355, 6.202755088082414, 5.685949302276117,
   5.248324348086186, 4.873204003649857, 4.5483236996146585,
   4.264459482747917, 4.014532790397335, 3.7930339560854938,
   3.5956130450569255, 3.4187907702727114, 3.259765623821998,
   3.116263124736558, 2.9864196062891506, 2.868715598188374,
   2.7618916064251007, 2.6649334010445216, 2.577006806373881,
   2.49744140180639, 2.4257380368287387, 2.36149961287304,
   2.30454310848976, 2.2546361711122795, 2.211949872156462,
   2.176486424050091, 2.1488285065012924, 2.1297593976297686,
   2.1197410864630797, 2.1200099384640447, 2.133318744013215,
   2.159801924117392, 2.208161283258402, 2.278358268171038,
   2.3859178077488314, 2.5221002614889048, 2.735787838200772,
   3.0589449017018593, 3.5337903678436287, 4.483072657771913,
   6.610083260474584
}};

const std::array<double, 50> Kink_profile_guesser::delta_grid_thin_3d = {{
   3333.65, 633.37, 350.034, 241.893, 184.827, 149.564, 125.612,
   108.282, 95.1601, 84.88, 76.6086, 69.8097, 64.122, 59.2939, 55.144,
   51.539, 48.3781, 45.584, 43.0964, 40.8675, 38.8589, 37.0396, 35.3839,
   33.8709, 32.4827, 31.2046, 30.024, 28.9302, 27.9139, 26.9672,
   26.0831, 25.2558, 24.4798, 23.7506, 23.064, 22.4165, 21.8048,
   21.2259, 20.6774, 20.157, 19.6624, 19.1919, 18.7438, 18.3164,
   17.9084, 17.5185, 17.1455, 16.7883, 16.4461, 16.1177
}};

const std::array<double,46> Kink_profile_guesser::lw_grid_3d = {{
   1.978459028942803, 1.9783391533653445, 1.980043620351959,
   1.9833183792496498, 1.9878550345245096, 1.9933922235764878,
   1.9997509461627188, 2.006889929263504, 2.014866291525423,
   2.0237981316444924, 2.033830355729957, 2.0451032141949237,
   2.0577449824900045, 2.071869407260343, 2.0875710013249194,
   2.104938421173479, 2.124057138088258, 2.145017377292984,
   2.1679122623900757, 2.1928592764450814, 2.219982766555708,
   2.2494508959415302, 2.2814351023846937, 2.3161624866512835,
   2.353904767566588, 2.394955773457391, 2.439715745121684,
   2.488559304923917, 2.542169733478082, 2.6010349012110594,
   2.666156356944669, 2.7382744146654523, 2.8184, 2.908544332472897,
   3.0106358023449724, 3.1262251590796586, 3.2611183199814784,
   3.414688887232135, 3.598807799913082, 3.8148923549129865,
   4.1005191039906865, 4.449432554257999, 4.910068487998784,
   5.623565898325162, 6.7095735786822, 9.468675012670683
}};

const std::array<double,50> Kink_profile_guesser::alpha_grid_thin_4d = {{
   0.5001, 0.500527, 0.500953, 0.50138, 0.501806, 0.502233, 0.502659,
   0.503086, 0.503512, 0.503939, 0.504365, 0.504792, 0.505218, 0.505645,
   0.506071, 0.506498, 0.506924, 0.507351, 0.507778, 0.508204, 0.508631,
   0.509057, 0.509484, 0.50991, 0.510337, 0.510763, 0.51119, 0.511616,
   0.512043, 0.512469, 0.512896, 0.513322, 0.513749, 0.514176, 0.514602,
   0.515029, 0.515455, 0.515882, 0.516308, 0.516735, 0.517161, 0.517588,
   0.518014, 0.518441, 0.518867, 0.519294, 0.51972, 0.520147, 0.520573,
   0.521
}};

const std::array<double,44> Kink_profile_guesser::alpha_grid_4d = {{
   0.53, 0.535, 0.54, 0.545, 0.55, 0.555, 0.56, 0.565, 0.57, 0.575,
   0.58, 0.585, 0.59, 0.595, 0.6, 0.605, 0.61, 0.615, 0.62, 0.625,
   0.63, 0.635, 0.64, 0.645, 0.65, 0.655, 0.66, 0.665, 0.67, 0.675,
   0.68, 0.685, 0.69, 0.695, 0.7, 0.705, 0.71, 0.715, 0.72, 0.725,
   0.73, 0.735, 0.74, 0.745
}};

const std::array<double,44> Kink_profile_guesser::phi0_grid_4d = {{
   0.9999141709240298, 0.9998414842279946, 0.9997210500020388,
   0.9995268952112457, 0.9992242029930557, 0.9987705244339796,
   0.9981181124984616, 0.9972175993269528, 0.9960185973241946,
   0.9944735499229372, 0.9925382119835573, 0.9901722726359306,
   0.9873385812085601, 0.9840030619353041, 0.9801330122384383,
   0.9756962485383627, 0.9706583978537124, 0.9649824145474486,
   0.9586249329012729, 0.9515359099943811, 0.9436552518345606,
   0.9349102241623806, 0.9252134935980744, 0.9144601969894771,
   0.9025250825209266, 0.889258109081031, 0.8744833042541158,
   0.8579917056433094, 0.8395440710509201, 0.8188595500222173,
   0.7956105375419259, 0.7694379924848905, 0.7399297435387302,
   0.7066287141400789, 0.6690428076412978, 0.6266415529455245,
   0.5788928016787105, 0.5252831220500885, 0.46536938556791574,
   0.39893460155848304, 0.32606063251956485, 0.24743515221041074,
   0.16466040885708588, 0.08063706717114749
}};

const std::array<double,44> Kink_profile_guesser::delta_grid_4d = {{
   17.014894262950175, 14.610028788855214, 12.800671189297544,
   11.388389548619989, 10.254164050435781, 9.322275602753079,
   8.542272010441318, 7.879224731566165, 7.308229590060747,
   6.811001429309676, 6.373827028154778, 5.986200667646987,
   5.639951332382715, 5.3286044388913085, 5.0469741633258645,
   4.790842036367612, 4.556758794482891, 4.341858143037507,
   3.9604595901218995, 3.790266654774057, 3.6317516751117074,
   3.4837041965355393, 3.3451012152698394, 3.2150778258445367,
   3.0929205912403455, 3.0929205912403455, 2.978040500905868,
   2.8699958171049054, 2.7684347628427344, 2.6731754668612733,
   2.5842339494473108, 2.5017100804638246, 2.4260077761703513,
   2.3578695501593203, 2.298451406267639, 2.2497129307648507,
   2.2146397850561472, 2.198121465069571, 2.208773532335836,
   2.2602870981624075, 2.3801778022550835, 2.6258827676267256,
   3.1429424551592477, 4.498054491401064
}};

const std::array<double, 50> Kink_profile_guesser::delta_grid_thin_4d = {{
   5000.47, 950.055, 525.051, 362.84, 277.241, 224.346, 188.419,
   162.423, 142.74, 127.32, 114.913, 104.714, 96.183, 88.9408, 82.7161,
   77.3085, 72.5671, 68.376, 64.6446, 61.3012, 58.2884, 55.5594,
   53.0759, 50.8063, 48.724, 46.8069, 45.036, 43.3953, 41.8708, 40.4508,
   39.1247, 37.8837, 36.7197, 35.6259, 34.596, 33.6247, 32.7071,
   31.8389, 31.0162, 30.2355, 29.4937, 28.7879, 28.1157, 27.4746,
   26.8626, 26.2777, 25.7182, 25.1825, 24.6691, 24.1766
}};

const std::array<double,44> Kink_profile_guesser::lw_grid_4d = {{
   1.9675754641180123, 1.9669686361875212, 1.9675684754392428,
   1.9692935533810652, 1.9720420676162898, 1.9757017740428529,
   1.9765326531478618, 1.9801619870000806, 1.9853477052970907,
   1.9911895601726401, 1.9976692121812891, 2.004802765595383,
   2.0126443146996595, 2.0212657580229862, 2.030762873784664,
   2.041244603971598, 2.052835375514101, 2.0656593514562167,
   2.0798588663907114, 2.1129636599979453, 2.1321986035241296,
   2.153465859759776, 2.1769749714123137, 2.202966433768833,
   2.2317241480939427, 2.2635750036552706, 2.2635750036552706,
   2.2989230797415345, 2.338229004268684, 2.3821155890530012,
   2.43130198377452, 2.4866304500518366, 2.5493444959232843,
   2.620939177741102, 2.703373914257843, 2.79937476709974,
   2.912552742082108, 3.0481985390720396, 3.2140493335989815,
   3.4214023205717288, 3.6909247111377415, 4.057576779696893,
   4.5943599048717525, 5.492260348908478
}};

void Kink_profile_guesser::set_max_domain_start(double r)
{
   if (r < 0.) {
      throw Domain_error("lower boundary of domain must be non-negative");
   }
   max_domain_start = r;
}

void Kink_profile_guesser::set_min_domain_end(double r)
{
   if (r < 0.) {
      throw Domain_error("upper boundary of domain must be non-negative");
   }
   min_domain_end = r;
}

void Kink_profile_guesser::compute_vacuum_distance(
   const Potential& potential, const Eigen::VectorXd& true_vacuum)
{
   dist_true_vacuum = true_vacuum.norm();

   assert(dist_true_vacuum > 0);

   if (dist_true_vacuum == 0) {
      throw Setup_error("Kink_profile_guesser::get_profile_guess: "
                        "True and false vacua are coincident");
   }

   num_fields = potential.get_number_of_fields();
   Eigen::VectorXd origin(Eigen::VectorXd::Zero(num_fields));
   if (Abs(potential(true_vacuum) - potential(origin)) < 1.e-12) {
      throw Setup_error("Kink_profile_guesser::get_profile_guess: "
                        "True and false vacua are degenerate");
   }
}

void Kink_profile_guesser::fit_ansatz_parameters(
   const Potential& potential, const Eigen::VectorXd& true_vacuum, int d)
{
   calculate_potential_parameters(potential, true_vacuum);

   // Look up / interpolate the three ansatz parameters
   if (d == 3) {
      phi0 = quadratic_lagrange_interpolation_at_point(
              alpha, alpha_grid_3d, phi0_grid_3d);
   } else {
      phi0 = quadratic_lagrange_interpolation_at_point(
              alpha, alpha_grid_4d, phi0_grid_4d);
   }

   if (d == 3) {
      if (alpha < 0.521) {
         delta = quadratic_lagrange_interpolation_at_point(
                 alpha, alpha_grid_thin_3d, delta_grid_thin_3d);
      }
      else {
         delta = quadratic_lagrange_interpolation_at_point(
                 alpha, alpha_grid_3d, delta_grid_3d);
      }
   } else {
      if (alpha < 0.521) {
         delta = quadratic_lagrange_interpolation_at_point(
                 alpha, alpha_grid_thin_4d, delta_grid_thin_4d);
      }
      else {
         delta = quadratic_lagrange_interpolation_at_point(
                 alpha, alpha_grid_4d, delta_grid_4d);
      }
   }
   if (d == 3) {
      lw = quadratic_lagrange_interpolation_at_point(
              alpha, alpha_grid_3d, lw_grid_3d);
   } else {
      lw = quadratic_lagrange_interpolation_at_point(
              alpha, alpha_grid_4d, lw_grid_4d);
   }

   logger.log_message(logging::Log_level::Trace,
                      "phi0: " + std::to_string(phi0));
   logger.log_message(logging::Log_level::Trace,
                      "Delta: " + std::to_string(delta));
   logger.log_message(logging::Log_level::Trace, "Lw: " + std::to_string(lw));
   logger.log_message(logging::Log_level::Trace, "global_min: "
                      + std::to_string(dist_true_vacuum));
}

/**
 * The initial guess is based on a kink ansatz corresponding to the fitted
 * solution of an approximate one-dimensional problem.  This is obtained
 * by first rotating the given field basis into one in which the vector
 * from the false to the true vacuum lies along the first axis in field
 * space.  The function then fits a quartic potential to the given potential
 * along the line between the two minima.  The returned guess is obtained
 * from the fitted solution to this one-dimensional problem, rotated back
 * into the initial field basis.
 *
 * If the start of the integration domain is to be guessed, the guessed value
 * is determined by finding the radial distance at which the derivative of
 * the one-dimensional fitted profile is equal to 1e-5.
 *
 * If the end of the integration domain is to be guessed, the guessed value
 * is determined by finding the radial distance at which the value of the
 * one-dimensional fitted profile is less than 1e-5.
 *
 * @note the false vacuum is assumed to be located at the origin
 *
 * @sa Profile_guesser::get_profile_guess
 *
 */
Field_profiles Kink_profile_guesser::get_profile_guess(
   const Potential& potential, const Eigen::VectorXd& true_vacuum, int d,
   double domain_start, double domain_end, double initial_step_size,
   double interpolation_points_fraction)
{
   compute_vacuum_distance(potential, true_vacuum);
   fit_ansatz_parameters(potential, true_vacuum, d);

   return calculate_field_profiles(d, domain_start, domain_end,
                                   initial_step_size,
                                   interpolation_points_fraction);
}

void Kink_profile_guesser::calculate_potential_parameters(
   const Potential& potential, const Eigen::VectorXd& true_vacuum)
{
   // Make a copy of the potential for use in ansatz construction
   auto ansatz_potential = std::unique_ptr<Potential>(potential.clone());

   // Get change of basis matrix that orients first component
   // in direction of true vacuum
   cob_matrix = calculate_rotation_to_target(true_vacuum).transpose();

   std::stringstream log_str;
   log_str << "COB: " << cob_matrix;
   logger.log_message(logging::Log_level::Trace, log_str.str());

   // Apply the change of basis matrix. Note that we add
   // the distance to the true vacuum as a scale factor so that the
   // true vacuum is located at (1,0,0,...,0).
   ansatz_potential->apply_basis_change(dist_true_vacuum * cob_matrix);

   // Add a constant term to enforce V(false_vacuum) = 0
   ansatz_potential->add_constant_term(-1*ansatz_potential->operator()(
         Eigen::VectorXd::Zero(num_fields)));

   // New coordinates of the true vacuum
   Eigen::VectorXd _true_vacuum = Eigen::VectorXd::Zero(num_fields);
   _true_vacuum(0) = 1;

   // Calculate ansatz parameters
   logger.log_message(logging::Log_level::Trace,
                      "d: " + std::to_string(trial_dist));

   Eigen::VectorXd trial_vec = Eigen::VectorXd::Zero(num_fields);
   trial_vec(0) = trial_dist;

   double vtrial = ansatz_potential->operator()(trial_vec);
   double vmin = ansatz_potential->operator()(_true_vacuum);
   double r = vmin / vtrial;

   aE = 2. * Abs(((2. - trial_dist * trial_dist) * vmin
                  - vtrial / (trial_dist * trial_dist))
                 / ((trial_dist - 1.) * (trial_dist - 1.)));
   alpha = (r * trial_dist * trial_dist * (trial_dist - 1.5) + 0.5)
      / (1. - r * trial_dist * trial_dist * (2. - trial_dist * trial_dist));

   // if ((Abs(alpha) < alpha_threshold)) {
   //    throw Thin_wall_error("Ansatz alpha = " + std::to_string(alpha) +
   //       " < " + std::to_string(alpha_threshold) + " indicates thin walled "
   //       "bubble; profiler will not converge");
   // }

   if ((Abs(alpha - 0.5) <= std::numeric_limits<double>::epsilon())
       || (Abs(alpha - 0.75) <= std::numeric_limits<double>::epsilon())) {
      throw Setup_error(
         "Kink_profile_guesser::calculate_potential_parameters: "
         "cannot construct ansatz for this potential");
   }

   logger.log_message(logging::Log_level::Trace,
                      "Alpha: " + std::to_string(alpha));
   logger.log_message(logging::Log_level::Trace, "|E|: " + std::to_string(aE));
}

Field_profiles Kink_profile_guesser::calculate_field_profiles(
   int d, double domain_start, double domain_end,
   double initial_step_size, double interpolation_points_fraction)
{
   double min_rho = domain_start >= 0. ? domain_start
      : guess_domain_start();
   double max_rho = domain_end >= 0. ? domain_end
      : guess_domain_end();

   if (min_rho > max_rho) {
      std::swap(min_rho, max_rho);
   }

   const int num_grid_points = 1 + static_cast<int>(
      std::ceil((max_rho - min_rho) / initial_step_size));

   const Eigen::VectorXd rho(
      Eigen::VectorXd::LinSpaced(num_grid_points, min_rho, max_rho));

   Eigen::VectorXd ansatz =
      rho.unaryExpr([this](double r) { return this->evaluate_ansatz_at(r); });

   // We still need to undo the change of basis
   Eigen::MatrixXd cob_reverse = cob_matrix.transpose();

   // Note - we're building all the profiles here, then
   // passing them over one at a time. This is inefficient,
   // but (I think) the path of least resistance right now.
   // Could be improved later.
   Eigen::MatrixXd m_profiles(num_grid_points, num_fields);

   Eigen::VectorXd temp_field_vec = Eigen::VectorXd::Zero(num_fields);

   for (int x = 0; x < num_grid_points; ++x) {
      // Ansatz basis field vec @ this radius
      temp_field_vec(0) = ansatz(x);

      // Transform it to the original field basis (and undo the
      // field coordinate scaling)
      Eigen::VectorXd orig_field_vec
         = cob_reverse * (dist_true_vacuum) * temp_field_vec;

      m_profiles.row(x) = orig_field_vec;
   }

   Field_profiles profiles(rho, m_profiles, interpolation_points_fraction);
   profiles.set_number_of_dimensions(d);

   return profiles;
}

double Kink_profile_guesser::evaluate_ansatz_at(double rho) const
{
   const double sE = Sqrt(aE);
   // The ansatz form assumes that coordinates are scaled
   // by this factor.
   const double coord_scaling = dist_true_vacuum / sE;

   // Note we divide by the scale factor to recover the
   // original spatial coordinate scale.
   const double kink_profile = 0.5 * phi0 *
      (1. - Tanh((rho / coord_scaling - delta) / lw));

   const double correction = -0.5 * phi0 * Exp(-rho * sE)
      * Sech(delta / lw) * Sech(delta / lw) / (lw * dist_true_vacuum);

   return kink_profile + correction;
}

double Kink_profile_guesser::evaluate_ansatz_deriv_at(double rho) const
{
   const double sE = Sqrt(aE);
   const double coord_scaling = dist_true_vacuum / sE;

   const double sechr = Sech((rho / coord_scaling - delta) / lw);
   const double sech0 = Sech(delta / lw);

   const double kink_deriv = -0.5 * phi0 * sechr * sechr / (lw * coord_scaling);

   const double correction_deriv = 0.5 * phi0 * sech0 * sech0
      * Exp(-rho * sE) / (lw * coord_scaling);

   return kink_deriv + correction_deriv;
}

double Kink_profile_guesser::guess_domain_start() const
{
   const double tolerance = 1.e-5;

   const double target_value = 1.e-5;


   if (Abs(evaluate_ansatz_deriv_at(max_domain_start)) <= target_value) {
      return max_domain_start;
   }

   const auto target_fn = [this, target_value](double r) {
      return Abs(this->evaluate_ansatz_deriv_at(r)) / target_value - 1.;
   };

   const double lower_value = target_fn(0.);
   const double upper_value = target_fn(max_domain_start);
   if (lower_value * upper_value >= 0.) {
      return max_domain_start;
   }

   double r_guess = max_domain_start;
   try {
      const auto convergence_test
         = [&target_fn, tolerance](double lower, double upper) {
         return (Abs(target_fn(lower)) <= tolerance &&
                 Abs(target_fn(upper)) <= tolerance);
      };

      const auto solution = boost::math::tools::bisect(
         target_fn, 0., max_domain_start, convergence_test);

      r_guess = 0.5 * (solution.first + solution.second);
   } catch (const boost::exception& e) {
      r_guess = max_domain_start;
   }

   logger.log_message(logging::Log_level::Trace, "guessed domain start: "
                      + std::to_string(r_guess));

   return r_guess;
}

double Kink_profile_guesser::guess_domain_end() const
{
   const double threshold = 1.e-5;

   const double sE = Sqrt(aE);
   const double coord_scaling = dist_true_vacuum / sE;

   double result = coord_scaling*(lw * ArcTanh(1 - ((2*threshold)/phi0))
                                  + delta);
   logger.log_message(logging::Log_level::Trace, "guessed domain end: "
                      + std::to_string(result));

   return std::max(min_domain_end, result);
}

} // namespace BubbleProfiler
