#include "mock_potentials.hpp"
#include "error.hpp"
#include "math_wrappers.hpp"

namespace BubbleProfiler {

Scaled_1D_potential::Scaled_1D_potential(double alpha_)
   : alpha(alpha_)
{
}

Scaled_1D_potential::Scaled_1D_potential(double alpha_, double E_)
   : alpha(alpha_)
   , E(E_)
{
}

double Scaled_1D_potential::get_global_minimum_location() const
{
   double global_min_loc;
   if (E == 0.) {
      global_min_loc = 0.;
   } else if (E > 0.) {
      if (alpha < 0.) {
         global_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else {
         // potential is unbounded from below
         throw Error("Scaled_1D_potential::get_global_minimum_location: "
                     "global minimum does not exist");
      }
   } else {
      if (alpha <= 0.) {
         // potential is unbounded from below
         throw Error("Scaled_1D_potential::get_global_minimum_location: "
                     "global minimum does not exist");
      } else if (alpha > 0. && alpha <= 0.25) {
         global_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else if (alpha > 0.25 && alpha <= 0.5) {
         global_min_loc = 0.;
      } else {
         global_min_loc = 1.;
      }
   }
   return global_min_loc - origin;
}

double Scaled_1D_potential::get_local_minimum_location() const
{
   double local_min_loc;
   if (E == 0.) {
      local_min_loc = 0.;
   } else if (E > 0.) {
      if (alpha <= 0.375) {
         local_min_loc = 1.;
      } else if (alpha > 0.375 && alpha <= 0.75) {
         local_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else {
         local_min_loc = 0.;
      }
   } else {
      if (alpha <= 0.25) {
         local_min_loc = 0.;
      } else if (alpha > 0.25 && alpha <= 0.375) {
         local_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else if (alpha > 0.375 && alpha <= 0.5) {
         local_min_loc = 1.;
      } else if (alpha > 0.5 && alpha <= 0.75) {
         local_min_loc = 0.;
      } else {
         local_min_loc = 0.25 * (3. - 4. * alpha) / alpha;
      }
   }
   return local_min_loc - origin;
}

double Scaled_1D_potential::get_local_maximum_location() const
{
   double local_max_loc;
   if (E == 0.) {
      local_max_loc = 0.;
   } else if (E > 0.) {
      if (alpha <= 0.25) {
         local_max_loc = 0.;
      } else if (alpha > 0.25 && alpha <= 0.375) {
         local_max_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else if (alpha > 0.375 && alpha <= 0.5) {
         local_max_loc = 1.;
      } else if (alpha > 0.5 && alpha <= 0.75) {
         local_max_loc = 0.;
      } else {
         local_max_loc = 0.25 * (3. - 4. * alpha) / alpha;
      }
   } else {
      if (alpha <= 0.375) {
         // potential has one or two local maxima, or no barrier,
         // returns the smaller barrier if multiple local maxima
         local_max_loc = 1.;
      } else if (alpha > 0.375 && alpha <= 0.75) {
         local_max_loc = 0.25 * (3. - 4. * alpha) / alpha;
      } else {
         local_max_loc = 0.;
      }
   }
   return local_max_loc - origin;
}

double Scaled_1D_potential::get_fitted_alpha() const
{
   const double phi_barrier =
      (get_local_maximum_location() - get_local_minimum_location())
      / (get_global_minimum_location() - get_local_minimum_location());
   return 0.75 / (1. + phi_barrier);
}

double Scaled_1D_potential::get_fitted_E() const
{
   const double alpha_fit = get_fitted_alpha();
   Eigen::VectorXd local_min_loc(1);
   Eigen::VectorXd global_min_loc(1);

   local_min_loc << get_local_minimum_location();
   global_min_loc << get_global_minimum_location();

   const double V_loc_min = this->operator()(local_min_loc);
   const double V_global_min = this->operator()(global_min_loc);

   return (V_global_min - V_loc_min) / (alpha_fit - 0.5);
}

double Scaled_1D_potential::operator()(const Eigen::VectorXd& coords) const
{
   if (coords.size() != 1) {
      throw Setup_error("Scaled_1D_potential::partial: "
                        "number of coordinates must be one");
   }

   const double phip = scale*coords(0) + origin;
   return 0.5 * (4. * alpha - 3.) * E * phip * phip
      + E * phip * phip * phip
      - alpha * E * phip * phip * phip * phip
      + offset;
}

double Scaled_1D_potential::partial(const Eigen::VectorXd& coords, int i) const
{
   if (coords.size() != 1) {
      throw Setup_error("Scaled_1D_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Scaled_1D_potential::partial: invalid field index "
         + std::to_string(i));
   }

   const double phip = scale*coords(0) + origin;
   return (4. * alpha - 3.) * E * phip + 3. * E * phip * phip
      - 4. * alpha * E * phip * phip * phip;
}

double Scaled_1D_potential::partial(const Eigen::VectorXd& coords,
                                    int i, int j) const
{
   if (coords.size() != 1) {
      throw Setup_error("Scaled_1D_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Scaled_1D_potential::partial: invalid field index "
         + std::to_string(i));
   }

   if (j != 0) {
      throw Out_of_bounds_error(
         j, "Scaled_1D_potential::partial: invalid field index "
         + std::to_string(j));
   }

   const double phip = scale*coords(0) + origin;
   return (4. * alpha - 3.) * E + 6. * E * phip
      - 12. * alpha * E * phip * phip;
}

Generalized_fubini_potential::Generalized_fubini_potential(double u_, double v_,
                                                           double n_)
   : u(u_), v(v_), n(n_)
{
   if (u <= 0.) {
      throw Setup_error("Generalized_fubini_potential: u must be positive");
   }
   if (v <= 0.) {
      throw Setup_error("Generalized_fubini_potential: v must be positive");
   }
   if (n <= 1.) {
      throw Setup_error("Generalized_fubini_potential: n must be greater than 1");
   }
}

double Generalized_fubini_potential::operator()(
   const Eigen::VectorXd& coords) const
{
   if (coords.size() != 1) {
      throw Setup_error("Generalized_fubini_potential::partial: "
                        "number of coordinates must be one");
   }

   const double phip = Abs(scale * coords(0) + origin);

   const auto pow1 = 2. + 1. / n;
   const auto pow2 = 2. + 2. / n;

   return (4. * u * n * n * (n - 1.) * std::pow(phip, pow1)) /
      (2. * n + 1) - 2. * u * v * n * n * std::pow(phip, pow2);
}

double Generalized_fubini_potential::partial(
   const Eigen::VectorXd& coords, int i) const
{
   if (coords.size() != 1) {
      throw Setup_error("Generalized_fubini_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Generalized_fubini_potential::partial: invalid field index "
         + std::to_string(i));
   }

   const double phip = Abs(scale * coords(0) + origin);

   const auto pow1 = 1. + 1. / n;
   const auto pow2 = 1. + 2. / n;

   return 4. * u * n * (n - 1.) * std::pow(phip, pow1)
      - 4. * u * v * n * (n + 1.) * std::pow(phip, pow2);
}

double Generalized_fubini_potential::partial(
   const Eigen::VectorXd& coords, int i, int j) const
{
   if (coords.size() != 1) {
      throw Setup_error("Generalized_fubini_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Generalized_fubini_potential::partial: invalid field index "
         + std::to_string(i));
   }

   if (j != 0) {
      throw Out_of_bounds_error(
         j, "Generalized_fubini_potential::partial: invalid field index "
         + std::to_string(j));
   }

   const double phip = Abs(scale * coords(0) + origin);

   const auto pow1 = 1. / n;
   const auto pow2 = 2. / n;

   return 4. * u * (n * n - 1.) * std::pow(phip, pow1)
      - 4. * u * v * (n * n + 3. * n + 2.) * std::pow(phip, pow2);
}

double Generalized_fubini_potential::get_local_minimum_location() const
{
   return -origin / scale;
}

double Generalized_fubini_potential::get_local_maximum_location() const
{
   const double phip = std::pow(n - 1., n) / std::pow(v * (n + 1.), n);

   return (phip - origin) / scale;
}

double Generalized_fubini_potential::get_bounce_solution_at(double r) const
{
   const double phip = 1. / std::pow(u * r * r + v, n);

   return (phip - origin) / scale;
}

Eigen::VectorXd Generalized_fubini_potential::get_bounce_solution_at(
   const Eigen::VectorXd& rho_values) const
{
   return rho_values.unaryExpr(
      [this](double r) {
         return this->get_bounce_solution_at(r);
      });
}

double Generalized_fubini_potential::get_action() const
{
   return n * Pi * Pi / ((4. * n * n - 1.) * u * std::pow(v, 2. * n - 1.));
}

Solvable_logarithmic_potential::Solvable_logarithmic_potential(double m_,
                                                               double w_)
   : m(m_), w(w_)
{
}

double Solvable_logarithmic_potential::operator()(
   const Eigen::VectorXd& coords) const
{
   if (coords.size() != 1) {
      throw Setup_error("Solvable_logarithmic_potential::partial: "
                        "number of coordinates must be one");
   }

   const double phip = scale * coords(0) + origin;

   if (Abs(phip) < std::numeric_limits<double>::min()) {
      return 0.;
   }

   return 0.5 * m * m * phip * phip * (1. - Log(phip * phip / (w * w)));
}

double Solvable_logarithmic_potential::partial(
   const Eigen::VectorXd& coords, int i) const
{
   if (coords.size() != 1) {
      throw Setup_error("Solvable_logarithmic_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Solvable_logarithmic_potential::partial: invalid field index "
         + std::to_string(i));
   }

   const double phip = scale * coords(0) + origin;

   if (Abs(phip) < std::numeric_limits<double>::min()) {
      return 0.;
   }

   return -m * m * phip * Log(phip * phip / (w * w));
}

double Solvable_logarithmic_potential::partial(
   const Eigen::VectorXd& coords, int i, int j) const
{
   if (coords.size() != 1) {
      throw Setup_error("Solvable_logarithmic_potential::partial: "
                        "number of coordinates must be one");
   }

   if (i != 0) {
      throw Out_of_bounds_error(
         i, "Solvable_logarithmic_potential::partial: invalid field index "
         + std::to_string(i));
   }

   if (j != 0) {
      throw Out_of_bounds_error(
         j, "Solvable_logarithmic_potential::partial: invalid field index "
         + std::to_string(j));
   }

   const double phip = scale * coords(0) + origin;

   if (Abs(phip) < std::numeric_limits<double>::min()) {
      return 0.;
   }

   return -m * m * (2. + Log(phip * phip / (w * w)));
}

double Solvable_logarithmic_potential::get_local_minimum_location() const
{
   return -origin / scale;
}

double Solvable_logarithmic_potential::get_local_maximum_location() const
{
   const double phip = w;

   return (phip - origin) / scale;
}

double Solvable_logarithmic_potential::get_bounce_solution_at(double r) const
{
   const double phip = w * Exp(-0.5 * m * m * r *r + 2.);
   return (phip - origin) / scale;
}

Eigen::VectorXd Solvable_logarithmic_potential::get_bounce_solution_at(
   const Eigen::VectorXd& rho_values) const
{
   return rho_values.unaryExpr(
      [this](double r) {
         return this->get_bounce_solution_at(r);
      });
}

double Solvable_logarithmic_potential::get_action() const
{
   return 0.5 * Pi * Pi * Exp(4.) * w * w / (m * m);
}

} // namespace BubbleProfiler
