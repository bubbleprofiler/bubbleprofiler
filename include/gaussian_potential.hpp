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

#ifndef BUBBLEPROFILER_GAUSSIAN_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_GAUSSIAN_POTENTIAL_HPP_INCLUDED

/**
 * @file gaussian_potential.hpp
 * @brief contains definitions of Gaussian_potential class
 */

#include "potential.hpp"
#include "error.hpp"

#include <cmath>
#include <ginac/ginac.h>

namespace BubbleProfiler {

/**
 * @class Gaussian_potential
 * @brief Potential composed of two multivariate Gaussian functions
 *
 * This class implements a potential made of two (multivariate)
 * Gaussian functions. The functional form of the potential is
 * \f[
 *    V(\phi) = -(N(\phi, 0) + \gamma N(\phi, \mu)),
 * \f]
 * where \f$N(\phi,\mu) = \frac{1}{(2\pi)^{n/2}} \exp(-\frac{1}{2}
 * | \phi - \mu|^2)\f$ is a unit Gaussian, \f$\gamma\f$ controls
 * the relative depth of the two minima, and \f$\mu = \frac{1}{\sqrt{n_f}}
 * (\lambda, \ldots, \lambda)\f$, where \f$n_f\f$ is the number of fields
 * and the parameter \f$\lambda\f$ controls the geometric distance between
 * the minima.
 *
 * Note that we enforce \f$\gamma > 1\f$, so that the Gaussian at the origin
 * is a false vacuum.
 */
class Gaussian_potential : public Potential {
public:

   /**
    * @brief instantiate a new Gaussian potential
    * @param gamma_ relative depth of the two minima
    * @param lambda_ geometric distance between minima
    *                (\f$\mu = \frac{1}{\sqrt{n_fields}}(\lambda, ..., \lambda)\f$)
    * @param n_fields_ number of fields
    */
   Gaussian_potential(double gamma_, double lambda_, std::size_t n_fields_);
   virtual ~Gaussian_potential() = default;

   /**
    * @brief return a copy of this potential
    * @return a copy of this object
    */
   virtual Gaussian_potential * clone() const override {
      return new Gaussian_potential(*this);
   };

   virtual double operator()(const Eigen::VectorXd& coords) const override;
   virtual double partial(const Eigen::VectorXd& coords, int i) const override;
   virtual double partial(const Eigen::VectorXd& coords, int i, int j) const override;
   virtual std::size_t get_number_of_fields() const override;

   virtual void translate_origin(const Eigen::VectorXd&) override;
   virtual void apply_basis_change(const Eigen::MatrixXd&) override;
   virtual void add_constant_term(double) override;


private:
   double gamma{2.};
   double lambda{1.};
   std::size_t n_fields{0};
   double constant_term{0.};
   Eigen::VectorXd origin{};
   Eigen::VectorXd origin_translation{};
   Eigen::VectorXd mu{};
   Eigen::MatrixXd basis_transform{};

   double gaussian(const Eigen::VectorXd&, const Eigen::VectorXd&) const;
   double gaussian_deriv(const Eigen::VectorXd&, const Eigen::VectorXd&,
                         int i) const;
   double gaussian_deriv2(const Eigen::VectorXd&, const Eigen::VectorXd&,
                          int i, int j) const;
};

}

#endif
