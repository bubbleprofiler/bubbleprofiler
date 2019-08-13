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

#ifndef BUBBLEPROFILER_POTENTIAL_HPP_INCLUDED
#define BUBBLEPROFILER_POTENTIAL_HPP_INCLUDED

#include <cstddef>
#include <Eigen/Core>
#include <tuple>
#include <vector>

#include "error.hpp"

namespace BubbleProfiler {

//! Abstract base class for a generic potential.
/*!
 * Implementations must supply methods for evaluation the potential at a
 * point in field space, and the computation of first and second partial
 * derivatives with respect to field space coordinates.
 *
 * Additionally, implementations must provide methods for translating the
 * origin of the potential, applying a basis change, and adding a constant term.
 */
class Potential {
public:
   virtual ~Potential() = default;

   //! Subclasses must implement a clone method
   virtual Potential * clone() const = 0;

   //! Evaluate potential at point
   /*!
    * @param coords Coordinates at which to evaluate
    * @return Value of potential at coordinates
    */
   virtual double operator()(const Eigen::VectorXd& coords) const = 0;

   //! Partial derivative WRT coordinate i at a point
   /*!
    * @param coords Coordinates at which to evaluate
    * @param i Index of coordinate to be differentiated
    * @return Value of specified partial at point
    */
   virtual double partial(const Eigen::VectorXd& coords, int i) const = 0;

   //! Partial derivative WRT coordinates i, j at a a point
   /*!
    * @param coords Coordinates at which to evaluate
    * @param i Index of first coordinate to be differentiated
    * @param j Index of second coordinate to be differentiated
    * @return Value of specified partial at point
    */
   virtual double partial(const Eigen::VectorXd& coords, int i, int j) const = 0;

   virtual std::size_t get_number_of_fields() const = 0;

   //! Shift the location of the origin by a specified vector
   /*!
    * @param translation shift of origin
    */
   virtual void translate_origin(const Eigen::VectorXd &translation) = 0;

   //! Apply a change of basis matrix
   /*!
    * @param cob_matrix
    */
   virtual void apply_basis_change(const Eigen::MatrixXd& cob_matrix) = 0;

   //! Add a constant offset to the potential
   /*!
    * @param offset constant term to add to potential
    */
   virtual void add_constant_term(double offset) = 0;

   //! Utility method for plotting 2d potentials.
   std::vector<std::tuple<double, double, double>> get_2d_potential_grid(
      unsigned int axis_size, double x_min, double x_max, double y_min, double y_max) {
         
      if (get_number_of_fields() != 2) {
         throw Not_implemented_error("Can only get potential grid for 2 field potentials");
      }
      assert(x_min < x_max);
      assert(y_min < y_max);

      std::vector<std::tuple<double, double, double>> grid;

      double x_step = (x_max - x_min) / axis_size;
      double y_step = (y_max - y_min) / axis_size;

      for (unsigned int ix = 0; ix < axis_size; ix++) {
         for (unsigned int iy = 0; iy < axis_size; iy++) {
            double x = x_min + ix*x_step;
            double y = y_min + iy*y_step;
            Eigen::Vector2d eval_coord(x, y);
            double z = (*this)(eval_coord);
            grid.push_back(std::make_tuple(x, y, z));
         }
      }
      
      return grid;
   }
};

} // namespace BubbleProfiler

#endif
