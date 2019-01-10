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
};

} // namespace BubbleProfiler

#endif
