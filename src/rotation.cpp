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

/*
Notes
- The method for obtaining the new coordinate system
  seems a bit inelegant. There is probably a cleaner way.
*/

#include "rotation.hpp"
#include "error.hpp"

#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>

#include <Eigen/Dense>

namespace BubbleProfiler {

Eigen::MatrixXd calculate_rotation_to_target(const Eigen::VectorXd& target)
{
   // Target must not be origin
   if (target.isZero()) {
      throw Setup_error("Rotation:calculate_rotation_to_target: "
                        "Invalid target: cannot rotate to origin.");
   }

   const int num_dims = target.rows();

   // Matrix with target as first *row*. We use LU decomp.
   // and find the null basis, which spans the complement
   // of the target vector.
   Eigen::MatrixXd mat_fr_target(
      Eigen::MatrixXd::Zero(num_dims, num_dims));
   mat_fr_target.row(0) = target;

   Eigen::FullPivLU<Eigen::MatrixXd> mat_t_lu(mat_fr_target);

   // Sanity check
   assert(mat_t_lu.dimensionOfKernel() == num_dims - 1);
   Eigen::MatrixXd mat_fr_kernel = mat_t_lu.kernel();

   // Adjoin target vector (as first column) to kernel matrix
   Eigen::MatrixXd mat_basis = Eigen::MatrixXd(num_dims, num_dims);
   mat_basis.col(0) = target;
   for (int i = 1; i != num_dims; ++i) {
      mat_basis.col(i) = mat_fr_kernel.col(i - 1);
   }

   // Now we use a QR decomp to obtain an orthonormal basis
   // for the new coordinate system. This is in fact the
   // change of basis matrix that maps new -> old
   Eigen::MatrixXd m_cob_new_to_old(
      mat_basis.fullPivHouseholderQr().matrixQ());

   // The QR decomp may flip the sign of the first column; to be
   // sure, we replace it with the normalized target vector
   m_cob_new_to_old.col(0) = target.normalized();

   return m_cob_new_to_old;
}

} // namespace BubbleProfiler
