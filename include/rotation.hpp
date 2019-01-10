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

#ifndef BUBBLEPROFILER_ROTATION_HPP_INCLUDED
#define BUBBLEPROFILER_ROTATION_HPP_INCLUDED

#include <Eigen/Core>

namespace BubbleProfiler {

//! Calculate the rotation matrix for a rotation to the target vector
/*!
 * Given a target vector, this method computes the rotation matrix
 * to transform the current coordinate system to one in which the
 * first basis vector is aligned with the target vector.
 *
 * Notes:
 * 1. this change of coordinates is not uniquely specified,
 * and may involve a reflection.
 *
 * 2. Right multiplication will take coordinates in the rotated system
 * to coordinates in the unrotated system. Take the transpose for the
 * opposite effect.
 *
 * @param target the target vector to align with
 * @return the rotation matrix
 */
Eigen::MatrixXd calculate_rotation_to_target(const Eigen::VectorXd& target);

} // namespace BubbleProfiler

#endif
