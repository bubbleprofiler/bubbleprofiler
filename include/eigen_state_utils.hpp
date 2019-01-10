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

#ifndef BUBBLEPROFILER_EIGEN_STATE_UTILS_HPP_INCLUDED
#define BUBBLEPROFILER_EIGEN_STATE_UTILS_HPP_INCLUDED

#include <Eigen/Core>

#include <boost/version.hpp>

#if BOOST_VERSION >= 105600
#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#else
#include <boost/numeric/odeint/algebra/range_algebra.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#endif

namespace Eigen {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject
operator+(const Eigen::MatrixBase<Derived>& m,
          const typename Eigen::MatrixBase<Derived>::Scalar& s)
{
   return m + s * Eigen::MatrixBase<Derived>::Ones(m.rows(), m.cols());
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject
operator+(const typename Eigen::MatrixBase<Derived>::Scalar& s,
          const Eigen::MatrixBase<Derived>& m)
{
   return m + s * Eigen::MatrixBase<Derived>::Ones(m.rows(), m.cols());
}

template <typename DerivedA, typename DerivedB>
auto operator/(const Eigen::MatrixBase<DerivedA>& x1,
               const Eigen::MatrixBase<DerivedB>& x2)
   -> decltype(x1.cwiseQuotient(x2))
{
   return x1.cwiseQuotient(x2);
}

template <typename Derived>
auto abs(const Eigen::MatrixBase<Derived>& m)
   -> decltype(m.cwiseAbs())
{
   return m.cwiseAbs();
}

} // namespace Eigen

#if BOOST_VERSION >= 105600

namespace boost {
namespace numeric {
namespace odeint {

template <typename Scalar, int Rows, int Cols, int Opts,
          int MaxRows, int MaxCols>
struct vector_space_norm_inf<
   Eigen::Matrix<Scalar,Rows,Cols,Opts,MaxRows,MaxCols> > {
   using result_type = typename Eigen::Matrix<Scalar,Rows,Cols,Opts,
                                              MaxRows,MaxCols>::RealScalar;

   result_type operator()(const Eigen::Matrix<Scalar,Rows,Cols,
                          Opts,MaxRows,MaxCols>& m) const {
      return m.template lpNorm<Eigen::Infinity>();
   }
};

template< class Derived >
struct algebra_dispatcher_sfinae< Derived ,
                                  typename boost::enable_if<
                                     typename boost::is_base_of<
                                        Eigen::MatrixBase< Derived > ,
                                        Derived >::type >::type >
{
   typedef vector_space_algebra algebra_type;
};


template < class Derived  >
struct algebra_dispatcher_sfinae< Derived ,
                      typename boost::enable_if<
                         typename boost::is_base_of<
                            Eigen::ArrayBase< Derived > ,
                            Derived >::type >::type >
{
   typedef vector_space_algebra algebra_type;
};

} // namespace odeint
} // namespace numeric
} // namespace boost

namespace BubbleProfiler {

template <class State>
struct State_algebra_dispatcher {
   using algebra_type =
      typename boost::numeric::odeint::algebra_dispatcher<State>::algebra_type;
};

} // namespace BubbleProfiler

#else

namespace BubbleProfiler {

namespace detail {

template <class State, class Enable = void>
struct State_algebra_dispatcher_helper {
   using algebra_type = boost::numeric::odeint::range_algebra;
};

template< class Derived >
struct State_algebra_dispatcher_helper< Derived ,
                                        typename boost::enable_if<
                                           typename boost::is_base_of<
                                              Eigen::MatrixBase< Derived > ,
                                              Derived >::type >::type >
{
   using algebra_type = boost::numeric::odeint::vector_space_algebra;
};


template < class Derived  >
struct State_algebra_dispatcher_helper< Derived ,
                                        typename boost::enable_if<
                                           typename boost::is_base_of<
                                              Eigen::ArrayBase< Derived > ,
                                              Derived >::type >::type >
{
   using algebra_type = boost::numeric::odeint::vector_space_algebra;
};

} // namespace detail

template <class State>
struct State_algebra_dispatcher :
      public detail::State_algebra_dispatcher_helper<State> {};


} // namespace BubbleProfiler

#endif

#endif
