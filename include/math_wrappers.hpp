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

#ifndef BUBBLEPROFILER_MATH_WRAPPERS_HPP_INCLUDED
#define BUBBLEPROFILER_MATH_WRAPPERS_HPP_INCLUDED

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/sign.hpp>

#include <cmath>
#include <complex>
#include <type_traits>

namespace BubbleProfiler {

static const double Pi = boost::math::double_constants::pi;

template <typename T>
T Abs(T x) noexcept
{
   return std::abs(x);
}

template <typename T>
T Abs(const std::complex<T>& z) noexcept
{
   return std::abs(z);
}

template <typename Order, typename Argument>
auto BesselI(Order v, Argument z)
   -> decltype(boost::math::cyl_bessel_i(v, z))
{
   return boost::math::cyl_bessel_i(v, z);
}

template <typename Order, typename Argument>
auto BesselJ(Order v, Argument z)
   -> decltype(boost::math::cyl_bessel_j(v, z))
{
   return boost::math::cyl_bessel_j(v, z);
}

template <typename Order, typename Argument>
auto BesselK(Order v, Argument z)
   -> decltype(boost::math::cyl_bessel_k(v, z))
{
   return boost::math::cyl_bessel_k(v, z);
}

template <typename Order, typename Argument>
auto BesselY(Order v, Argument z)
   -> decltype(boost::math::cyl_neumann(v, z))
{
   return boost::math::cyl_neumann(v, z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Exp(T x) noexcept
{
   return std::exp(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Exp(T n) noexcept
{
   return std::exp(n);
}

template <typename T>
std::complex<T> Exp(const std::complex<T>& z) noexcept
{
   return std::exp(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Log(T x) noexcept
{
   return std::log(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Log(T n) noexcept
{
   return std::log(n);
}

template <typename T>
std::complex<T> Log(const std::complex<T>& z) noexcept
{
   return std::log(z);
}

template <typename T>
int Sign(const T& x)
{
   return boost::math::sign(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Cos(T x) noexcept
{
   return std::cos(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Cos(T n) noexcept
{
   return std::cos(n);
}

template <typename T>
std::complex<T> Cos(const std::complex<T>& z) noexcept
{
   return std::cos(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Sin(T x) noexcept
{
   return std::sin(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Sin(T n) noexcept
{
   return std::sin(n);
}

template <typename T>
std::complex<T> Sin(const std::complex<T>& z) noexcept
{
   return std::sin(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Tan(T x) noexcept
{
   return std::tan(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Tan(T n) noexcept
{
   return std::tan(n);
}

template <typename T>
std::complex<T> Tan(const std::complex<T>& z) noexcept
{
   return std::tan(z);
}

template <typename T>
auto Cot(T x) noexcept -> decltype(Tan(x))
{
   return 1. / Tan(x);
}

template <typename T>
auto Csc(T x) noexcept -> decltype(Sin(x))
{
   return 1. / Sin(x);
}

template <typename T>
auto Sec(T x) noexcept -> decltype(Cos(x))
{
   return 1. / Cos(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T ArcCos(T x) noexcept
{
   return std::acos(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double ArcCos(T n) noexcept
{
   return std::acos(n);
}

template <typename T>
std::complex<T> ArcCos(const std::complex<T>& z) noexcept
{
   return std::acos(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T ArcSin(T x) noexcept
{
   return std::asin(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double ArcSin(T n) noexcept
{
   return std::asin(n);
}

template <typename T>
std::complex<T> ArcSin(const std::complex<T>& z) noexcept
{
   return std::asin(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T ArcTan(T x) noexcept
{
   return std::atan(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double ArcTan(T n) noexcept
{
   return std::atan(n);
}

template <typename T>
std::complex<T> ArcTan(const std::complex<T>& z) noexcept
{
   return std::atan(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Cosh(T x) noexcept
{
   return std::cosh(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Cosh(T n) noexcept
{
   return std::cosh(n);
}

template <typename T>
std::complex<T> Cosh(const std::complex<T>& z) noexcept
{
   return std::cosh(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Sinh(T x) noexcept
{
   return std::sinh(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Sinh(T n) noexcept
{
   return std::sinh(n);
}

template <typename T>
std::complex<T> Sinh(const std::complex<T>& z) noexcept
{
   return std::sinh(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Tanh(T x) noexcept
{
   return std::tanh(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Tanh(T n) noexcept
{
   return std::tanh(n);
}

template <typename T>
std::complex<T> Tanh(const std::complex<T>& z) noexcept
{
   return std::tanh(z);
}

template <typename T>
auto Coth(T x) noexcept -> decltype(Tanh(x))
{
   return 1. / Tanh(x);
}

template <typename T>
auto Csch(T x) noexcept -> decltype(Sinh(x))
{
   return 1. / Sinh(x);
}

template <typename T>
auto Sech(T x) noexcept -> decltype(Cosh(x))
{
   return 1. / Cosh(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T ArcCosh(T x) noexcept
{
   return std::acosh(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double ArcCosh(T n) noexcept
{
   return std::acosh(n);
}

template <typename T>
std::complex<T> ArcCosh(const std::complex<T>& z) noexcept
{
   return std::acosh(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T ArcSinh(T x) noexcept
{
   return std::asinh(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double ArcSinh(T n) noexcept
{
   return std::asinh(n);
}

template <typename T>
std::complex<T> ArcSinh(const std::complex<T>& z) noexcept
{
   return std::asinh(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T ArcTanh(T x) noexcept
{
   return std::atanh(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double ArcTanh(T n) noexcept
{
   return std::atanh(n);
}

template <typename T>
std::complex<T> ArcTanh(const std::complex<T>& z) noexcept
{
   return std::atanh(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_arithmetic<T>::value,T>::type>
T Im(T /* x */) noexcept
{
   return 0.;
}

template <typename T>
T Im(const std::complex<T>& z) noexcept
{
   return std::imag(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_arithmetic<T>::value,T>::type>
T Re(T x) noexcept
{
   return x;
}

template <typename T>
T Re(const std::complex<T>& z) noexcept
{
   return std::real(z);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_floating_point<T>::value,T>::type>
T Sqrt(T x) noexcept
{
   return std::sqrt(x);
}

template <typename T,
          typename = typename std::enable_if<
             std::is_integral<T>::value,T>::type>
double Sqrt(T n) noexcept
{
   return std::sqrt(n);
}

} // namespace BubbleProfiler

#endif
