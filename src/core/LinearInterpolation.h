//  Copyright 2016 National Renewable Energy Laboratory
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//

#ifndef LINEARINTERPOLATION_H
#define LINEARINTERPOLATION_H

/** \file LinearInterpolation.h
 *  1-D linear interpolation utilities
 *
 *  This file implements a single public API method `linear_interp` that is used
 *  to perform piecewise linear interpolations.
 */

#include <iostream>
#include <vector>
#include <stdexcept>
#include <utility>

namespace sierra {
namespace nalu {
namespace utils {

/** Flags and actions for out-of-bounds operation
 */
struct OutOfBounds
{
  //! Out of bounds limit types
  enum boundLimits {
    LOWLIM = -2, //!< xtgt < xarray[0]
    UPLIM = -1,  //!< xtgt > xarray[N]
    VALID = 0    //!< xarray[0] <= xtgt <= xarray[N]
  };

  //! Flags indicating action to perform on Out of Bounds situation
  enum OobAction {
    ERROR = 0,  //!< Raise runtime error
    WARN,       //!< Warn and then CLAMP
    CLAMP,      //!< Clamp values to the end points
    EXTRAPOLATE //!< Extrapolate linearly based on end point
  };
};

template <typename T>
using Array1D = std::vector<T>;

/** Useful type definitions for tracking interpolation data
 */
template <typename T>
struct InterpTraits
{
  typedef typename Array1D<T>::size_type size_type;

  /** A pair indicating whether the point to be interpolate is within bounds
   * or out of bounds and a corresponding index. If the point is within bounds
   * then the index returned is such that `xarray[i] <= xtarget <
   * xarray[i+1]`. For out of bounds, the index is either 0 or IMAX-1.
   */
  typedef typename std::pair<OutOfBounds::boundLimits, size_type> index_type;
};

/**
 * Determine whether the given value is within the limits of the interpolation
 * table
 *
 * \param xinp 1-D array of monotonically increasing values
 * \param x The value to check for
 *
 * \result A std::pair containing the OutOfBounds flag and the index (0 or MAX)
 */
template <typename T>
inline typename InterpTraits<T>::index_type
check_bounds(const Array1D<T>& xinp, const T& x)
{
  auto sz = xinp.size();

  if (sz < 2) {
    throw std::runtime_error(
      "Interpolation table contains less than 2 entries.");
  }

  if (x < xinp[0]) {
    return std::make_pair(OutOfBounds::LOWLIM, 0);
  } else if (x > xinp[sz - 1]) {
    return std::make_pair(OutOfBounds::UPLIM, sz - 1);
  } else {
    return std::make_pair(OutOfBounds::VALID, 0);
  }
}

/**
 * Return an index object corresponding to the x-value based on interpolation
 * table.
 *
 * \param xinp 1-D array of monotonically increasing values
 * \param x The value to check for
 *
 * \return The `std::pair` returned contains two values: the bounds indicator and the
 * index of the element in the interpolation table such that `xarray[i] <= x <
 * xarray[i+1]`
 */
template <typename T>
inline typename InterpTraits<T>::index_type
find_index(const Array1D<T>& xinp, const T& x)
{
  auto idx = check_bounds(xinp, x);
  if (idx.first == OutOfBounds::UPLIM || idx.first == OutOfBounds::LOWLIM)
    return idx;

  auto sz = xinp.size();
  for (size_t i = 1; i < sz; i++) {
    if (x <= xinp[i]) {
      idx.second = i - 1;
      break;
    }
  }
  return idx;
}

/**
 * Perform a 1-D linear interpolation
 *
 * \param xinp A 1-d vector of monotonically increasing x-values
 * \param yinp Corresponding 1-d vector of y-values
 * \param xout Target x-value for interpolation
 * \param yout Interpolated value at `xout`
 * \param oob  (Optional) Out-of-bounds handling (default: CLAMP)
 */
template <typename T>
void
linear_interp(
  const Array1D<T>& xinp,
  const Array1D<T>& yinp,
  const T& xout,
  T& yout,
  OutOfBounds::OobAction oob = OutOfBounds::CLAMP)
{
  auto idx = find_index(xinp, xout);

  switch (idx.first) {
  case OutOfBounds::LOWLIM:
  case OutOfBounds::UPLIM: {
    switch (oob) {
    case OutOfBounds::ERROR:
      throw std::runtime_error("Out of bounds error in interpolation");
      break;

    case OutOfBounds::WARN:
      std::cout
        << std::endl
        << "WARNING: Out of bound values encountered during interpolation"
        << std::endl;
    // no break here... allow fallthrough

    case OutOfBounds::CLAMP:
      yout = yinp[idx.second];
      break;

    case OutOfBounds::EXTRAPOLATE: {
      auto ii = idx.second;
      if ((ii + 1) == xinp.size())
        --ii;

      yout = yinp[ii] +
             (yinp[ii + 1] - yinp[ii]) / (xinp[ii + 1] - xinp[ii]) *
               (xout - xinp[ii]);
      break;
    }
    }
    break;
  }
  case OutOfBounds::VALID:
    auto j = idx.second;
    T fac = (xout - xinp[j]) / (xinp[j + 1] - xinp[j]);
    yout = (static_cast<T>(1.0) - fac) * yinp[j] + fac * yinp[j + 1];
    break;
  }
}

} // namespace utils
} // namespace nalu
} // namespace sierra

#endif /* LINEARINTERPOLATION_H */
