/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file weight_utils.h
 *
 *  Copyright (C) 2016-2019 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_WEIGHT_UTILS_H
#define PLANCK_WEIGHT_UTILS_H

#include <vector>
#include <memory>
#include "healpix_map.h"

/*! Applies the vector \a wgt containing compressed full weights to \a map. */
template<typename T>void apply_fullweights (Healpix_Map<T> &map,
  const std::vector<double> &wgt);

/*! Computes full weights for a map of the given \a nside at a maximum
    multipole \a lmax. The solution is obtained via CGNE iteration, which stops
    once the norm of the residual falls below \a epsilon times the norm of the
    initial residual, or once \a itmax iterations have been performed.
    \returns a vector containing the compressed full weights.
    \note \a lmax must be even; odd \a l do not contribute to the weights. */
std::vector<double> get_fullweights(int nside, int lmax, double epsilon,
  int itmax, double &epsilon_out);

/*! Computes ring weights for a map of the given \a nside at a maximum
    multipole \a lmax. The solution is obtained via CGNE iteration, which stops
    once the norm of the residual falls below \a epsilon times the norm of the
    initial residual, or once \a itmax iterations have been performed.
    \returns a vector of size \c 2*nside containing the ring weights.
    \note \a lmax must be even; odd \a l do not contribute to the weights. */
std::vector<double> get_ringweights(int nside, int lmax, double epsilon,
  int itmax, double &epsilon_out);

namespace weight_utils_detail {

class FullWeightImpl;

}

class FullWeightComputer
  {
  private:
    std::unique_ptr<weight_utils_detail::FullWeightImpl> impl;
  public:
    FullWeightComputer(int nside, int lmax);
    ~FullWeightComputer();
    void iterate(int niter);
    std::vector<double> current_alm() const;
    std::vector<double> alm2wgt(const std::vector<double> &alm) const;
    double current_epsilon() const;
    int current_iter() const;
  };

#endif
