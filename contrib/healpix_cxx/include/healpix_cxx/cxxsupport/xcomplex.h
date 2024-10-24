/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file xcomplex.h
 *  convenience additions for std::complex
 *
 *  Copyright (C) 2003-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_XCOMPLEX_H
#define PLANCK_XCOMPLEX_H

#include <complex>

#define xcomplex std::complex
typedef xcomplex<float>  fcomplex;
typedef xcomplex<double> dcomplex;

template<typename T> inline xcomplex<T> times_i(xcomplex<T> inp)
  { return xcomplex<T>(-inp.imag(),inp.real()); }

#endif
