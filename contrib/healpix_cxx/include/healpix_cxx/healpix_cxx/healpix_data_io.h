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

/*
 *  Copyright (C) 2003-2016 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef HEALPIX_DATA_IO_H
#define HEALPIX_DATA_IO_H

#include <string>
#include <vector>

class paramfile;
template<typename T> class arr;

/*! Reads a set of full pixel weights for a map of the given \a nside from the
    FITS file \a weightfile. The code checks that the number of weights in the
    file is compatible with the provided \a nside parameter.
    \returns a vector containing the compressed full weights.
 */
std::vector<double> read_fullweights_from_fits(const std::string &weightfile,
  int nside);

void read_weight_ring (const std::string &dir, int nside, arr<double> &weight);

void get_ring_weights (paramfile &params, int nside, arr<double> &weight);

void read_pixwin (const std::string &file, arr<double> &temp);
void read_pixwin (const std::string &file, arr<double> &temp, arr<double> &pol);

void get_pixwin (paramfile &params, int lmax, arr<double> &pixwin);
void get_pixwin (paramfile &params, int lmax, arr<double> &pixwin,
  arr<double> &pixwin_pol);

#endif
