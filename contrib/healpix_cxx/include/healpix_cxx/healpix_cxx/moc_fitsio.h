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

/*! \file moc_fitsio.h
 *  FITS I/O for multi-order coverage information
 *
 *  Copyright (C) 2014-2017 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_MOC_FITSIO_H
#define PLANCK_MOC_FITSIO_H

#include <string>
#include "moc.h"

/*! \defgroup moc_fitsio_group FITS-based I/O of MOC objects */
/*! \{ */

template<typename T> Moc<T> read_Moc_from_fits
  (const std::string &filename, bool peano=false);

template<typename T> void write_Moc_to_fits
  (const std::string &outfile, const Moc<T> &moc, bool peano=false);

/*! \} */

#endif
