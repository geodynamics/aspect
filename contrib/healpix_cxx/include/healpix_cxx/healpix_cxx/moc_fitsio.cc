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
 *  Copyright (C) 2014-2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "moc_fitsio.h"
#include "fitshandle.h"

using namespace std;

template<typename T> Moc<T> read_Moc_from_fits
  (const std::string &filename, bool peano)
  {
  fitshandle inp;
  inp.open (filename);
  inp.goto_hdu (2);
  vector<T> data;
  inp.read_entire_column(1,data);
  if (peano)
    Moc<T>::uniq_nest2peano(data);
  return Moc<T>::fromUniq(data);
  }

template Moc<int> read_Moc_from_fits
  (const std::string &filename, bool peano);
template Moc<int64> read_Moc_from_fits
  (const std::string &filename, bool peano);

template<typename T> void write_Moc_to_fits
  (const std::string &outfile, const Moc<T> &moc, bool peano)
  {
  PDT outtype=PLANCK_INT16;
  vector<T> data=moc.toUniq();
  if (peano)
    Moc<T>::uniq_peano2nest(data);
  if (data.size()>0)
    {
    if (data.back()>0x7fff) outtype=PLANCK_INT32;
    if (data.back()>0x7fffffffLL) outtype=PLANCK_INT64;
    }
  fitshandle out;
  out.create(outfile);
  vector<fitscolumn> cols;
  cols.push_back (fitscolumn ("PIXEL","",1,outtype));
  out.insert_bintab(cols);
  out.set_key("PIXTYPE", string("HEALPIX"), "HEALPix magic value");
  out.set_key("ORDERING", string("NUNIQ"), "NUNIQ coding method");
  out.set_key("COORDSYS", string("C"), "mandated by MOC standard");
  out.set_key<int>("MOCORDER", moc.maxOrder(), "MOC resolution (best order)");
  out.write_column(1,data);
  }

template void write_Moc_to_fits
  (const std::string &outfile, const Moc<int> &moc, bool peano);
template void write_Moc_to_fits
  (const std::string &outfile, const Moc<int64> &moc, bool peano);
