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

/*! \file colour.h
 *  Simple class for representing colours. Adapted from RAY++ source code.
 *
 *  Copyright (C) 2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_COLOUR_H
#define PLANCK_COLOUR_H

#include <iostream>

template<typename T> class RGB_tuple
  {
  public:
    T r, g, b;

    RGB_tuple () {}
    RGB_tuple (T rv, T gv, T bv)
      : r (rv), g (gv), b (bv) {}
    template<typename T2> explicit RGB_tuple (const RGB_tuple<T2> &orig)
      : r(orig.r), g(orig.g), b(orig.b) {}

    const RGB_tuple &operator= (const RGB_tuple &Col2)
      { r=Col2.r; g=Col2.g; b=Col2.b; return *this; }
    const RGB_tuple &operator+= (const RGB_tuple &Col2)
      { r+=Col2.r; g+=Col2.g; b+=Col2.b; return *this; }
    const RGB_tuple &operator*= (T fac)
      { r*=fac; g*=fac; b*=fac; return *this; }
    RGB_tuple operator+ (const RGB_tuple &Col2) const
      { return RGB_tuple (r+Col2.r, g+Col2.g, b+Col2.b); }
    RGB_tuple operator- (const RGB_tuple &Col2) const
      { return RGB_tuple (r-Col2.r, g-Col2.g, b-Col2.b); }
    template<typename T2> RGB_tuple operator* (T2 factor) const
      { return RGB_tuple (r*factor, g*factor, b*factor); }
    template<typename T2> friend inline RGB_tuple operator* (T2 factor,
      const RGB_tuple &Col)
      { return RGB_tuple (Col.r*factor, Col.g*factor, Col.b*factor); }

    void Set (T r2, T g2, T b2)
      { r=r2; g=g2; b=b2; }

    friend std::ostream &operator<< (std::ostream &os, const RGB_tuple &c)
      {
      os << "(" << c.r << ", " << c.g << ", " << c.b << ")";
      return os;
      }
  };

typedef RGB_tuple<float> Colour;

#endif
