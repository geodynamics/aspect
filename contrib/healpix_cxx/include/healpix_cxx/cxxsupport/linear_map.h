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

/*! \file linear_map.h
 *  Simple class for 1D mapping with linear interpolation.
 *  Adapted from RAY++ source code.
 *
 *  Copyright (C) 2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_LINEAR_MAP_H
#define PLANCK_LINEAR_MAP_H

#include <vector>
#include "colour.h"
#include "sort_utils.h"
#include "math_utils.h"

template<typename T> class linearMap
  {
  private:
    bool sorted;
    std::vector<double> x;
    std::vector<T> y;

  public:
    void addVal (double x_, const T &val)
      {
      sorted=x.empty()||(x_>x.back());
      x.push_back(x_);
      y.push_back(val);
      }

    void sortMap()
      {
      if (sorted) return;
      std::vector<size_t> idx;
      buildIndex(x.begin(),x.end(),idx);
      sortByIndex(x.begin(),x.end(),idx);
      sortByIndex(y.begin(),y.end(),idx);
      sorted = true;
      }

    void clear()
      {
      x.clear(); y.clear(); sorted=true;
      }

    T getVal_const (double x_) const
      {
      planck_assert(x.size()>0,"trying to access an empty map");
      planck_assert(sorted,"map must be sorted");
      if (x.size()==1) return y[0];
      if (x_>=x.back())
        return y.back();
      if (x_<=x[0])
        return y[0];
      tsize index;
      double frac;
      interpol_helper (x.begin(), x.end(), x_, index, frac);
      return (1.-frac)*y[index]+frac*y[index+1];
      }
    T getVal (double x_)
      {
      if (!sorted) sortMap();
      return getVal_const(x_);
      }

    size_t size() const { return x.size(); }
    double getX (size_t idx) { if (!sorted) sortMap(); return x[idx]; }
    T getY (size_t idx) { if (!sorted) sortMap(); return y[idx]; }
  };

#endif
