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

/*! \file moc_query.h
 *  Copyright (C) 2014-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef MOC_QUERY_H
#define MOC_QUERY_H

#include "vec3.h"
#include "moc.h"

enum MocQueryOp { AND,OR,XOR,NOT,NONE };

class MocQueryComponent
  {
  public:
    MocQueryOp op;
    int nops;
    vec3 center;
    double radius;
    MocQueryComponent(MocQueryOp op_)
      : op(op_)
      {
      planck_assert(op_!=NONE,"bad operator");
      switch (op)
        {
        case AND:
        case OR:
        case XOR:
          nops=2;
          break;
        case NOT:
          nops=1;
          break;
        case NONE:
          nops=0;
          break;
        }
      }
    MocQueryComponent(MocQueryOp op_, int nops_)
      {
      op= op_;
      nops=nops_;
      switch (op)
        {
        case AND:
        case OR:
          planck_assert(nops>=2,"bad nops");
          break;
        case XOR:
          planck_assert(nops==2,"bad nops");
          break;
        case NOT:
          planck_assert(nops==1,"bad nops");
          break;
        case NONE:
          planck_fail("bad operator");
          break;
        }
      }
    MocQueryComponent(const vec3 &cnt, double rad)
      : op (NONE), nops(0), center(cnt.Norm()), radius(rad) {}
  };

template<typename I> Moc<I> mocQuery (int order,
  const std::vector<MocQueryComponent> &comp);

template<typename I> Moc<I> mocQueryInclusive (int order, int omax,
  const std::vector<MocQueryComponent> &comp);

std::vector<MocQueryComponent> prepPolygon (const std::vector<vec3> &vertex);

template<typename I> inline Moc<I> mocQuery (int order,
  const std::vector<vec3> &vertex)
  { return mocQuery<I>(order, prepPolygon(vertex)); }

template<typename I> inline Moc<I> mocQueryInclusive (int order, int omax,
  const std::vector<vec3> &vertex)
  { return mocQueryInclusive<I>(order, omax, prepPolygon(vertex)); }

#endif
