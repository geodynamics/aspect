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

/*! \file moc.h
 *  Copyright (C) 2014-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef HEALPIX_MOC_H
#define HEALPIX_MOC_H

#include "healpix_base.h"
#include "compress_utils.h"

template<typename I> class Moc
  {
  public:
    enum{maxorder=T_Healpix_Base<I>::order_max};

  private:
    rangeset<I> rs;

    static Moc fromNewRangeSet(const rangeset<I> &rngset)
      {
      Moc res;
      res.rs = rngset;
      return res;
      }

  public:
    const rangeset<I> &Rs() { return rs; }
    tsize maxOrder() const
      {
      I combo=0;
      for (tsize i=0; i<rs.nranges(); ++i)
        combo|=rs.ivbegin(i)|rs.ivend(i);
      return maxorder-(trailingZeros(combo)>>1);
      }
    Moc degradedToOrder (int order, bool keepPartialCells) const
      {
      int shift=2*(maxorder-order);
      I ofs=(I(1)<<shift)-1;
      I mask = ~ofs;
      I adda = keepPartialCells ? I(0) : ofs,
        addb = keepPartialCells ? ofs : I(0);
      rangeset<I> rs2;
      for (tsize i=0; i<rs.nranges(); ++i)
        {
        I a=(rs.ivbegin(i)+adda)&mask;
        I b=(rs.ivend  (i)+addb)&mask;
        if (b>a) rs2.append(a,b);
        }
      return fromNewRangeSet(rs2);
      }
    void addPixelRange (int order, I p1, I p2)
      {
      int shift=2*(maxorder-order);
      rs.add(p1<<shift,p2<<shift);
      }
    void appendPixelRange (int order, I p1, I p2)
      {
      int shift=2*(maxorder-order);
      rs.append(p1<<shift,p2<<shift);
      }
    void appendPixel (int order, I p)
      { appendPixelRange(order,p,p+1); }
    /*! Returns a new Moc that contains the union of this Moc and \a other. */
    Moc op_or (const Moc &other) const
      { return fromNewRangeSet(rs.op_or(other.rs)); }
    /*! Returns a new Moc that contains the intersection of this Moc and
        \a other. */
    Moc op_and (const Moc &other) const
      { return fromNewRangeSet(rs.op_and(other.rs)); }
    Moc op_xor (const Moc &other) const
      { return fromNewRangeSet(rs.op_xor(other.rs)); }
    /*! Returns a new Moc that contains all parts of this Moc that are not
        contained in \a other. */
    Moc op_andnot (const Moc &other) const
      { return fromNewRangeSet(rs.op_andnot(other.rs)); }
    /*! Returns the complement of this Moc. */
    Moc complement() const
      {
      rangeset<I> full; full.append(I(0),I(12)*(I(1)<<(2*maxorder)));
      return fromNewRangeSet(full.op_andnot(rs));
      }
    /** Returns \c true if \a other is a subset of this Moc, else \c false. */
    bool contains(const Moc &other) const
      { return rs.contains(other.rs); }
    /** Returns \c true if the intersection of this Moc and \a other is not
        empty. */
    bool overlaps(const Moc &other) const
      { return rs.overlaps(other.rs); }

    /** Returns a vector containing all HEALPix pixels (in ascending NUNIQ
        order) covered by this Moc. The result is well-formed in the sense that
        every pixel is given at its lowest possible HEALPix order. */
    std::vector<I> toUniq() const
      {
      std::vector<I> res;
      std::vector<std::vector<I> > buf(maxorder+1);
      for (tsize i=0; i<rs.nranges(); ++i)
        {
        I start=rs.ivbegin(i), end=rs.ivend(i);
        while(start!=end)
          {
          int logstep=std::min<int>(maxorder,trailingZeros(start)>>1);
          logstep=std::min(logstep,ilog2(end-start)>>1);
          buf[maxorder-logstep].push_back(start);
          start+=I(1)<<(2*logstep);
          }
        }
      for (int o=0; o<=maxorder; ++o)
        {
        I ofs=I(1)<<(2*o+2);
        int shift=2*(maxorder-o);
        for (tsize j=0; j<buf[o].size(); ++j)
          res.push_back((buf[o][j]>>shift)+ofs);
        }
      return res;
      }

    static Moc fromUniq (const std::vector<I> &vu)
      {
      rangeset<I> r, rtmp;
      int lastorder=0;
      int shift=2*maxorder;
      for (tsize i=0; i<vu.size(); ++i)
        {
        int order = ilog2(vu[i]>>2)>>1;
        if (order!=lastorder)
          {
          r=r.op_or(rtmp);
          rtmp.clear();
          lastorder=order;
          shift=2*(maxorder-order);
          }
        I pix = vu[i]-(I(1)<<(2*order+2));
        rtmp.append (pix<<shift,(pix+1)<<shift);
        }
      r=r.op_or(rtmp);
      return fromNewRangeSet(r);
      }

    static void uniq_nest2peano(std::vector<I> &vu)
      {
      using namespace std;
      if (vu.empty()) return;

      tsize start=0;
      int order=-1;
      T_Healpix_Base<I> base;
      I offset=0;
      for (tsize j=0; j<vu.size(); ++j)
        {
        int neworder=ilog2(vu[j]>>2)>>1;
        if (neworder>order)
          {
          sort(vu.begin()+start,vu.begin()+j);
          order=neworder;
          start=j;
          base.Set(order,NEST);
          offset=I(1)<<(2*order+2);
          }
        vu[j]=base.nest2peano(vu[j]-offset)+offset;
        }
      sort(vu.begin()+start,vu.end());
      }
    static void uniq_peano2nest(std::vector<I> &vu)
      {
      using namespace std;
      if (vu.empty()) return;

      tsize start=0;
      int order=-1;
      T_Healpix_Base<I> base;
      I offset=0;
      for (tsize j=0; j<vu.size(); ++j)
        {
        int neworder=ilog2(vu[j]>>2)>>1;
        if (neworder>order)
          {
          sort(vu.begin()+start,vu.begin()+j);
          order=neworder;
          start=j;
          base.Set(order,NEST);
          offset=I(1)<<(2*order+2);
          }
        vu[j]=base.peano2nest(vu[j]-offset)+offset;
        }
      sort(vu.begin()+start,vu.end());
      }

    std::vector<uint8> toCompressed() const
      {
      obitstream obs;
      interpol_encode(rs.data().begin(),rs.data().end(),obs);
      return obs.state();
      }
    static Moc fromCompressed(const std::vector<uint8> &data)
      {
      ibitstream ibs(data);
      std::vector<I> v;
      interpol_decode(v,ibs);
      Moc out;
      out.rs.setData(v);
      return out;
      }

    bool operator==(const Moc &other) const
      {
      if (this == &other)
        return true;
      return rs==other.rs;
      }

    tsize nranges() const
      { return rs.nranges(); }
    I nval() const
      { return rs.nval(); }
  };

#endif
