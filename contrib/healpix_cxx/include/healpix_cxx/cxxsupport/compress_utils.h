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

/*! \file compress_utils.h
 *  Support for compression of integer arrays.
 *
 *  Copyright (C) 2013-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_COMPRESS_UTILS_H
#define PLANCK_COMPRESS_UTILS_H

#include <limits>
#include <vector>
#include <iterator>

#include "datatypes.h"
#include "error_handling.h"
#include "math_utils.h"

class obitstream
  {
  private:
    std::vector<uint8> data;
    tsize bitpos;

    template<typename T> void put_internal (const T &val, uint8 bits)
      {
      int bitsleft = 8-(bitpos&7);
      if (bitsleft==8) data.push_back(0);
      if (bits<=bitsleft) // wbits==bits
        {
        data.back() |= ((val&((T(1)<<bits)-1))<<(bitsleft-bits));
        bitpos+=bits;
        }
      else // wbits==bitsleft
        {
        data.back() |= ((val>>(bits-bitsleft))&((1<<bitsleft)-1));
        bitpos+=bitsleft;
        put_internal(val,bits-bitsleft);
        }
      }

  public:
    obitstream() : bitpos(0) {}

    template<typename T> void put (const T &val, uint8 bits)
      {
      if (bits==0) return;
      planck_assert(bits<=(sizeof(T)<<3),"too many bits for this data type");
//      planck_assert(val-(val&((T(1)<<bits)-1))==0,"value has too many bits");
      put_internal(val,bits);
      }

    const std::vector<uint8> &state() const
      { return data; }
  };

class ibitstream
  {
  private:
    const std::vector<uint8> data;
    tsize bitpos;

    template<typename T> void get_internal (T &val, int bits)
      {
      int bitsleft = 8-(bitpos&7);
      if (bits<=bitsleft)
        {
        val |= ((data[bitpos>>3]>>(bitsleft-bits))&((1<<bits)-1));
        bitpos+=bits;
        }
      else
        {
        val |= T((data[bitpos>>3])&((1<<bitsleft)-1))<<(bits-bitsleft);
        bitpos+=bitsleft;
        get_internal(val,bits-bitsleft);
        }
      }
  public:
    ibitstream(const std::vector<uint8> &indata) : data(indata), bitpos(0) {}

    template<typename T> T get (uint8 bits)
      {
      if (bits==0) return T(0);
      planck_assert(bits<=(sizeof(T)<<3),"too many bits for this data type");
      planck_assert((bitpos+bits)<=8*data.size(),"reading past end of stream");
      T res=T(0);
      get_internal(res,bits);
      return res;
      }
  };

template<typename T,bool> struct notNegativeHelper__
  { notNegativeHelper__(const T &) {} };
template<typename T> struct notNegativeHelper__<T,true>
  {
  notNegativeHelper__(const T &val)
    { planck_assert(val>=T(0),"numbers must be nonnegative");}
  };
template<typename T> void assertNotNegative(const T &val)
  { notNegativeHelper__<T,std::numeric_limits<T>::is_signed> dummy(val); }

template<typename Iter> void interpol_encode2 (Iter l, Iter r, obitstream &obs,
  int shift)
  {
  if (r-l<=1) return;
  typedef std::iterator_traits<Iter> traits;
  typedef typename traits::value_type T;

  Iter m=l+(r-l)/2;
  T nval = ((*r-*l)>>shift) - (r-l) + 1;
  if (nval<=1) return;

  uint8 nb = 1+ilog2(nval-1);
  T val = ((*m)>>shift)-(((*l)>>shift)+(m-l));
  T nshort=(T(1)<<nb)-nval;
#if 0
  // optional rotation
  T nrot=nval-(nshort>>1);
  if (val>=nrot)
    val-=nrot;
  else
    val+=(nval-nrot);
#endif
  if (val<nshort)
    obs.put(val,nb-1);
  else
    obs.put(val+nshort,nb);
  interpol_encode2(l,m,obs,shift);
  interpol_encode2(m,r,obs,shift);
  }

template<typename Iter> void interpol_encode (Iter l, Iter r, obitstream &obs)
  {
  typedef std::iterator_traits<Iter> traits;
  typedef typename traits::value_type T;

  if (l==r) // empty range
    { obs.put(0,8); return; }

  assertNotNegative(*l);

  T combo=*l;
  for (Iter i=l+1; i!=r; ++i)
    {
    planck_assert(*i>*(i-1),"numbers not strictly increasing");
    combo|=*i;
    }
  int shift = trailingZeros(combo);
  T maxnum=(*(r-1))>>shift;
  if (T(r-l)>maxnum) maxnum=T(r-l);
  uint8 maxbits=1+ilog2(maxnum);
  obs.put(maxbits,8);
  obs.put(shift,8);
  obs.put(r-l,maxbits);
  obs.put((*l)>>shift,maxbits);
  if (r-l==1) return;
  obs.put((*(r-1))>>shift,maxbits);
  interpol_encode2(l,r-1,obs,shift);
  }

template<typename Iter> void interpol_decode2 (Iter l, Iter r, ibitstream &ibs,
  int shift)
  {
  if (r-l<=1) return;
  typedef std::iterator_traits<Iter> traits;
  typedef typename traits::value_type T;

  Iter m=l+(r-l)/2;
  T nval = ((*r-*l)>>shift) - (r-l) + 1;
  T val=0;

  if (nval>1)
    {
    uint8 nb = 1+ilog2(nval-1);
    T nshort=(T(1)<<nb)-nval;
    val=ibs.get<T>(nb-1);
    if (val>=nshort)
      val=(val<<1)+ ibs.get<T>(1) - nshort;
#if 0
    // optional rotation
    T nrot=nval-(nshort>>1);
    if (val<(nval-nrot))
      val+=nrot;
    else
      val-=(nval-nrot);
#endif
    }
  *m=*l+(((m-l)+val)<<shift);

  interpol_decode2(l,m,ibs,shift);
  interpol_decode2(m,r,ibs,shift);
  }

template<typename T> void interpol_decode (std::vector<T> &v, ibitstream &ibs)
  {
  uint8 maxbits=ibs.get<uint8>(8);
  if (maxbits==0) { v.clear(); return; }
  int shift = ibs.get<int>(8);
  v.resize(ibs.get<tsize>(maxbits));
  v[0]=ibs.get<T>(maxbits)<<shift;
  if (v.size()==1) return;
  v[v.size()-1]=ibs.get<T>(maxbits)<<shift;
  interpol_decode2(v.begin(),v.end()-1,ibs,shift);
  }

#endif
