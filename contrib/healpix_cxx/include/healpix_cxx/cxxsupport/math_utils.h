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

/*! \file math_utils.h
 *  Various convenience mathematical functions.
 *
 *  Copyright (C) 2002-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_MATH_UTILS_H
#define PLANCK_MATH_UTILS_H

#include <cmath>
#include <vector>
#include <algorithm>
#include "datatypes.h"

/*! \defgroup mathutilsgroup Mathematical helper functions */
/*! \{ */

/*! Returns \e true if | \a a-b | <= \a epsilon * | \a b |, else \e false. */
template<typename F> inline bool approx (F a, F b, F epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) <= (epsilon*abs(b));
  }

/*! Returns \e true if | \a a-b | <= \a epsilon, else \e false. */
template<typename F> inline bool abs_approx (F a, F b, F epsilon=1e-5)
  {
  using namespace std;
  return abs(a-b) <= epsilon;
  }

/*! Returns the largest integer which is smaller than (or equal to) \a arg. */
template<typename I, typename F> inline I ifloor (F arg)
  {
  using namespace std;
  return I(floor(arg));
  }

/*! Returns the integer which is nearest to \a arg. */
template<typename I, typename F> inline I nearest (F arg)
  { return ifloor<I>(arg+0.5); }

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
inline double fmodulo (double v1, double v2)
  {
  using namespace std;
  if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
  double tmp=fmod(v1,v2)+v2;
  return (tmp==v2) ? 0. : tmp;
//  return (v1>=0) ? ((v1<v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2);
  }

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
template<typename I> inline I imodulo (I v1, I v2)
  { I v=v1%v2; return (v>=0) ? v : v+v2; }

/*! Returns -1 if \a signvalue is negative, else +1. */
template<typename T> inline T sign (const T& signvalue)
  { return (signvalue>=0) ? 1 : -1; }

/*! Returns \a val*pow(-1,m) */
template<typename T, typename I> inline T xpow (I m, T val)
  { return (m&1) ? -val : val; }

template<typename I, bool g4> struct isqrt_helper__
  {};
template<typename I> struct isqrt_helper__ <I, false>
  {
  static uint32 isqrt (I arg)
    {
    using namespace std;
    return uint32 (sqrt(arg+0.5));
    }
  };
template<typename I> struct isqrt_helper__ <I, true>
  {
  static uint32 isqrt (I arg)
    {
    using namespace std;
    I res = sqrt(double(arg)+0.5);
    if (arg<(int64(1)<<50)) return uint32(res);
    if (res*res>arg)
      --res;
    else if ((res+1)*(res+1)<=arg)
      ++res;
    return uint32(res);
    }
  };

/*! Returns the integer \a n, which fulfills \a n*n<=arg<(n+1)*(n+1). */
template<typename I> inline uint32 isqrt (I arg)
  { return isqrt_helper__<I,(sizeof(I)>4)>::isqrt(arg); }

/*! Returns the largest integer \a n that fulfills \a 2^n<=arg. */
template<typename I> inline int ilog2 (I arg)
  {
#ifdef __GNUC__
  if (arg==0) return 0;
  if (sizeof(I)==sizeof(int))
    return 8*sizeof(int)-1-__builtin_clz(arg);
  if (sizeof(I)==sizeof(long))
    return 8*sizeof(long)-1-__builtin_clzl(arg);
  if (sizeof(I)==sizeof(long long))
    return 8*sizeof(long long)-1-__builtin_clzll(arg);
#endif
  int res=0;
  while (arg > 0xFFFF) { res+=16; arg>>=16; }
  if (arg > 0x00FF) { res|=8; arg>>=8; }
  if (arg > 0x000F) { res|=4; arg>>=4; }
  if (arg > 0x0003) { res|=2; arg>>=2; }
  if (arg > 0x0001) { res|=1; }
  return res;
  }

template<typename I> inline int ilog2_nonnull (I arg)
  {
#ifdef __GNUC__
  if (sizeof(I)<=sizeof(int))
    return 8*sizeof(int)-1-__builtin_clz(arg);
  if (sizeof(I)==sizeof(long))
    return 8*sizeof(long)-1-__builtin_clzl(arg);
  if (sizeof(I)==sizeof(long long))
    return 8*sizeof(long long)-1-__builtin_clzll(arg);
#endif
  return ilog2 (arg);
  }

template<typename I> inline int trailingZeros(I arg)
  {
  if (arg==0) return sizeof(I)<<3;
#ifdef __GNUC__
  if (sizeof(I)<=sizeof(int))
    return __builtin_ctz(arg);
  if (sizeof(I)==sizeof(long))
    return __builtin_ctzl(arg);
  if (sizeof(I)==sizeof(long long))
    return __builtin_ctzll(arg);
#endif
  int res=0;
  while ((arg&0xFFFF)==0) { res+=16; arg>>=16; }
  if ((arg&0x00FF)==0) { res|=8; arg>>=8; }
  if ((arg&0x000F)==0) { res|=4; arg>>=4; }
  if ((arg&0x0003)==0) { res|=2; arg>>=2; }
  if ((arg&0x0001)==0) { res|=1; }
  return res;
  }

/*! Returns \a atan2(y,x) if \a x!=0 or \a y!=0; else returns 0. */
inline double safe_atan2 (double y, double x)
  {
  using namespace std;
  return ((x==0.) && (y==0.)) ? 0.0 : atan2(y,x);
  }

/*! Helper function for linear interpolation (or extrapolation).
    The array must be ordered in ascending order; no two values may be equal. */
template<typename T, typename Iter, typename Comp> inline void interpol_helper
  (const Iter &begin, const Iter &end, const T &val, Comp comp, tsize &idx,
  T &frac)
  {
  using namespace std;
  planck_assert((end-begin)>1,"sequence too small for interpolation");
  idx = lower_bound(begin,end,val,comp)-begin;
  if (idx>0) --idx;
  idx = min(tsize(end-begin-2),idx);
  frac = (val-begin[idx])/(begin[idx+1]-begin[idx]);
  }

/*! Helper function for linear interpolation (or extrapolation).
    The array must be ordered in ascending order; no two values may be equal. */
template<typename T, typename Iter> inline void interpol_helper
  (const Iter &begin, const Iter &end, const T &val, tsize &idx, T &frac)
  { interpol_helper (begin,end,val,std::less<T>(),idx,frac); }

/*! \} */

#if (__cplusplus>=201103L)

template<typename T1, typename T2>
inline bool multiequal (const T1 &a, const T2 &b)
  { return (a==b); }

template<typename T1, typename T2, typename... Args>
inline bool multiequal (const T1 &a, const T2 &b, Args... args)
  { return (a==b) && multiequal (a, args...); }

#else

template<typename T> inline bool multiequal (const T &a, const T &b)
  { return (a==b); }

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c)
  { return (a==b) && (a==c); }

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c,
  const T &d)
  { return (a==b) && (a==c) && (a==d); }

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c,
  const T &d, const T &e)
  { return (a==b) && (a==c) && (a==d) && (a==e); }

template<typename T> inline bool multiequal (const T &a, const T &b, const T &c,
  const T &d, const T &e, const T &f)
  { return (a==b) && (a==c) && (a==d) && (a==e) && (a==f); }

#endif

template<typename T> class kahan_adder
  {
  private:
    T sum, c;
  public:
    kahan_adder(): sum(0), c(0) {}

    void add (const T &val)
      {
      volatile T tc=c; // volatile to disable over-eager optimizers
      volatile T y=val-tc;
      volatile T t=sum+y;
      tc=(t-sum)-y;
      sum=t;
      c=tc;
      }
    T result() const { return sum; }
  };

template<typename T> class tree_adder
  {
  private:
    std::vector<T> state;
    tsize n;

  public:
    tree_adder(): state(1,T(0)), n(0) {}

    void add (const T &val)
      {
      state[0]+=val; ++n;
      if (n>(tsize(1)<<(state.size()-1)))
        state.push_back(T(0));
      int shift=0;
      while (((n>>shift)&1)==0)
        {
        state[shift+1]+=state[shift];
        state[shift]=T(0);
        ++shift;
        }
      }
    T result() const
      {
      T sum(0);
      for (tsize i=0; i<state.size(); ++i)
        sum+=state[i];
      return sum;
      }
  };

template<typename Iter> bool checkNan (Iter begin, Iter end)
  {
  while (begin!=end)
    {
    if (*begin != *begin) return true;
    ++begin;
    }
  return false;
  }

#endif
