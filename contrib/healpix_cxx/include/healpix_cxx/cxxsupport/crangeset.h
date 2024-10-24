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

/*! \file crangeset.h
 *  Class for storing sets of ranges of integer numbers
 *
 *  Copyright (C) 2015-2021 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_CRANGESET_H
#define PLANCK_CRANGESET_H

#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include "datatypes.h"
#include "error_handling.h"
#include "math_utils.h"

/*
compact rangeset (CRS)

The underlying type T must be signed.
BUT: All "on" ranges must lie entirely in N_0.

A nonnegative number indicates the start of an "on" range at this number.
A negative number indicates the start of an "off" range at the absolute value
of this number.

All numbers r are sorted in a way that |r_i| < |r_{i+1}|.

If two consecutive numbers are both nonnegative or negative, it means that the
interval between them contains in fact _two_ ranges: one of length 1 and the
other filling the remainder.

Consequences:

- The first number in CRS must be nonnegative
- numbers r_i and r_{i+1} in the CRS have the property
  |r_{i+1}| > |r_i| + 1 (because of the special treatment of short intervals)

Example:
The range set
[5;7[; [10;15[; [16;22[; [23;24[; [25;30[; [35;36[; [40;45[

would be encoded as
5 -7 10 -15 -22 -24 -30 35 40 -45
*/

/*! Class for storing sets of ranges of integer numbers
    T must be a signed integer type, but all numbers entered into the range set
    must be nonnegative! */
template<typename T> class crangeset
  {
  public:
    typedef std::vector<T> rtype;

  private:
    rtype r;

    struct abscomp
      {
      bool operator()(const T &a, const T &b)
        {
        using namespace std;
        return abs(a)<abs(b);
        }
      };

    class RsIter
      {
      private:
        tsize idx;
        bool extra;
        const crangeset &ref;
        T val;

        bool singleVal() const
          {
          if (idx==ref.r.size()-1)
            return (ref.r[idx]>=0);
          else
            return ((ref.r[idx]^ref.r[idx+1])>=0);
          }

      public:
      RsIter(const crangeset &ref_) : idx(0), extra(false), ref(ref_), val(0)
          {
          using namespace std;
          if (!atEnd())
            val=abs(ref.r[0]);
          }
        bool atEnd() const { return idx>=ref.r.size(); }
        T operator*() const
          {
          return val;
          }
        RsIter &operator++ ()
          {
          using namespace std;
          if (extra)
            {
            ++idx;
            extra=false;
            if (!atEnd()) val=abs(ref.r[idx]);
            }
          else
            {
            if (singleVal())
              {
              extra=true;
              ++val;
              }
            else
              {
              ++idx;
              if (!atEnd()) val=abs(ref.r[idx]);
              }
            }
          return *this;
          }
        void advance_up_to(T goal)
          {
          using namespace std;
          tdiff idx2=ref.iiv(goal-1);
          if (idx2>tdiff(idx))
            {
            idx=idx2;
            extra=false;
            val=abs(ref.r[idx]);
            }
          }
        bool onoff() const
          { return (ref.r[idx]>=0)^extra; }
      };
  public:
    class IvIter
      {
      private:
        RsIter rsi;
        T b,e;

      public:
        IvIter(const crangeset &ref_) : rsi(ref_)
          {
          if (rsi.atEnd()) return;
          b=*rsi;
          ++rsi;
          e=*rsi;
          }
        bool atEnd() const { return rsi.atEnd(); }
        T ivbegin() const { return b; }
        T ivend() const { return e; }
        tdiff ivlen() const { return e-b; }
        IvIter &operator++ ()
          {
          ++rsi;
          if (rsi.atEnd()) return *this;
          b=*rsi;
          ++rsi;
          e=*rsi;
          return *this;
          }
      };
    class RsOutputIter
      {
      private:
        T val;
        bool val_set;
        crangeset &ref;

      public:
        RsOutputIter(crangeset &ref_) : val(T(0)), val_set(false), ref(ref_)
          { ref.clear(); }
        RsOutputIter &operator*() { return *this; }
        RsOutputIter &operator++ () { return *this; }
        RsOutputIter &operator++ (int) { return *this; }
        RsOutputIter &operator= (const T &v)
          {
          if (val_set)
            ref.append(val,v);
          else
            val=v;
          val_set=!val_set;
          return *this;
          }
      };
    /*! Returns the index of the last number in \c r whose absolute value
        is <= \a val
        If \a val is smaller than all absolute values in \c r, returns -1. */
    tdiff iiv (const T &val) const
      {
      return tdiff(std::upper_bound(r.begin(),r.end(),val,abscomp())
                   -r.begin())-1;
      }

    /*! Estimate a good strategy for set operations involving two rangesets. */
    static int strategy (tsize sza, tsize szb)
      {
      const double fct1 = 1.;
      const double fct2 = 1.;
      tsize slo = sza<szb ? sza : szb,
            shi = sza<szb ? szb : sza;
      double cost1 = fct1 * (sza+szb);
      double cost2 = fct2 * slo * std::max(1,ilog2(shi));
      return (cost1<=cost2) ? 1 : (slo==sza) ? 2 : 3;
      }

    static crangeset generalUnion1 (const crangeset &a, const crangeset &b,
      bool flip_a, bool flip_b)
      {
      crangeset res;
      bool state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
      RsIter ia(a), ib(b);
      RsOutputIter io(res);
      bool runa = !ia.atEnd(), runb = !ib.atEnd();
      while(runa||runb)
        {
        T va = runa ? *ia : T(0),
          vb = runb ? *ib : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = !ia.atEnd(); }
        if (adv_b) { state_b=!state_b; ++ib; runb = !ib.atEnd(); }
        if ((state_a||state_b)!=state_res)
          { *io++ = (adv_a ? va : vb); state_res = !state_res; }
        }
      return res;
      }
    static crangeset generalUnion2 (const crangeset &a, const crangeset &b,
      bool flip_a, bool flip_b)
      {
      crangeset res;
      bool state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
      RsIter ia(a), ib(b);
      RsOutputIter io(res);
      bool runa = !ia.atEnd(), runb = !ib.atEnd();
      while(runa||runb)
        {
        T va = runa ? *ia : T(0);
        if (state_a && runb) // changes in b are irrelevant while state_a is true
          {
          if (!runa) break; // we are at the end
          if (*ib<va)
            {
            ib.advance_up_to (va);
            state_b=!(flip_b^ib.onoff());
            }
          }
        T vb = runb ? *ib : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = !ia.atEnd(); }
        if (adv_b) { state_b=!state_b; ++ib; runb = !ib.atEnd(); }
        if ((state_a||state_b)!=state_res)
          { *io++ = (adv_a ? va : vb); state_res = !state_res; }
        }
      return res;
      }
    static crangeset generalUnion (const crangeset &a, const crangeset &b,
      bool flip_a, bool flip_b)
      {
      if (a.r.empty())
        return flip_a ? crangeset() : b;
      if (b.r.empty())
        return flip_b ? crangeset() : a;

      int strat = strategy (a.r.size(), b.r.size());
      return (strat==1) ? generalUnion1(a,b,flip_a,flip_b) :
               ((strat==2) ? generalUnion2(a,b,flip_a,flip_b)
                           : generalUnion2(b,a,flip_b,flip_a));
      }
    static crangeset generalXor (const crangeset &a, const crangeset &b)
      {
      crangeset res;
      bool state_a=false, state_b=false, state_res=state_a||state_b;
      RsIter ia(a), ib(b);
      RsOutputIter io(res);
      bool runa = !ia.atEnd(), runb = !ib.atEnd();
      while(runa||runb)
        {
        T va = runa ? *ia : T(0),
          vb = runb ? *ib : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = !ia.atEnd(); }
        if (adv_b) { state_b=!state_b; ++ib; runb = !ib.atEnd(); }
        if ((state_a^state_b)!=state_res)
          { *io++ = (adv_a ? va : vb); state_res = !state_res; }
        }
      return res;
      }

    static bool generalAllOrNothing1 (const crangeset &a, const crangeset &b,
      bool flip_a, bool flip_b)
      {
      bool state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
      RsIter ia(a), ib(b);
      bool runa = !ia.atEnd(), runb = !ib.atEnd();
      while(runa||runb)
        {
          using namespace std;
        T va = runa ? *ia : T(0),
          vb = runb ? *ib : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = !ia.atEnd(); }
        if (adv_b) { state_b=!state_b; ++ib; runb = !ib.atEnd(); }
        if ((state_a||state_b)!=state_res)
          return false;
        }
      return true;
      }
    static bool generalAllOrNothing2 (const crangeset &a, const crangeset &b,
      bool flip_a, bool flip_b)
      {
      bool state_a=flip_a, state_b=flip_b, state_res=state_a||state_b;
      RsIter ia(a), ib(b);
      bool runa = !ia.atEnd(), runb = !ib.atEnd();
      while(runa||runb)
        {
        T va = runa ? *ia : T(0);
        if (state_a && runb) // changes in b are irrelevant while state_a is true
          {
          if (!runa) break; // we are at the end
          if (*ib<va)
            {
            ib.advance_up_to (va);
            state_b=!(flip_b^ib.onoff());
            }
          }
        T vb = runb ? *ib : T(0);
        bool adv_a = runa && (!runb || (va<=vb)),
             adv_b = runb && (!runa || (vb<=va));
        if (adv_a) { state_a=!state_a; ++ia; runa = !ia.atEnd(); }
        if (adv_b) { state_b=!state_b; ++ib; runb = !ib.atEnd(); }
        if ((state_a||state_b)!=state_res)
          return false;
        }
      return true;
      }

    static bool generalAllOrNothing (const crangeset &a, const crangeset &b,
      bool flip_a, bool flip_b)
      {
      if (a.r.empty())
        return flip_a ? true : b.r.empty();
      if (b.r.empty())
        return flip_b ? true : a.r.empty();
      int strat = strategy (a.r.size(), b.r.size());
      return (strat==1) ? generalAllOrNothing1(a,b,flip_a,flip_b) :
               ((strat==2) ? generalAllOrNothing2(a,b,flip_a,flip_b)
                           : generalAllOrNothing2(b,a,flip_b,flip_a));
      }
  public:
    /*! Removes all rangeset entries. */
    void clear() { r.clear(); }
    bool empty() const { return r.empty(); }
    const rtype &data() const { return r; }
    void checkConsistency() const
      {
      using namespace std;
      if (r.size()==0) return;
      planck_assert(r[0]>=0,"incorrect first element in range set");
      for (tsize i=1; i<r.size(); ++i)
        planck_assert(abs(r[i])>abs(r[i-1]+1),"inconsistent entries");
      }
    void setData (const rtype &inp)
      {
      r=inp;
      checkConsistency();
      }
    /*! Appends \a [v1;v2[ to the rangeset. \a v1 must be larger
        than the minimum of the last range in the rangeset. */
    void append(const T &v1, const T &v2)
      {
      using namespace std;
      if (v2<=v1) return;
      T le = -100;
      if (!empty())
        le=(r.back()<0) ? -r.back() : r.back()+1;

      if (v1>le) // clean append
        {
        if ((v1>le+1)||(r.back()>0))
          {
          r.push_back(v1);
          if (v2-v1>1) r.push_back(-v2);
          }
        else // short off interval
          r.push_back(-v2);
        return;
        }
        T lastbegin=-200;
        if (!r.empty())
          {
          if (r.back()>=0) lastbegin= r.back();
          else if (r[r.size()-2]<0) lastbegin = -r[r.size()-2]+1;
          else lastbegin = r[r.size()-2];
          }
      planck_assert(v1>=lastbegin,"bad append operation");
      // merge with the last interval
      T endval=max(abs(r.back()),v2);
      if (r.back()<0)
        r.back()=-endval;
      else
        if (endval>r.back()+1) r.push_back(-endval);
      }
    /*! Appends \a [v;v+1[ to the rangeset. \a v must be larger
        than the minimum of the last range in the rangeset. */
    void append(const T &v)
      { append(v,v+1); }

    /*! Appends \a other to the rangeset. All values in \a other must be larger
        than the minimum of the last range in the rangeset. */
    void append (const crangeset &other)
      {
      typename crangeset<T>::IvIter iter(other);
      while (!iter.atEnd())
        {
        append(iter.ivbegin(),iter.ivend());
        ++iter;
        }
      }
    /*! Returns the total number of elements in the rangeset. */
    T nval() const
      {
      T res=0;
      typename crangeset<T>::IvIter iter(*this);
      while (!iter.atEnd())
        {res+=iter.ivlen(); ++iter;}
      return res;
      }

    crangeset op_or (const crangeset &other) const
      { return generalUnion (*this,other,false,false); }
    crangeset op_and (const crangeset &other) const
      { return generalUnion (*this,other,true,true); }
    crangeset op_andnot (const crangeset &other) const
      { return generalUnion (*this,other,true,false); }
    crangeset op_xor (const crangeset &other) const
      { return generalXor (*this,other); }

      /*! Returns \a true if the rangeset is identical to \a other, else \a false.
        */
    bool operator== (const crangeset &other) const
      { return r==other.r; }

    /*! Returns \a true if the rangeset contains all values in the range
        \a [a;b[, else \a false. */
    bool contains (T a,T b) const
      {
      tdiff res=iiv(a);
      if (res<0) return false;
      if (r[res]<0)
        {
        if (res==tdiff(r.size())-1) return false; // beyond end
        if (r[res+1]>=0) return false; // long "off" range
        return ((a>-r[res])&&(b<=-r[res+1])); // mixed range
        }
      // r[res]>=0
      if ((res==tdiff(r.size())-1) || (r[res+1]>=0)) // short interval
        return ((a==r[res])&&(b==a+1));
      return b<=-r[res+1];
      }
    /*! Returns \a true if the rangeset contains the value \a v,
        else \a false. */
    bool contains (T v) const
      {
      tdiff res=iiv(v);
      if (res<0) return false;
      if (r[res]<0)
        {
        if (res==tdiff(r.size())-1) return false; // beyond end
        if (r[res+1]>=0) return false; // long "off" range
        return (v>-r[res]); // mixed range
        }
      if ((res<tdiff(r.size())-1) && (r[res+1]<0))
        return true;
      return (r[res]==v);
      }
    /*! Returns \a true if the rangeset contains all values stored in \a other,
        else \a false. */
    bool contains (const crangeset &other) const
      { return generalAllOrNothing(*this,other,false,true); }

    /** Returns true if any of the numbers [a;b[ are contained in the set,
        else false. */
    bool overlaps (T a,T b) const
      {
      using namespace std;
      if (empty()) return false;
      tdiff res=iiv(a);
      if (res<tdiff(r.size())-1)
        if (b>abs(r[res+1])) return true; // at least one switch
      // now we know that [a;b[ lies entirely inside an interval,
      // but it may still contain a short sub-interval
      if (res==-1) return false; // no sub-intervals in that one
      if (res==tdiff(r.size())-1)
        return a==r[res]; // only overlaps if r[res]>0 and a==abs(r[res])
      // we are somewhere in the middle
      if (r[res]>=0)
        return (a==r[res]) || (r[res+1]<0);
      return (r[res+1]<0) && (b>-r[res]+1);
      }

    /** Returns true if there is overlap between the set and "other",
        else false. */
    bool overlaps (const crangeset &other) const
      { return !generalAllOrNothing(*this,other,true,true); }
  };

template<typename T> inline std::ostream &operator<< (std::ostream &os,
  const crangeset<T> &rs)
  {
  os << "{ ";
  typename crangeset<T>::IvIter iter(rs);
  while (!iter.atEnd())
    {
    os << "["<<iter.ivbegin()<<";"<<iter.ivend()<<"[ ";
    ++iter;
    }
  return os << "}";
  }

#endif
