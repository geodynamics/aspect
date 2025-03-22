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

/*! \file sort_utils.h
 *  Various helper functions for sorting sequences.
 *
 *  Copyright (C) 2002-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_SORT_UTILS_H
#define PLANCK_SORT_UTILS_H

#include <algorithm>
#include <functional>
#include <vector>

template<typename It, typename Comp> class IdxComp__
  {
  private:
    It begin;
    Comp comp;
  public:
    IdxComp__ (It begin_, Comp comp_): begin(begin_), comp(comp_) {}
    bool operator() (std::size_t a, std::size_t b) const
      { return comp(*(begin+a),*(begin+b)); }
  };

/*! Performs an indirect sort on the supplied iterator range and returns in
    \a idx a \a vector containing the indices of the smallest, second smallest,
    third smallest, etc. element, according to \a comp. */
template<typename It, typename T2, typename Comp>
  inline void buildIndex (It begin, It end, std::vector<T2> &idx, Comp comp)
  {
  using namespace std;
  T2 num=end-begin;
  idx.resize(num);
  for (T2 i=0; i<num; ++i) idx[i] = i;
  sort (idx.begin(),idx.end(),IdxComp__<It,Comp>(begin,comp));
  }

/*! Performs an indirect sort on the supplied iterator range and returns in
    \a idx a \a vector containing the indices of the smallest, second smallest,
    third smallest, etc. element. */
template<typename It, typename T2> inline void buildIndex (It begin, It end,
  std::vector<T2> &idx)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  buildIndex(begin,end,idx,less<T>());
  }

/*! Sorts the supplied iterator range according to the order given by \a idx.
    The operation is done out of place and requires temporary extra storage. */
template<typename It, typename T2> inline void sortByIndex (It begin, It end,
  const std::vector<T2> &idx)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  T2 num=end-begin;
  T *tmp= new T[num];
  for (T2 i=0; i<num; ++i) swap(tmp[i],*(begin+i));
  for (T2 i=0; i<num; ++i) swap(*(begin+i),tmp[idx[i]]);
  delete[] tmp;
  }

/*! Sorts the supplied iterator range according to the order given by \a idx.
    The operation is done in place. */
template<typename It, typename T2> inline void sortByIndex_inplace
  (It begin, It end, const std::vector<T2> &idx)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  T2 num=end-begin;
  vector<bool> done(num,false);
  T2 cnt=0;
  while (cnt<num)
    {
    if (!done[cnt]) // new cycle
      {
      T tmp(*(begin+cnt));
      T2 cnt2 = cnt;
      T2 cnt3 = idx[cnt];
      while (cnt3!=cnt)
        {
        done[cnt2]=true;
        *(begin+cnt2)=*(begin+cnt3);
        cnt2=cnt3;
        cnt3=idx[cnt3];
        }
      *(begin+cnt2) = tmp;
      }
    ++cnt;
    }
  }

template<typename It, typename Comp> inline void indirectSort (It begin, It end,
  Comp comp)
  {
  using namespace std;
  vector<std::size_t> idx;
  buildIndex (begin,end,idx,comp);
  sortByIndex (begin,end,idx);
  }

template<typename It> inline void indirectSort (It begin, It end)
  {
  using namespace std;
  typedef typename iterator_traits<It>::value_type T;
  indirectSort(begin,end,less<T>());
  }

#endif
