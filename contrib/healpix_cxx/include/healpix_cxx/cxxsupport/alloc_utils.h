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

/*! \file alloc_utils.h
 *  Classes providing raw memory allocation and deallocation support.
 *
 *  Copyright (C) 2011-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef PLANCK_ALLOC_UTILS_H
#define PLANCK_ALLOC_UTILS_H

#include <cstdlib>
#include <cstddef>
#include "datatypes.h"

template <typename T> class normalAlloc__
  {
  public:
    static T *alloc(tsize sz) { return (sz>0) ? new T[sz] : 0; }
    static void dealloc (T *ptr) { delete[] ptr; }
  };

template <typename T, int align> class alignAlloc__
  {
  private:
//#if (__cplusplus>=201103L)
//    enum { max_nat_align = alignof(std::max_align_t) };
//#else
    enum { max_nat_align = sizeof(void *) };
//#endif

  public:
    static T *alloc(tsize sz)
      {
      using namespace std;
      if (sz==0) return 0;
      planck_assert((align&(align-1))==0,"alignment must be power of 2");
      void *res;
      if (align<=max_nat_align)
        {
        res=malloc(sz*sizeof(T));
        planck_assert(res!=0,"error in malloc()");
        }
      else
        {
        tsize overhead=align-1+sizeof(void*);
        void *ptr=malloc(sz*sizeof(T)+overhead);
        planck_assert(ptr!=0,"error in malloc()");
        tsize sptr=reinterpret_cast<tsize>(ptr);
        sptr = (sptr+overhead) & ~(align-1);
        void **ptr2 = reinterpret_cast<void **>(sptr);
        ptr2[-1]=ptr;
        res=ptr2;
        }
      return static_cast<T *>(res);
      }
    static void dealloc(T *ptr)
      {
      using namespace std;
      if (align<=max_nat_align)
        free(ptr);
      else
        {
        if (ptr==0) return;
        void **ptr2 = reinterpret_cast<void **>(ptr);
        free (ptr2[-1]);
        }
      }
  };

#endif
