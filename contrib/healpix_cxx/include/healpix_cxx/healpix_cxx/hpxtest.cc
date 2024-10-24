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
 *  Copyright (C) 2004-2017 Max-Planck-Society
 *  Author: Martin Reinecke
 */

/*

Candidates for testing the validity of the Healpix routines:

- done: ang2pix(pix2ang(i)) = i for all Healpix_Bases
- done: pix2ang(ang2pix(ptg)) dot ptg > 1-delta for all Healpix_Bases
- done: ring2nest(nest2ring(i)) = i for all hierarchical Healpix_Bases
- done: downgrade(upgrade(map)) = map for all maps
- done: map and downgraded map should have same average
- done: alm2map(map2alm(map)) approx map (same for pol)
- partly done: neighbor tests
- powspec -> alm -> powspec (should produce similar powspecs, also for pol)
- done: two swap_schemes() should produce original map
- done: query_disc tests (dot products etc.)
- a_lms: test Set(), Scale(), Add(), alm(l,m) = alm.mstart(m)[l], etc.

*/

#include <iostream>
#include "healpix_base.h"
#include "healpix_map.h"
#include "arr.h"
#include "planck_rng.h"
#include "lsconstants.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "geom_utils.h"
#include "walltimer.h"
#include "announce.h"
#include "compress_utils.h"
#include "moc.h"
#include "moc_fitsio.h"
#include "crangeset.h"
#include "weight_utils.h"
#include "powspec.h"

using namespace std;

#define UNITTESTS

#ifdef UNITTESTS
int errcount=0;
#define FAIL(a) {a; if (++errcount>100) planck_fail("unit test errors");}
#else
#define FAIL(a) a;
#endif

const int nsamples = 1000000;

planck_rng rng;

namespace {

void random_dir (pointing &ptg)
  {
  ptg.theta = acos(rng.rand_uni()*2-1);
  ptg.phi = rng.rand_uni()*twopi;
  }
void random_zphi (double &z, double &phi)
  {
  z = rng.rand_uni()*2-1;
  phi = rng.rand_uni()*twopi;
  }

template<typename I> string bname()
  { return string("(basetype: ")+type2typename<I>()+")"; }


template<typename I> rangeset<I> randomRangeSet(int num, I start, int dist)
  {
  rangeset<I> rs;
  I curval=start;
  for (int i=0; i<num; ++i)
    {
    I v1=curval+1+(rng.int_rand_uni()%dist);
    I v2=v1+1+(rng.int_rand_uni()%dist);
    rs.append(v1,v2);
    curval=v2;
    }
  return rs;
  }
template<typename I> Moc<I> randomMoc(int num, I start, int dist)
  {
  Moc<I> moc;
  I curval=start+(I(1)<<(2*Moc<I>::maxorder));
  for (int i=0; i<num; ++i)
    {
    I v1=curval+1+(rng.int_rand_uni()%dist);
    I v2=v1+1+(rng.int_rand_uni()%dist);
    moc.addPixelRange(Moc<I>::maxorder,v1,v2);
    curval=v2;
    }
  return moc;
  }

template<typename I> rangeset<I> makeRS(const string &input)
  {
  istringstream is(input);
  rangeset<I> res;
  I a,b;
  while (is)
    {
    is>>a>>b;
    planck_assert (is||is.eof(),"aargh");
    if (is) res.append(a,b);
    }
  return res;
  }
template<typename I> void rsOps(const rangeset<I> &a, const rangeset<I> &b)
  {
  planck_assert(!b.overlaps(a.op_andnot(b)),"error");
  planck_assert(a.op_or(b).nval()>=a.nval(),"error");
  planck_assert(a.op_or(b).nval()>=b.nval(),"error");
  planck_assert(a.op_and(b).nval()<=a.nval(),"error");
  planck_assert(a.op_and(b).nval()<=b.nval(),"error");
  planck_assert(a.op_or(b).contains(a),"error");
  planck_assert(a.op_or(b).contains(b),"error");
  planck_assert(!a.op_andnot(b).overlaps(b),"error");
  planck_assert(!b.op_andnot(a).overlaps(a),"error");
  planck_assert(a.op_or(b).nval()==a.nval()+b.nval()-a.op_and(b).nval(),
    "error");
  planck_assert(a.op_or(b).op_andnot(a.op_and(b)).nval()==
                a.op_or(b).nval()-a.op_and(b).nval(),"error");
  planck_assert(a.op_xor(b)==a.op_andnot(b).op_or(b.op_andnot(a)),"error");
  }
template<typename I> void check_rangeset()
  {
  cout << "testing rangeset " << bname<I>() << endl;
  {
  rangeset<I> b;
  b.append(1,11);
  planck_assert(b==makeRS<I>("1 11"),"mismatch");
  b.append(10,15);
  planck_assert(b==makeRS<I>("1 15"),"mismatch");
  b.append(1,15);
  planck_assert(b==makeRS<I>("1 15"),"mismatch");
  b.append(7,15);
  planck_assert(b==makeRS<I>("1 15"),"mismatch");
  b.append(30,41);
#if 0
  planck_assert(b==makeRS<I>("1 15 30 41"),"mismatch");
  try
    {
    b.append(29,31);
    planck_fail("Should have raised an IllegalArgumentException");
    }
  catch (PlanckError E) {}
#endif
  }
  {
  rangeset<I> b=makeRS<I>("1 11 30 41");
  planck_assert(!b.contains(0),"error");
  planck_assert(b.contains(1),"error");
  planck_assert(b.contains(5),"error");
  planck_assert(b.contains(10),"error");
  planck_assert(!b.contains(11),"error");
  planck_assert(!b.contains(29),"error");
  planck_assert(b.contains(30),"error");
  planck_assert(b.contains(35),"error");
  planck_assert(b.contains(40),"error");
  planck_assert(!b.contains(41),"errror");
  }
  {
  rangeset<I> b;
  b.add(5, 11);
  planck_assert(b==makeRS<I>("5 11"),"error");
  b.add(1, 7);
  planck_assert(b==makeRS<I>("1 11"),"error");
  b.add(1, 11);
  planck_assert(b==makeRS<I>("1 11"),"error");
  b.add(30, 41);
  planck_assert(b==makeRS<I>("1 11 30 41"),"error");
  b.add(1, 11);
  planck_assert(b==makeRS<I>("1 11 30 41"),"error");
  b.add(-1,0);
  planck_assert(b==makeRS<I>("-1 0 1 11 30 41"),"error");
  b.add(-2,-1);
  planck_assert(b==makeRS<I>("-2 0 1 11 30 41"),"error");
  b.add(-2,-1);
  planck_assert(b==makeRS<I>("-2 0 1 11 30 41"),"error");
  b.add(2, 11);
  planck_assert(b==makeRS<I>("-2 0 1 11 30 41"),"error");
  b.add(1, 10);
  planck_assert(b==makeRS<I>("-2 0 1 11 30 41"),"error");
  b.add(15, 21);
  planck_assert(b==makeRS<I>("-2 0 1 11 15 21 30 41"),"error");
  }
  {
  rangeset<I> b = makeRS<I>("0 11 20 31");
  b.remove(5,25);
  planck_assert(b==makeRS<I>("0 5 25 31"),"error");
  b.remove(31,32);
  planck_assert(b==makeRS<I>("0 5 25 31"),"error");
  b.remove(35,38);
  planck_assert(b==makeRS<I>("0 5 25 31"),"error");
  b.remove(-90,-80);
  planck_assert(b==makeRS<I>("0 5 25 31"),"error");
  b.remove(27,29);
  planck_assert(b==makeRS<I>("0 5 25 27 29 31"),"error");
  b.remove(25,26);
  planck_assert(b==makeRS<I>("0 5 26 27 29 31"),"error");
  b.remove(4,6);
  planck_assert(b==makeRS<I>("0 4 26 27 29 31"),"error");
  b.remove(-20,40);
  planck_assert(b==makeRS<I>(""),"error");
  b.remove(-20,40);
  planck_assert(b==makeRS<I>(""),"error");
  }
  {
  rangeset<I> b = makeRS<I>("0 11 20 31");
  b.intersect(2,29);
  planck_assert(b==makeRS<I>("2 11 20 29"),"error");
  b.intersect(-8,50);
  planck_assert(b==makeRS<I>("2 11 20 29"),"error");
  b.intersect(2,50);
  planck_assert(b==makeRS<I>("2 11 20 29"),"error");
  b.intersect(2,29);
  planck_assert(b==makeRS<I>("2 11 20 29"),"error");
  b.intersect(-18,29);
  planck_assert(b==makeRS<I>("2 11 20 29"),"error");
  b.intersect(3,11);
  planck_assert(b==makeRS<I>("3 11"),"error");
  b = makeRS<I>("0 11 20 31");
  b.intersect(3,15);
  planck_assert(b==makeRS<I>("3 11"),"error");
  b = makeRS<I>("0 11 20 31");
  b.intersect(17,30);
  planck_assert(b==makeRS<I>("20 30"),"error");
  b = makeRS<I>("0 11 20 31");
  b.intersect(11,20);
  planck_assert(b==makeRS<I>(""),"error");
  b = makeRS<I>("0 11 20 31");
  b.intersect(-8,-7);
  planck_assert(b==makeRS<I>(""),"error");
  b = makeRS<I>("0 11 20 31");
  b.intersect(31,35);
  planck_assert(b==makeRS<I>(""),"error");
  }
  {
  planck_assert(makeRS<I>("1 11 20 31 40 56")==
                makeRS<I>("20 31 40 51").op_or
                (makeRS<I>("1 11 45 56")),"error");
  planck_assert(makeRS<I>("1 11 45 56")==
                makeRS<I>("").op_or
                (makeRS<I>("1 11 45 56")),"error");
  planck_assert(makeRS<I>("1 11 45 56")==
                makeRS<I>("1 11 45 56").op_or
                (makeRS<I>("")),"error");
  }
  {
  planck_assert(makeRS<I>("22 24 45 51")==
                makeRS<I>("20 31 40 51").op_and
                (makeRS<I>("1 11 22 24 45 56")),"error");
  planck_assert(makeRS<I>("20 31 40 51 90 101 110 121 200 201")==
                makeRS<I>("10 101 110 121 200 221").op_and
                (makeRS<I>("20 31 40 51 90 201")),"error");
  planck_assert(makeRS<I>("")==
                makeRS<I>("20 31 40 51").op_and
                (makeRS<I>("")),"error");
  planck_assert(makeRS<I>("")==
                makeRS<I>("").op_and
                (makeRS<I>("20 31 40 51")),"error");
  }
  {
  planck_assert(makeRS<I>("20 31 40 45")==
                makeRS<I>("20 31 40 51").op_andnot
                (makeRS<I>("1 11 45 56")),"error");
  planck_assert(makeRS<I>("")==
                makeRS<I>("").op_andnot
                (makeRS<I>("1 11 45 56")),"error");
  planck_assert(makeRS<I>("1 11 45 56")==
                makeRS<I>("1 11 45 56").op_andnot
                (makeRS<I>("")),"error");
  }
  {
  rangeset<I> b = makeRS<I>("20 31 40 51");

  planck_assert(!b.contains(0,11),"error");
  planck_assert(!b.contains(10,21),"error");
  planck_assert(!b.contains(19,20),"error");
  planck_assert(b.contains(20,21),"error");
  planck_assert(b.contains(21,22),"error");
  planck_assert(b.contains(20,31),"error");
  planck_assert(!b.contains(25,36),"error");
  planck_assert(b.contains(30,31),"error");
  planck_assert(!b.contains(31,32),"error");
  planck_assert(!b.contains(35,38),"error");
  planck_assert(!b.contains(35,46),"error");
  planck_assert(b.contains(40,41),"error");
  planck_assert(!b.contains(45,56),"error");
  planck_assert(!b.contains(60,71),"error");
  }
  {
  rangeset<I> b = makeRS<I>("20 31 40 51");

  planck_assert(b.contains(makeRS<I>("20 31 40 51")),"error");
  planck_assert(b.contains(makeRS<I>("20 21")),"error");
  planck_assert(b.contains(makeRS<I>("50 51")),"error");
  planck_assert(!b.contains(makeRS<I>("19 31 40 51")),"error");
  planck_assert(!b.contains(makeRS<I>("20 31 40 52")),"error");
  planck_assert(!b.contains(makeRS<I>("20 51")),"error");
  planck_assert(!b.contains(makeRS<I>("0 1")),"error");
  planck_assert(!b.contains(makeRS<I>("0 20 31 40 51 100")),"error");
  }
  {
  rangeset<I> b = makeRS<I>("20 31 40 51");

  planck_assert(!b.overlaps(0,11),"error");
  planck_assert(b.overlaps(10,21),"error");
  planck_assert(!b.overlaps(19,20),"error");
  planck_assert(b.overlaps(20,21),"error");
  planck_assert(b.overlaps(21,22),"error");
  planck_assert(b.overlaps(20,31),"error");
  planck_assert(b.overlaps(25,36),"error");
  planck_assert(b.overlaps(30,37),"error");
  planck_assert(!b.overlaps(31,32),"error");
  planck_assert(!b.overlaps(35,38),"error");
  planck_assert(b.overlaps(35,46),"error");
  planck_assert(b.overlaps(40,41),"error");
  planck_assert(b.overlaps(45,56),"error");
  planck_assert(!b.overlaps(60,71),"error");
  }
  {
  rangeset<I> b = makeRS<I>("20 31 40 51");

  planck_assert(b.overlaps(makeRS<I>("20 31 40 51")),"error");
  planck_assert(b.overlaps(makeRS<I>("20 21")),"error");
  planck_assert(b.overlaps(makeRS<I>("50 51")),"error");
  planck_assert(b.overlaps(makeRS<I>("19 31 40 51")),"error");
  planck_assert(b.overlaps(makeRS<I>("20 31 40 52")),"error");
  planck_assert(b.overlaps(makeRS<I>("20 51")),"error");
  planck_assert(!b.overlaps(makeRS<I>("0 1")),"error");
  planck_assert(!b.overlaps(makeRS<I>("0 20 31 40 51 100")),"error");
  }
  {
  const int niter = 1000;
  for (int iter=0; iter<niter; ++iter)
    {
    rangeset<I> a = randomRangeSet<I>(1000, 0, 100);
    rangeset<I> b = randomRangeSet<I>(1000, 0, 100);
    rangeset<I> c = randomRangeSet<I>(10, 10000, 100);
    rsOps(a,b);
    rsOps(a,c);
    rsOps(c,a);
    }
  }
  }
template<typename I> crangeset<I> randomCRangeSet(int num, I start, int dist)
  {
  crangeset<I> rs;
  I curval=start;
  for (int i=0; i<num; ++i)
    {
    I v1=curval+1+(rng.int_rand_uni()%dist);
    I v2=v1+1+(rng.int_rand_uni()%dist);
    rs.append(v1,v2);
    curval=v2;
    }
  return rs;
  }
template<typename I> crangeset<I> makeCRS(const string &input)
  {
  istringstream is(input);
  crangeset<I> res;
  I a,b;
  while (is)
    {
    is>>a>>b;
    planck_assert (is||is.eof(),"aargh");
    if (is) res.append(a,b);
    }
  return res;
  }
template<typename I> void crsOps(const crangeset<I> &a, const crangeset<I> &b)
  {
  planck_assert(!b.overlaps(a.op_andnot(b)),"error");
  planck_assert(a.op_or(b).nval()>=a.nval(),"error");
  planck_assert(a.op_or(b).nval()>=b.nval(),"error");
  planck_assert(a.op_and(b).nval()<=a.nval(),"error");
  planck_assert(a.op_and(b).nval()<=b.nval(),"error");
  planck_assert(a.op_or(b).contains(a),"error");
  planck_assert(a.op_or(b).contains(b),"error");
  planck_assert(!a.op_andnot(b).overlaps(b),"error");
  planck_assert(!b.op_andnot(a).overlaps(a),"error");
  planck_assert(a.op_or(b).nval()==a.nval()+b.nval()-a.op_and(b).nval(),
    "error");
  planck_assert(a.op_or(b).op_andnot(a.op_and(b)).nval()==
    a.op_or(b).nval()-a.op_and(b).nval(),"error");
  planck_assert(a.op_xor(b)==a.op_andnot(b).op_or(b.op_andnot(a)),"error");
  }
template<typename I> void check_crangeset()
  {
  cout << "testing crangeset " << bname<I>() << endl;
  {
  crangeset<I> b;
  b.append(1,11);
  planck_assert(b==makeCRS<I>("1 11"),"mismatch");
  b.append(10,15);
  planck_assert(b==makeCRS<I>("1 15"),"mismatch");
  b.append(1,15);
  planck_assert(b==makeCRS<I>("1 15"),"mismatch");
  b.append(7,15);
  planck_assert(b==makeCRS<I>("1 15"),"mismatch");
  b.append(30,41);
#if 0
  planck_assert(b==makeRS<I>("1 15 30 41"),"mismatch");
  try
    {
    b.append(29,31);
    planck_fail("Should have raised an IllegalArgumentException");
    }
  catch (PlanckError E) {}
#endif
  }
  {
  crangeset<I> b=makeCRS<I>("1 11 30 41");
  planck_assert(!b.contains(0),"error");
  planck_assert(b.contains(1),"error");
  planck_assert(b.contains(5),"error");
  planck_assert(b.contains(10),"error");
  planck_assert(!b.contains(11),"error");
  planck_assert(!b.contains(29),"error");
  planck_assert(b.contains(30),"error");
  planck_assert(b.contains(35),"error");
  planck_assert(b.contains(40),"error");
  planck_assert(!b.contains(41),"errror");
  }
  {
  planck_assert(makeCRS<I>("1 11 20 31 40 56")==
                makeCRS<I>("20 31 40 51").op_or
                (makeCRS<I>("1 11 45 56")),"error");
  planck_assert(makeCRS<I>("1 11 45 56")==
                makeCRS<I>("").op_or
                (makeCRS<I>("1 11 45 56")),"error");
  planck_assert(makeCRS<I>("1 11 45 56")==
                makeCRS<I>("1 11 45 56").op_or
                (makeCRS<I>("")),"error");
  }
  {
  planck_assert(makeCRS<I>("22 24 45 51")==
                makeCRS<I>("20 31 40 51").op_and
                (makeCRS<I>("1 11 22 24 45 56")),"error");
  planck_assert(makeCRS<I>("20 31 40 51 90 101 110 121 200 201")==
                makeCRS<I>("10 101 110 121 200 221").op_and
                (makeCRS<I>("20 31 40 51 90 201")),"error");
  planck_assert(makeCRS<I>("")==
                makeCRS<I>("20 31 40 51").op_and
                (makeCRS<I>("")),"error");
  planck_assert(makeCRS<I>("")==
                makeCRS<I>("").op_and
                (makeCRS<I>("20 31 40 51")),"error");
  }
  {
  planck_assert(makeCRS<I>("20 31 40 45")==
                makeCRS<I>("20 31 40 51").op_andnot
                (makeCRS<I>("1 11 45 56")),"error");
  planck_assert(makeCRS<I>("")==
                makeCRS<I>("").op_andnot
                (makeCRS<I>("1 11 45 56")),"error");
  planck_assert(makeCRS<I>("1 11 45 56")==
                makeCRS<I>("1 11 45 56").op_andnot
                (makeCRS<I>("")),"error");
  }
  {
  crangeset<I> b = makeCRS<I>("20 31 40 51");

  planck_assert(!b.contains(0,11),"error");
  planck_assert(!b.contains(10,21),"error");
  planck_assert(!b.contains(19,20),"error");
  planck_assert(b.contains(20,21),"error");
  planck_assert(b.contains(21,22),"error");
  planck_assert(b.contains(20,31),"error");
  planck_assert(!b.contains(25,36),"error");
  planck_assert(b.contains(30,31),"error");
  planck_assert(!b.contains(31,32),"error");
  planck_assert(!b.contains(35,38),"error");
  planck_assert(!b.contains(35,46),"error");
  planck_assert(b.contains(40,41),"error");
  planck_assert(!b.contains(45,56),"error");
  planck_assert(!b.contains(60,71),"error");
  }
  {
  crangeset<I> b = makeCRS<I>("20 31 40 51");

  planck_assert(b.contains(makeCRS<I>("20 31 40 51")),"error");
  planck_assert(b.contains(makeCRS<I>("20 21")),"error");
  planck_assert(b.contains(makeCRS<I>("50 51")),"error");
  planck_assert(!b.contains(makeCRS<I>("19 31 40 51")),"error");
  planck_assert(!b.contains(makeCRS<I>("20 31 40 52")),"error");
  planck_assert(!b.contains(makeCRS<I>("20 51")),"error");
  planck_assert(!b.contains(makeCRS<I>("0 1")),"error");
  planck_assert(!b.contains(makeCRS<I>("0 20 31 40 51 100")),"error");
  }
  {
  crangeset<I> b = makeCRS<I>("20 31 40 51");

  planck_assert(!b.overlaps(0,11),"error");
  planck_assert(b.overlaps(10,21),"error");
  planck_assert(!b.overlaps(19,20),"error");
  planck_assert(b.overlaps(20,21),"error");
  planck_assert(b.overlaps(21,22),"error");
  planck_assert(b.overlaps(20,31),"error");
  planck_assert(b.overlaps(25,36),"error");
  planck_assert(b.overlaps(30,37),"error");
  planck_assert(!b.overlaps(31,32),"error");
  planck_assert(!b.overlaps(35,38),"error");
  planck_assert(b.overlaps(35,46),"error");
  planck_assert(b.overlaps(40,41),"error");
  planck_assert(b.overlaps(45,56),"error");
  planck_assert(!b.overlaps(60,71),"error");
  }
  {
  crangeset<I> b = makeCRS<I>("20 31 40 51");

  planck_assert(b.overlaps(makeCRS<I>("20 31 40 51")),"error");
  planck_assert(b.overlaps(makeCRS<I>("20 21")),"error");
  planck_assert(b.overlaps(makeCRS<I>("50 51")),"error");
  planck_assert(b.overlaps(makeCRS<I>("19 31 40 51")),"error");
  planck_assert(b.overlaps(makeCRS<I>("20 31 40 52")),"error");
  planck_assert(b.overlaps(makeCRS<I>("20 51")),"error");
  planck_assert(!b.overlaps(makeCRS<I>("0 1")),"error");
  planck_assert(!b.overlaps(makeCRS<I>("0 20 31 40 51 100")),"error");
  }
  {
  const int niter = 1000;
  for (int iter=0; iter<niter; ++iter)
    {
    crangeset<I> a = randomCRangeSet<I>(1000, 0, 100);
    crangeset<I> b = randomCRangeSet<I>(1000, 0, 100);
    crangeset<I> c = randomCRangeSet<I>(10, 10000, 100);
    crsOps(a,b);
    crsOps(a,c);
    crsOps(c,a);
    }
  }
  }
template<typename I> void check_Moc()
  {
  cout << "testing MOC " << bname<I>() << endl;
  Moc<I> moc;
  moc.addPixelRange(0,4,5);
  moc.addPixelRange(0,6,7);
  moc.addPixelRange(2,4,17);
  moc.addPixelRange(10,3000000,3000001);

  planck_assert(moc==moc.complement().complement(),"error");
  planck_assert(moc==Moc<I>::fromUniq(moc.toUniq()),"error");
  planck_assert(moc.maxOrder()==10,"error");
  Moc<I> xtmp = moc.degradedToOrder(8,false);
  planck_assert(moc.contains(xtmp),"error");
  planck_assert(!xtmp.contains(moc),"error");
  planck_assert(xtmp.overlaps(moc),"error");
  xtmp=moc.degradedToOrder(8,true);
  planck_assert(!moc.contains(xtmp),"error");
  planck_assert(xtmp.contains(moc),"error");
  planck_assert(xtmp.overlaps(moc),"error");
  planck_assert(moc==Moc<I>::fromCompressed(moc.toCompressed()),"error");
#if 0
  assertEquals("inconsistency",moc,MocUtil.mocFromString(" 0/4, 6 2/ \t 4 -16 10/3000000 \t\n "));
  assertEquals("inconsistency",moc,MocUtil.mocFromString("0/6 2/ 5 2/4 2/6- 16 0/4  10/3000000"));
  assertEquals("inconsistency",moc,MocUtil.mocFromString
    ("{\"0\":[6] , \"2\": [5 ], \"2\":[  4,6,7,8,9,10,11,12,13,14,15,16], \"0\":[4],  \"10\":[3000000]}"));
  assertEquals("inconsistency",moc,MocUtil.mocFromString(MocUtil.mocToStringASCII(moc)));
  assertEquals("inconsistency",moc,MocUtil.mocFromString(MocUtil.mocToStringJSON(moc)));
  ByteArrayOutputStream out= new ByteArrayOutputStream();
  MocUtil.mocToFits(moc,out);
  ByteArrayInputStream inp = new ByteArrayInputStream(out.toByteArray());
  assertEquals("inconsistency",moc,MocUtil.mocFromFits(inp));
#endif
  {
  tsize niter = 100;
  Moc<I> full; full.addPixelRange(0,0,12);
  Moc<I> empty;
  for (tsize iter=0; iter<niter; ++iter)
    {
    Moc<I> a = randomMoc<I>(1000, 0, 100);
    planck_assert(a.complement().complement()==a,"error");
    planck_assert(!a.overlaps(a.complement()),"error");
    planck_assert(a.op_or(a.complement())==full,"error");
    planck_assert(a.op_and(a.complement())==empty,"error");
#if 0
    write_Moc_to_fits("!healpixtestmoctmp",a);
    planck_assert(a==read_Moc_from_fits<I>("healpixtestmoctmp"),"FITS problem");
#endif
    }
  }
  }
template<typename I> void check_compress()
  {
  cout << "testing interpolation coding " << bname<I>() << endl;
  planck_assert(trailingZeros(4)==2,"error");
  planck_assert(trailingZeros(5)==0,"error");
  planck_assert(trailingZeros(int64(1)<<48)==48,"error");
  for (tsize x=0; x<100; ++x)
    {
    rangeset<I> a = randomRangeSet<I>(1000, 0, 100);
    rangeset<I> b;
    for (tsize i=0; i<a.nranges(); ++i)
      b.append(a.ivbegin(i)<<6,a.ivend(i)<<6);
    const vector<I> &v=b.data();
    obitstream obs;
    interpol_encode(v.begin(),v.end(),obs);
    vector<uint8> comp=obs.state();
    vector<I> v2;
    ibitstream ibs(comp);
    interpol_decode(v2,ibs);
    planck_assert(v==v2,"data mismatch");
    }
  }

template<typename I> void check_ringnestring()
  {
  cout << "testing ring2nest(nest2ring(m))==m " << bname<I>() << endl;
  for (int order=0; order<=T_Healpix_Base<I>::order_max; ++order)
    {
    T_Healpix_Base<I> base (order,RING);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      if (base.ring2nest(base.nest2ring(pix))!=pix)
        FAIL(cout<<"  PROBLEM: order = "<<order<<", pixel = "<<pix<<endl)
      }
    }
  }

template<typename I> void check_nestpeanonest()
  {
  cout << "testing peano2nest(nest2peano(m))==m " << bname<I>() << endl;
  for (int order=0; order<=T_Healpix_Base<I>::order_max; ++order)
    {
    T_Healpix_Base<I> base (order,NEST);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      if (base.peano2nest(base.nest2peano(pix))!=pix)
        FAIL(cout<<"  PROBLEM: order = "<<order<<", pixel = "<<pix<<endl)
      }
    }
  }

template<typename I> void check_pixzphipix()
  {
  cout << "testing zphi2pix(pix2zphi(m))==m " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base1 (order,RING), base2 (order,NEST);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi;
      I pix = I(rng.rand_uni()*base1.Npix());
      base1.pix2zphi(pix,z,phi);
      if (base1.zphi2pix(z,phi)!=pix)
        FAIL(cout<<"  PROBLEM: order = "<<order<<", pixel = "<<pix<<endl)
      base2.pix2zphi(pix,z,phi);
      if (base2.zphi2pix(z,phi)!=pix)
        FAIL(cout<<"  PROBLEM: order = "<<order<<", pixel = "<<pix<<endl)
      }
    }
  for (I nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi;
      I pix = I(rng.rand_uni()*base.Npix());
      base.pix2zphi(pix,z,phi);
      if (base.zphi2pix(z,phi)!=pix)
        FAIL(cout<<"  PROBLEM: nside = "<<nside<<", pixel = "<<pix<<endl)
      }
    }
  }

template<typename I> void check_zphipixzphi()
  {
  cout << "testing pix2zphi(zphi2pix(ptg)) approx zphi " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base1 (order,NEST), base2 (order,RING);
    double mincos = min (cos(base1.max_pixrad()),0.999999999999999);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi,z2,phi2;
      random_zphi (z,phi);
      base1.pix2zphi(base1.zphi2pix(z,phi),z2,phi2);
      if (cosdist_zphi(z,phi,z2,phi2)<mincos)
        FAIL(cout << "  PROBLEM: order = " << order
                  << ", zphi = " << z << ", " << phi << endl)
      base2.pix2zphi(base2.zphi2pix(z,phi),z2,phi2);
      if (cosdist_zphi(z,phi,z2,phi2)<mincos)
        FAIL(cout << "  PROBLEM: order = " << order
                  << ", zphi = " << z << ", " << phi << endl)
      }
    }
  for (int nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    double mincos = min (cos(base.max_pixrad()),0.999999999999999);
    for (int m=0; m<nsamples; ++m)
      {
      double z,phi,z2,phi2;
      random_zphi (z,phi);
      base.pix2zphi(base.zphi2pix(z,phi),z2,phi2);
      if (cosdist_zphi(z,phi,z2,phi2)<mincos)
        FAIL(cout << "  PROBLEM: nside = " << nside
                  << ", zphi = " << z << ", " << phi << endl)
      }
    }
  }

template<typename I> void check_pixangpix()
  {
  cout << "testing ang2pix(pix2ang(m))==m " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base1 (order,RING), base2 (order,NEST);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base1.Npix());
      if (base1.ang2pix(base1.pix2ang(pix))!=pix)
        FAIL(cout<<"  PROBLEM: order = "<<order<<", pixel = "<<pix<<endl)
      if (base2.ang2pix(base2.pix2ang(pix))!=pix)
        FAIL(cout<<"  PROBLEM: order = "<<order<<", pixel = "<<pix<<endl)
      }
    }
  for (I nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    for (int m=0; m<nsamples; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      if (base.ang2pix(base.pix2ang(pix))!=pix)
        FAIL(cout<<"  PROBLEM: nside = "<<nside<<", pixel = "<<pix<<endl)
      }
    }
  }

template<typename I> void check_neighbors()
  {
  cout << "testing neighbor function " << bname<I>() << endl;
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,NEST), base2(order,RING);
    double maxang = 2.01*base.max_pixrad();
    for (int m=0; m<nsamples/10; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      fix_arr<I,8> nb,nb2;
      vec3 pixpt = base.pix2vec(pix);
      base.neighbors(pix,nb);
      base2.neighbors(base.nest2ring(pix),nb2);
      for (int n=0; n<8; ++n)
        if (nb[n]<0)
          planck_assert(nb2[n]<0,"neighbor inconsistency");
        else
          planck_assert(base.nest2ring(nb[n])==nb2[n],"neighbor inconsistency");
      sort(&nb[0],&nb[0]+8);
      int check=0;
      for (int n=0; n<8; ++n)
        {
        if (nb[n]<0)
          ++check;
        else
          {
          if (v_angle(base.pix2vec(nb[n]),pixpt)>maxang)
            FAIL(cout<<"  PROBLEM: order = "<<order<<", pix = "<<pix<<endl)
          if ((n>0) && (nb[n]==nb[n-1]))
            FAIL(cout<<"  PROBLEM: order = "<<order<<", pix = "<<pix<<endl)
          }
        }
      planck_assert((check<=1)||((order==0)&&(check<=2)),"too few neighbors");
      }
    }
  for (I nside=3; nside<(I(1)<<omax); nside+=nside/2+1)
    {
    T_Healpix_Base<I> base (nside,RING,SET_NSIDE);
    double maxang = 2.01*base.max_pixrad();
    for (int m=0; m<nsamples/10; ++m)
      {
      I pix = I(rng.rand_uni()*base.Npix());
      fix_arr<I,8> nb;
      vec3 pixpt = base.pix2vec(pix);
      base.neighbors(pix,nb);
      for (int n=0; n<8; ++n)
        if ((nb[n]>=0) && (v_angle(base.pix2vec(nb[n]),pixpt)>maxang))
          FAIL(cout<<"  PROBLEM: nside = "<<nside<<", pix = "<<pix<<endl)
      }
    }
  }

void check_swap_scheme()
  {
  cout << "testing whether double swap_scheme() returns the original map"
       << endl << "(for orders 0 to 10)." << endl;
  for (int order=0; order<=10; ++order)
    {
    Healpix_Map<uint8> map(order,NEST);
    for (int m=0; m<map.Npix(); ++m) map[m]=uint8(m&0xFF);
    map.swap_scheme();
    map.swap_scheme();
    for (int m=0; m<map.Npix(); ++m)
      if (map[m]!=(m&0xFF))
        FAIL(cout<<"  PROBLEM: order = "<<order<<", pix = "<<m<<endl)
    }
  }

void check_issue_229 (Healpix_Ordering_Scheme scheme)
  {
  cout << "checking issue #229" << endl;
  int order=8;
  Healpix_Map<bool> map (order,scheme);
  map.fill(false);
  Healpix_Map<vec3> vmap(order,scheme);
  for (int m=0; m<vmap.Npix(); ++m)
    vmap[m]=vmap.pix2vec(m);
  pointing ptg(halfpi-0.1,0);
  double rad=0.1;
  auto pixset = map.query_disc(ptg,rad);
  vec3 vptg=ptg;
  double cosrad=cos(rad);
  for (tsize j=0; j<pixset.nranges(); ++j)
    for (int i=pixset.ivbegin(j); i<pixset.ivend(j); ++i)
      map[i] = true;
  for (int i=0; i<map.Npix(); ++i)
    {
    bool inside = dotprod(vmap[i],vptg)>cosrad;
    if (inside^map[i])
      FAIL(cout << "  PROBLEM: issue 229" << endl)
    }
  }

void check_query_disc_strict (Healpix_Ordering_Scheme scheme)
  {
  cout << "testing whether all pixels found by query_disc() really" << endl
       << "lie inside the disk (and vice versa)" << endl;
  cout << "Ordering scheme: " << (scheme==RING ? "RING" : "NEST") << endl;
  rng.seed(42);
  for (int order=0; order<=5; ++order)
    {
    Healpix_Map<bool> map (order,scheme);
    map.fill(false);
    Healpix_Map<vec3> vmap(order,scheme);
    for (int m=0; m<vmap.Npix(); ++m)
      vmap[m]=vmap.pix2vec(m);
    for (int m=0; m<100000; ++m)
      {
      pointing ptg;
      random_dir (ptg);
      double rad = pi/1 * rng.rand_uni();
      auto pixset = map.query_disc(ptg,rad);
      vec3 vptg=ptg;
      double cosrad=cos(rad);
      for (tsize j=0; j<pixset.nranges(); ++j)
        for (int i=pixset.ivbegin(j); i<pixset.ivend(j); ++i)
          map[i] = true;
      for (int i=0; i<map.Npix(); ++i)
        {
        bool inside = dotprod(vmap[i],vptg)>cosrad;
        if (inside^map[i])
          FAIL(cout<<"  PROBLEM: order = "<<order<<", ptg = "<<ptg<<endl)
        }
      for (tsize j=0; j<pixset.nranges(); ++j)
        for (int i=pixset.ivbegin(j); i<pixset.ivend(j); ++i)
          map[i] = false;
      }
    }
  }

template<typename I>void check_query_disc()
  {
  cout << "checking query_disc() " << bname<I>() << endl;
  rng.seed(48);
  int omax=min<int>(20,T_Healpix_Base<I>::order_max);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> rbase (order,RING), nbase (order,NEST);
    int niter=max(1,min(1000,100000>>order));
    for (int m=0; m<niter; ++m)
      {
      pointing ptg;
      random_dir (ptg);
      double rad = pi/1 * rng.rand_uni();
      auto pixset = rbase.query_disc(ptg,rad);
      rangeset<I> pslast=pixset;
      for (tsize fct=5; fct>0; --fct)
        {
        auto psi = rbase.query_disc_inclusive(ptg,rad,fct);
        if (!psi.contains(pslast))
          cout << "  Potential problem: RING pixel sets inconsistent" << endl;
        swap(pslast,psi);
        }
      I nval = pixset.nval();
      pixset = nbase.query_disc(ptg,rad);
      pslast=pixset;
      for (tsize fct=8; fct>0; fct>>=1)
        {
        auto psi = nbase.query_disc_inclusive(ptg,rad,fct);
        if (!psi.contains(pslast))
          FAIL(cout << "  PROBLEM: NEST pixel sets inconsistent" << endl)
        swap(pslast,psi);
        }
      if (nval!=pixset.nval())
        FAIL(cout << "  PROBLEM: number of pixels different: "
                  << nval << " vs. " << pixset.nval() << endl)
      }
    }
  }
template<typename I>void check_query_polygon()
  {
  cout << "checking query_polygon() " << bname<I>() << endl;
  rng.seed(42);
  int omax=min<int>(20,T_Healpix_Base<I>::order_max);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> rbase (order,RING), nbase (order,NEST);
    int niter=max(1,min(1000,100000>>order));
    for (int m=0; m<niter; ++m)
      {
      vector<pointing> corner(3);
      random_dir(corner[0]); random_dir(corner[1]); random_dir(corner[2]);
      auto pixset = rbase.query_polygon(corner);
      I nval = pixset.nval();
      pixset = nbase.query_polygon(corner);
      if (nval!=pixset.nval())
        FAIL(cout << "  PROBLEM: number of pixels different: "
                  << nval << " vs. " << pixset.nval() << endl)
      pixset = rbase.query_polygon_inclusive(corner,4);
      I nv1=pixset.nval();
      pixset = nbase.query_polygon_inclusive(corner,4);
      I nv2=pixset.nval();
      if (nv1<nv2)
        FAIL(cout << "  PROBLEM: inclusive(RING)<inclusive(NEST): "
                  << nv1 << " vs. " << nv2 << endl)
      if (nv2<nval)
        FAIL(cout << "  PROBLEM: inclusive(NEST)<non-inclusive: "
                  << nv2 << " vs. " << nval << endl)
      }
    }
  }

void helper_oop (int order)
  {
  Healpix_Map<double> map (order,RING), map2 (order,NEST), map3 (order,RING);
  for (int m=0; m<map.Npix(); ++m) map[m] = rng.rand_uni()+0.01;
  map2.Import(map);
  map3.Import(map2);
  for (int m=0; m<map.Npix(); ++m)
    if (!approx(map[m],map3[m],1e-12))
      FAIL(cout << "PROBLEM: order = " << order << endl)
  }
void helper_udgrade (int order, Healpix_Ordering_Scheme s1,
  Healpix_Ordering_Scheme s2)
  {
  Healpix_Map<double> map (order,s1), map2 (order+2,s2), map3 (order, s1);
  for (int m=0; m<map.Npix(); ++m) map[m] = rng.rand_uni()+0.01;
  map2.Import(map);
  map3.Import(map2);
  for (int m=0; m<map.Npix(); ++m)
    if (!approx(map[m],map3[m],1e-12))
      FAIL(cout << "PROBLEM: order = " << order << endl)
  }
void helper_udgrade2 (int nside)
  {
  Healpix_Map<double> map (nside,RING,SET_NSIDE), map2 (nside*3,RING,SET_NSIDE),
    map3 (nside, RING,SET_NSIDE);
  for (int m=0; m<map.Npix(); ++m) map[m] = rng.rand_uni()+0.01;
  map2.Import(map);
  map3.Import(map2);
  for (int m=0; m<map.Npix(); ++m)
    if (!approx(map[m],map3[m],1e-12))
      FAIL(cout << "PROBLEM: nside = " << nside << endl)
  }

void check_import()
  {
  cout << "testing out-of-place swapping" << endl;
  for (int order=0; order<=7; ++order)
    helper_oop(order);
  cout << "testing downgrade(upgrade(map)) == map" << endl;
  for (int order=0; order<=7; ++order)
    {
    helper_udgrade(order,RING,RING);
    helper_udgrade(order,RING,NEST);
    helper_udgrade(order,NEST,NEST);
    helper_udgrade(order,NEST,RING);
    }
  for (int nside=3; nside<500; nside+=nside/2+1)
    helper_udgrade2(nside);
  }

void check_average()
  {
  cout << "testing whether average(map) == average(downgraded map)" << endl;
  for (int order=1; order<=10; ++order)
    {
    Healpix_Map<double> map (order,RING), map2(1,RING);
    for (int m=0; m<map.Npix(); ++m)
      map[m] = rng.rand_uni()+0.01;
    map2.Import(map);
    double avg=map.average(), avg2=map2.average();
    if (!approx(avg,avg2,1e-13))
      FAIL(cout << "PROBLEM: order = " << order << " " << avg/avg2-1 << endl)
    }
  for (int nside=3; nside<1000; nside += nside/2+1)
    {
    Healpix_Map<double> map (nside,RING,SET_NSIDE), map2(1,RING,SET_NSIDE);
    for (int m=0; m<map.Npix(); ++m)
      map[m] = rng.rand_uni()+0.01;
    map2.Import(map);
    double avg=map.average(), avg2=map2.average();
    if (!approx(avg,avg2,1e-13))
      FAIL(cout << "PROBLEM: nside = " << nside << " " << avg/avg2-1 << endl)
    }
  }

void random_alm (Alm<xcomplex<double> >&almT, Alm<xcomplex<double> >&almG,
  Alm<xcomplex<double> >&almC, int lmax, int mmax)
  {
  almT.Set(lmax,mmax); almG.Set(lmax,mmax); almC.Set(lmax,mmax);

  for (int l=0; l<=lmax; ++l)
    {
    almT(l,0)=dcomplex(rng.rand_gauss(),0.);
    almG(l,0)=dcomplex(rng.rand_gauss(),0.);
    almC(l,0)=dcomplex(rng.rand_gauss(),0.);
    }

  for (int m=1; m<=mmax; ++m)
    for (int l=m; l<=lmax; ++l)
      {
      almT(l,m).real(rng.rand_gauss()); almT(l,m).imag(rng.rand_gauss());
      almG(l,m).real(rng.rand_gauss()); almG(l,m).imag(rng.rand_gauss());
      almC(l,m).real(rng.rand_gauss()); almC(l,m).imag(rng.rand_gauss());
      }
  almG(0,0)=almC(0,0)=almG(1,0)=almC(1,0)=almG(1,1)=almC(1,1)=0;
  }

void random_alm (Alm<xcomplex<double> >&alm, int lmax, int mmax)
  {
  alm.Set(lmax,mmax);

  for (int l=0; l<=lmax; ++l)
    { alm(l,0)=dcomplex(rng.rand_gauss(),0.); }

  for (int m=1; m<=mmax; ++m)
    for (int l=m; l<=lmax; ++l)
      { alm(l,m).real(rng.rand_gauss()); alm(l,m).imag(rng.rand_gauss()); }
  }

void check_alm (const Alm<xcomplex<double> >&oalm,
  const Alm<xcomplex<double> >&alm, double epsilon)
  {
  for (int m=0; m<=alm.Mmax(); ++m)
    for (int l=m; l<=alm.Lmax(); ++l)
      {
      if (!abs_approx(oalm(l,m).real(),alm(l,m).real(),epsilon))
        FAIL(cout << "Problemr " << l << " " << m << endl)
      if (!abs_approx(oalm(l,m).imag(),alm(l,m).imag(),epsilon))
        FAIL(cout << "Problemi " << l << " " << m <<  endl)
      }
  }

void check_alm2map2alm (int lmax, int mmax, int nside)
  {
  cout << "testing whether a_lm->map->a_lm returns original a_lm" << endl;
  cout << "lmax=" << lmax <<", mmax=" << mmax << ", nside=" << nside << endl;
  const double epsilon = 1e-8;
  Alm<xcomplex<double> > oalmT(lmax,mmax),almT(lmax,mmax),
    oalmG(lmax,mmax),almG(lmax,mmax),oalmC(lmax,mmax),almC(lmax,mmax);
  Healpix_Map<double> mapT(nside,RING,SET_NSIDE), mapQ(nside,RING,SET_NSIDE),
    mapU(nside,RING,SET_NSIDE);

  const double eps0=1e-11;
  const double eps=1e-12;

  random_alm(oalmT,oalmG,oalmC,lmax,mmax);
  alm2map(oalmT,mapT);
  map2alm_iter2(mapT,almT,eps0,eps0);
  check_alm (oalmT, almT, epsilon);

  alm2map_spin(oalmG,oalmC,mapQ,mapU,1);
  map2alm_spin_iter2(mapQ,mapU,almG,almC,1,eps,eps);
  check_alm (oalmG, almG, epsilon);
  check_alm (oalmC, almC, epsilon);

  alm2map_pol(oalmT,oalmG,oalmC,mapT,mapQ,mapU);
  map2alm_pol_iter2(mapT,mapQ,mapU,almT,almG,almC,max(eps,eps0),max(eps,eps0));
  check_alm (oalmT, almT, epsilon);
  check_alm (oalmG, almG, epsilon);
  check_alm (oalmC, almC, epsilon);
  }

void check_smooth_alm ()
  {
  cout << "testing whether unsmooth(smooth(a_lm)) returns a_lm" << endl;
  const double epsilon = 1e-14;
  const double fwhm = 100.*arcmin2rad;
  const int lmax=300, mmax=300;
  Alm<xcomplex<double> > oalmT(lmax,mmax),almT(lmax,mmax),
    oalmG(lmax,mmax),almG(lmax,mmax),oalmC(lmax,mmax),almC(lmax,mmax);

  random_alm(oalmT,oalmG,oalmC,lmax,mmax);

  almT=oalmT; almG=oalmG; almC=oalmC;
  smoothWithGauss (almT, fwhm);
  smoothWithGauss (almT, -fwhm);
  check_alm (oalmT, almT, epsilon);
  almT=oalmT;
  smoothWithGauss (almT, almG, almC, fwhm);
  smoothWithGauss (almT, almG, almC, -fwhm);
  check_alm (oalmT, almT, epsilon);
  check_alm (oalmG, almG, epsilon);
  check_alm (oalmC, almC, epsilon);
  }

void check_rot_alm ()
  {
  cout << "testing whether rot^-1(rot(a_lm)) returns a_lm" << endl;
  const double epsilon = 2e-13;
  const int lmax=300;
  Alm<xcomplex<double> > oalm(lmax,lmax),alm(lmax,lmax);

  random_alm(oalm,lmax,lmax);

  alm=oalm;
  rotate_alm (alm,3,4,5);
  rotate_alm (alm,-5,-4,-3);
  check_alm (oalm, alm, epsilon);
  }

double map_residual (const Healpix_Map<double> &a, const Healpix_Map<double> &b)
  {
  planck_assert(a.conformable(b), "maps are not conformable");
  double res=0;
  for (int i=0; i<a.Npix(); ++i)
    {
    double diff=a[i]-b[i];
    res+=diff*diff;
    }
  return res;
  }

void check_ringweights ()
  {
  cout << "testing the accuracy of ring weights" << endl;
  int nside=128;
  double dummy;
  auto rwgt0=get_ringweights(nside,3*nside,1e-6,3000,dummy);
  arr<double> rwgt(2*nside);
  for (tsize i=0; i<rwgt.size(); ++i) rwgt[i]=rwgt0[i]+1.;
  Alm<xcomplex<double> > alm; random_alm(alm,1.5*nside,0);
  Healpix_Map<double> omap(nside,RING,SET_NSIDE),map2(omap);
  alm2map(alm,omap);
  map2alm(omap,alm,arr<double>(2*nside,1.));
  alm2map(alm,map2);
  double ressimple = map_residual(omap,map2);
  map2alm(omap,alm,rwgt);
  alm2map(alm,map2);
  double reswgt = map_residual(omap,map2);
  double resquot=sqrt(reswgt/ressimple);
  if (resquot>1e-5)
    FAIL(cout << "  PROBLEM with ringweights: resquot =  " << resquot << endl)
  }

void check_fullweights ()
  {
  cout << "testing the accuracy of pixel weights" << endl;
  int nside=128;
  double dummy;
  auto fwgt=get_fullweights(nside,3*nside,1e-6,3000,dummy);
  Alm<xcomplex<double> > alm; random_alm(alm,1.5*nside,1.5*nside);
  Healpix_Map<double> omap(nside,RING,SET_NSIDE),map2(omap);
  alm2map(alm,omap);
  map2alm(omap,alm,arr<double>(2*nside,1.));
  alm2map(alm,map2);
  double ressimple = map_residual(omap,map2);
  map2=omap;
  apply_fullweights(map2,fwgt);
  map2alm(map2,alm,arr<double>(2*nside,1.));
  alm2map(alm,map2);
  double reswgt = map_residual(omap,map2);
  double resquot=sqrt(reswgt/ressimple);
  if (resquot>1e-5)
    FAIL(cout << "  PROBLEM with full weights: resquot =  " << resquot << endl)
  }

void check_isqrt()
  {
  cout << "testing whether isqrt() works reliably" << endl;
  uint64 val=uint64(0xF234)<<16, valsq=val*val;
  if (isqrt(valsq)!=val) FAIL(cout << "PROBLEM1" << endl)
  if (isqrt(valsq-1)!=val-1) FAIL(cout << "PROBLEM2" << endl)
  }

void check_pix2ang_acc()
  {
  cout << "testing accuracy of pix2ang at the poles" << endl;
  for (int m=0; m<=29;++m)
    {
    Healpix_Base2 base(m,RING);
    if (base.pix2ang(1).theta==0.)
      FAIL(cout << "PROBLEM: order " << m << endl)
    }
  }

const int nsteps=1000000;

template<typename I>void perf_neighbors(const string &name,
  Healpix_Ordering_Scheme scheme)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    fix_arr<I,8> nres;
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      base.neighbors(pix,nres);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_pix2ang(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      pointing p(base.pix2ang(pix));
      dummy+=p.theta+p.phi;
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_pix2vec(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      vec3 v(base.pix2vec(pix));
      dummy+=v.x+v.y+v.z;
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_pix2zphi(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    I dpix=max(base.Npix()/nsteps,I(1));
    double z,phi;
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      base.pix2zphi(pix,z,phi);
      dummy+=z+phi;
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_zphi2pix(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  double dz=2./sqrt(nsteps);
  double dph=twopi/sqrt(nsteps);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    for (double z=-1; z<1; z+=dz)
      for (double phi=0; phi<twopi; phi+=dph)
        {
        dummy+=base.zphi2pix(z,phi);
        ++cnt;
        }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_ang2pix(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  double dth=pi/sqrt(nsteps);
  double dph=twopi/sqrt(nsteps);
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,scheme);
    for (double theta=0; theta<pi; theta+=dth)
      for (double phi=0; phi<twopi; phi+=dph)
        {
        dummy+=base.ang2pix(pointing(theta+1e-15*phi,phi));
        ++cnt;
        }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_ring2nest(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,RING);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.ring2nest(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_nest2ring(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,RING);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.nest2ring(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_peano2nest(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,NEST);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.peano2nest(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_nest2peano(const string &name,double &dummy)
  {
  tsize cnt=0;
  wallTimers.start(name);
  int omax=T_Healpix_Base<I>::order_max;
  for (int order=0; order<=omax; ++order)
    {
    T_Healpix_Base<I> base (order,NEST);
    I dpix=max(base.Npix()/nsteps,I(1));
    for (I pix=0; pix<base.Npix(); pix+=dpix)
      {
      dummy+=base.nest2peano(pix);
      ++cnt;
      }
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-6 << "MOps/s" << endl;
  }
template<typename I>void perf_query_disc(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  T_Healpix_Base<I> base(1024,scheme,SET_NSIDE);
  wallTimers.start(name);
  for (int m=0; m<1000; ++m)
    {
    dummy += base.query_disc(vec3(1,0,0),halfpi/9).nranges();
    ++cnt;
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-3 << "kOps/s" << endl;
  }
template<typename I>void perf_query_triangle(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  T_Healpix_Base<I> base(1024,scheme,SET_NSIDE);
  vector<pointing> corner;
  corner.push_back(vec3(1,0.01,0.01));
  corner.push_back(vec3(0.01,1,0.01));
  corner.push_back(vec3(0.01,0.01,1));
  wallTimers.start(name);
  for (int m=0; m<1000; ++m)
    {
    dummy += base.query_polygon(corner).nranges();
    ++cnt;
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-3 << "kOps/s" << endl;
  }
template<typename I>void perf_query_polygon(const string &name,
  Healpix_Ordering_Scheme scheme, double &dummy)
  {
  tsize cnt=0;
  T_Healpix_Base<I> base(1024,scheme,SET_NSIDE);
  vector<pointing> corner;
  corner.push_back(vec3(1,0.01,0.01));
  corner.push_back(vec3(1,1,-0.3));
  corner.push_back(vec3(0.01,1,0.01));
  corner.push_back(vec3(0.01,0.01,1));
  wallTimers.start(name);
  for (int m=0; m<1000; ++m)
    {
    dummy += base.query_polygon(corner).nranges();
    ++cnt;
    }
  wallTimers.stop(name);
  cout << name << ": " << cnt/wallTimers.acc(name)*1e-3 << "kOps/s" << endl;
  }

void perftest()
  {
  double dummy=0;
  cout << "Measuring performance of Healpix_Base methods." << endl;
  perf_pix2zphi<int>   ("pix2zphi (RING):int  ",RING,dummy);
  perf_pix2zphi<int>   ("pix2zphi (NEST):int  ",NEST,dummy);
  perf_pix2zphi<int64> ("pix2zphi (RING):int64",RING,dummy);
  perf_pix2zphi<int64> ("pix2zphi (NEST):int64",NEST,dummy);
  perf_zphi2pix<int>   ("zphi2pix (RING):int  ",RING,dummy);
  perf_zphi2pix<int>   ("zphi2pix (NEST):int  ",NEST,dummy);
  perf_zphi2pix<int64> ("zphi2pix (RING):int64",RING,dummy);
  perf_zphi2pix<int64> ("zphi2pix (NEST):int64",NEST,dummy);
  perf_pix2ang<int>    ("pix2ang  (RING):int  ",RING,dummy);
  perf_pix2ang<int>    ("pix2ang  (NEST):int  ",NEST,dummy);
  perf_pix2ang<int64>  ("pix2ang  (RING):int64",RING,dummy);
  perf_pix2ang<int64>  ("pix2ang  (NEST):int64",NEST,dummy);
  perf_ang2pix<int>    ("ang2pix  (RING):int  ",RING,dummy);
  perf_ang2pix<int>    ("ang2pix  (NEST):int  ",NEST,dummy);
  perf_ang2pix<int64>  ("ang2pix  (RING):int64",RING,dummy);
  perf_ang2pix<int64>  ("ang2pix  (NEST):int64",NEST,dummy);
  perf_pix2vec<int>    ("pix2vec  (RING):int  ",RING,dummy);
  perf_pix2vec<int>    ("pix2vec  (NEST):int  ",NEST,dummy);
  perf_pix2vec<int64>  ("pix2vec  (RING):int64",RING,dummy);
  perf_pix2vec<int64>  ("pix2vec  (NEST):int64",NEST,dummy);
  perf_neighbors<int>  ("neighbors(NEST):int  ",NEST);
  perf_neighbors<int>  ("neighbors(RING):int  ",RING);
  perf_neighbors<int64>("neighbors(NEST):int64",NEST);
  perf_neighbors<int64>("neighbors(RING):int64",RING);
  perf_ring2nest<int>  ("ring2nest      :int  ",dummy);
  perf_ring2nest<int64>("ring2nest      :int64",dummy);
  perf_nest2ring<int>  ("nest2ring      :int  ",dummy);
  perf_nest2ring<int64>("nest2ring      :int64",dummy);
  perf_peano2nest<int>  ("peano2nest     :int  ",dummy);
  perf_peano2nest<int64>("peano2nest     :int64",dummy);
  perf_nest2peano<int>  ("nest2peano     :int  ",dummy);
  perf_nest2peano<int64>("nest2peano     :int64",dummy);
  perf_query_disc<int>      ("query_disc    (RING):int  ",RING,dummy);
  perf_query_disc<int>      ("query_disc    (NEST):int  ",NEST,dummy);
  perf_query_disc<int64>    ("query_disc    (RING):int64",RING,dummy);
  perf_query_disc<int64>    ("query_disc    (NEST):int64",NEST,dummy);
  perf_query_triangle<int>  ("query_triangle(RING):int  ",RING,dummy);
  perf_query_triangle<int>  ("query_triangle(NEST):int  ",NEST,dummy);
  perf_query_triangle<int64>("query_triangle(RING):int64",RING,dummy);
  perf_query_triangle<int64>("query_triangle(NEST):int64",NEST,dummy);
  perf_query_polygon<int>   ("query_polygon (RING):int  ",RING,dummy);
  perf_query_polygon<int>   ("query_polygon (NEST):int  ",NEST,dummy);
  perf_query_polygon<int64> ("query_polygon (RING):int64",RING,dummy);
  perf_query_polygon<int64> ("query_polygon (NEST):int64",NEST,dummy);

  if (dummy<0) cout << dummy << endl;
  }

} // unnamed namespace

int main(int argc, const char **argv)
  {
  module_startup ("hpxtest",argc,argv,1,"");
  perftest();
  check_compress<int>();
  check_compress<unsigned>();
  check_compress<int64>();
  check_compress<uint64>();
  check_Moc<int>();
  check_Moc<int64>();
  check_rangeset<int>();
  check_rangeset<int64>();
  check_crangeset<int>();
  check_crangeset<int64>();
  check_isqrt();
  check_pix2ang_acc();
  check_smooth_alm();
  check_rot_alm();
  check_alm2map2alm(620,620,256);
  check_alm2map2alm(620,2,256);
  check_average();
  check_import();
  check_ringnestring<int>();
  check_ringnestring<int64>();
  check_nestpeanonest<int>();
  check_nestpeanonest<int64>();
  check_pixzphipix<int>();
  check_pixzphipix<int64>();
  check_zphipixzphi<int>();
  check_zphipixzphi<int64>();
  check_pixangpix<int>();
  check_pixangpix<int64>();
  check_neighbors<int>();
  check_neighbors<int64>();
  check_swap_scheme();
  check_query_disc_strict(RING);
  check_query_disc_strict(NEST);
  check_issue_229(RING);
  check_issue_229(NEST);
  check_query_disc<int>();
  check_query_disc<int64>();
  check_query_polygon<int>();
  check_query_polygon<int64>();
  check_ringweights();
  check_fullweights();
#ifdef UNITTESTS
  if (errcount>0) planck_fail("unit test errors");
#endif
  return 0;
  }
