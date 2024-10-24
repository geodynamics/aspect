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

/*! \file moc_query.cc
 *  Copyright (C) 2014-2015 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include "moc_query.h"
#include "geom_utils.h"
#include "lsconstants.h"

using namespace std;

namespace {

template<typename I> class queryHelper
  {
  private:
    int order, omax;
    bool inclusive;
    vector<MocQueryComponent> comp;
    vector<T_Healpix_Base<I> > base;
    vector<int> shortcut;
    arr<double> cr;
    arr2<double> crmin;
    arr2<double> crmax;

    vector<pair<I,int> > stk; // stack for pixel numbers and their orders
    I pix;
    int o;
    int stacktop; // a place to save a stack position
    vec3 pv;

    void check_pixel (int zone, Moc<I> &pixset)
      {
      if (zone==0) return;
      if (o<order)
        {
        if (zone>=3)
          pixset.appendPixel(o,pix); // output all subpixels
        else // (zone>=1)
          for (int i=0; i<4; ++i)
            stk.push_back(make_pair(4*pix+3-i,o+1)); // add children
        }
      else if (o>order) // this implies that inclusive==true
        {
        if (zone>=2) // pixel center in shape
          {
          pixset.appendPixel(order,pix>>(2*(o-order))); // output parent pixel
          stk.resize(stacktop); // unwind the stack
          }
        else // (zone>=1): pixel center in safety range
          {
          if (o<omax) // check sublevels
            for (int i=0; i<4; ++i) // add children in reverse order
              stk.push_back(make_pair(4*pix+3-i,o+1));
          else // at resolution limit
            {
            pixset.appendPixel(order,pix>>(2*(o-order))); // output parent pixel
            stk.resize(stacktop); // unwind the stack
            }
          }
        }
      else // o==order
        {
        if (zone>=2)
          pixset.appendPixel(order,pix);
        else if (inclusive) // and (zone>=1)
          {
          if (order<omax) // check sublevels
            {
            stacktop=stk.size(); // remember current stack position
            for (int i=0; i<4; ++i) // add children in reverse order
              stk.push_back(make_pair(4*pix+3-i,o+1));
            }
          else // at resolution limit
            pixset.appendPixel(order,pix); // output the pixel
          }
        }
      }

    void correctLoc(int &loc) const
      {
      int myloc=loc--;
      planck_assert((myloc>=0)&&(myloc<int(comp.size())),"inconsistency");
      for (int i=0; i<comp[myloc].nops; ++i)
        correctLoc(loc);
      }

    int getZone (int &loc, int zmin, int zmax) const
      {
      if (zmin==zmax) { loc=shortcut[loc]; return zmin; } // short-circuit
      int myloc=loc--;
//      planck_assert((myloc>=0)&&(myloc<int(comp.size())),"inconsistency");
      switch (comp[myloc].op)
        {
        case AND:
          {
          int z1=zmax;
          for (int i=0; i<comp[myloc].nops; ++i)
            z1 = getZone(loc,zmin,z1);
          return z1;
          }
        case OR:
          {
          int z1=zmin;
          for (int i=0; i<comp[myloc].nops; ++i)
            z1 = getZone(loc,z1,zmax);
          return z1;
          }
        case XOR:
          {
          int z1=getZone(loc,0,3);
          int z2=getZone(loc,0,3);
          return max(zmin,min(zmax,max(min(z1,3-z2),min(3-z1,z2))));
          }
        case NOT:
          return 3-getZone(loc,3-zmax,3-zmin);
        case NONE:
          {
          int res=zmax;
          double crad=dotprod(pv,comp[myloc].center);
          if (crad<=crmax(o,myloc)) res=0;
          else if (crad<=cr[myloc]) res=1;
          else if (crad<=crmin(o,myloc)) res=2;
          return max(zmin,min(zmax,res));
          }
        }
      planck_fail("must not get here");
      }

  public:
    queryHelper (int order_, int omax_, bool inclusive_,
      const vector<MocQueryComponent> &comp_)
      : order(order_), omax(omax_), inclusive(inclusive_), comp(comp_),
        base(omax+1), shortcut(comp.size()), cr(comp.size()),
        crmin(omax+1,comp.size()), crmax(omax+1,comp.size())
      {
      planck_assert(comp.size()>=1,"bad query component vector");
      planck_assert(order<=omax,"order>omax");
      if (!inclusive) planck_assert(order==omax,"inconsistency");
      planck_assert(omax<=T_Healpix_Base<I>::order_max,"omax too high");

      for (tsize i=0; i<comp.size(); ++i)
        if (comp[i].op==NONE) // it's a cap
          cr[i]=cos(comp[i].radius);
      for (o=0; o<=omax; ++o) // prepare data at the required orders
        {
        base[o].Set(o,NEST);
        double dr=base[o].max_pixrad(); // safety distance
        for (tsize i=0; i<comp.size(); ++i)
          if (comp[i].op==NONE) // it's a cap
            {
            crmax(o,i) = (comp[i].radius+dr>=pi) ? -1.01:cos(comp[i].radius+dr);
            crmin(o,i) = (comp[i].radius-dr<=0.) ?  1.01:cos(comp[i].radius-dr);
            }
        }

      for (tsize i=0; i<comp.size(); ++i)
        {
        int loc=int(i);
        correctLoc(loc);
        shortcut[i]=loc;
        }
      }
    Moc<I> result()
      {
      Moc<I> pixset;
      stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
      for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
        stk.push_back(make_pair(I(11-i),0));

      stacktop=0; // a place to save a stack position

      while (!stk.empty()) // as long as there are pixels on the stack
        {
        // pop current pixel number and order from the stack
        pix=stk.back().first;
        o=stk.back().second;
        stk.pop_back();
        pv = base[o].pix2vec(pix);

        int loc=comp.size()-1;
        tsize zone=getZone(loc,0,3);
        check_pixel (zone, pixset);
        planck_assert(loc==-1,"stack not used up");
        }
      return pixset;
      }
  };

double isLeft (const vec3 &a, const vec3 &b, const vec3 &c)
  { return dotprod(crossprod(a,b),c); }

// adapted from code available at http://geomalgorithms.com/a12-_hull-3.html
// Original copyright notice follows:
// Copyright 2001 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
vector<int> getHull (const vector<vec3> &vert, const vector<int> &P)
  {
  // initialize a deque D[] from bottom to top so that the
  // 1st three vertices of P[] are a ccw triangle
  int n = P.size();
  arr<int> D(2*n+1);
  int bot = n-2, top = bot+3;    // initial bottom and top deque indices
  D[bot] = D[top] = P[2];      // 3rd vertex is at both bot and top
  if (isLeft(vert[P[0]], vert[P[1]], vert[P[2]]) > 0)
    {
    D[bot+1] = P[0];
    D[bot+2] = P[1];           // ccw vertices are: 2,0,1,2
    }
  else
    {
    D[bot+1] = P[1];
    D[bot+2] = P[0];           // ccw vertices are: 2,1,0,2
    }

  // compute the hull on the deque D[]
  for (int i=3; i<n; i++)
    {   // process the rest of vertices
    // test if next vertex is inside the deque hull
    if ((isLeft(vert[D[bot]], vert[D[bot+1]], vert[P[i]]) > 0) &&
        (isLeft(vert[D[top-1]], vert[D[top]], vert[P[i]]) > 0) )
      continue;         // skip an interior vertex

    // incrementally add an exterior vertex to the deque hull
    // get the rightmost tangent at the deque bot
    while (isLeft(vert[D[bot]], vert[D[bot+1]], vert[P[i]]) <= 0)
      ++bot;                 // remove bot of deque
    D[--bot] = P[i];         // insert P[i] at bot of deque

    // get the leftmost tangent at the deque top
    while (isLeft(vert[D[top-1]], vert[D[top]], vert[P[i]]) <= 0)
      --top;                 // pop top of deque
    D[++top] = P[i];         // push P[i] onto top of deque
    }

  // transcribe deque D[] to the output hull array H[]
  int nout = top-bot;
  vector<int> res(nout);
  for (int h=0; h<nout; h++)
    res[h] = D[bot + h +1];

  return res;
  }

void prepPolyHelper (const vector<vec3> &vv, const vector<int> &P,
  vector<MocQueryComponent> &comp, bool doLast)
  {
  vector<int> hull=getHull(vv,P);
  vector<bool> addHull(hull.size());

  // sync both sequences at the first point of the convex hull
  int ihull=0, ipoly=0, nhull=hull.size(), npoly=P.size();
  while (hull[ihull]!=P[ipoly]) ++ipoly;

  // iterate over the pockets between the polygon and its convex hull
  int npockets=0;
  if (P.size()==3) // really P.size(), not vv.size()?
    for (int i=0; i<3; i++) addHull[i]=true;
  else
    {
    do
      {
      int ihull_next = (ihull+1)%nhull,
          ipoly_next = (ipoly+1)%npoly;
      if (hull[ihull_next]==P[ipoly_next]) // no pocket found
        { addHull[ihull]=true; ihull=ihull_next; ipoly=ipoly_next; }
      else // query pocket
        {
        int nvpocket=2; // number of vertices for this pocket
        while (P[ipoly_next]!=hull[ihull_next])
          {
          ipoly_next = (ipoly_next+1)%npoly;
          ++nvpocket;
          }
        vector<int> ppocket(nvpocket);
        int idx=0;
        int ipoly_bw=ipoly_next;
        while (P[ipoly_bw]!=hull[ihull])
          {
          ppocket[idx++]=P[ipoly_bw];
          ipoly_bw=(ipoly_bw+npoly-1)%npoly;
          }
        ppocket[idx]=hull[ihull];
        // process pocket recursively
        ++npockets;
        prepPolyHelper (vv, ppocket, comp, false);
        ihull=ihull_next;
        ipoly=ipoly_next;
        }
      } while (ihull!=0);
    }
  if (npockets>1) 
    comp.push_back(MocQueryComponent(OR,npockets));
  if (npockets>0) 
    comp.push_back(MocQueryComponent(NOT));

  if (!doLast)
    addHull.back()=false;

  // add convex hull
  for (tsize i=0; i<hull.size(); ++i)
    if (addHull[i])
      comp.push_back(MocQueryComponent
        (crossprod(vv[hull[i]],vv[hull[(i+1)%hull.size()]]).Norm(),0.5*pi));

  int num_and = 0;
  for (tsize i=0; i<hull.size(); ++i)
    if (addHull[i]) ++num_and;
  if (npockets>0) ++num_and;
  if (num_and>1) 
    comp.push_back(MocQueryComponent(AND,num_and));
  }

} // unnamed namespace

template<typename I> Moc<I> mocQuery (int order,
  const vector<MocQueryComponent> &comp)
  { return queryHelper<I>(order,order,false,comp).result(); }
template Moc<int> mocQuery (int order, const vector<MocQueryComponent> &comp);
template Moc<int64> mocQuery (int order, const vector<MocQueryComponent> &comp);

template<typename I> Moc<I> mocQueryInclusive (int order, int omax,
  const vector<MocQueryComponent> &comp)
  { return queryHelper<I>(order,omax,true,comp).result(); }
template Moc<int> mocQueryInclusive (int order, int omax,
  const vector<MocQueryComponent> &comp);
template Moc<int64> mocQueryInclusive (int order, int omax,
  const vector<MocQueryComponent> &comp);

vector<MocQueryComponent> prepPolygon (const vector<vec3> &vertex)
  {
  planck_assert(vertex.size()>=3,"not enough vertices in polygon");
  vector<vec3> vv(vertex.size());
  for (tsize i=0; i<vertex.size(); ++i)
    vv[i]=vertex[i].Norm();

  vector<int> P(vv.size());
  for (tsize i=0; i<P.size(); ++i)
    P[i]=i;
  vector<MocQueryComponent> comp;
  prepPolyHelper(vv,P,comp,true);
  return comp;
  }
