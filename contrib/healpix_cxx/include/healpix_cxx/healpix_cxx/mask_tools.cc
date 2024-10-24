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

/*! \file mask_tools.cxx
 *  Copyright (C) 2003-2019 Max-Planck-Society
 *  \author Martin Reinecke
 */

#include <vector>
#include <functional>
#include "mask_tools.h"
#include "lsconstants.h"

using namespace std;

namespace {

double pseudodist(const vec3 &v1, const vec3 &v2)
  { return (v1-v2).SquaredLength(); }
double dist2pseudodist(double pdist)
  { double tmp=2*sin(0.5*pdist); return tmp*tmp; }
double pseudodist2dist(double dist)
  { return 2*asin(sqrt(dist)*0.5); }

} //unnamed namespace

Healpix_Map<double> dist2holes(const Healpix_Map<double> &mask, double maxdist)
  {
  constexpr int FULLY_IN_MASK=1;
  constexpr int BORDER=2;
  int maxord = mask.Order();
  planck_assert(maxord>=0, "Nside must be a power of 2");
  vector<int> orders;
  if (maxord<6)
    orders={maxord};
  else
    orders={(maxord+1)/2, maxord};
  vector<Healpix_Map<uint8>> omask; // masks at requested orders
  vector<double> buf; // safety distance at requested orders
  for (size_t i=0; i<orders.size(); ++i)
    {
    omask.emplace_back(orders[i], NEST);
    buf.push_back(2*omask[i].max_pixrad());
    }
  auto &maxmask(omask.back());
  if (mask.Scheme()==RING)
#pragma omp parallel for schedule(static)
    for (int i=0; i<mask.Npix(); ++i)
      maxmask[i] = (mask[mask.nest2ring(i)]==0) ? FULLY_IN_MASK : 0;
  else
#pragma omp parallel for schedule(static)
    for (int i=0; i<mask.Npix(); ++i)
      maxmask[i] = (mask[i]==0) ? FULLY_IN_MASK : 0;
  // find border pixels
#pragma omp parallel for schedule(dynamic,10000)
  for (int i=0; i<mask.Npix(); ++i)
    if (maxmask[i])
      {
      fix_arr<int, 8> nb;
      maxmask.neighbors(i, nb);
      for (size_t j=0; j<8; ++j)
        if ((nb[j]!=-1) && (maxmask[nb[j]]==0))
          {
          maxmask[i] |= BORDER;
          break;
          }
      }

  // compute coarsened masks
  for (int j=omask.size()-2; j>=0; --j)
    {
    int fct = 1<<(2*(omask[j+1].Order()-omask[j].Order()));
#pragma omp parallel for schedule(static)
    for (int i=0; i<omask[j].Npix(); ++i)
      {
      uint8 And=255, Or=0;
      for (int k=0; k<fct; ++k)
        {
        And &= omask[j+1][fct*i+k];
        Or  |= omask[j+1][fct*i+k];
        }
      omask[j][i] = (And&FULLY_IN_MASK) | (Or&BORDER);
      }
    }

  Healpix_Map<double> tdist(maxord, NEST);
#pragma omp parallel for schedule(static)
  for (int i=0; i<mask.Npix(); ++i)
    tdist[i] = (maxmask[i]>0) ? 0. : maxdist;

  function<void(int, int, const vector<int> &, const vector<vec3> &)> process;
  process = [&omask,&tdist,maxord,&process,&buf,maxdist]
    (int oo, int pix, const vector<int> &submask, const vector<vec3> &subvec)->void
    {
    int order=omask[oo].Order();
    if (submask.empty()) return;
    if (omask[oo][pix]&FULLY_IN_MASK) return; // pixel lies fully in the mask
    vec3 v=omask[oo].pix2vec(pix);
    if (order==maxord)
      {
      double pd=10;
      for (auto subv:subvec)
        pd=min(pd, pseudodist(subv, v));
      tdist[pix]=min(maxdist,pseudodist2dist(pd));
      return;
      }
    vector<double> psubdist(submask.size());
    double pdlim=10;
    for (size_t i=0; i<submask.size(); ++i)
      {
      psubdist[i] = pseudodist(v, subvec[i]);
      pdlim = min(pdlim,psubdist[i]);
      }

    double lim0=dist2pseudodist(min(pi,maxdist+buf[oo]));
    if (pdlim>lim0) return;
    double lim1=dist2pseudodist(min(pi,pseudodist2dist(pdlim)+2*buf[oo]));
    // not yet full order;
    vector<int> submask2;
    vector<vec3> subvec2;
    int fct = 1<<(2*(omask[oo+1].Order()-omask[oo].Order()));
    for (size_t i=0; i<submask.size(); ++i)
      if (psubdist[i]<lim1)
        if (psubdist[i]<lim0)
          for (int j=submask[i]*fct; j<(submask[i]+1)*fct; ++j)
            if (omask[oo+1][j]&BORDER)
              {
              submask2.push_back(j);
              subvec2.push_back(omask[oo+1].pix2vec(j));
              }
    for (int pix2=pix*fct; pix2<(pix+1)*fct; ++pix2)
      process(oo+1, pix2, submask2, subvec2);
    };

  vector<int> submask;
  vector<vec3> subvec;
  for (int i=0; i<omask[0].Npix(); ++i)
    if (omask[0][i]&BORDER)
      {
      submask.push_back(i);
      subvec.push_back(omask[0].pix2vec(i));
      }
#pragma omp parallel for schedule(dynamic)
  for (int i=0; i<omask[0].Npix(); ++i)
    process(0, i, submask, subvec);

  if (mask.Scheme()==RING) tdist.swap_scheme();
  return tdist;
  }
