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

/*! \file alice3.cc
 *  Copyright (C) 2005-2015 David Larson, Max-Planck-Society
 *  \author David Larson \author Martin Reinecke
 */

#include <iostream>
#include "announce.h"
#include "paramfile.h"
#include "healpix_map_fitsio.h"
#include "lsconstants.h"
#include "arr.h"
#include "fitshandle.h"
#include "vec3.h"
#include "string_utils.h"
#include "alm.h"
#include "alm_healpix_tools.h"
#include "planck_rng.h"
#include "healpix_map.h"

using namespace std;

class PolarizationHolder
  {
  public:
    Healpix_Map<double> Q, U;

    void getQU(const pointing &p, double &q, double &u) const
      {
      fix_arr<int,4> pix;
      fix_arr<double,4> wgt;
      Q.get_interpol(p,pix,wgt);
      q=u=0.f;
      for (tsize i=0;i<4;++i) {q+=Q[pix[i]]*wgt[i];u+=U[pix[i]]*wgt[i]; }
      }

    vec3 getQUDir(const vec3 &loc) const
      {
      double q,u;
      getQU(loc,q,u);
      vec3 east(1,0,0);
      if (abs(loc.x)+abs(loc.y) > 0.0)
        east = vec3(-loc.y,loc.x,0).Norm();
      vec3 north = crossprod(loc, east);
      double angle = 0.5*safe_atan2(u,q);
      return north*(-cos(angle)) + east*sin(angle);
      }

    // Return the magnitude of the polarization at some pointing.
    double getQUMagnitude(const pointing& p) const
      {
      double q,u;
      getQU(p,q,u);
      return sqrt(q*q + u*u);
      }
  };

/*! Steps from loc in direction dir for an angle theta and updates loc and dir. */
void get_step(const PolarizationHolder &ph, vec3 &loc, vec3 &dir, double theta)
  {
  loc=(loc+dir*theta).Norm();
  vec3 tdir=ph.getQUDir(loc);
  dir = (dotprod(dir,tdir)<0) ? -tdir : tdir;
  }

/*! Performs one Runge-Kutta second order step. Updates loc and dir. */
void runge_kutta_step(vec3 &loc, vec3 &dir, const PolarizationHolder &ph,
  double theta)
  {
  // Take a half-theta step
  vec3 tloc=loc;
  get_step(ph, tloc, dir, theta/2.0);

  // Then take a full step with the new direction
  get_step(ph, loc, dir, theta);
  }

/*! Second order Runge-Kutta integration on the sphere.  Given a
  starting location, a qu map of the sky, and a step size theta, this
  subroutine returns an array of vectors extending in both
  directions from the starting location.  */
void runge_kutta_2(const vec3 &location, const PolarizationHolder &ph,
  double theta, arr<vec3> &locs)
  {
  vec3 first_dir=ph.getQUDir(location);
  vec3 dir = first_dir;
  vec3 loc = location;

  locs[locs.size()/2] = loc;

  for(int i = 1 + locs.size()/2; i<int(locs.size()); i++)
    {
    runge_kutta_step(loc, dir, ph, theta);
    locs[i] = loc;
    }

  dir = -first_dir;
  loc = location;
  for(int i = -1 + locs.size()/2; i>=0; i--)
    {
    runge_kutta_step(loc, dir, ph, theta);
    locs[i] = loc;
    }
  }

/*! Create a sinusoidal kernel. */
void make_kernel(arr<double> &kernel)
  {
  for(tsize i=0; i<kernel.size(); i++)
    {
    double sinx = sin(pi*(i+1) / (kernel.size()+1));
    kernel[i] = sinx*sinx;
    }
  }

/*! Convolve an array with a kernel. */
void convolve(const arr<double> &kernel, const arr<double> &raw, arr<double> &convolution)
  {
  convolution.alloc(raw.size()-kernel.size()+1);
  for(tsize i=0; i<convolution.size(); i++)
    {
    double total=0;
    for (tsize j=0; j<kernel.size(); j++)
      total += kernel[j] * raw[i+j];
    convolution[i] = total;
    }
  }

/*! Perform line integral convolution on sphere. */
int lic_function(Healpix_Map<double> &hitcount, Healpix_Map<double> &texture,
  const PolarizationHolder &ph, const Healpix_Map<double> &th, int steps,
  int kernel_steps, double step_radian)
  {
  arr<double> kernel(kernel_steps), convolution, rawtexture;
  make_kernel(kernel);
  arr<vec3> curve(steps);

  texture.fill(0.);
  int num_curves=0;

  for(int i=0; i<texture.Npix(); i++)
    {
    if (hitcount[i]<1.0)
      {
      num_curves++;
      runge_kutta_2(texture.pix2vec(i), ph, step_radian, curve);
      rawtexture.alloc(curve.size());
      for (tsize i2=0; i2<curve.size(); i2++)
        rawtexture[i2] = th.interpolated_value(curve[i2]);
      convolve(kernel, rawtexture, convolution);
      for (tsize j=0; j<convolution.size(); j++)
        {
        int k = texture.vec2pix(curve[j+kernel.size()/2]);
        texture[k] += convolution[j];
        hitcount[k] += 1.;
        }
      }
    }
  return num_curves;
  }

/*! Expose line integral convolution for external use. */
void lic_main(const Healpix_Map<double> &Q, const Healpix_Map<double> &U, const Healpix_Map<double> &th,
  Healpix_Map<double> &hit, Healpix_Map<double> &tex, Healpix_Map<double> &mag,
  int steps, int kernel_steps, double step_radian, double polmin, double polmax)
  {
  PolarizationHolder ph;
  ph.Q = Q;
  ph.U = U;

  hit.fill(0.);

  for (int i=0; i<mag.Npix(); i++)
    {
    pointing p = mag.pix2ang(i);

    mag[i] = min(polmax,max(polmin,ph.getQUMagnitude(p)));
    tex[i] = th.interpolated_value(p);
    }

  lic_function(hit, tex, ph, th, steps, kernel_steps, step_radian);

  for (int i=0; i<tex.Npix(); ++i)
    tex[i]/=hit[i];
  double tmin,tmax,mmin,mmax;
  tex.minmax(tmin,tmax);
  mag.minmax(mmin,mmax);
  for (int i=0; i<tex.Npix(); ++i)
    {
    mag[i]*=(tex[i]-tmin);
    tex[i]=1.0-(tex[i]-tmin)/(tmax-tmin);
    }
  mag.minmax(mmin,mmax);
  for (int i=0; i<mag.Npix(); ++i)
    mag[i]=1.0-(mag[i]-mmin)/(mmax-mmin);
  }

int main(int argc, const char** argv)
  {
  module_startup ("alice3", argc, argv);
  paramfile params (getParamsFromCmdline(argc,argv));

  Healpix_Map<double> Q, U;
  // Load a polarized fits file, with Q and U as the second
  // and third columns (the standard form).
  read_Healpix_map_from_fits(params.find<string>("in"), Q, 2, 2);
  read_Healpix_map_from_fits(params.find<string>("in"), U, 3, 2);

  int nside = params.find<int>("nside");
  int steps = params.find<int>("steps",100);
  double step_radian=arcmin2rad*params.find<double>("step_arcmin",10.);
  int kernel_steps = params.find<int>("kernel_steps",50);
  double polmin = params.find<double>("polmin",-1e30),
        polmax = params.find<double>("polmax",1e30);
  string out = params.find<string>("out");

  Healpix_Map<double> th(nside,RING,SET_NSIDE);
  if (params.param_present("texture"))
    read_Healpix_map_from_fits(params.find<string>("texture"),th);
  else
    {
    planck_rng rng(params.find<int>("rand_seed",42));
    if (params.param_present("ell"))
      {
      int ell = params.find<int>("ell");
      if (ell<0) ell=2*nside;
      Alm<xcomplex<double> > a(ell,ell);
      a.SetToZero();
      cout << "Background texture using ell = " << ell << endl;

      a(ell,0)=fcomplex(rng.rand_gauss(),0.);
      for (int m=0; m<=ell; m++)
        { a(ell,m).real(rng.rand_gauss()); a(ell,m).imag(rng.rand_gauss()); }
      alm2map(a, th);
      }
    else
      {
      for (int i=0; i<th.Npix(); i++)
        th[i] = rng.rand_uni() - 0.5;
      }
    }

  Healpix_Map<double> hit(nside,RING,SET_NSIDE),
                     tex(nside,RING,SET_NSIDE),
                     mag(nside,RING,SET_NSIDE);

  {
  PolarizationHolder ph;
  ph.Q = Q;
  ph.U = U;

  for (int i=0; i<mag.Npix(); i++)
    tex[i] = th.interpolated_value(mag.pix2ang(i));
  write_Healpix_map_to_fits(out+"_background.fits",tex,PLANCK_FLOAT32);
  }

  lic_main(Q, U, th, hit, tex, mag, steps, kernel_steps, step_radian, polmin, polmax);

  write_Healpix_map_to_fits(out+"_texture.fits",tex,PLANCK_FLOAT32);
  write_Healpix_map_to_fits(out+"_mod_texture.fits",mag,PLANCK_FLOAT32);
  }
