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
 *  Copyright (C) 2003-2014 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "xcomplex.h"
#include "paramfile.h"
#include "healpix_data_io.h"
#include "alm.h"
#include "alm_fitsio.h"
#include "healpix_map.h"
#include "healpix_map_fitsio.h"
#include "alm_healpix_tools.h"
#include "alm_powspec_tools.h"
#include "fitshandle.h"
#include "levels_facilities.h"
#include "lsconstants.h"
#include "announce.h"
#include "planck_rng.h"

using namespace std;

namespace {

template<typename T> void alm2map_cxx (paramfile &params)
  {
  int nlmax = params.template find<int>("nlmax");
  int nmmax = params.template find<int>("nmmax",nlmax);
  planck_assert(nmmax<=nlmax,"nmmax must not be larger than nlmax");
  string infile = params.template find<string>("infile");
  string outfile = params.template find<string>("outfile");
  int nside = params.template find<int>("nside");
  double fwhm = arcmin2rad*params.template find<double>("fwhm_arcmin",0);
  int cw_lmin=-1, cw_lmax=-1;
  if (params.param_present("cw_lmin"))
    {
    cw_lmin = params.template find<int>("cw_lmin");
    cw_lmax = params.template find<int>("cw_lmax");
    }

  arr<double> temp, pol;
  get_pixwin (params,nlmax,temp,pol);

  bool deriv = params.template find<bool>("derivatives",false);
  if (deriv)
    {
    Alm<xcomplex<T> > alm;
    read_Alm_from_fits(infile,alm,nlmax,nmmax,2);
    if (fwhm>0) smoothWithGauss (alm, fwhm);
    if (cw_lmin>=0) applyCosineWindow(alm, cw_lmin, cw_lmax);
    Healpix_Map<T> map(nside,RING,SET_NSIDE),
                   mapdth(nside,RING,SET_NSIDE),
                   mapdph(nside,RING,SET_NSIDE);
    alm.ScaleL(temp);

    double offset = alm(0,0).real()/sqrt(fourpi);
    alm(0,0) = 0;
    alm2map_der1(alm,map,mapdth,mapdph);
    map.Add(T(offset));
    write_Healpix_map_to_fits (outfile,map,mapdth,mapdph,planckType<T>());
    return;
    }

  bool polarisation = params.template find<bool>("polarisation");
  bool do_regnoise = params.param_present("regnoiseT");
  if (!polarisation)
    {
    Alm<xcomplex<T> > alm;
    read_Alm_from_fits(infile,alm,nlmax,nmmax,2);
    if (fwhm>0) smoothWithGauss (alm, fwhm);
    if (cw_lmin>=0) applyCosineWindow(alm, cw_lmin, cw_lmax);
    Healpix_Map<T> map(nside,RING,SET_NSIDE);
    alm.ScaleL(temp);

    double offset = alm(0,0).real()/sqrt(fourpi);
    alm(0,0) = 0;
    alm2map(alm,map);
    map.Add(T(offset));
    if (do_regnoise)
      {
      planck_rng rng(params.template find<int>("rand_seed",42));
      double rms = params.template find<double>("regnoiseT");
      for (int i=0; i<map.Npix(); ++i)
        map[i] += rms*rng.rand_gauss();
      }
    write_Healpix_map_to_fits (outfile,map,planckType<T>());
    }
  else
    {
    Alm<xcomplex<T> > almT, almG, almC;
    read_Alm_from_fits(infile,almT,almG,almC,nlmax,nmmax,2);
    if (fwhm>0) smoothWithGauss (almT, almG, almC, fwhm);
    if (cw_lmin>=0) applyCosineWindow(almT, almG, almC, cw_lmin, cw_lmax);
    Healpix_Map<T> mapT(nside,RING,SET_NSIDE), mapQ(nside,RING,SET_NSIDE),
                   mapU(nside,RING,SET_NSIDE);
    almT.ScaleL(temp);
    almG.ScaleL(pol); almC.ScaleL(pol);

    double offset = almT(0,0).real()/sqrt(fourpi);
    almT(0,0) = 0;
    alm2map_pol(almT,almG,almC,mapT,mapQ,mapU);
    mapT.Add(T(offset));
    if (do_regnoise)
      {
      planck_rng rng(params.template find<int>("rand_seed",42));
      double rmsT  = params.template find<double>("regnoiseT"),
             rmsQU = params.template find<double>("regnoiseQU");
      for (int i=0; i<mapT.Npix(); ++i)
        {
        mapT[i] += rmsT *rng.rand_gauss();
        mapQ[i] += rmsQU*rng.rand_gauss();
        mapU[i] += rmsQU*rng.rand_gauss();
        }
      }
    write_Healpix_map_to_fits (outfile,mapT,mapQ,mapU,planckType<T>());
    }
  }

} // unnamed namespace

int alm2map_cxx_module (int argc, const char **argv)
  {
  module_startup ("alm2map_cxx", argc, argv);
  paramfile params (getParamsFromCmdline(argc,argv));

  bool dp = params.find<bool> ("double_precision",false);
  dp ? alm2map_cxx<double>(params) : alm2map_cxx<float>(params);
  return 0;
  }
