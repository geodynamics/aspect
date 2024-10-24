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
 *  Copyright (C) 2003, 2005, 2009 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_data_io.h"
#include "arr.h"
#include "fitshandle.h"
#include "paramfile.h"
#include "string_utils.h"

using namespace std;

namespace {

void read_wring (const string &weightfile, int nside, arr<double> &weight)
  {
  fitshandle inp;
  inp.open(weightfile);
  inp.goto_hdu(2);
  planck_assert(nside==inp.get_key<int>("NSIDE"),"incorrect Nside parameter");
  inp.read_entire_column(1,weight);
  planck_assert(weight.size()==size_t(2*nside),
    "incorrect number of weights in ring weight file");
  }

} // unnamed namespace

void read_weight_ring (const string &dir, int nside, arr<double> &weight)
  {
  read_wring(dir+"/weight_ring_n"+intToString(nside,5)+".fits", nside, weight);
  }

void get_ring_weights (paramfile &params, int nside, arr<double> &weight)
  {
  string weightfile = params.find<string>("ringweights","");
  weight.alloc (2*nside);
  if (weightfile!="")
    {
    read_wring (weightfile, nside, weight);
    for (tsize m=0; m<weight.size(); ++m) weight[m]+=1;
    }
  else
    weight.fill(1);
  }

vector<double> read_fullweights_from_fits(const std::string &weightfile,
  int nside)
  {
  fitshandle inp;
  inp.open(weightfile);
  inp.goto_hdu(2);
  planck_assert(inp.colname(1)=="COMPRESSED PIXEL WEIGHTS","wrong column name");
  planck_assert(inp.get_key<int>("NSIDE")==nside,"incorrect NSIDE parameter");
  vector<double> res;
  inp.read_entire_column(1,res);
  return res;
  }

void read_pixwin (const string &file, arr<double> &temp)
  {
  fitshandle inp;
  inp.open(file);
  inp.goto_hdu(2);
  if (temp.size()==0)
    inp.read_entire_column(1,temp);
  else
    inp.read_column(1,temp);
  }
void read_pixwin (const string &file, arr<double> &temp, arr<double> &pol)
  {
  fitshandle inp;
  inp.open(file);
  inp.goto_hdu(2);
  if (temp.size()==0)
    inp.read_entire_column(1,temp);
  else
    inp.read_column(1,temp);
  if (pol.size()==0)
    inp.read_entire_column(2,pol);
  else
    inp.read_column(2,pol);
  }

void get_pixwin (paramfile &params, int lmax, arr<double> &pixwin)
  {
  string windowfile = params.find<string>("windowfile","");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  if (windowfile!="")
    read_pixwin (windowfile,pixwin);
  }
void get_pixwin (paramfile &params, int lmax, arr<double> &pixwin,
  arr<double> &pixwin_pol)
  {
  string windowfile = params.find<string>("windowfile","");
  pixwin.alloc(lmax+1);
  pixwin.fill(1);
  pixwin_pol.alloc(lmax+1);
  pixwin_pol.fill(1);
  if (windowfile!="")
    read_pixwin (windowfile,pixwin,pixwin_pol);
  }
