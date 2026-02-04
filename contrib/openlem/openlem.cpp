// OpenLEM Version 44 experimental 2025-07-10
//
// Copyright (C) 2012-2025 Stefan Hergarten
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details. To view a copy of this license, visit
// https://www.gnu.org/licenses/gpl-3.0.txt or write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// For any questions concerning the codes and bug reports, please contact the
// developer at stefan.hergarten@geologie.uni-freiburg.de.
//
// Changelog
//
// Version 44 2025-07-10
// Reimplemented the method fillLakes for better efficiency.
// Fixed problem with flow directions if boundaries are changed during the
// simulation.
// Changed status of flooded areas (water level > topography) from hillslope
// to channel for option CHANNEL.
//
// Version 43 2023-08-22
// Introduced faster computation for neighbors.
// Changed some default parameter values.
// Renamed method computeDischarge to computeFluxes.
// Removed option INTOLAKES.
// Introduced new concept for computing the flow pattern in lakes (option
// LAKES without LAKEFLOWDIR, highly experimental).
// Introduced option REDREC for reducing the level of recursion (highly
// experimental)
//
// Version 42 2023-06-29
// Fixed problem with stream power in lakes (option LAKES).
// Introduced option BORROW for treating the situation that the erosion
// of deposits exceeds the available thickness of the deposit layer
// (in combination with the option LAYERS). This option removes the ability
// to keep track of the time at which the deposit layer was removed.
// Introduced option NOCHANNELDIFFUS to switch off diffusion at channelized
// sites.
// Extended exponent uexp for the U-shape of glacial valleys towards elliptical
// shapes (uexp < 0).
// Fixed problem with combination of options LAKES and DEPOS.
//
// Version 41 2023-02-21
// Introduced stream-power exponent spexph for hillslopes used in combination
// of the options CHANNEL, SHAREDSP, and SPEXP.
// Introduced exponent uexp for the U-shape of glacial valleys.
// Reorganized default definitions to fix a problem with the combination
// CHANNEL and SHAREDSP.
//
// Version 40 2022-07-06
// Fixed problem in constructor Grid() in combination with option PRECIP.
// Extended methods for reading so that files written with a different value
// of LAYERS can be read.
//
// Version 39 2022-06-15
// Changed member spfac to public in order to allow custom stream-power laws.
// Changed sites that contain ice to channelized in order to achive a
// compatibility of the options ICE and CHANNEL.
// Introduced erodibilities kdh and kth for hillslopes used in combination
// of the options CHANNEL and SHAREDSP.
// Removed option RESCALE.
//
// Version 38 2022-04-15
// Introduced option CHANNEL for the distinction between channel and hillslope
// sites by the flow pattern.
//
// Version 37 2022-01-18
// Modified implementation of diffusion in such a way that the option
// DIFFSEMIIMPL without DIFFUS and EIGHTNEIGHBORS implements diffusion
// in D8 flow direction.
// Added option DEPOS for switching to fully transport-limited conditons
// in domains of sediment deposition.
// Added option MELTBOUNDARY for melting all ice at the boundaries.
// Fixed a problem in constructing a grid of size 1x1.
// Fixed a problem in reading unrecognized keys.
// Fixed a problem with option MELTRATE.
//
// Version 36 2021-09-08
// Removed decprecated options SEDIM, CONV and fixed-point iteration for
// method erode().
// Improved treatment of boundary points for the LFPM.
//
// Version 35 2021-08-19
// Included LFPM for orographic precipitation.
//
// Version 34 2021-08-03
// Added option DIFFUSIVITY for spatially variable diffusivity.
// Added option ICEFRAC for using own expressions for the elevation-dependent
// ice production.
// Added option MELTRATE for defining a melting rate instead of the default
// mixing model.
//
// Version 33 2021-03-29
// Added support for stream-power exponents (SPEXP) in all modes.
// Improved option LAKES for linear decline model.
// Added option LAKEFLOWDIR for improved flow direction in lakes.
// Removed option BOUNDARY.
// Fixed problem in reading array data.
//
// Version 32 2021-03-02
// Changed interpretation of variable bottom, relative to surface now.
// Added option for writing individual layers of array-valued variables.
// Fixed error in the definition of the keys for erodibilities.
// Added continued reading after finding an unrecognized key in a file.
//
// Version 31 2021-01-19
// Fixed problem with glacial erodibilities under option NOMELTWATER
// Added file handling for array-valued variables
// Added support for layers with different erodibilities
//
// Version 30 2021-01-08
// Fixed problems with unrecognized keys.
// Added initializer for variable excav.
//
// Version 29 2020-12-22
// Fixed problem with option ICE if glacial erosion is not detachment-limited.
// Changed mean of erodibilities for mixed ice/water conditions.
// Shared stream-power model is now default for option ICE.
// Added option NOMELTWATER.
//
// Version 28 2020-12-17
// Fixed problem with option SHAREDSP.
// Proceeded with work on option ICE.
//
// Version 27 2020-11-27
// Introduced some technical simplificiations that should not affect the results.
// Proceeded with work on option ICE.
//
// Version 26 2020-10-21
// Introduced updated methods for obtaining the neighbors.
// Simplified the computation for the transport-limited model and the linear
// decline model.
// Added option LAKES for reducing stream power in lakes.
// Improved information about unrecognized keys when reading files.
//
// Version 25 2020-10-12
// Fixed problem with lowercase letters in keys.
//
// Version 24 2020-10-06
// Fixed problem with resizing.
// Introduced option EIGHTNEIGHBORS for a 9-point diffusion stencil.
// Introduced option DIFFSEMIIMPL for semi-implicit diffusion scheme.
// Renamed option SEDIMIMPL TO LINDECL and added option SHAREDSP.
// Replaced option REMOVE by a parameter excav.
// Proceeded with work on option ICE.
//
// Version 23 2020-09-06
// Added partly implicit scheme for hillslope diffusion for transport-limited
// erosion. Cannot be combined with option REMOVE.
// Removed method diffuse().
// Renamed variable ac to athr (fluvial threshold).
// Introduced new methods for reading and writing data.
// Reduced writing precision for the discharge with option PRECIP to 4 Bytes.
// Changed meaning of the third argument in method write().
// Added option ICE (experimental).
//
// Version 22 2020-07-29
// Removed test part in method erode as this altered the computed sediment
// fluxes
//
// Version 21 2020-06-03
// Introduced implicit scheme for the transport model of Davy and Lague (2009)
//
// Version 20 2020-04-14
// Introduced transport-limited erosion TRANSPORT.
// This option is not compatible with the option PRECIP and with a fluvial
// threshold ac > 0.
// Introduced method setErosionFactors(), has to be called explicitly if
// the vales A^theta shall be preset for higher performance.
// Method computeDischarge() is now called automatically by the method erode().
// Changed definition of class to class template. Use Grid<> instead of Grid.
//
// Version 19 2020-03-09
// Introduced threshold catchment size ac
// Added new functions for setting parameters
//
// Version 18 2020-02-13
// Introduced option PRECIP
//
// Version 17 2020-02-03
// Changed definition of alpha according to Turcotte and Schubert

#include <stdlib.h>
#include <stdio.h>
#include <cfloat>
#include <math.h>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>

namespace openlem
{
#define SHAREDSP
#define LAYERS 2
#define DEFORM
#ifdef ICE
#define SHAREDSP
#define BUFFERS
#endif

#ifdef DEPOS
#define SHAREDSP
#define BUFFERS
#endif

#ifdef SHAREDSP
#define LINDECL
#ifdef CHANNEL
#define BUFFERS
#endif
#endif

#ifdef LINDECL
#define TRANSPORT
#ifdef LAKES
#define BUFFERS
#endif
#endif

#ifdef RIGIDITY
#define FLEXURE
#endif

#ifdef EXTPRECIP
#define PRECIP
#endif

using namespace std;

class Keymapentry
{
  public:
  char  keystring[4], varname[16];
  int   keyint, offset, size, fsize, grid, readonly, array;
};
  
class Point
{
  public:
  short int  i, j;

  Point ( int i = 0, int j = 0 )
  {
    this->i = i;
    this->j = j;
  }

  int equals ( int i, int j )
  { 
    return i == this->i && j == this->j;
  }  

  int isDiag ( Point p )
  { 
    return i != p.i && j != p.j;
  }  

  int isDiag ( int i, int j )
  { 
    return i != this->i && j != this->j;
  }  

  int operator== ( Point p )
  {
    return i == p.i && j == p.j;
  }  

  int operator!= ( Point p )
  {
    return i != p.i || j != p.j;
  }  
public:
  bool operator== (const Point p) const
  {
    return i == p.i && j == p.j;
  }

};

class Node
{
  public:	
  double        h,   // surface elevation
                l;   // water level, >= h
  float         u;   // uplift rate
#ifdef TRANSPORT		
  double        qs;  // sediment flux 
  double        qp;  // derivative of sediment flux
#endif		     
#ifdef ICE		
  double        qi;  // ice flux
  double        th;  // ice thickness
  Point         trunk; // donor with the largest catchment size
#endif		     
#ifdef LAKES 
#ifndef ICE
  Point         trunk; // additional donor in lakes 
#endif		     
#endif		     
#ifdef ERODIBILITY
  float	        e;   // erodibility
#endif		     
#ifdef DIFFUS		
  float	        r;   // diffusive erosion rate, does not include fluxes
                     // in flow direction for transport-limited erosion
#endif
#ifdef DIFFUSIVITY		
  float	        diff; // diffusivity
#endif
#ifdef LAYERS		
  float	        bottom[LAYERS-1]; // lower boundary of the respective layer
                                  // relative to the surface
#endif
#ifdef FLEXURE		
  double        w;   // vertical deflection
  float         g;   // load
#endif
#ifdef RIGIDITY
  float         alpha; // flexural parameter (4D/(rhom g))^0.25
#endif
#ifdef PRECIP
  float         p;   // precipitation rate
  double        q;   // discharge
#ifdef EXTPRECIP
  float         qv;  // vapor flux
  float         qc;  // cloud water flux
  float         ptot; // total precipitation
#endif
#else  
#ifdef DEFORM
  double        q;   // discharge (catchment size)
#else
  unsigned int  q;   // discharge (catchment size)
#endif  
#endif  
#ifdef DEFORM
  double        x;   // x-position (row)
  double        y;   // y-position (column)
  double        cs;   // cell size
#endif
#ifdef CHANNEL
  char          channel; // channel marker (1 for channelized sites, 0 else)
#endif  
  char          b,   // boundary marker (1 for boundary nodes, 0 else)
                m;   // marker (only internal use)
  Point         d;   // flow target
#ifdef BUFFERS
  double        lambdaq;
#endif  
#ifdef ICE
  double        qeff;
#endif  

  int drainsTo ( Point p )
  {
    return p.i == d.i && p.j == d.j;
  }  

  int drainsTo ( int i, int j )
  {
    return i == d.i && j == d.j;
  }  

#ifdef LAYERS
  void adjustLayers ( double dh, double t = FLT_MAX )
  {
#ifdef BORROW
    if ( dh > 0 || bottom[0] < 0 )  bottom[0] -= dh;
    float borrowed = bottom[0] >= 0 ? bottom[0] : 0;    
    for ( int k = 1; k < LAYERS-1; ++k )
      if ( bottom[k] <= 0 && (bottom[k]-=dh) >= 0 && t > 0 )  bottom[k] = t;
#else
    if ( bottom[0] > 0 && dh > 0 )  bottom[0] = 0;
    for ( int k = 0; k < LAYERS-1; ++k )
      if ( bottom[k] <= 0 && (bottom[k]-=dh) >= 0 && t > 0 )  bottom[k] = t;
#endif
  }

  int getLayer()
  {
    for ( int i = 0; i < LAYERS-1; ++i )
      if ( bottom[i] < 0 )  return i;
    return LAYERS-1;
  }
#endif 

#ifdef DEFORM
  double distSquare ( Node *pn )
  {
    return  (x-pn->x)*(x-pn->x)+(y-pn->y)*(y-pn->y);
  }

  double distance ( Node *pn )
  {
    return  sqrt(distSquare(pn));
  }
#endif 
};

class Vect
{
  public:
  double  u1, u2;

  Vect ( double u1 = 0, double u2 = 0 )
  {
    this->u1 = u1;
    this->u2 = u2;
  }

  Vect operator -= ( Vect u )
  {
    u1 -= u.u1;
    u2 -= u.u2;
    return *this;
  }
};

class Matrix
{
  public:
  double  a11, a12, a21, a22;  

  Matrix ( double a11 = 0, double a12 = 0, double a21 = 0, double a22 = 0 )
  {
    this->a11 = a11;
    this->a12 = a12;
    this->a21 = a21;
    this->a22 = a22;
  } 

  Matrix operator -= ( Matrix a )
  {
    a11 -= a.a11;
    a12 -= a.a12;
    a21 -= a.a21;
    a22 -= a.a22;
    return *this;
  }

  Vect operator * ( Vect b )
  {
    return Vect ( a11*b.u1+a12*b.u2, a21*b.u1+a22*b.u2 );
  }

  Matrix operator * ( Matrix b )
  {
    return Matrix ( a11*b.a11+a12*b.a21, a11*b.a12+a12*b.a22, 
                    a21*b.a11+a22*b.a21, a21*b.a12+a22*b.a22 );
  }

  Matrix inv ()
  {
    double det = a11*a22-a12*a21;
    return  Matrix ( a22/det, -a12/det, -a21/det, a11/det );  
  }
};

template <typename valuetype>
class PointValue
{
  public:
  Point      p;
  valuetype  d;

  PointValue ( Point p = Point(0,0), valuetype d = 0 )
  {
    this->p = p;
    this->d = d;
  }
};

template <class Node = Node>
class Grid
{
  vector<vector<Node> >    u;         // lattice of nodes 
  Point                    shift[8];  // relative position of the 8 neighbors
  vector<int>              imap;      // mapping for coarsening  
  				      // (only internal use)
  vector<int>              jmap;      // mapping for coarsening  
  			              // (only internal use)
  unsigned int		   c;         // coarsening level for flexure
                                      // (only internal use)
  vector<vector<vector<double> > >  w, res, frhs; // for flexure
#ifdef RIGIDITY
  vector<vector<vector<double> > >  rig; 
#else
  vector<double>           rig;
#endif  
  double                   at;         // a^theta (only internal use)

  public:
  short int                m, n;      // lattice size m times n
  double  t,         // time
          theta,     // concavity index
          spexp,     // exponent in the stream power law
          spexph,    // exponent in the stream power law at hillslopes
	  lambda,    // sediment deposition parameter in the linear decline model
#ifdef LAYERS
          kd[LAYERS],   // erodibility for detachment-limited conditions
          kt[LAYERS],   // erodibility for transport-limited conditions
          kdg[LAYERS],  // erodibility for detachment-limited conditions
                        // for glacial erosion
          ktg[LAYERS],  // erodibility for transport-limited conditions
                        // for glacial erosion
          kdh[LAYERS],  // erodibility for detachment-limited conditions
                        // at hillslopes
          kth[LAYERS],  // erodibility for transport-limited conditions
                        // at hillslopes
#else
          kd,        // erodibility for detachment-limited conditions
          kt,        // erodibility for transport-limited conditions
          kdg,       // erodibility for detachment-limited conditions
                     // for glacial erosion
          ktg,       // erodibility for transport-limited conditions
                     // for glacial erosion
          kdh,       // erodibility for detachment-limited conditions
                     // at hillslopes
          kth,       // erodibility for transport-limited conditions
                     // at hillslopes
#endif
          diff,      // diffusivity
          excav,     // factor defining which fraction of the material
                     // moving from the hillslopes into rivers is
                     // excavated for detachment-limited erosion
// Option LAKES
          lff,       // fraction of the actual lake depth attempted to be
                     // filled with sediments
          mld,       // minimum lake depth
// Option ICE
          ela,       // glacial equilibrium line altitude
          fsa,       // altitude where all precipitation is snow
          gwf, gwe,  // parameters in the relation
                     // width = gwf*qi^gwe
          ghf, ghe,  // not used
          ghw,       // thickness to width ratio of glaciers
          uexp,      // exponent for the U-shape of glacial valley profiles
          ge,        // not used
          eps,	     // weight factor in the geometric mean q^eps*qi^(1-eps)
	             // used in the stream-power term for qi > q
// Option FLEXURE
          alpha,     // flexural parameter (4*D/(rhom g))^0.25
          tau,       // flexural relaxation time
          rhoc,      // density ratio crust/mantle

// Option PRECIP
          lc,        // length scale of condensation
          lf,        // length scale of fallout
          ll,        // length scale of long-range transport
          ld,        // length scale of dispersion
          refheight, // reference elevation  
          refslope,  // reference slope for evaporation
          evap,      // evaporation fraction at sea level
          qin,       // influx
// Option CHANNEL
          ah;        // apparent catchment size for non-channel nodes
  unsigned int  athr,  // threshold catchment size where fluvial erosion starts
                a,     // shift in catchment size used in the stream-power law
                dir;   // edge used for influx
  vector<float>            spfac;     // factor in stream power law
  map<int,Keymapentry>     km;
  vector<Matrix>  diag, upper, lower, right, bottom;
  vector<Vect>    rhs;
  vector<double>  psi;
// Option REDREC
  vector<Point>  seq; 
  int  nrec;         
  vector<int>  prevrow, nextrow, prevcol, nextcol;
  class Lake;
  vector<Lake>             lakes;     // all detected lakes
//  double  qdev;
// Option LAKES
  vector<PointValue<float>>  offset;

#ifdef DEPOS
  int etrue, efalse, dtrue, dfalse, oppos;
  int etruesum, efalsesum, dtruesum, dfalsesum, oppossum;
#endif

  Grid ( int m = 1, int n = 1, double theta = 0.5, unsigned int a = 0,   
         double spexp = 1., double lambda = 1., double diff = 0.,
         double alpha = 500., double rhoc = 0.85, double tau = 0. )
  {
    this->m = this->n = 0;
    t = 0.;
    setErosionParams(theta,spexp,0,a);   
    setSharedStreamPowerParams();   
    this->lambda = lambda;
    this->diff = diff;
    excav = 0.;
    lff = 1;
    mld = 0;
// Default parameters for orographic precipitation
// useful for 100 m grid spacing
    lc = lf = ld = 250;
    ll = 10000;
    refheight = 10;
    evap = dir = 0;
    qin = 5000;
// Default parameters for glacial erosion    
    ela = 10;
    fsa = 11;
    gwf = 1.5;
    gwe = 0.3;
    ghw = 0.1;
    uexp = 0;
    eps = 0.25;

    setFlexureParams(alpha,rhoc,tau);
    resize ( m, n );
    Point  p(0,0);
    int    k = 0;
    for ( p.i = -1; p.i <= 1; ++p.i )
      for ( p.j = -1; p.j <= 1; ++p.j )
	if ( p.i || p.j  )  shift[k++] = p;       
    setKeys();
  }    

  void addKey ( const char *key, const char *varname, void *p, unsigned int size, unsigned int fsize, unsigned int grid, int readonly = 0, int array = 1 )
  {
    Keymapentry e;  
    strncpy(e.keystring,key,4);
    if ( grid ) 
      for ( int i = 0; i < 4; ++i )  e.keystring[i] = toupper(e.keystring[i]);
    memcpy(&e.keyint,e.keystring,4);    
    strncpy(e.varname,varname,15);    
    e.offset = (char*)p - ( grid ? (char*)getNode(0,0) : (char*)&t ); 
    e.size = size;
    e.fsize = fsize;
    e.grid = grid;
    e.readonly = readonly;
    e.array = array;
    km[e.keyint] = e;
  }

  Keymapentry getKey ( unsigned int keyint )
  {
    return  km[keyint];
  }

  Keymapentry getKey ( const char *keystring )
  {
    unsigned int  keyint;
    strncpy((char*)&keyint,keystring,4);    
    return  km[keyint];
  }

  void setKeys()
  {
    addKey ( "t", "t", &t, sizeof(t), 8, 0 );
    addKey ( "c", "theta", &theta, sizeof(theta), 8, 0, 1 );
    addKey ( "thf", "theta", &theta, sizeof(theta), 8, 0 );
    addKey ( "e", "spexp", &spexp, sizeof(spexp), 8, 0, 1 );
    addKey ( "exf", "spexp", &spexp, sizeof(spexp), 8, 0 );
    addKey ( "exh", "spexph", &spexph, sizeof(spexph), 8, 0 );
    addKey ( "l", "lambda", &lambda, sizeof(lambda), 8, 0, 1 );
    addKey ( "sdc", "lambda", &lambda, sizeof(lambda), 8, 0 );
    addKey ( "d", "diff", &diff, sizeof(diff), 8, 0 );
    addKey ( "exc", "excav", &excav, sizeof(excav), 8, 0 );
    addKey ( "lff", "lff", &lff, sizeof(lff), 8, 0 );
    addKey ( "mld", "mld", &mld, sizeof(mld), 8, 0 );
#ifdef LAYERS
    addKey ( "kd", "kd", kd, sizeof(kd[0]), 8, 0, 0, LAYERS );
    addKey ( "kt", "kt", kt, sizeof(kt[0]), 8, 0, 0, LAYERS );
    addKey ( "kdg", "kdg", kdg, sizeof(kdg[0]), 8, 0, 0, LAYERS );
    addKey ( "ktg", "ktg", ktg, sizeof(ktg[0]), 8, 0, 0, LAYERS );
    addKey ( "kdh", "kdh", kdh, sizeof(kdh[0]), 8, 0, 0, LAYERS );
    addKey ( "kth", "kth", kth, sizeof(kth[0]), 8, 0, 0, LAYERS );
#else
    addKey ( "kd", "kd", &kd, sizeof(kd), 8, 0 );
    addKey ( "kt", "kt", &kt, sizeof(kt), 8, 0 );
    addKey ( "kdg", "kdg", &kdg, sizeof(kdg), 8, 0 );
    addKey ( "ktg", "ktg", &ktg, sizeof(ktg), 8, 0 );
    addKey ( "kdh", "kdh", &kdh, sizeof(kdh), 8, 0 );
    addKey ( "kth", "kth", &kth, sizeof(kth), 8, 0 );
#endif
    addKey ( "ela", "ela", &ela, sizeof(ela), 8, 0 );
    addKey ( "fsa", "fsa", &fsa, sizeof(fsa), 8, 0 );
    addKey ( "gwf", "gwf", &gwf, sizeof(gwf), 8, 0 );
    addKey ( "gwe", "gwe", &gwe, sizeof(gwe), 8, 0 );
    addKey ( "ghf", "ghf", &ghf, sizeof(ghf), 8, 0, 1 );
    addKey ( "ghe", "ghe", &ghe, sizeof(ghe), 8, 0, 1 );
    addKey ( "ghw", "ghw", &ghw, sizeof(ghw), 8, 0 );
    addKey ( "uex", "uexp", &uexp, sizeof(uexp), 8, 0 );
    addKey ( "ge", "ge", &ge, sizeof(ge), 8, 0, 1 );
    addKey ( "eps", "eps", &eps, sizeof(eps), 8, 0 );
    addKey ( "p", "alpha", &alpha, sizeof(alpha), 8, 0, 1 );
    addKey ( "alp", "alpha", &alpha, sizeof(alpha), 8, 0 );
    addKey ( "u", "tau", &tau, sizeof(tau), 8, 0, 1 );
    addKey ( "tau", "tau", &tau, sizeof(tau), 8, 0 );
    addKey ( "r", "rhoc", &rhoc, sizeof(rhoc), 8, 0, 1 );
    addKey ( "rho", "rhoc", &rhoc, sizeof(rhoc), 8, 0 );
    addKey ( "lc", "lc", &lc, sizeof(lc), 8, 0 );
    addKey ( "lf", "lf", &lf, sizeof(lf), 8, 0 );
    addKey ( "ll", "ll", &ll, sizeof(ll), 8, 0 );
    addKey ( "ld", "ld", &ld, sizeof(ld), 8, 0 );
    addKey ( "rh", "refheight", &refheight, sizeof(refheight), 8, 0 );
    addKey ( "rs", "refslope", &refslope, sizeof(refslope), 8, 0 );
    addKey ( "ev", "evap", &evap, sizeof(evap), 8, 0 );
    addKey ( "qin", "qin", &qin, sizeof(qin), 8, 0 );
    addKey ( "ah", "ah", &ah, sizeof(ah), 8, 0 );
    addKey ( "dir", "dir", &dir, sizeof(dir), 4, 0 );
    addKey ( "f", "athr", &athr, sizeof(athr), 4, 0, 1 );
    addKey ( "ath", "athr", &athr, sizeof(athr), 4, 0 );
    addKey ( "a", "a", &a, sizeof(a), 4, 0, 1 );
    addKey ( "acr", "a", &a, sizeof(a), 4, 0 );

    Node  *p = getNode(0,0);
    addKey ( "H", "h", &p->h, sizeof(p->h), 8, 1 );
    addKey ( "H8", "h", &p->h, sizeof(p->h), 8, 1 );
    addKey ( "H4", "h", &p->h, sizeof(p->h), 4, 1 );
    addKey ( "L", "l", &p->l, sizeof(p->l), 8, 1 );
    addKey ( "L8", "l", &p->l, sizeof(p->l), 8, 1 );
    addKey ( "L4", "l", &p->l, sizeof(p->l), 4, 1 );
    addKey ( "WL", "l", &p->l, sizeof(p->l), 8, 1, 1 );
    addKey ( "WL8", "l", &p->l, sizeof(p->l), 8, 1, 1 );
    addKey ( "WL4", "l", &p->l, sizeof(p->l), 4, 1, 1 );
    addKey ( "U", "u", &p->u, sizeof(p->u), 4, 1 );
#ifdef TRANSPORT		
    addKey ( "QS", "qs", &p->qs, sizeof(p->qs), 4, 1 );
#endif		     
#ifdef ICE		
    addKey ( "QI", "qi", &p->qi, sizeof(p->qi), 4, 1 );
    addKey ( "TH", "th", &p->th, sizeof(p->th), 4, 1 );
#endif		     
#ifdef ERODIBILITY
    addKey ( "E", "e", &p->e, sizeof(p->e), 4, 1 );
#endif		     
#ifdef DIFFUS	
    addKey ( "R", "r", &p->r, sizeof(p->r), 4, 1 );
#endif
#ifdef DIFFUSIVITY		
    addKey ( "DIF", "diff", &p->diff, sizeof(p->diff), 4, 1 );
#endif
#ifdef LAYERS		
    addKey ( "BOT", "bottom", p->bottom, sizeof(p->bottom[0]), 4, 1, 0, LAYERS-1 );
#endif
#ifdef FLEXURE		
    addKey ( "W", "w", &p->w, sizeof(p->w), 8, 1 );
    addKey ( "W8", "w", &p->w, sizeof(p->w), 8, 1 );
    addKey ( "W4", "w", &p->w, sizeof(p->w), 4, 1 );
    addKey ( "G", "g", &p->g, sizeof(p->g), 4, 1 );
    addKey ( "LOA", "g", &p->g, sizeof(p->g), 4, 1 );
#endif
#ifdef RIGIDITY
//    addKey ( "X", "alpha", &p->alpha, sizeof(p->alpha), 4, 1, 1 );
    addKey ( "ALP", "alpha", &p->alpha, sizeof(p->alpha), 4, 1 );
#endif
#ifdef PRECIP
    addKey ( "P", "p", &p->p, sizeof(p->p), 4, 1 );
#endif  
#ifdef EXTPRECIP
    addKey ( "QV", "qv", &p->qv, sizeof(p->qv), 4, 1 );
    addKey ( "QC", "qc", &p->qc, sizeof(p->qc), 4, 1 );
    addKey ( "PT", "ptot", &p->ptot, sizeof(p->ptot), 4, 1 );
#endif
    addKey ( "Q", "q", &p->q, sizeof(p->q), 4, 1 );
#ifdef DEFORM
    addKey ( "X", "x", &p->x, sizeof(p->x), 8, 1 );
    addKey ( "X8", "x", &p->x, sizeof(p->x), 8, 1 );
    addKey ( "X4", "x", &p->x, sizeof(p->x), 4, 1 );
    addKey ( "Y", "y", &p->y, sizeof(p->y), 8, 1 );
    addKey ( "Y8", "y", &p->y, sizeof(p->y), 8, 1 );
    addKey ( "Y4", "y", &p->y, sizeof(p->y), 4, 1 );
#endif
#ifdef CHANNEL
    addKey ( "CHA", "channel", &p->channel, sizeof(p->channel), 1, 1 );
#endif  
    addKey ( "B", "b", &p->b, sizeof(p->b), 1, 1 );
    addKey ( "D", "d", &p->d, sizeof(p->d), sizeof(p->d), 1 );
    addKey ( "S", "s", &p, 0, sizeof(float), 1 );
  }

  void setSharedStreamPowerParams ( int layer = 0, double kd = 2, double kt = 2,
                                                   double kdg = 10, double ktg = 1e100,
                                                   double kdh = 6, double kth = 6 )
  {
#ifdef LAYERS
    this->kd[layer] = kd;
    this->kt[layer] = kt;
    this->kdg[layer] = kdg;
    this->ktg[layer] = ktg;
    this->kdh[layer] = kdh;
    this->kth[layer] = kth;
#else
    this->kd = kd;
    this->kt = kt;
    this->kdg = kdg;
    this->ktg = ktg;
    this->kdh = kdh;
    this->kth = kth;
#endif
  }

  void setErosionParams ( double theta = 0.5, double spexp = 1.,
                          unsigned int athr = 0, unsigned int a = 0 )   
  {
    this->theta = theta;
    this->spexp = spexp;
    this->spexph = 1;
    this->athr = athr;
    this->a = a;
    this->ah = 9;
  }

  void setFlexureParams ( double alpha = 500., double rhoc = 0.85, double tau = 0. )
  {
    this->alpha = alpha;
    this->rhoc = rhoc;
    this->tau = tau;
  }

  void setPrecipParams ( double lc, double lf, double ll,
                         double ld, double refheight, double refslope,
                         double evap, double qin, int dir = 0 )
  {
    this->lc = lc;
    this->lf = lf;
    this->ll = ll;
    this->ld = ld;
    this->refheight = refheight;
    this->refslope = refslope; 
    this->evap = evap; 
    setInflux ( qin, dir );
  }

  void setInflux ( double qin, int dir = 0 )
  {
    this->qin = qin;
    this->dir = dir;
  }

  void resize ( int m, int n )
  {
    if ( this->m == m && this->n == n )  return;
    this->m = m;
    this->n = n;
    u.resize(m);
    prevrow.resize(m);
    prevrow[0] = m-1;
    for ( int i = 1; i < m; ++i )  prevrow[i] = i-1;
    nextrow.resize(m);
    for ( int i = 0; i < m-1; ++i )  nextrow[i] = i+1;
    nextrow[m-1] = 0;
    prevcol.resize(n);
    prevcol[0] = n-1;
    for ( int j = 1; j < n; ++j )  prevcol[j] = j-1;
    nextcol.resize(n);
    for ( int j = 0; j < n-1; ++j )  nextcol[j] = j+1;
    nextcol[n-1] = 0;
#ifdef DEFORM
    prevrow[0] = prevcol[0] = 0;
    nextrow[m-1] = m-1;
    nextcol[n-1] = n-1;
#endif
#ifdef REDREC
    seq.clear();
    seq.reserve(m*n);
#endif
    for ( int i = 0; i < m; ++i )
    {
      u[i].resize(n);
      for ( int j = 0; j < n; ++j )
      {
        Node *pn = getNode(i,j);
        pn->h = pn->l = pn->u = pn->b = 0;	      
        pn->d = Point(i,j);	      
#ifdef TRANSPORT
        pn->qs = pn->qp = 0;
#endif
#ifdef CHANNEL
        pn->channel = 0;
#endif
#ifdef PRECIP
        pn->p = 1;
#endif
#ifdef DEFORM
        pn->x = i;
        pn->y = j;
        pn->cs = 1.;
        if ( i == 0 || i == m-1 || j == 0 || j == n-1 )  pn->b = 1; 
#endif
#ifdef ERODIBILITY
	pn->e = 1.;
#endif
#ifdef DIFFUSIVITY
	pn->diff = diff;
#endif
#ifdef RIGIDITY
	pn->alpha = alpha;
#endif
#ifdef ICE
	pn->qi = pn->th = 0;
#endif
#ifdef FLEXURE
        pn->w = pn->g = 0;
#endif
#ifdef REDREC
        seq.push_back(Point(i,j));
#endif
      }
    }
#ifdef FLEXURE
    w.clear();
    res.clear();
    frhs.clear();
    rig.clear();
    do		
    {
      w.push_back(*new vector<vector<double> >);
      w.back().resize(m);
      for ( int i = 0; i < m; ++i )  w.back()[i].resize(n);     
      res.push_back(*new vector<vector<double> >);
      res.back().resize(m);
      for ( int i = 0; i < m; ++i )  res.back()[i].resize(n);     
      frhs.push_back(*new vector<vector<double> >);
      frhs.back().resize(m);
      for ( int i = 0; i < m; ++i )  frhs.back()[i].resize(n);     
#ifdef RIGIDITY
      rig.push_back(*new vector<vector<double> >);
      rig.back().resize(m);
      for ( int i = 0; i < m; ++i )  rig.back()[i].resize(n);     
#else  // NOT RIGIDITY
      rig.push_back(0.);
#endif  // NOT RIGIDITY
    }
    while ( m%2 == 0 && n%2 == 0 && (m/=2)>=2  && (n/=2)>=2 );
#endif  // FLEXURE
#ifdef LAKES
    offset.resize(m*n);
    int  k = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        offset[k].p = Point(i,j);
        int di = i <= m/2 ? i : m-i;
        int dj = j <= n/2 ? j : n-j;
        offset[k++].d = di*di+dj*dj;
      }
    offset[0].d = -1;
    sort(offset);
#endif
#ifdef PRECIP
    n = this->n >= 2 ? this->n : 2;
    if ( this->m > n )  n = this->m;
    diag.resize(n);
    upper.resize(n-1);
    lower.resize(n-1);
    right.resize(n-2);
    bottom.resize(n-2);
    rhs.resize(n);
    rhs.resize(n);
    psi.resize(n);
#endif  // PRECIP

  }    

#ifndef PRECIP
  void setErosionFactors()
  {
    spfac.resize(m*n);
    at = pow(a,theta);
    for ( int q = 0; q < spfac.size(); ++q ) 
      spfac[q] = q >= athr ? pow(q,theta)+at : 0.; 
  }
#endif

  Node * getNode ( int i, int j )
  {
    return  &u[i][j];  	  
  }

  Node * getNode ( Point p )
  {
    return  &u[p.i][p.j];  	  
  }

  Node * getNode ( Point *p )
  {
    return  &u[p->i][p->j];  	  
  }

  Node * getNode ( vector<Point>::iterator p )
  {
    return  &u[p->i][p->j];  	  
  }

  Node * getNodeP ( int i, int j )
  {
    return  &u[(i+m)%m][(j+n)%n];  	  
  }

  vector<Node> & operator [] ( int i )
  {
    return u[i];
  }  

  Node & operator [] ( Point p )
  {
    return u[p.i][p.j];
  }  

  Node & operator [] ( Point *p )
  {
    return u[p->i][p->j];
  }  

  Node & operator [] ( vector<Point>::iterator p )
  {
    return u[p->i][p->j];
  }  

  Point getNeighbor ( Point p, int i )
  {
    return  Point ( (p.i+shift[i].i+m)%m, (p.j+shift[i].j+n)%n );
  }

  vector<Point> getNeighbors ( Point p )
// Returns a vector of the 8 nearest and second nearest neighbors
  {
    vector<Point>  neigh = { Point(nextrow[p.i],p.j), 
                             Point(p.i,nextcol[p.j]),
                             Point(prevrow[p.i],p.j),
                             Point(p.i,prevcol[p.j]), 
                             Point(nextrow[p.i],prevcol[p.j]), 
                             Point(nextrow[p.i],nextcol[p.j]),
                             Point(prevrow[p.i],nextcol[p.j]),
                             Point(prevrow[p.i],prevcol[p.j]) };
//    vector<Point>  neigh = { Point(prevrow[p.i],prevcol[p.j]), 
//                             Point(prevrow[p.i],p.j),
//                             Point(prevrow[p.i],nextcol[p.j]),
//                             Point(p.i,prevcol[p.j]), 
//                             Point(p.i,nextcol[p.j]),
//                             Point(nextrow[p.i],prevcol[p.j]), 
//                             Point(nextrow[p.i],p.j),
//                             Point(nextrow[p.i],nextcol[p.j]) };
    return neigh;
  }

  vector<Node*> getNeighborsN ( Point p )
// Returns a vector of the 8 nearest and second nearest neighbors
  {
    vector<Node*>  neigh = { getNode(nextrow[p.i],p.j), 
                             getNode(p.i,nextcol[p.j]),
                             getNode(prevrow[p.i],p.j),
                             getNode(p.i,prevcol[p.j]), 
                             getNode(nextrow[p.i],prevcol[p.j]), 
                             getNode(nextrow[p.i],nextcol[p.j]),
                             getNode(prevrow[p.i],nextcol[p.j]),
                             getNode(prevrow[p.i],prevcol[p.j]) };
    return neigh;
  }

  vector<Point> getNearestNeighbors ( Point p )
// Returns a vector of the 4 nearest neighbors
  {
    vector<Point>  neigh = { Point(prevrow[p.i],p.j),
                             Point(p.i,prevcol[p.j]), 
                             Point(p.i,nextcol[p.j]),
                             Point(nextrow[p.i],p.j) };
    return neigh;
  }

  int computeDistSquare ( Point p, Point q )
  {
    int  di = abs(p.i-q.i);
    if ( 2*di > m )  di = m-di;
    int  dj = abs(p.j-q.j);
    if ( 2*dj > n )  dj = n-dj;
    return  di*di+dj*dj;
  }

#ifdef DEFORM
  double  det ( Node *pn1, Node *pn2, Node *qn1, Node *qn2 )
  {
    return  (pn1->x-pn2->x)*(qn1->y-qn2->y)-(qn1->x-qn2->x)*(pn1->y-pn2->y);
  }

  void computeCellSizes()
  {
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        Node  *pn = getNode(i,j); 
        vector<Node*>  neigh = getNeighborsN(Point(i,j));	    
        pn->cs = ( det(neigh[0],pn,neigh[5],neigh[4])
                 + det(neigh[1],pn,neigh[6],neigh[5])
                 + det(neigh[2],pn,neigh[7],neigh[6])
                 + det(neigh[3],pn,neigh[0],neigh[7]) ) / 8;
      }
  }
#endif

#ifdef REDREC
#ifdef DOWN
  void setSeq()
// Defines a sequence of points in downstream order. Boundary points are
// not included and marked with a positive value.  
  {
    int  nb = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  nb += getNode(i,j)->m = getNode(i,j)->b;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  ++getNode(getNode(i,j)->d)->m;
    seq.clear();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        Point  p(i,j);
        Node   *pn = getNode(p);
        while ( pn->m == 0 )
        {
          seq.push_back(p);
          --pn->m;
/*
          if ( pn->d == p )
          {
            printf ( "Fatal error\n" );
            exit(0);
          }
*/
          p = pn->d;  
          --(pn=getNode(p))->m;
        }
      }
    for ( vector<Point>::iterator p = seq.begin(); p != seq.end(); ++p )
      getNode(p)->m = 0;
    if ( seq.size()+nb != m*n )
    {
      printf ( "Fatal error: Wrong size of sequence %i\n", seq.size() );
      exit ( 0 );
    }
/*
    printf ( "%i\n", seq.size() );
    vector<int>  num(20,0);
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  ++num[getNode(i,j)->m+10];
    for ( int i = 0; i < num.size(); ++i )
      printf ( "%i %i\n", i-10, num[i] );
*/    
  }

#else

  void setSeq()
  {
    seq.clear();
    vector<Point>::iterator s = seq.begin();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( getNode(i,j)->b )
        {
          seq.push_back(Point(i,j));
          for ( ; s != seq.end(); ++s )
          {
            vector<Point>  neigh = getNeighbors(*s);	    
            for ( vector<Point>::iterator p = neigh.begin(); p != neigh.end(); ++p )
              if ( getNode(p)->drainsTo(*s) )
                seq.push_back(*p);
#ifdef CHANNEL
              else
                getNode(s)->channel -= ( getNode(s)->h >= getNode(p)->h ); 
#endif
#ifdef LAKES
#ifndef LAKEFLOWDIR
            if ( getNode(*s)->trunk != *s )  seq.push_back(getNode(*s)->trunk);
#endif
#endif
          }
        }
    if ( seq.size() != m*n )
    {
      printf ( "Fatal error: Wrong size of sequence %i\n", seq.size() );
      exit ( 0 );
    }
  }

#endif
#endif

  private:
  int computeWaterLevel ( Point p )
  {
    int     nc = 0;	  
    double  l = 0.;	  
    if ( getNode(p)->m )  return nc;
    getNode(p)->m = 1;
    l = getNode(p)->h;
    Point  dest = getNode(p)->d;
    if ( dest != p )  nc += computeWaterLevel(dest);
    if ( getNode(dest)->l > l )  l = getNode(dest)->l;
    if ( l != getNode(p)->l )  ++nc;
    getNode(p)->l = l;    
    return nc;
  }

  public:
  int computeWaterLevel ( )
// Recomputes the water level using given flow directions and returns the number
// of sites where the water level has changed. Calling this function makes sense
// after changing surface elevations by any other process than fluvial erosion
// before recomputing the flow directions. If the surface has been changed only
// by fluvial erosion, calling this method has no effect.
  {
    int nc = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        nc += computeWaterLevel(Point(i,j));
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
    return nc;
  }

#ifdef LAKEFLOWDIR
  unsigned int computeFlowDirection ( Point p )
  {
// Computes the flow direction of a node. Leaves it unchanged if no downward
// direction was found or if the steepest downward direction results in an
// invalid flow pattern.	  
    Node           *pn = getNode(p); 
    if ( pn->b )   return 0;
    Point          pd = pn->d;
    double         s, maxs = 0.; 
    vector<Point>  neigh = getNeighbors(p);	    
    for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
    {
#ifdef DEFORM
      s = (pn->h-getNode(u)->h)/pn->distance(getNode(u));
#else      
      s = pn->h - getNode(u)->h;
      if ( p.isDiag(*u) )  s *= sqrt(0.5);
#endif
      if ( s > maxs )
      {
        pd = *u;
	maxs = s;
      }
    }
    if ( pd == pn->d )  return 0;
    Point d = pd;
    Node *dn;
    while ( !(dn=getNode(d))->b )
    { 
      if ( d == p )
      {
//        pn->q = 1;
//        printf ( "Loop %i %i\n", p.i, p.j );
        return 0;
      }
      Node *dn = getNode(d);
      if ( dn->l > pn->l )
      {
//        pn->q = 2;
//        printf ( "Increasing water level %i %i %g %i %i %g\n", p.i, p.j, pn->l, d.i, d.j, dn->l );
        return 0;
      }
      if ( dn->l < pn->l )
      {
//        if ( pn->l > pn->h )
//          printf ( "Success %i %i %g %i %i %g\n", p.i, p.j, pn->l, d.i, d.j, dn->l );
        pn->d = pd;
        return 1;
      }
      d = dn->d;
    }
//    if ( pd != d )
//      printf ( "Boundary %i %i %g %i %i %g\n", p.i, p.j, pn->l, d.i, d.j, dn->l );
    pn->d = pd;
    return 1;
  }

#else  // NOT LAKEFLOWDIR
  unsigned int computeFlowDirection ( Point p )
  {
// Computes the flow direction of a node. Leaves it unchanged if no downward
// direction was found.	  
    Node           *pn = getNode(p); 
    if ( pn->b )  return 0;
    Point          pd = pn->d;
    double         s, maxs = 0.; 
    vector<Point>  neigh = getNeighbors(p);	    
    for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
    {
#ifdef DEFORM
      s = (pn->l-getNode(u)->l)/pn->distance(getNode(u));
#else      
      s = pn->l - getNode(u)->l;
      if ( p.isDiag(*u) )  s *= sqrt(0.5);
#endif
      if ( s > maxs )
      {
        pn->d = *u;
	maxs = s;
      }	
    }
    return pd != pn->d;
  }  
#endif  // NOT LAKEFLOWDIR

  unsigned int computeFlowDirection()
// Computes the flow direction of all nodes based on the water level (not
// surface elevation). Leaves those nodes unchanged where no downward
// direction was found. Returns the number of nodes where the flow direction
// has changed.	  
  {
    int  nchanges = 0; 
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        nchanges += computeFlowDirection(Point(i,j));
    return nchanges;
  }

//  private:
#ifdef ICE
  void propagateUpstream ( Point p, Point pb )
  {
    Point  pr = pb;
    for ( int i = 0; ; ++i )
    {	    
      int  di = 2*abs(p.i-pr.i);
      if ( di > m )  di = 2*m-di;
      int  dj = 2*abs(p.j-pr.j);
      if ( dj > n )  dj = 2*n-dj;
      double dist;
      if ( ( dist = (di*di+dj*dj)/(getNode(pr)->th*getNode(pr)->th) ) <= 1 )
      {	      
        if ( getNode(pb)->qi > getNode(p)->qi )
        {
          getNode(p)->qi = getNode(pb)->qi;
          getNode(p)->qp = dist;
        }
        vector<Point>  neigh = getNeighbors(p);	    
        for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
          if ( getNode(*u)->drainsTo(p) )  
            propagateUpstream(*u,pb);
	return;
      }
      Point pnew = getNode(pr)->trunk;
      if ( pnew == pr )  return;
      pr = pnew;
    }
  }

  double computeIceFlux ( Point p )
  {
    Node  *pn = getNode(p);
    if ( pn->m )  return pn->qi;
    pn->m = 1;
#ifdef ICEFRAC
    double  f = (ICEFRAC);
#else
    double  f = (pn->h+pn->th-ela)/(fsa-ela);
#endif
#ifdef MELTRATE
    if ( f < 0 )  f = 0;
#endif
#ifdef PRECIP
    pn->q = pn->p;
    pn->qi = f*pn->p;
    double  qin, qmax = 0;
#else
    ++pn->q;
    pn->qi += f;
    int  qin, qmax = 0;
#endif
#ifdef MELTBOUNDARY
    if ( pn->b )  return pn->qi = pn->th = 0;
#endif
    vector<Point>  neigh = getNeighbors(p);
    pn->trunk = p;    
    for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
      if ( getNode(u)->drainsTo(p) )  
      {
        pn->qi += computeIceFlux(*u);
        pn->q += qin = getNode(u)->q;
        if ( qin > qmax )
	{
	  qmax = qin;
	  pn->trunk = *u;
	}  
      }    
    if ( pn->qi > 0. )
    {
      if ( pn->qi > pn->q )  pn->qi = pn->q;	      
#ifdef MELTRATE
      f = (MELTRATE);
//printf ( "th %i %i %e\n", p.i, p.j, pn->qi );
      if ( f > 0 && ( pn->qi -= f*gwf*pow(pn->qi,gwe) ) < 0 )  pn->qi = 0;
#endif
// Set thickness to glacier width for the moment
      pn->th = gwf*pow(pn->qi,gwe);
//printf ( "th %i %i %e\n", p.i, p.j, pn->qi );
      for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
        if ( getNode(u)->drainsTo(p) && *u != pn->trunk )  
//        if ( getNode(u)->drainsTo(p) && getNode(u)->q < qmax )  
          propagateUpstream ( *u, p );
    }
#ifdef MELTRATE
    else  pn->qi = pn->th = 0;
#else
    else  pn->th = 0;
#endif
    return  pn->qi;
  }      

//  public:
  void computeIceFlux()
  {
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        getNode(i,j)->qp = getNode(i,j)->qi = getNode(i,j)->q = getNode(i,j)->m = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        computeIceFlux ( Point(i,j) );
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  
      {
        if ( getNode(i,j)->qi < 0 )  getNode(i,j)->qi = 0;
	getNode(i,j)->m = 0;
      }
  }  

  void computeIceThickness ( Point p )
  {
    Node  *pn = getNode(p);
    if ( pn->m )  return;
    pn->m = 1;
    pn->th *= ghw;
    Point  dest = pn->d;
    Node  *pd = getNode(dest);
    if ( dest != p )
    {
      computeIceThickness(dest);
// not sure whether this is completely safe
      if ( pd->qi == pn->qi && pd->th > pn->th )  pn->th = pd->th; 
    }
  }

  double addIceLayer ( Point p )
  {
    double  v = 0., dh;	  
    Node  *pn = getNode(p);
    if ( pn->m )  return 0;
    pn->m = 1;
    if ( uexp > 0 && pn->qp > 0 )  pn->th *= 1-pow(pn->qp,uexp/2);
    if ( uexp < 0 && pn->qp > 0 )  pn->th *= sqrt(1-pn->qp);
    Point  dest = pn->d;
    Node  *pd = getNode(dest);
    if ( dest != p )
    {
      v += addIceLayer(dest);
      if ( pd->th > 0 && ( dh = pd->h-pn->h-pn->th ) > 0 )  pn->th += dh;
    }
    v += pn->th; 
    pn->h += pn->th; 
    return v;
  }

  double addIceLayer ( )
  {
    double  v = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        computeIceThickness ( Point(i,j) );
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        v += addIceLayer ( Point(i,j) );
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
    return v;
  }

  void removeIceLayer ( Point p )
  {
    Node  *pn = getNode(p);
    if ( pn->m )  return;
    pn->m = 1;
    pn->h -= pn->th;
    pn->l = pn->h;
    Point  dest = pn->d;
    if ( dest != p )
    {
      removeIceLayer(dest);
      if ( getNode(dest)->l > pn->l )  pn->l = getNode(dest)->l;
    }
  }

  void removeIceLayer ( )
  {
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        removeIceLayer ( Point(i,j) );
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
  }
#endif  // ICE

#if defined(PRECIP) || defined(DEFORM)
  double computeDischarge ( Point p )
#else
  int computeDischarge ( Point p )
#endif
  {
// Computes the discharge of a node
    Node  *pn = getNode(p);
    if ( pn->m )  return pn->q;
    pn->m = 1;
#ifdef PRECIP
#ifdef DEFORM
    pn->q += pn->cs*pn->p;
#else
    pn->q += pn->p;
#endif
#else
#ifdef DEFORM
    pn->q += pn->cs;
#else
    ++(pn->q);
#endif
#endif
    vector<Point>  neigh = getNeighbors(p);	    
    for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
      if ( getNode(u)->drainsTo(p) )  
        pn->q += computeDischarge(*u);
    return  pn->q;
  }  

//  public:
  void computeDischarge()
  {
// Computes all discharges
#ifdef REDREC
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
#ifdef PRECIP
#ifdef DEFORM
        getNode(i,j)->q = getNode(i,j)->cs*getNode(i,j)->p;
#else
        getNode(i,j)->q = getNode(i,j)->p;
#endif
#else
#ifdef DEFORM
        getNode(i,j)->q = getNode(i,j)->cs;
#else
        getNode(i,j)->q = 1;
#endif
#endif
    setSeq();
#ifdef DOWN
    for ( vector<Point>::iterator p = seq.begin(); p != seq.end(); ++p )
      getNode(getNode(*p)->d)->q += getNode(*p)->q;
#else
    for ( vector<Point>::reverse_iterator p = seq.rbegin(); p != seq.rend(); ++p )
      if ( !getNode(*p)->b ) 
        getNode(getNode(*p)->d)->q += getNode(*p)->q;
#endif
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
#if !defined(PRECIP) && !defined(DEFORM)
    int q = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( getNode(i,j)->b )  q += getNode(i,j)->q;
    if ( q != m*n )
    {
      printf ( "Error: Total discharge is %i, should be %i\n", q, m*n );
      exit ( -1 );
    } 
#endif
#else
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        Node  *pn = getNode(i,j);
        pn->q = pn->m = 0;
      }
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        computeDischarge ( Point(i,j) );
#endif
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
  }  

  double computeFluxes ( Point p, double dt = 0. )
  {
// Computes the discharge of a node
    Node  *pn = getNode(p);
    if ( pn->m )  return pn->q;
    pn->m = 1;
    Point   dest = pn->d;
    Node    *destn = getNode(dest);
#ifdef TRANSPORT
    double  qssum = 0., qpsum = 0.;
#endif
    vector<Point>  neigh = getNeighbors(p);	    
#ifdef LAKES
#ifndef LAKEFLOWDIR
    if ( pn->trunk != p )  neigh.push_back(pn->trunk);
#endif
#endif
    for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
    {
      if ( getNode(u)->drainsTo(p) )  
      {
        if ( !getNode(u)->m )  ++nrec;
#ifdef REDREC
#ifdef DOWN
        pn->q += getNode(u)->q;
#endif
#else
//        pn->q += computeFluxes(*u,dt);
        computeFluxes(*u,dt);
#endif
#ifdef TRANSPORT
        qssum += getNode(u)->qs; 
        qpsum += getNode(u)->qp; 
#endif
#ifdef CHANNEL
#if !defined(REDREC) || defined(DOWN) 
//        pn->channel += 10*getNode(u)->channel;
#endif
#endif
      }
#ifdef CHANNEL
#if !defined(REDREC) || defined(DOWN) 
      else
        pn->channel -= ( pn->h >= getNode(u)->h ); 
#endif
#endif
    }      
#ifdef CHANNEL
    pn->channel = pn->channel > 0 && !pn->b;
#endif
#ifdef TRANSPORT
//    if ( !pn->b )
//    {
//      if ( fabs((qssum-pn->qs)/qssum) > qdev )  qdev = fabs((qssum-pn->qs)/qssum);
//      if ( fabs((qpsum-pn->qp)/qpsum) > qdev )  qdev = fabs((qpsum-pn->qp)/qpsum);
//    }
    qssum = pn->qs;
    qpsum = pn->qp;
    if ( !pn->b && dt > 0. )
    {
#ifdef DEFORM
      double  cs = pn->cs;
#else
      double  cs = 1;
#endif
      double  alpha = cs-dt*qpsum;
      double  u = pn->u;
      double  q = pn->q;
      double  spf = 1e-10;
#ifdef DIFFUSIVITY		
      double  diff = pn->diff;
#else
      double  diff = this->diff;
#endif
#ifdef DIFFUS
// pn->r must still be thickness per time
      u += pn->r;
#ifdef EIGHTNEIGHBORS
      diff *= sqrt(2)-1;
#endif  // EIGHTNEIGHBORS
#endif  // DIFFUS

#ifdef CHANNEL
#ifdef NOCHANNELDIFFUS
    if ( pn->channel )  diff = 0;
#endif
#endif

#ifdef ICE
      double  ice = pn->qi/pn->q;
      if ( ice > 1 )
      {	      
	q = pn->qi*pow(ice,-eps);      
        ice = 1;
      }
#ifdef NOMELTWATER
      else
        if ( ice > 0 )
	{
	  q = pn->qi;
          ice = 1;
        }	  
#endif  // NOMELTWATER
//      else  if ( ice < 0 )  ice = 0;
      pn->qeff = q;

#endif  // ICE

#ifdef LAKES
#ifndef LINDECL
// Not sure whether this is correct
      if ( pn->l <= pn->h || dt*(cs*u+qssum) >= lff*(pn->l-pn->h)*alpha )   
//      if ( pn->l <= pn->h || dt*(u+qssum) >= lff*(pn->l-pn->h)*(1-dt*qpsum) )   
#endif  // NOT LINDECL
#endif  // LAKES

#ifdef CHANNEL
#ifdef ICE
      if ( !pn->channel && ice == 0 )  spf = 1;
      else
#else
#ifdef SHAREDSP
      if ( !pn->channel )  spf = 1;
      else
#else
      if ( !pn->channel )  spf = pow(ah,theta);
      else
#endif
#endif
#endif

#if defined(ICE) || defined(PRECIP) || defined(DEFORM) 
      spf = q >= athr ? pow(q,theta)+at : 0.;
#else  
      spf = spfac.size() ? spfac[pn->q] : ( q >= athr ? pow(q,theta)+at : 0. );
#endif 

#ifdef LINDECL
      double  lambda = this->lambda;
#ifdef SHAREDSP

#ifdef LAYERS
      int layer = pn->getLayer(); 
      double  kd = this->kd[layer];
      double  kt = this->kt[layer];
      double  kdg = this->kdg[layer];
      double  ktg = this->ktg[layer];
      double  kdh = this->kdh[layer];
      double  kth = this->kth[layer];
#else  // NOT LAYERS
      double  kd = this->kd;
      double  kt = this->kt;
#endif  // NOT LAYERS

#ifdef CHANNEL
#ifdef ICE
      if ( !pn->channel && ice == 0 )
      {
        kd = kdh;
        kt = kth;
      }
#else
      if ( !pn->channel )
      {
        kd = kdh;
        kt = kth;
        double  spexp = pn->channel ? spexph : this->spexp;       
      }
#endif
#endif  // CHANNEL

      lambda = kd/kt;

#ifdef ICE
      if ( ice > 0 )
	if ( ice >= 1 )
	{
          kd = kdg;
          lambda = kdg/ktg;
        }
	else
	{	
          double  tmp1 = (1-ice)*kt/kd;
          double  tmp2 = ice*ktg/kdg;
          lambda = 1 / ( tmp1 > tmp2 ? tmp1 : tmp2 );
	  kd = pow(1-ice,theta*spexp)*kd + pow(ice,theta*spexp)*kdg;
	}  
#endif  // ICE
#endif  // SHAREDSP
#endif  // LINDECL

      double  f = spf;
#ifdef DEFORM 
      double  df = 1/pn->distance(destn);
#else
      double  df = p.isDiag(dest) ? sqrt(0.5) : 1.;
#endif

#ifdef ERODIBILITY
      f *= pn->e;
#endif  // ERODIBILITY 

#ifdef SPEXP
      f *= pow(spf*fabs(pn->h-destn->h)*df,spexp-1.);
#endif  // SPEXP

#ifdef LINDECL
#ifdef LAKES
      if ( pn->l > pn->h )
      {
        double  me = lff*(pn->l-pn->h)/dt-u;
        if ( me > 0 )
        {
          double  qstot = qssum+lff*(pn->l-pn->h)*qpsum-cs*me;
//          double  lambdanew = qstot > 0 ? q*me/qstot : 1e100;        
//          if ( lambdanew > lambda )  lambda = lambdanew;
          lambda = qstot > 0 ? q*me/qstot : 1e100;        
          f = 0;
        }
      }
#endif  // LAKES

#ifdef SHAREDSP
      f *= kd;
#else  // NOT SHAREDSP
      f *= 1+lambda;
#endif  // NOT SHAREDSP

#ifdef DIFFSEMIIMPL
// Not compatible with DEFORM
#ifdef EIGHTNEIGHBORS
      f += lambda*diff/q;
#else
      f += lambda*diff/q*df;
#endif
#endif  // DIFFSEMIIMPL

      f *= df;
      double  den = alpha*lambda/q+1+dt*f;
      pn->qs = (alpha*(f*(pn->h-destn->h)-u)+(cs*u+qssum)*(1+dt*f))/den;

#ifdef DEPOS 
#ifdef LAKES
      if ( pn->l <= pn->h )
#endif
      if ( lambda/q*pn->qs > f*(pn->h+dt*u-destn->h) )
      {
        f *= 1e100;
        lambda *= 1e100;
        den = alpha*lambda/q+1+dt*f;
        pn->qs = (alpha*(f*(pn->h-destn->h)-u)+(cs*u+qssum)*(1+dt*f))/den;
      }
#endif
      pn->qp = -alpha*f/den;
#ifdef BUFFERS
      pn->lambdaq = lambda/q;
#endif  // BUFFERS    

#else  // NO LINDECL

      f *= q;
  
#ifdef DIFFSEMIIMPL
// Not compatible with DEFORM
#ifdef EIGHTNEIGHBORS
      f += diff;
#else
      f += diff*df;
#endif
#endif  // DIFFSEMIIMPL

      f *= df;
      double  den = dt+alpha/f;
      pn->qp = -alpha/den;
      pn->qs = (dt*(cs*u+qssum)+alpha*(pn->h-destn->h))/den;
#endif  // NO LINDECL
      pn->l = f;
#ifndef DOWN
      destn->qs += pn->qs;
      if ( !destn->b )  destn->qp += pn->qp;
#endif
    }
#endif  // TRANSPORT
    if ( !pn->b )
    {
      destn->q += pn->q;
#ifdef CHANNEL
      destn->channel += 10*pn->channel;
#endif
    }
    return  pn->q;
  }  

//  public:
  void computeFluxes ( double dt = 0. )
  {
// Computes all discharges and checks for a consistent drainage pattern
#ifdef ICE
    computeIceFlux();
    addIceLayer();
#endif
#ifdef LAKES
#ifndef LAKEFLOWDIR
    findLakes();
#endif
#endif
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        Node  *pn = getNode(i,j);
        pn->m = 0;
#ifdef PRECIP
#ifdef DEFORM
    pn->q = pn->cs*pn->p;
#else
    pn->q = pn->p;
#endif
#else
#ifdef DEFORM
    pn->q = pn->cs;
#else
    pn->q = 1;
#endif
#endif
#ifdef TRANSPORT
        pn->qs = pn->qp = 0.;
#endif
#ifdef CHANNEL
        pn->channel = pn->l == pn->h ? 2 : 10;
#endif
#ifdef DIFFUS
        pn->r = 0.;
#endif
      }
#ifdef DIFFUS  
// Not yet compatible with DEFORM
    double  diff = this->diff;
#ifdef EIGHTNEIGHBORS
    diff *= sqrt(2)-1;
#endif
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        Point  p(i,j);
        Node  *pn = getNode(p);
#ifdef DIFFUSIVITY		
        diff = pn->diff;
#ifdef EIGHTNEIGHBORS
        diff *= sqrt(2)-1;
#endif  // EIGHTNEIGHBORS
#endif  // DIFFUSIVITY

#ifdef CHANNEL
#ifdef NOCHANNELDIFFUS
    if ( pn->channel )  diff = 0;
#endif
#endif

#ifdef EIGHTNEIGHBORS
        vector<Point>  neigh = getNeighbors(p);	    
        for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
#ifdef TRANSPORT
#ifdef DIFFSEMIIMPL
          if ( !pn->drainsTo(*u) )  
#endif
#endif
          {
            double  dh = pn->h-getNode(u)->h;
            if ( dh > 0. )
            {
              double  flux = diff*dh;
              if ( u->isDiag(p) )  flux *= sqrt(0.5);
              pn->r -= flux;
              getNode(u)->r += flux;
            }
          }
#else  // NOT EIGHTNEIGHBORS
        for ( int k = i-1; k <= i+1; k += 2 )
        {
          int kper = (k+m)%m;
          double  dh = pn->h-getNode(kper,j)->h;
          if ( dh > 0. )
#ifdef TRANSPORT
#ifdef DIFFSEMIIMPL
          if ( kper != pn->d.i )
#endif
#endif
          {
            double  flux = diff*dh;
            pn->r -= flux;
            getNode(kper,j)->r += flux;
          }
        }   
        for ( int k = j-1; k <= j+1; k += 2 )
        {
          int kper = (k+n)%n;
          Node *q = getNode(i,kper);
          double  dh = pn->h-getNode(i,kper)->h;
          if ( dh > 0. )
#ifdef TRANSPORT
#ifdef DIFFSEMIIMPL
          if ( kper != pn->d.j )
#endif
#endif
          {
            double  flux = diff*dh;
            pn->r -= flux;
            getNode(i,kper)->r += flux;
          }
        }   
#endif  // NOT EIGHTNEIGHBORS
      }
#endif  // DIFFUS

    nrec = 0;
//    qdev = 0;
#ifdef REDREC
    setSeq();

#ifdef DOWN
    for ( vector<Point>::iterator p = seq.begin(); p != seq.end(); ++p )
#else
    for ( vector<Point>::reverse_iterator p = seq.rbegin(); p != seq.rend(); ++p )
#endif
    {
      computeFluxes ( *p, dt );
/*
      Node  *pn = getNode(*p);
      if ( pn->m ) 
      {
        printf ( "Fatal error\n" );
        exit(0);
      }
      int  tmp = computeFluxes ( *p, dt );
      if ( pn->b )  q += tmp;
*/
    }
//    printf ( "qdev = %e\n", qdev );
#else        
    int  q = 0;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        auto point = Point(i,j);
        int  tmp = computeFluxes ( point, dt );
	if ( getNode(i,j)->b )  q += tmp;
      }	
#if !defined(PRECIP) && !defined(DEFORM)
    if ( q != m*n )
    {
      printf ( "Error: Total discharge is %i, should be %i\n", q, m*n );
      exit ( -1 );
    } 
#endif
#endif
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
  }  

  private:
  double erode ( double dt, Point p )
// Performs a timestep of length dt for a given point and returns
// the absolute change in elevation
  {
    double  e, maxc = 0.; 
    Node  *pn = getNode(p);
    Point  dest = pn->d;
#ifndef REDREC
    if ( pn->m )  return 0;
    pn->m = 1;
    if ( pn->b )  return 0;
    if ( dest == p )
    {
      printf ( "Error: Detected node without valid flow target: %i %i\n",
	       p.i, p.j );
      exit ( -1 );
    } 
#ifndef PRECIP
    if ( pn->q == 0 )
    {
      printf ( "Error: Detected node without discharge: %i %i\n",
	       p.i, p.j );
      exit ( -1 );
    } 
#endif    
    if ( ( e = erode(dt,dest) ) > maxc )  maxc = e;
#endif
    double  hold = pn->h;
    double  q = pn->q;
    double  spf;
    double  u = pn->u;

#ifdef TRANSPORT

#ifdef DIFFUS
    u += pn->r;
#endif    

    pn->qs += pn->qp*getNode(dest)->qp;

#ifdef LINDECL

#ifdef SHAREDSP
   
#ifdef BUFFERS
    double  lambdaq = pn->lambdaq;
#ifdef ICE
    q = pn->qeff;
#endif  // ICE

#else  // NO BUFFERS

#ifdef LAYERS
    int layer = pn->getLayer(); 
    double  kd = this->kd[layer];
    double  kt = this->kt[layer];
    double  kdg = this->kdg[layer];
    double  ktg = this->ktg[layer];
#else  // NO LAYERS
    double  kd = this->kd;
#endif  // NO LAYERS
    double  lambda = kd/kt;
    double  lambdaq = lambda/q;
#endif  // NO BUFFERS

#else  // NO SHAREDSP

#ifdef BUFFERS
    double  lambdaq = pn->lambdaq;
#else  // NO BUFFERS
    double  lambdaq = lambda/q;
#endif

#endif  // NO SHAREDSP
    double  hnew = (pn->h+dt*(pn->l*getNode(dest)->h+lambdaq*pn->qs+u))/(1+dt*pn->l);

#ifdef DEPOS
    if ( pn->qs < 0 )  ++oppos;
    if ( lambdaq > 1e10 )
      if ( hnew > pn->h+u*dt )  ++dtrue;
      else                      
      {
        ++efalse;
//        printf ( "%i %e\n", pn->q, (pn->h-hnew)/dt+u );
      }
    else
      if ( hnew > pn->h+u*dt )  ++dfalse;
      else                      ++etrue;
#endif

#else  // NO LINDECL
    double  hnew = getNode(dest)->h+pn->qs/pn->l;
#endif  // NO LINDECL

    pn->qp = hnew-pn->h;
    pn->h = hnew;
    pn->l = pn->h >= getNode(dest)->l ? pn->h : getNode(dest)->l;

#else  // Detachment-limited 

#ifdef PRECIP
    spf = q >= athr ? pow(q,theta)+at : 0.;
#else
    spf = spfac.size() ? spfac[pn->q] : ( q >= athr ? pow(q,theta)+at : 0. );
#endif
    double  f = spf;
    double  df = p.isDiag(dest) ? sqrt(0.5) : 1.;
#ifdef ERODIBILITY
    f *= pn->e;
#endif    

#ifdef SPEXP
    if ( pn->l <= getNode(dest)->l )  f = 0;
    else  f *= pow(spf*(pn->l-getNode(dest)->l)*df,spexp-1.);
#endif

//    e = 0.;

#ifdef DIFFUS
    u += f > 0 && pn->r > 0  ? (1-excav)*pn->r : pn->r;
#endif    

    f *= df;
    pn->l = (pn->h+dt*(u+f*getNode(dest)->l))/(1.+dt*f);
    if ( pn->l < getNode(dest)->l )
    {
      pn->l = getNode(dest)->l;
      pn->h += dt*u;
    }
    else
      pn->h = pn->l;  
#endif
    double dh = pn->h-hold;
    if ( ( e = fabs(dh) ) > maxc )  maxc = e;
    return maxc;
  }

  public:
  double erode ( double dt )
// Performs a timestep of length dt for the entire grid and returns the
// maximum absolute change in elevation
  {
#ifdef DEPOS
    etrue = efalse = dtrue = dfalse = oppos = 0;
#endif
    t += dt;
    at = pow(a,theta);
    computeFluxes(dt);
    double maxc = 0.;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        getNode(i,j)->m = 0;

#ifdef REDREC
#ifdef DOWN
    for ( vector<Point>::reverse_iterator p = seq.rbegin(); p != seq.rend(); ++p )
#else
    for ( vector<Point>::iterator p = seq.begin(); p != seq.end(); ++p )
#endif
    {
      double  e = erode ( dt, *p );
      if ( e > maxc)  maxc = e;
    }
#else
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        double  e = erode ( dt, Point(i,j) );
        if ( e > maxc)  maxc = e;
      }  
#endif
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  getNode(i,j)->m = 0;
#ifdef TRANSPORT
#ifdef ICE
    removeIceLayer();
#endif
#ifdef LAKES
#ifndef LAKEFLOWDIR
    for ( int i = 0; i < lakes.size(); ++i )
      lakes[i].computeFlowDirection();
    computeWaterLevel();
#endif
#endif
#ifdef LAYERS
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( getNode(i,j)->b == 0 )
          getNode(i,j)->adjustLayers(getNode(i,j)->qp-getNode(i,j)->u*dt,t);
#endif
#endif
    return maxc;
  }

  private:
  void adjustWaterLevel ( Point p )
  {
    double  e, maxc = 0.;

    if ( getNode(p)->m )  return;
    getNode(p)->m = 1;
    getNode(p)->l = getNode(p)->h;
    if ( !getNode(p)->b )
    {
      Point  dest = getNode(p)->d;
      if ( dest == p )
      {
        printf ( "Error: Detected node without valid flow target: %i %i\n",
                 p.i, p.j );
        exit ( -1 );
      }
#ifndef PRECIP
      if ( getNode(p)->q == 0 )
      {
        printf ( "Error: Detected node without valid flow target: %i %i\n",
                 p.i, p.j );
        exit ( -1 );
      }
#endif
      adjustWaterLevel(dest);
      if ( getNode(p)->l < getNode(dest)->l )  getNode(p)->l = getNode(dest)->l;
    }
  }

  public:
  double computeLakeVolume()
// Computes to total volume of all filled lakes.
  {
    double  v = 0.;	  
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  v += getNode(i,j)->l - getNode(i,j)->h;
    return v;
  }

  void sort ( vector<Point>& p )
  {
    int     i, j, l = p.size()>>1, ir = p.size()-1;
    double  h;
    Point   pxt;

    if ( ir == 0 )  return;
    for ( ; ; )
    {
      if ( l )  h = getNode(pxt=p[--l])->l;
      else
      {
        h = getNode(pxt=p[ir])->l;
        p[ir] = p[0];
        if ( !(--ir) )
        {
          p[0] = pxt;
          return;
        }
      }
      j = ( (i=l) << 1 ) + 1;
      while ( j <= ir )
      {
        if ( j < ir && getNode(p[j])->l > getNode(p[j+1])->l )  j++;
        if ( h > getNode(p[j])->l )
        {
          p[i] = p[j];
          j += (i=j)+1;
        }
        else  j = ir+1;
      }
      p[i] = pxt;
    }
  }

  template <typename valuetype>
  void sort ( vector<PointValue<valuetype>>& p )
  {
    int         i, j, l = p.size()>>1, ir = p.size()-1;
    valuetype   h;
    PointValue<valuetype>  pxt;

    if ( ir == 0 )  return;
    for ( ; ; )
    {
      if ( l )  h = (pxt=p[--l]).d;
      else
      {
        h = (pxt=p[ir]).d;
        p[ir] = p[0];
        if ( !(--ir) )
        {
          p[0] = pxt;
          return;
        }
      }
      j = ( (i=l) << 1 ) + 1;
      while ( j <= ir )
      {
        if ( j < ir && p[j].d < p[j+1].d )  j++;
        if ( h < p[j].d )
        {
          p[i] = p[j];
          j += (i=j)+1;
        }
        else  j = ir+1;
      }
      p[i] = pxt;
    }
  }

#ifdef FLEXURE 
  void restrict ( vector<vector<double> > &src, vector<vector<double> > &dest,
                  vector<vector<double> > &tmp )
  {
    int  m = src.size();
    int  n = src[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; j += 2 )
        getNode(i,j/2)->w = src[i][(j+n-1)%n]+2*src[i][j]+src[i][(j+1)%n];
    for ( int i = 0; i < m; i += 2 )
      for ( int j = 0; j < n/2; ++j )
        dest[i/2][j] = (getNode((i+m-1)%m,j)->w+2*getNode(i,j)->w+getNode((i+1)%m,j)->w)/16;
  }

  void prolongate ( vector<vector<double> > &src, vector<vector<double> > &dest )
  {
    int  m = src.size();
    int  n = src[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  dest[2*i][2*j] = src[i][j];
    m = 2*m;
    n = 2*n;
    for ( int i = 0; i < m; i += 2 )
      for ( int j = 1; j < n; j += 2 )
        dest[i][j] = 0.5*(dest[i][(j+n-1)%n]+dest[i][(j+1)%n]); 
    for ( int i = 1; i < m; i += 2 )
      for ( int j = 0; j < n; ++j )
        dest[i][j] = 0.5*(dest[(i+m-1)%m][j]+dest[(i+1)%m][j]);
  }

#ifdef RIGIDITY
  void gaussseidel ( vector<vector<double> > &u, vector<vector<double> > &rhs,
                     double tau, vector<vector<double> > &rig )
  {
    double  maxchange;
    double  den = 1 + tau + 20*rig[0][0];
    int  m = u.size();
    int  n = u[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( rhs[i][j] < 1e50 )
          u[i][j] = ( rhs[i][j]
                      + 4*rig[i][j] * ( u[(i+m-1)%m][j] + u[(i+1)%m][j]
                                      + u[i][(j+n-1)%n] + u[i][(j+1)%n] )
                      + rig[(i+m-1)%m][j] * ( 4*u[(i+m-1)%m][j] - u[(i+m-2)%m][j] 
                                            - u[(i+m-1)%m][(j+n-1)%n] - u[(i+m-1)%m][(j+1)%n] )
                      + rig[(i+1)%m][j] * ( 4*u[(i+1)%m][j] - u[(i+2)%m][j] 
                                            - u[(i+1)%m][(j+n-1)%n] - u[(i+1)%m][(j+1)%n] )
                      + rig[i][(j+n-1)%n] * ( 4*u[i][(j+n-1)%n] - u[i][(j+n-2)%n]
                                            - u[(i+m-1)%m][(j+n-1)%n] - u[(i+1)%m][(j+n-1)%n] )
                      + rig[i][(j+1)%n] * ( 4*u[i][(j+1)%n] - u[i][(j+2)%n]
                                          - u[(i+m-1)%m][(j+1)%n] - u[(i+1)%m][(j+1)%n] ) )
                  / ( 1 + tau + 16*rig[i][j] 
                      + ( rig[(i+m-1)%m][j] + rig[(i+1)%m][j]
                            + rig[i][(j+n-1)%n] + rig[i][(j+1)%n] ) );
  }

  void residuum ( vector<vector<double> > &u, vector<vector<double> > &rhs,
                  vector<vector<double> > &res, double tau,
                  vector<vector<double> > &rig )
  {
    int  m = u.size();
    int  n = u[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        getNode(i,j)->w = rig[i][j] * ( u[(i+m-1)%m][j] + u[(i+1)%m][j]
                                      + u[i][(j+n-1)%n] + u[i][(j+1)%n]
                                      - 4*u[i][j] );
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        res[i][j] = rhs[i][j] - getNode((i+m-1)%m,j)->w - getNode((i+1)%m,j)->w
                              - getNode(i,(j+n-1)%n)->w - getNode(i,(j+1)%n)->w
                              + 4*getNode(i,j)->w - (1+tau)*u[i][j];
  }

#else
  void gaussseidel ( vector<vector<double> > &u, vector<vector<double> > &rhs,
                     double tau, double rig )
  {
    double  maxchange;
    double  den = 1 + tau + 20*rig;
    int  m = u.size();
    int  n = u[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( rhs[i][j] < 1e50 )
          u[i][j] = ( rhs[i][j]
                      + 8*rig * ( u[(i+m-1)%m][j] + u[(i+1)%m][j]
                                + u[i][(j+n-1)%n] + u[i][(j+1)%n] )
                      - rig * ( u[(i+m-2)%m][j] + u[(i+2)%m][j]
                                + u[i][(j+n-2)%n] + u[i][(j+2)%n] 
                                + 2 * ( u[(i+m-1)%m][(j+n-1)%n] 
                                      + u[(i+m-1)%m][(j+1)%n]
                                      + u[(i+1)%m][(j+n-1)%n]
                                      + u[(i+1)%m][(j+1)%n] ) ) ) / den;
  }

  void residuum ( vector<vector<double> > &u, vector<vector<double> > &rhs,
                  vector<vector<double> > &res, double tau, double rig )
  {
    double  diag = 1 + tau + 20*rig;
    int  m = u.size();
    int  n = u[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        getNode(i,j)->w = rig * ( u[(i+m-1)%m][j] + u[(i+1)%m][j]
                                + u[i][(j+n-1)%n] + u[i][(j+1)%n]
                                - 4*u[i][j] );
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        res[i][j] = rhs[i][j] - getNode((i+m-1)%m,j)->w - getNode((i+1)%m,j)->w
                              - getNode(i,(j+n-1)%n)->w - getNode(i,(j+1)%n)->w
                              + 4*getNode(i,j)->w - (1+tau)*u[i][j];
  }
#endif

  void setzero ( vector<vector<double> > &u )
  {
    int  m = u.size();
    int  n = u[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  u[i][j] = 0;
  }

  void add ( vector<vector<double> > &u, vector<vector<double> > &v )
  {
    int  m = u.size();
    int  n = u[0].size();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )  u[i][j] += v[i][j];
  }

  void cycle ( double tau, int l, int ct, int ng, int ngl )
  {
    for ( int i = 0; i < ng; ++i )
      gaussseidel ( w[l], frhs[l], tau, rig[l] );
    if ( l+1 < w.size() )
    {
      residuum ( w[l], frhs[l], res[l], tau, rig[l] );
      restrict ( res[l], frhs[l+1], res[l] );
      setzero ( w[l+1] );
      for ( int i = 0; i < ct; ++i )
        cycle ( tau, l+1, ct, ng, ngl );
      prolongate ( w[l+1], res[l] );
      add ( w[l], res[l] );
    }
    else ng = ngl;
    for ( int i = 0; i < ng; ++i )
      gaussseidel ( w[l], frhs[l], tau, rig[l] );
  }  

  void applyFlexure ( double dt = 1e99, int ct = 2, int ng = 1, int ngl = 100, int iter = 1 )
  {
    double  tau = this->tau/dt; 

    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        Node *p = getNode(i,j);
        p->h -= p->w;
        w[0][i][j] = p->w;
        frhs[0][i][j] = //p->b ? 1e99 : 
                              -rhoc*(p->h+p->g) + tau*p->w;
      }
#ifdef RIGIDITY
    if ( rig[0][0][0] != pow(getNode(0,0)->alpha,4)/4 )
    {
      for ( int i = 0; i < m; ++i )
        for ( int j = 0; j < n; ++j )  rig[0][i][j] = pow(getNode(i,j)->alpha,4)/4;
      for ( int l = 1; l < rig.size(); ++l )
      {
        restrict ( rig[l-1], rig[l], rig[l-1] );
        for ( int i = 0; i < rig[l].size(); ++i )
          for ( int j = 0; j < rig[l][i].size(); ++j )  rig[l][i][j] /= 16;
      }
    }
#else
    rig[0] = pow(alpha,4)/4; 
    for ( int l = 1; l < rig.size(); ++l )  rig[l] = rig[l-1]/16;
#endif
    for ( int it = 0; it < iter; ++it )
      cycle ( tau, 0, ct, ng, ngl );
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        Node *p = getNode(i,j);
        p->w = w[0][i][j];
        p->h += p->w;
      }
    computeWaterLevel();
  }
#endif

#ifdef PRECIP
  void solveRowCol ( int mode, int n, int i )
// mode = 0 -> row i; mode = 1 -> col i  
  {
    double  beta = (1-lc/ll)*(ll/lf-1);
    Matrix  dmat(-ld,0,0,-ld);
    for ( int j = 0; j < n; ++j )
    {
      Node *p = mode ? getNode(j,i) : getNode(i,j);
      double  f = exp(-p->l/refheight);
      Matrix  dmat(-ld,0,0,-ld);
      if ( !p->b )
      {
        if ( refslope > 0 )
        {
          double  s = p->l - getNode(p->d)->l;
          if ( p->d.isDiag(i,j) )  s *= sqrt(0.5);
          psi[j] = s > 0. ? exp(-refslope/s) : 0.;
          psi[j] = 1-f*evap*(1-psi[j]);
        }
        else  psi[j] = 1-f*evap;
        f *= beta/lc;
        diag[j] = Matrix(2*ld+1+1/lc,-f-(1-psi[j])/lf,-1/lc,2*ld+1+f+1/lf);
      }
      else
      {
        diag[j] = Matrix(1,0,0,1);
        dmat = Matrix(0,0,0,0);
        f *= lf/(ll-lf);
        rhs[j] = Vect(qin/(1+f),qin*f/(1+f));
      }
      if ( j > 0 )
        lower[j-1] = dmat;
      else
        right[0] = dmat;
      if ( j < n-1 )
        upper[j] = dmat;
      else
        bottom[0] = dmat;
    }
//    for ( int j = 0; j < n-1; ++j )  lower[j] = upper[j] = Matrix(-ld,0,0,-ld);
//    right[0] = bottom[0] = Matrix(-ld,0,0,-ld);
    for ( int j = 1; j < n-2; ++j )  right[j] = bottom[j] = Matrix(0,0,0,0);

    for ( int j = 0; j < n-1; ++j )
    {
      Matrix  f = diag[j].inv()*lower[j];
      diag[j+1] -= f*upper[j]; 
      rhs[j+1] -= f*rhs[j]; 

      if ( j < n-2 )
      {
        Matrix  g = diag[j].inv()*bottom[j];
        if ( j < n-3)
        {
          right[j+1] -= f*right[j]; 
          bottom[j+1] -= g*upper[j];
        }
        else
        {
          upper[j+1] -= f*right[j];
          lower[j+1] -= g*upper[j];
        }
        diag[n-1] -= g*right[j];
        rhs[n-1] -= g*rhs[j];
      }
    }
    rhs[n-1] = diag[n-1].inv()*rhs[n-1];    
    rhs[n-2] -= upper[n-2]*rhs[n-1];    
    rhs[n-2] = diag[n-2].inv()*rhs[n-2];    
    for ( int j = n-3; j >= 0; --j )
    {
      rhs[j] -= upper[j]*rhs[j+1];
      rhs[j] -= right[j]*rhs[n-1];
      rhs[j] = diag[j].inv()*rhs[j];
    }
    for ( int j = 0; j < n; ++j )
    {
      if ( rhs[j].u1 < 0 || rhs[j].u2 < 0 )
      {
        printf ( "At least one of the fluxes is negative.\n" 
                 "So far I thought that this cannot happen.\n"
                 "i = %i, j = %i, qv = %e, qc = %e\n",
                 i, j, rhs[j].u1, rhs[j].u2 );
        exit(1);
      }
      Node *p = mode ? getNode(j,i) : getNode(i,j);
      p->p = psi[j]*rhs[j].u2/lf;
#ifdef EXTPRECIP
      p->qv = rhs[j].u1;
      p->qc = rhs[j].u2;
      p->ptot = rhs[j].u2/lf;
#endif
    }
  }

  void computePrecipitation()
  {
    int  m = dir < 2 ? this->m : this->n; 
    int  n = dir < 2 ? this->n : this->m; 
    for ( int i = 0; i < n; ++i )
    {
      Node *p = dir < 2 ? getNode ( dir ? m-1 : 0, i ) : getNode ( i, dir > 2 ? m-1 : 0 );
      double  f = lf/(ll-lf)*exp(p->l/refheight);
      rhs[i] = Vect(qin/(1+f),qin*f/(1+f));
#ifdef EXTPRECIP
      p->qv = rhs[i].u1;
      p->qc = rhs[i].u2;
#endif
    }
    if ( dir%2 )
      for ( int i = m-2; i >= 0; --i )  solveRowCol ( dir > 2, n, i );
    else
      for ( int i = 1; i < m; ++i )  solveRowCol ( dir > 1, n, i );
  }
#endif // PRECIP

  void write ( string filename, string mask, int writeheaders = 1, int layer = -1 )
// Writes the lattice to a file. The mask defines which variables are written.
// For writeheaders = 1, all non-gridded variables are automatically included.
// For writeheaders = 0, only the raw grid data are written.
  {
    for ( int i = 0; i < mask.size(); ++i )  mask[i] = toupper(mask[i]);
    FILE *fp = fopen ( filename.c_str(), "wb" );
    if ( writeheaders )
    {
      fwrite ( &m, sizeof(short int), 1, fp );
      fwrite ( &n, sizeof(short int), 1, fp );
    }
    map<int,Keymapentry>::reverse_iterator it;
    for ( it = km.rbegin(); it != km.rend(); ++it )
    {
      Keymapentry e = it->second;
      if ( e.grid )
      {
        int  pos = mask.find(e.keystring);        
        if ( pos != std::string::npos )
        {
          for ( int i = pos; i < pos+strlen(e.keystring); ++i )  mask[i] = ' ';
//          printf ( "%s\n", mask.c_str() );
          if ( writeheaders )
          {
            if ( e.array > 1 && layer < 0 )
            {
              fwrite ( "ARR", sizeof(int), 1, fp );
              fwrite ( &e.array, sizeof(int), 1, fp );
            }
            fwrite ( &e.keyint, sizeof(int), 1, fp );
          }
          if ( strcmp(e.keystring,"S") )
            if ( e.size == 8 && e.fsize == 4 )
              for ( int i = 0; i < m; ++i )
                for ( int j = 0; j < n; ++j ) 
                {
                  float  buf = *(double*)((char*)&getNode(i,j)->h+e.offset);
                  if ( layer < 0 )   
                    fwrite ( &buf, e.fsize, e.array, fp );
                  else
                    fwrite ( &buf+layer, e.fsize, 1, fp );
                }
            else
              for ( int i = 0; i < m; ++i )
                for ( int j = 0; j < n; ++j ) 
                  if ( layer < 0 )   
                    fwrite ( (char*)&getNode(i,j)->h+e.offset, e.fsize, e.array, fp );
                  else
                    fwrite ( (char*)&getNode(i,j)->h+e.offset+layer*e.fsize, e.fsize, 1, fp );
          else
            for ( int i = 0; i < m; ++i )
              for ( int j = 0; j < n; ++j )
              {
	        float  s = getNode(i,j)->l - getNode(getNode(i,j)->d)->l;
                if ( getNode(i,j)->d.isDiag(i,j) )  s *= sqrt(0.5);
                fwrite ( &s, 4, 1, fp );
	  }
        }
      }
      else
        if ( writeheaders && !e.readonly )
        {
          if ( e.array > 1 )
          {
            fwrite ( "ARR", sizeof(int), 1, fp );
            fwrite ( &e.array, sizeof(int), 1, fp );
          }
          fwrite ( &e.keyint, sizeof(int), 1, fp );
          fwrite ( (char*)&t+e.offset, e.size, e.array, fp );
        }
    }
    fclose(fp);
  }

  void read ( FILE *fp, Keymapentry e, int arr = 1 )
  {
    if ( e.size == 8 && e.fsize == 4 )
      for ( int i = 0; i < m; ++i )
        for ( int j = 0; j < n; ++j ) 
          for ( int k = 0; k < arr; ++k )
          {
            float  buf;
            fread ( &buf, 4, arr, fp );
            double  bufd = buf;
            if ( k < e.array )
              memcpy ( (char*)&getNode(i,j)->h+e.offset+8*k, &bufd, 8 );   
          }
    else 
      if ( e.size > 0 )
        for ( int i = 0; i < m; ++i )
          for ( int j = 0; j < n; ++j )
            if ( arr > e.array )
            {
              fread ( (char*)&getNode(i,j)->h+e.offset, e.fsize, e.array, fp );
              fseek ( fp, e.fsize*(arr-e.array), SEEK_CUR );
            }
            else
              fread ( (char*)&getNode(i,j)->h+e.offset, e.fsize, arr, fp );
      else  
        fseek ( fp, m*n*e.fsize*arr, SEEK_CUR );    
  }

  void read ( string filename, string mask, int arr = 0 )
  {
    Keymapentry  e = getKey(mask.c_str());
    if ( e.grid )
    {
      FILE *fp = fopen ( filename.c_str(), "rb" );
      read ( fp, e, arr > 0 ? arr : e.array );
      fclose(fp);
    }
    else
      printf ( "The given key %s could not be assigned to any variable.\n", mask.c_str() ); 
  }

  double read ( string filename )
// Reads data from a file and returns time.
  {
    short int     m, n;
    int           buf;
    unsigned int  arr = 1, unrec = 0;
    FILE *fp = fopen ( filename.c_str(), "rb" );
    fread ( &m, sizeof(short int), 1, fp );
    fread ( &n, sizeof(short int), 1, fp );
    resize(m,n);
    while ( fread ( &buf, sizeof(int), 1, fp ) == 1 )
    {
      if ( km.find(buf) == km.end() )
      {
        char  cbuf[5];
        memcpy(cbuf,&buf,4);
        cbuf[4] = '\0';
        if ( strcmp(cbuf,"ARR") )
        {
          if ( !unrec )
          {
            unrec = ftell(fp)-4;
            printf ( "Unrecognized key %i at file position %i, looks like '%s'.\n",
                     buf, unrec, cbuf );
          } 
        }
        else  fread ( &arr, sizeof(int), 1, fp );
      }
      else
      {
        Keymapentry e = km[buf];
        if ( unrec )
        {
          unrec = 0;
          printf ( "Continued reading at file position %i with key '%s'.\n", 
                   ftell(fp)-4, e.keystring );
        }
        if ( e.grid )  read(fp,e,arr);
        else
          if ( arr > e.array )
          {
            fread ( (char*)&t+e.offset, e.size, e.array, fp );
            fseek ( fp, e.size*(arr-e.array), SEEK_CUR );
          }
          else
            fread ( (char*)&t+e.offset, e.size, arr, fp );
        arr = 1;
      }
    }
    fclose ( fp ); 
    return t;
  }

  void writeMatlabClass ( string classname = "", string path = "." )
  {
    int  arrkey;
    memcpy(&arrkey,"ARR",4);
    if ( classname.length() == 0 )  classname = "openlem"; 
    FILE *fp = fopen ( (path+"/"+classname+".m").c_str(), "w" );
    fprintf ( fp, "classdef openlem < handle\n\n  properties\n    m;\n    n;\n" );
    char    varname[32];
    string  varlist;
    map<int,Keymapentry>::reverse_iterator it;
    for ( it = km.rbegin(); it != km.rend(); ++it )
    {
      Keymapentry e = it->second;
      if ( strcmp(e.keystring,"D") )
        sprintf(varname,"    %s;\n",e.varname);
      else
        sprintf(varname,"    di;\n    dj;\n");
      if ( varlist.find(varname) == string::npos )
        varlist.append(varname);
    }
    fprintf ( fp, "%s", varlist.c_str() ); 
    fprintf ( fp, "  end\n\n  methods\n" 
                  "    function o = openlem(filename)\n" 
                  "      fid = fopen(filename);\n"
                  "      if fid < 0\n" 
                  "        o.m = 0;\n" 
                  "        o.n = 0;\n" 
                  "        return\n" 
                  "      end\n" 
                  "      arr = 1;\n"  
                  "      o.m = fread(fid,1,'int16=>int32');\n" 
                  "      o.n = fread(fid,1,'int16=>int32');\n" 
                  "      while true\n" 
                  "        if arr > 1\n" 
                  "          arrsize = [arr,o.n,o.m];\n" 
                  "          perm = [ 3 2 1 ];\n" 
                  "        else\n"
                  "          arrsize = [o.n,o.m];\n" 
                  "          perm = [ 2 1 ];\n" 
                  "        end\n"
                  "        key = fread(fid,1,'uint32=>uint32');\n" 
                  "        if numel(key) == 0\n" 
                  "          break\n" 
                  "        end\n" 
                  "        switch key\n" );
    for ( it = km.rbegin(); it != km.rend(); ++it )
    {
      Keymapentry e = it->second;
      fprintf ( fp, "          case %u\n", e.keyint ); 
      if ( !strcmp(e.keystring,"D" ) )
      {
        fprintf ( fp, "            tmp = fread(fid,[2,o.m*o.n],'uint16')+1;\n" 
                      "            o.di = reshape(tmp(1,:),[o.n,o.m])';\n" 
                      "            o.dj = reshape(tmp(2,:),[o.n,o.m])';\n" );
        continue;
      }
#if !defined(PRECIP) && !defined(DEFORM)
      if ( !strcmp(e.keystring,"Q" ) )
      {
        fprintf ( fp, "            o.q = fread(fid,[o.n,o.m],'uint32=>uint32')';\n" ); 
        continue;
      }
#endif
      if ( e.grid )
        fprintf ( fp, "            o.%s = permute(reshape(fread(fid,prod(arrsize),'", e.varname ); 
      else
        fprintf ( fp, "            o.%s = fread(fid,arr,'", e.varname ); 
      switch ( e.fsize )
      {
        case 1:  
          fprintf ( fp, "int8" );
          break;
        case 4:  
          fprintf ( fp, e.grid ? "float" : "uint32" );
          break;
        case 8:  
          fprintf ( fp, "double" );
          break;
      }
      fprintf ( fp, "=>" );
      switch ( e.size )
      {
        case 1:  
          fprintf ( fp, "int8" );
          break;
        case 8:  
          fprintf ( fp, "double" );
          break;
        default:  
          fprintf ( fp, e.grid ? "float" : "uint32" );
      }
      fprintf ( fp, "'" );
      if ( e.grid )  fprintf ( fp, "),arrsize),perm" );
      fprintf ( fp, ");\n" 
                    "            arr = 1;\n" );
    }
    fprintf ( fp, "          case %u\n", arrkey ); 
    fprintf ( fp, "            arr = fread(fid,1,'uint32');\n" 
                  "        end\n" 
                  "      end\n" 
                  "      fclose(fid);\n"
                  "    end\n\n"
                  "    function write(o,filename)\n" 
                  "      fid = fopen(filename,'w');\n"
                  "      fwrite(fid,o.m,'uint16');\n" 
                  "      fwrite(fid,o.n,'uint16');\n" );
    for ( it = km.rbegin(); it != km.rend(); ++it )
    {
      Keymapentry e = it->second;
      if ( e.readonly || strchr(e.keystring,'4') != NULL || strchr(e.keystring,'8') )
        continue;
      if ( !strcmp(e.keystring,"D" ) )
      {
        fprintf ( fp, "      if numel(o.di)\n" ); 
        fprintf ( fp, "        fwrite(fid,%u,'uint32');\n", e.keyint ); 
        fprintf ( fp, "        tmp = [ reshape(o.di',[1,o.m*o.n]); reshape(o.dj',[1,o.m*o.n]) ]-1;\n" 
                      "        fwrite(fid,tmp,'uint16')';\n"
                      "      end\n" );
        continue;
      }
      fprintf ( fp, "      if numel(o.%s)\n", e.varname );
      if ( e.grid )
        fprintf ( fp, "        arr = 1;\n" 
                      "        if numel(size(o.%s)) == 3\n"
                      "          arr = size(o.%s,3);\n"
                      "        end\n", e.varname, e.varname );
      else
        fprintf ( fp, "        arr = numel(o.%s);\n", e.varname );
      fprintf ( fp, "        if arr > 1\n"  
                    "          fwrite(fid,%u,'uint32');\n"
                    "          fwrite(fid,arr,'uint32');\n"
                    "        end\n", arrkey );
      fprintf ( fp, "        fwrite(fid,%u,'uint32');\n", e.keyint ); 
#ifndef PRECIP
      if ( !strcmp(e.keystring,"Q" ) )
      {
        fprintf ( fp, "        fwrite(fid,o.q','uint32')';\n"
                      "      end\n" ); 
        continue;
      }
#endif
//      fprintf ( fp, "        fwrite(fid,o.%s','", e.varname ); 
      if ( e.grid )
        fprintf ( fp, "        fwrite(fid,permute(o.%s,numel(size(o.%s)):-1:1),'" , e.varname, e.varname ); 
      else
        fprintf ( fp, "        fwrite(fid,o.%s,'", e.varname ); 
      switch ( e.fsize )
      {
        case 1:  
          fprintf ( fp, "int8" );
          break;
        case 4:  
          fprintf ( fp, e.grid ? "float" : "uint32" );
          break;
        case 8:  
          fprintf ( fp, "double" );
          break;
      }
      fprintf ( fp, "');\n" 
                    "      end\n" );
    }                  
    fprintf ( fp, "      fclose(fid);\n"
                  "    end\n\n"
                  "  end\n\n"
                  "end\n" );
    fclose ( fp );
  }

class Delta
{
  public:
  Point       first;
  Point       last;
  int         i;
  double      q;
  double      d;
  Grid<Node>  *g;

  Delta ( Grid<Node> *g, Point p )
  {
    this->g = g;
    first = last = p;
    i = 0;
    q = g->getNode(p)->q;
    d = g->offset[i].d/q;
  }

  void add ( Point p )
  {
    g->getNode(last)->d = p;    
  }

  Point point()
  {
    return Point((first.i+g->offset[i].p.i)%g->m,(first.j+g->offset[i].p.j)%g->n);
  }

  void step()
  {
    d = g->offset[++i].d/q;
  }
};

class Lake
{
  public:
  double         l;
  vector<Point>  points;
  vector<Point>  boundary;
  vector<Delta>  delta;
  map<double,set<unsigned int>>  rim;
  Point          outlet;
  Grid<Node>     *g;

  Node * getNode ( Point p )
  {
    return  g->getNode(p);  	  
  }

  Node * getNode ( vector<Point>::iterator p )
  {
    return  g->getNode(p);  	  
  }

  int add ( Point p )
  {
    Node  *pn = getNode(p);
    if ( pn->b || pn->m )  return 0;
    pn->m = 1;    
    points.push_back(p);
    int  n = 1;
    vector<Point>  neigh = g->getNeighbors(p);	    
    for ( vector<Point>::iterator it = neigh.begin(); it != neigh.end(); ++it )
    {
      pn = getNode(it);
      if ( !pn->m )
        rim[pn->l].insert(*(unsigned int*)&*it);
    }
    return n;
  }

  void addDonors ( Point p )
  {
// printf ( "Add %i %i\n", p.i, p.j );
    vector<Point>  neigh = g->getNeighbors(p);	    
    for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
      if ( g->getNode(u)->drainsTo(p) )  
      {
// printf ( "Neigh %i %i %f %f\n", u->i, u->j, g->getNode(u)->l, l );
        if ( g->getNode(u)->l == l )
        {
          points.push_back(*u);
          addDonors(*u);
        }
        else
          boundary.push_back(*u);
      }
  }

  public:
  Lake ( Grid<Node> *g )
  {
    this->g = g;
  }

  Lake ( Grid<Node> *g, Point p )
  {
    this->g = g;
    int     lake = 0;
    double  q = 0;
    Node    *pn = getNode(p);
    l = pn->l;
    while ( pn->l == l )
    {
      outlet = p;
      if ( pn->b )  break;
      pn = getNode(p=pn->d);
    }
    vector<Point>  neigh = g->getNeighbors(outlet);	    
    for ( vector<Point>::iterator p = neigh.begin(); p != neigh.end(); ++p )
      if ( (pn=getNode(p))->drainsTo(outlet) )  
        if ( pn->l == l )
        {
          points.push_back(*p);
          if ( pn->q > q )
          {
            q = pn->q;
            lake = 1;
          }
          addDonors(*p);
        }
        else
          if ( pn->q > q )
          {
            q = pn->q;
            lake = 0;
          }
    double  depth = 0;
    for ( vector<Point>::iterator p = points.begin(); p != points.end(); ++p )
      if ( getNode(p)->l-getNode(p)->h >= depth )  depth = getNode(p)->l-getNode(p)->h;
    if ( lake )
    {
      for ( vector<Point>::iterator p = points.begin(); p != points.end(); ++p )
      {
        getNode(p)->l = getNode(p)->q = 0;
        getNode(p)->m = 2;
      }
      int  nb = 0;
      for ( vector<Point>::iterator p = boundary.begin(); p != boundary.end(); ++p )
      {
        double         h = (pn=getNode(p))->h, maxs = 0.; 
        vector<Point>  neigh = g->getNeighbors(*p);	    
        for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
          if ( getNode(u)->m == 2 )
          {
            double  s = h-getNode(u)->h;
            if ( p->isDiag(*u) )  s *= sqrt(0.5);
            if ( s > maxs )
            {
              pn->d = *u;
	      maxs = s;
            }	
          }
        if ( getNode(pn->d)->q == 0 )  boundary[nb++] = pn->d;
        getNode(pn->d)->q += pn->q;
      }
      boundary.resize(nb);
// Now boundary should consist of all points that receive discharge from outside.
// new part
//      for ( int i = 0; i < nb; ++i )  delta.push_back(Delta(g,boundary[i]));
     
// Find biggest delta 
//      int  ilowest = 0;
//      for ( int i = 0; i < delta.size(); ++i )
//        if ( delta[i].d < delta[ilowest].d )  ilowest = i;
// 

// end of new part
      for ( vector<Point>::iterator p = points.begin(); p != points.end(); ++p )
      {
        pn = getNode(p);
        for ( vector<Point>::iterator b = boundary.begin(); b != boundary.end(); ++b )
          pn->l += getNode(b)->q/(g->computeDistSquare(*p,*b)+1.);
      }
      g->sort(points);
      points.push_back(outlet); 
      for ( vector<Point>::iterator p = points.begin()+1; p != points.end(); ++p )
      {
        getNode(p-1)->d = *p;
        if ( g->computeDistSquare(*(p-1),*p) > 2 )
          getNode(p)->trunk = *(p-1);
//        printf ( "(%i,%i) -> (%i,%i) %i (%i,%i)\n", (p-1)->i, (p-1)->j, p->i, p->j, g->computeDistSquare(*(p-1),*p), getNode(p)->trunk.i, getNode(p)->trunk.j );
      }  
//        printf ( "%i %i %f %f\n", p->j+1, p->i+1, g->getNode(p)->q, g->getNode(p)->l );
//      getNode(points.back())->d = outlet;
      for ( vector<Point>::iterator p = points.begin(); p != points.end(); ++p )
      {
        getNode(p)->l = l;
        getNode(p)->m = 3;
      }
      printf ( "Lake " );
    }
    else
    {
      for ( vector<Point>::iterator p = points.begin(); p != points.end(); ++p )
      {
        getNode(p)->l = getNode(p)->h;
        getNode(p)->m = 1;
      }
      printf ( "Floodplain " );
    }
    printf ( "outlet = (%i,%i), area = %i, depth = %f\n", outlet.i, outlet.j, points.size(), depth );
  }

/*
  void printRim()
  {
    for ( map<double,set<unsigned int>>::iterator it = rim.begin(); it != rim.end(); ++it )
      for ( set<unsigned int>::iterator q = it->second.begin(); q != it->second.end(); ++q )
      {
        Point xp = *(Point*)&*q;
        printf ( "%e -> (%i,%i) %i\n", it->first, xp.i, xp.j, *q );
      }
  }  
*/

  int fill ( Point p )
  {
    int    n;

    points.clear();
    rim.clear();
    if ( add(p) == 0 )  return 0;    
    n = 1;
    l = getNode(p)->l;
    set<unsigned int>  lowest;
    while ( 1 )
    {    
      double  minelev = rim.begin()->first;
      lowest = rim.begin()->second;
      if ( minelev < l )  break;
      rim.erase(rim.begin());
      for ( set<unsigned int>::iterator q = lowest.begin(); q != lowest.end(); ++q )
        n += add(*(Point*)&*q);
      if ( minelev == l )
        while ( rim.count(l) )
        {
          lowest = rim[l];
          rim.erase(l);
          for ( set<unsigned int>::iterator q = lowest.begin(); q != lowest.end(); ++q )
            n += add(*(Point*)&*q);
        }  
      l = minelev;
    }
    outlet = *(Point*)&*(lowest.begin());
    for ( int i = 0; i < points.size(); ++i )
      getNode(points[i])->l = l;
    return points.size();
  }

  void computeFlowDirection()
  {
    Node   *p;
    int  osize = points.size();
    for ( int i = 0; i < points.size(); ++i )  getNode(points[i])->m = 1;
    points.clear();
    points.push_back(outlet);
    getNode(outlet)->m = 0;
    boundary.clear();
    for ( int i = 0; i < points.size(); ++i )
    {
      vector<Point>  neigh = g->getNeighbors(points[i]);	    
      for ( vector<Point>::iterator p = neigh.begin(); p != neigh.end(); ++p )
        if ( getNode(p)->m )
	{
	  getNode(p)->m = 0;
          if ( !getNode(p)->b)  getNode(p)->d = points[i];
          getNode(p)->l = getNode(p)->h;
          if ( getNode(p)->l <= getNode(points[i])->l )
          {
            getNode(p)->l = getNode(points[i])->l;
            points.push_back(*p);
          }
          else
          {
//            printf ( "Silted up: (%i,%i)\n", p->i, p->j );	  
            boundary.push_back(*p);
          }
        } 
      if ( i == points.size()-1 && boundary.size() > 0 )
      {
        int  jmin = 0;
        for ( int j = 1; j < boundary.size(); ++j )
          if ( getNode(boundary[j])->h < getNode(boundary[jmin])->h )  jmin = j;
        points.push_back(boundary[jmin]);
//        printf ( "Added silted-up point (%i,%i) %i\n", boundary[jmin].i, boundary[jmin].j, points.size() ); 
        boundary[jmin] = boundary.back();
        boundary.pop_back();
      }
    }
    for ( int i = 0; i < points.size(); ++i )  getNode(points[i])->m = 0;
    if ( points.size() < osize )
      printf ( "Lost lake points %i -> %i\n", osize, points.size() );
  }

  void clear()
  {
    points.clear();
    boundary.clear();
  }
};

  void fillLakes()
// Fills all local depressions by determining a local water level.
// Flow directions of all points are computed.
  {
    Lake  lake(this);
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        getNode(i,j)->l = getNode(i,j)->h;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        getNode(i,j)->d = Point(i,j);	      
        if ( getNode(i,j)->b == 0 )
	  computeFlowDirection( Point(i,j) );
      }
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( getNode(i,j)->b == 0 && getNode(i,j)->d.equals(i,j) )
        {
          lake.fill(Point(i,j));
          lake.computeFlowDirection();
          lake.clear();
        }
    printf ( "Finished\n" );
    computeFlowDirection();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( getNode(i,j)->m )
        {
          printf ( "Unresolved marker: %i %i\n", i, j );	
          exit ( 0 );
        }
  }

  void findLakes()
  {
    Node  *pn;
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        getNode(i,j)->trunk = Point(i,j);
#ifndef ICE
    computeDischarge();
// Otherwise it has already been computed
#endif
    printf ( "WL2 %i\n", computeWaterLevel() );
    lakes.clear();
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
        if ( !(pn=getNode(i,j))->m && !pn->b && pn->l-pn->h > mld )
        {
          Lake  tmp(this,Point(i,j));
          if ( getNode(tmp.outlet)->m == 3 )  lakes.push_back(tmp);
        }
    for ( int i = 0; i < m; ++i )
      for ( int j = 0; j < n; ++j )
      {
        if ( (pn=getNode(i,j))->m != 3 )  pn->l = pn->h;
        pn->m = 0;
      }
//    for ( int i = 0; i < lakes.size(); ++i )
//      printf ( "Lake %i, outlet = (%i,%i), area = %i\n", i, lakes[i].outlet.i, lakes[i].outlet.j, lakes[i].points.size() );
  }

};

class Delta
{
  public:
  Point       first;
  Point       last;
  int         i;
  double      q;
  double      d;
  Grid<Node>  *g;

  Delta ( Grid<Node> *g, Point p )
  {
    this->g = g;
    first = last = p;
    i = 0;
    q = g->getNode(p)->q;
    d = g->offset[i].d/q;
  }

  void add ( Point p )
  {
    g->getNode(last)->d = p;    
  }

  Point point()
  {
    return Point((first.i+g->offset[i].p.i)%g->m,(first.j+g->offset[i].p.j)%g->n);
  }

  void step()
  {
    d = g->offset[++i].d/q;
  }
};

}
