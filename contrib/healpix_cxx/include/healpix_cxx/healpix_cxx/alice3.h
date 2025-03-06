#ifndef ALICE3_H
#define ALICE3_H

#include "healpix_map.h"

void lic_main(const Healpix_Map<double> &Q, const Healpix_Map<double> &U, const Healpix_Map<double> &th,
  Healpix_Map<double> &hit, Healpix_Map<double> &tex, Healpix_Map<double> &mag,
  int steps, int kernel_steps, double step_radian, double polmin, double polmax);

#endif
