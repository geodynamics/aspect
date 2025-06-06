# This script scales the prefactors terms from viscous flow laws commonly used
# in geodynamic modeling in a manner that is internally consistent with ASPECT's
# conventions following Dannberg et al. (2017, https://doi.org/10.1002/2017GC006944).

import numpy as np


# Gleason and Tullis (1995), Tectonophysics, v.247 p.1-23.
# "A flow law for dislocation creep of quartz aggregates determined with the molten salt cell"
#  https://doi.org/10.1016/0040-1951(95)00011-B

pfct_gt95_dqtz_disl_rprt = 1.1e-4 # Original reported prefactor in units of MPa^-n micrometers^m

sexp_gt95_dqtz_disl_rprt = 4.0    # Original reported stress exponent

gexp_gt95_dqtz_disl_rprt = 0      # Original reported grain size exponent

# Convert prefactor to units of Pa^-n meters^m s^-1

pfct_gt95_dqtz_disl_unit = (pfct_gt95_dqtz_disl_rprt) * \
                           (1.e-6**sexp_gt95_dqtz_disl_rprt) * \
                           (1.e-6**gexp_gt95_dqtz_disl_rprt)

# Calculate prefactor scaling term

pfct_gt95_dqtz_disl_sfac = 2.**(sexp_gt95_dqtz_disl_rprt-1.) * 3.**((sexp_gt95_dqtz_disl_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_gt95_dqtz_disl_aspt = pfct_gt95_dqtz_disl_unit * pfct_gt95_dqtz_disl_sfac

print('')

print('Gleason and Tullis (1995) Quarzite Dislocation Prefactor Scaling = ', pfct_gt95_dqtz_disl_sfac)

print('Gleason and Tullis (1995) Quarzite Dislocation Prefactor ASPECT  = ', pfct_gt95_dqtz_disl_aspt)

print('')


# Rybacki and Dresen (2000), J. Geophys. Res., v.111(B3).
# "Dislocation and diffusion creep of synthetic anorthite aggregates".
# https://doi.org/10.1029/2000JB900223

pfct_ry00_want_disl_rprt = 10.**(2.6) # Original reported log(prefactor) in units of MPa^(-n-r) micrometers^m

sexp_ry00_want_disl_rprt = 3.0        # Original reported stress exponent

gexp_ry00_want_disl_rprt = 0          # Original reported grain size exponent

wfug_ry00_want_disl_rprt = 0          # Original reported water fugacity

# Convert prefactor to units of Pa^(-n-r) meters^m s^-1

pfct_ry00_want_disl_unit = (pfct_ry00_want_disl_rprt) * \
                           (1.e-6**(sexp_ry00_want_disl_rprt + wfug_ry00_want_disl_rprt)) * \
                           (1.e-6**gexp_ry00_want_disl_rprt)

# Calculate prefactor scaling term

pfct_ry00_want_disl_sfac = 2.**(sexp_ry00_want_disl_rprt-1.) * 3.**((sexp_ry00_want_disl_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_ry00_want_disl_aspt = pfct_ry00_want_disl_unit * pfct_ry00_want_disl_sfac

print('')

print('Rybacki et al. (2000) Wet Anorthite Dislocation Prefactor Scaling = ', pfct_ry00_want_disl_sfac)

print('Rybacki et al. (2000) Wet Anorthite Dislocation Prefactor ASPECT  = ', pfct_ry00_want_disl_aspt)

print('')


# Hirth & Kohlstedt (2004),  Geophys. Monogr. Am. Geophys. Soc., v.138, p.83-105.
# "Rheology of the upper mantle and the mantle wedge:a view from the experimentalists".
# https://doi.org/10.1029/138GM06

pfct_hk04_doli_disl_rprt = 1.10e5 # Original reported prefactor in units of MPa^-n micrometers^m

sexp_hk04_doli_disl_rprt = 3.5    # Original reported stress exponent

gexp_hk04_doli_disl_rprt = 0      # Original reported grain size exponent

# Convert prefactor to units of Pa^-n meters^m s^-1

pfct_hk04_doli_disl_unit = (pfct_hk04_doli_disl_rprt) * \
                           (1.e-6**sexp_hk04_doli_disl_rprt) * \
                           (1.e-6**gexp_hk04_doli_disl_rprt)

# Calculate prefactor scaling term

pfct_hk04_doli_disl_sfac = 2.**(sexp_hk04_doli_disl_rprt-1.) * 3.**((sexp_hk04_doli_disl_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_hk04_doli_disl_aspt = pfct_hk04_doli_disl_unit * pfct_hk04_doli_disl_sfac

print('')

print('Hirth Kohlsted (2004) Dry Olivine Dislocation Prefactor Scaling = ', pfct_hk04_doli_disl_sfac)

print('Hirth Kohlsted (2004) Dry Olivine Dislocation Prefactor ASPECT  = ', pfct_hk04_doli_disl_aspt)

print('')
