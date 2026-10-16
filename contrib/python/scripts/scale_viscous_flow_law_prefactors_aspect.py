# This script scales the prefactors terms from viscous flow laws commonly used
# in geodynamic modeling in a manner that is internally consistent with ASPECT's
# conventions following Dannberg et al. (2017, https://doi.org/10.1002/2017GC006944).

import numpy as np

# Gleason and Tullis (1995) - Wet Quartzite - Dislocation Creep

pfct_gt95_wqtz_disl_rprt = 1.8e-8 # Original reported prefactor in units of MPa^-n micrometers^m

sexp_gt95_wqtz_disl_rprt = 4.0    # Original reported stress exponent

gexp_gt95_wqtz_disl_rprt = 0      # Original reported grain size exponent

# Convert prefactor to units of Pa^-n meters^m s^-1

pfct_gt95_wqtz_disl_unit = (pfct_gt95_wqtz_disl_rprt) * \
                           (1.e-6**sexp_gt95_wqtz_disl_rprt) * \
                           (1.e-6**gexp_gt95_wqtz_disl_rprt)

# Calculate prefactor scaling term

pfct_gt95_wqtz_disl_sfac = 2.**(sexp_gt95_wqtz_disl_rprt-1.) * 3.**((sexp_gt95_wqtz_disl_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_gt95_wqtz_disl_aspt = pfct_gt95_wqtz_disl_unit * pfct_gt95_wqtz_disl_sfac

print('')

print('Gleason and Tullis (1995) Quarzite w/ Melt Dislocation Prefactor Scaling = ', pfct_gt95_wqtz_disl_sfac)

print('Gleason and Tullis (1995) Quarzite w/ Melt Dislocation Prefactor ASPECT  = ', pfct_gt95_wqtz_disl_aspt)

print('')


# Gleason and Tullis (1995) - Dry Quartzite - Dislocation Creep

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

print('Gleason and Tullis (1995) Dry Quarzite Dislocation Prefactor Scaling = ', pfct_gt95_dqtz_disl_sfac)

print('Gleason and Tullis (1995) Dry Quarzite Dislocation Prefactor ASPECT  = ', pfct_gt95_dqtz_disl_aspt)

print('')


# Rybacki et al. (2000) - Wet Anorthite - Dislocation Creep

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

# Rybacki et al. (2006) - Wet Anorthite - Dislocation Creep

pfct_ry06_want_disl_rprt = 10.**(0.2) # Original reported log(prefactor) in units of MPa^(-n-r) micrometers^m

sexp_ry06_want_disl_rprt = 3.0        # Original reported stress exponent

gexp_ry06_want_disl_rprt = 0          # Original reported grain size exponent

wfug_ry06_want_disl_rprt = 1          # Original reported water fugacity

# Convert prefactor to units of Pa^(-n-r) meters^m s^-1

pfct_ry06_want_disl_unit = (pfct_ry06_want_disl_rprt) * \
                           (1.e-6**(sexp_ry06_want_disl_rprt + wfug_ry06_want_disl_rprt)) * \
                           (1.e-6**gexp_ry06_want_disl_rprt)

# Calculate prefactor scaling term

pfct_ry06_want_disl_sfac = 2.**(sexp_ry06_want_disl_rprt-1.) * 3.**((sexp_ry06_want_disl_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_ry06_want_disl_aspt = pfct_ry06_want_disl_unit * pfct_ry06_want_disl_sfac

print('')

print('Rybacki et al. (2006) Wet Anorthite Dislocation Prefactor Scaling = ', pfct_ry06_want_disl_sfac)

print('Rybacki et al. (2006) Wet Anorthite Dislocation Prefactor ASPECT  = ', pfct_ry06_want_disl_aspt)

print('')

# Hirth and Kohlstedt (2004) - Dry Olivine - Dislocation Creep

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

# Hirth and Kohlstedt (2004) - Dry Olivine - Diffusion Creep

pfct_hk04_doli_diff_rprt = 1.5e9  # Original reported prefactor in units of MPa^-n micrometers^m

sexp_hk04_doli_diff_rprt = 1.0    # Original reported stress exponent

gexp_hk04_doli_diff_rprt = 3.0    # Original reported grain size exponent

# Convert prefactor to units of Pa^-n meters^m s^-1

pfct_hk04_doli_diff_unit = (pfct_hk04_doli_diff_rprt) * \
                           (1.e-6**sexp_hk04_doli_diff_rprt) * \
                           (1.e-6**gexp_hk04_doli_diff_rprt)

# Calculate prefactor scaling term

pfct_hk04_doli_diff_sfac = 2.**(sexp_hk04_doli_diff_rprt-1.) * 3.**((sexp_hk04_doli_diff_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_hk04_doli_diff_aspt = pfct_hk04_doli_diff_unit * pfct_hk04_doli_diff_sfac

print('')

print('Hirth Kohlsted (2004) Dry Olivine Diffusion Prefactor Scaling = ', pfct_hk04_doli_diff_sfac)

print('Hirth Kohlsted (2004) Dry Olivine Diffusion Prefactor ASPECT  = ', pfct_hk04_doli_diff_aspt)

print('')

# Hirth and Kohlstedt (2004) - Wet Olivine - Dislocation Creep

pfct_hk04_woli_disl_rprt = 90. # Original reported prefactor in units of MPa^-n micrometers^m

sexp_hk04_woli_disl_rprt = 3.5    # Original reported stress exponent

gexp_hk04_woli_disl_rprt = 0      # Original reported grain size exponent

wfug_hk04_woli_disl_rprt = 0    # Original reported water fugacity

# Convert prefactor to units of Pa^(-n-r) meters^m s^-1

pfct_hk04_woli_disl_unit = (pfct_hk04_woli_disl_rprt) * \
                           (1.e-6**(sexp_hk04_woli_disl_rprt + wfug_hk04_woli_disl_rprt)) * \
                           (1.e-6**gexp_hk04_woli_disl_rprt)

# Calculate prefactor scaling term

pfct_hk04_woli_disl_sfac = 2.**(sexp_hk04_woli_disl_rprt-1.) * 3.**((sexp_hk04_woli_disl_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_hk04_woli_disl_aspt = pfct_hk04_woli_disl_unit * pfct_hk04_woli_disl_sfac

print('')

print('Hirth Kohlsted (2004) Wet Olivine Dislocation Prefactor Scaling = ', pfct_hk04_woli_disl_sfac)

print('Hirth Kohlsted (2004) Wet Olivine Dislocation Prefactor ASPECT  = ', pfct_hk04_woli_disl_aspt)

print('')

# Hirth and Kohlstedt (2004) - Wet Olivine - Diffusion Creep

pfct_hk04_woli_diff_rprt = 1.e6 # Original reported prefactor in units of MPa^-n micrometers^m

sexp_hk04_woli_diff_rprt = 1.0  # Original reported stress exponent

gexp_hk04_woli_diff_rprt = 3    # Original reported grain size exponent

wfug_hk04_woli_diff_rprt = 0.   # Original reported water fugacity

# Convert prefactor to units of Pa^(-n-r) meters^m s^-1

pfct_hk04_woli_diff_unit = (pfct_hk04_woli_diff_rprt) * \
                           (1.e-6**(sexp_hk04_woli_diff_rprt + wfug_hk04_woli_diff_rprt)) * \
                           (1.e-6**gexp_hk04_woli_diff_rprt)

# Calculate prefactor scaling term

pfct_hk04_woli_diff_sfac = 2.**(sexp_hk04_woli_diff_rprt-1.) * 3.**((sexp_hk04_woli_diff_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_hk04_woli_diff_aspt = pfct_hk04_woli_diff_unit * pfct_hk04_woli_diff_sfac

print('')

print('Hirth Kohlsted (2004) Wet Olivine Diffusion Prefactor Scaling = ', pfct_hk04_woli_diff_sfac)

print('Hirth Kohlsted (2004) Wet Olivine Diffusion Prefactor ASPECT  = ', pfct_hk04_woli_diff_aspt)

print('')



# Wilks and Carter (1990) - Gabbro - Dislocation Creep

pfct_wc90_gabbro_disl_rprt = 8.e-3 # Original reported prefactor in units of MPa^-n micrometers^m

sexp_wc90_gabbro_disl_rprt = 3.1  # Original reported stress exponent

gexp_wc90_gabbro_disl_rprt = 0    # Original reported grain size exponent

# Convert prefactor to units of Pa^(-n-r) meters^m s^-1

pfct_hk04_woli_diff_unit = (pfct_wc90_gabbro_disl_rprt) * \
                           (1.e-6**(sexp_wc90_gabbro_disl_rprt)) * \
                           (1.e-6**gexp_wc90_gabbro_disl_rprt)

# Calculate prefactor scaling term

pfct_wc90_gabbro_disl_sfac = 2.**(sexp_wc90_gabbro_disl_rprt-1.) * 3.**((sexp_wc90_gabbro_disl_rprt+1.)/2.)

# Calculated modified prefactor term for ASPECT

pfct_wc90_gabbro_disl_aspt = pfct_hk04_woli_diff_unit * pfct_wc90_gabbro_disl_sfac

print('')

print('Wilson Carter (1990) Gabbro Dislocation Prefactor Scaling = ', pfct_wc90_gabbro_disl_sfac)

print('Wilson Carter (1990) Gabbro Dislocation Prefactor ASPECT  = ', pfct_wc90_gabbro_disl_aspt)

print('')
