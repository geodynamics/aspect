import numpy as np

# Conversion of Continental Material Flow-Law Prefactors from Experiments to ASPECT
# Upper crust: Wet Quartzite from Gleason and Tullis (1995)
# The total conversion factor consists of two parts.
#
# (1) Unit conversion from the experimental formulation
#     (micrometers -> meters and MPa -> Pa):
#
#     F_unit = 10^(-6*n - 6*p)
#
# (2) Conversion between the strain-rate definitions used in the
#     experiments and in ASPECT.
#
#     This factor depends on the form of the flow law implemented in
#     the code. For example, the dislocation creep module in ASPECT
#     computes the viscosity as
#
#         viscosity = 0.5 * A^(-1/n) * edot_ii^((1-n)/n)
#                     * exp((E + P*V) / (n*R*T))
#
#     where edot_ii is the second invariant of the strain-rate tensor.
#     Combining this formulation with the experimental definition based
#     on the shear strain rate leads to the conversion factor
#
#         F_s = 3^((n + 1) / 2) / 2
#
# where:
#   n = stress exponent
#   p = grain-size exponent
#
# The prefactor used by ASPECT is related to the experimental prefactor by:
#
#     A_ASPECT = F_unit * F_s * A_exp
#
# or equivalently:
#
#     A_ASPECT = 10^(-6*n - 6*p) * 3^((n + 1) / 2) / 2 * A_exp
For additional details, see supporting information section S1.3 of Dannberg et al., 2016, G3,
The importance of grain size to mantle dynamics and seismological observations, 
https://doi.org/10.1002/2017GC006944.
# Gleason and Tullis (1995) - Wet Quartzite - Dislocation Creep
pfct_gt95_wqtz_disl_rprt = 1.1e-4 # Original reported prefactor in units of MPa^-n micrometers^m
sexp_gt95_wqtz_disl_rprt = 4.0    # Original reported stress exponent
gexp_gt95_wqtz_disl_rprt = 0      # Original reported grain size exponent

# Convert prefactor to units of Pa^-n meters^m s^-1
pfct_gt95_wqtz_disl_unit = (pfct_gt95_wqtz_disl_rprt) * \
                           (1.e-6**sexp_gt95_wqtz_disl_rprt) * \
                           (1.e-6**gexp_gt95_wqtz_disl_rprt)

# Calculate prefactor scaling term
pfct_gt95_wqtz_disl_sfac = 3.**((sexp_gt95_wqtz_disl_rprt+1.)/2.) / 2


# Calculated modified prefactor term for ASPECT
pfct_gt95_wqtz_disl_aspt = pfct_gt95_wqtz_disl_unit * pfct_gt95_wqtz_disl_sfac

print('')
print('Gleason and Tullis (1995) Wet Quarzite Dislocation Prefactor Scaling = ', pfct_gt95_wqtz_disl_sfac)
print('Gleason and Tullis (1995) Wet Quarzite Dislocation Prefactor ASPECT  = ', pfct_gt95_wqtz_disl_aspt)
print('')


# Lower crust: Wet Anorthite from Rybacki et al. (2006)
# Water fugacity conversion
#
# This study uses water fugacity explicitly in the original flow-law
# formulation. Therefore, an additional conversion factor associated
# with the water fugacity term must be included.
#
# Assuming a constant water fugacity:
#
#     F_H2O = f_H2O^r
#
# where:
#   f_H2O = water fugacity, assumed to be a constant value of 1 MPa
#   r     = water fugacity exponent reported in the experimental study
#
# The unit-conversion factor (F_unit) and the strain-rate conversion
# factor (F_s) are identical to those described above.
#
# The ASPECT prefactor is therefore related to the experimental
# prefactor by:
#
#     A_ASPECT = F_H2O * F_unit * F_s * A_exp

# Rybacki et al. (2006) - Wet Anorthite - Dislocation Creep
pfct_ry06_want_disl_rprt = 10.**(0.2) # Original reported log(prefactor) in units of MPa^(-n-r) micrometers^m
sexp_ry06_want_disl_rprt = 3.0        # Original reported stress exponent
gexp_ry06_want_disl_rprt = 0          # Original reported grain size exponent
wfug_ry06_want_disl_rprt = 1          # Original reported water fugacity
constant_fug = 1 # 1MPa

# Convert prefactor to units of Pa^(-n-r) meters^m s^-1
pfct_ry06_want_disl_unit = pfct_ry06_want_disl_rprt * \
                           (1.e-6**(sexp_ry06_want_disl_rprt)) * \
                           (1.e-6**gexp_ry06_want_disl_rprt)

# Calculate water fugacity term
pfct_ry06_want_disl_ffac = constant_fug**wfug_ry06_want_disl_rprt

# Calculate prefactor scaling term
pfct_ry06_want_disl_sfac = 3.**((sexp_ry06_want_disl_rprt+1.)/2.) / 2.0

# Calculated modified prefactor term for ASPECT
pfct_ry06_want_disl_aspt = pfct_ry06_want_disl_ffac * pfct_ry06_want_disl_unit * pfct_ry06_want_disl_sfac

print('')
print('Rybacki et al. (2006) Wet Anorthite Dislocation Prefactor Scaling = ', pfct_ry06_want_disl_sfac)
print('Rybacki et al. (2006) Wet Anorthite Dislocation Prefactor ASPECT  = ', pfct_ry06_want_disl_aspt)
print('')