/*
  Copyright (C) 2014 - 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/

#ifndef _aspect_material_model_viscoelastic_plastic_h
#define _aspect_material_model_viscoelastic_plastic_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
    * Additional output fields for the elastic shear modulus to be added to
    * the MaterialModel::MaterialModelOutputs structure and filled in the
    * MaterialModel::Interface::evaluate() function.
    */
    template <int dim>
    class ElasticAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        ElasticAdditionalOutputs(const unsigned int n_points);

        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * Elastic shear moduli at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> elastic_shear_moduli;
    };

    /**
     * A material model that combines non-linear plasticity with a simple linear
     * viscoelastic material behavior. The model is incompressible. Note that
     * this material model is based heavily on and combines functionality from
     * the following material models: DiffusionDislocation, DruckerPrager,
     * ViscoPlastic and Viscoelastic.
     *
     * Viscous stress is limited by plastic deformation, which follows
     * a Drucker Prager yield criterion:
     *  $\sigma_y = C*\cos(\phi) + P*\sin(\phi)$  (2D)
     * or in 3D
     *  $\sigma_y = \frac{6*C*\cos(\phi) + 2*P*\sin(\phi)}{\sqrt(3)*(3+\sin(\phi))}$
     * where
     *  $\sigma_y$ is the yield stress, $C$ is cohesion, $phi$ is the angle
     *  of internal friction and $P$ is pressure.
     * If the viscous stress ($2*v*{\varepsilon}_{ii})$) exceeds the yield
     * stress ($\sigma_{y}$), the viscosity is rescaled back to the yield
     * surface: $v_{y}=\sigma_{y}/(2.*{\varepsilon}_{ii})$
     * This form of plasticity is commonly used in geodynamic models.
     * See, for example, Thieulot, C. (2011), PEPI 188, pp. 47-68. "
     *
     * The viscoelastic rheology takes into account the elastic shear
     * strength (e.g., shear modulus), while the tensile and volumetric
     * strength (e.g., Young's and bulk modulus) are not considered. The model
     * is incompressible and allows specifying an arbitrary number of
     * compositional fields, where each field represents a different rock type
     * or component of the viscoelastic stress tensor. The stress tensor in 2D
     * and 3D, respectively, contains 3 or 6 components. The compositional fields
     * representing these components must be named and listed in a very specific
     * format, which is designed to minimize mislabeling stress tensor components
     * as distinct 'compositional rock types' (or vice versa). For 2D models,
     * the first three compositional fields must be labeled stress_xx, stress_yy
     * and stress_xy. In 3D, the first six compositional fields must be labeled
     * stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz.
     *
     * Combining this viscoelasticity implementation with non-linear viscous flow
     * and plasticity produces a constitutive relationship commonly referred to
     * as partial elastoviscoplastic (e.g., pEVP) in the geodynamics community.
     * While extensively discussed and applied within the geodynamics literature,
     * notable references include:
     * Moresi et al. (2003), J. Comp. Phys., v. 184, p. 476-497.
     * Gerya and Yuen (2007), Phys. Earth. Planet. Inter., v. 163, p. 83-105.
     * Gerya (2010), Introduction to Numerical Geodynamic Modeling.
     * Kaus (2010), Tectonophysics, v. 484, p. 36-47.
     * Choi et al. (2013), J. Geophys. Res., v. 118, p. 2429-2444.
     * Keller et al. (2013), Geophys. J. Int., v. 195, p. 1406-1442.
     *
     * The overview below directly follows Moresi et al. (2003) eqns. 23-32.
     * However, an important distinction between this material model and
     * the studies above is the use of compositional fields, rather than
     * tracers, to track individual components of the viscoelastic stress
     * tensor. The material model will be udpated when an option to track
     * and calculate viscoelastic stresses with tracers is implemented.
     *
     * Moresi et al. (2003) begins (eqn. 23) by writing the deviatoric
     * rate of deformation ($\hat{D}$) as the sum of elastic
     * ($\hat{D_{e}}$) and viscous ($\hat{D_{v}}$) components:
     * $\hat{D} = \hat{D_{e}} + \hat{D_{v}}$.
     * These terms further decompose into
     * $\hat{D_{v}} = \frac{\tau}{2\eta}$ and
     * $\hat{D_{e}} = \frac{\overset{\triangledown}{\tau}}{2\mu}$, where
     * $\tau$ is the viscous deviatoric stress, $\eta$ is the shear viscosity,
     * $\mu$ is the shear modulus and $\overset{\triangledown}{\tau}$ is the
     * Jaumann corotational stress rate. This later term (eqn. 24) contains the
     * time derivative of the deviatoric stress ($\dot{\tau}$) and terms that
     * account for material spin (e.g., rotation) due to advection:
     * $\overset{\triangledown}{\tau} = \dot{\tau} + {\tau}W -W\tau$.
     * Above, $W$ is the material spin tensor (eqn. 25):
     * $W_{ij} = \frac{1}{2} \left (\frac{\partial V_{i}}{\partial x_{j}} -
     * \frac{\partial V_{j}}{\partial x_{i}} \right )$.
     *
     * The Jaumann stress-rate can also be approximated using terms from the time
     * at the previous time step ($t$) and current time step ($t + \Delta t_^{e}$):
     * $\smash[t]{\overset{\triangledown}{\tau}}^{t + \Delta t^{e}} \approx
     * \frac{\tau^{t + \Delta t^{e} - \tau^{t}}}{\Delta t^{e}} -
     * W^{t}\tau^{t} + \tau^{t}W^{t}$.
     * In this material model, the size of the time step above ($\\Delta t^{e}$)
     * can be specified as the numerical time step size or an independent fixed time
     * step. If the latter case is a selected, the user has an option to apply a
     * stress averaging scheme to account for the differences between the numerical
     * and fixed elastic time step (eqn. 32). If one selects to use a fixed elastic time
     * step throughout the model run, this can still be achieved by using CFL and
     * maximum time step values that restrict the numerical time step to a specific time.
     *
     * The formulation above allows rewriting the total rate of deformation (eqn. 29) as
     * $\tau^{t + \Delta t^{e}} = \eta_{eff} \left (
     * 2\\hat{D}^{t + \\triangle t^{e}} + \\frac{\\tau^{t}}{\\mu \\Delta t^{e}} +
     * \\frac{W^{t}\\tau^{t} - \\tau^{t}W^{t}}{\\mu}  \\right )$.
     *
     * The effective viscosity (eqn. 28) is a function of the viscosity ($\eta$),
     * elastic time step size ($\Delta t^{e}$) and shear relaxation time
     * ($ \alpha = \frac{\eta}{\mu} $):
     * $\eta_{eff} = \eta \frac{\Delta t^{e}}{\Delta t^{e} + \alpha}$
     * The magnitude of the shear modulus thus controls how much the effective
     * viscosity is reduced relative to the initial viscosity.
     *
     * Elastic effects are introduced into the governing Stokes equations through
     * an elastic force term (eqn. 30) using stresses from the previous time step:
     * $F^{e,t} = -\frac{\eta_{eff}}{\mu \Delta t^{e}} \tau^{t}$.
     * This force term is added onto the right-hand side force vector in the
     * system of equations.
     *
     * The value of each compositional field representing distinct
     * rock types at a point is interpreted to be a volume fraction of that
     * rock type. If the sum of the compositional field volume fractions is
     * less than one, then the remainder of the volume is assumed to be
     * 'background material'.
     *
     * Several model parameters (densities, elastic shear moduli,
     * thermal expansivities, plasticity parameters, viscosity terms, etc) can
     * be defined per-compositional field. For each material parameter the
     * user supplies a comma delimited list of length N+1, where N is the
     * number of compositional fields. The additional field corresponds to
     * the value for background material. They should be ordered
     * ``background, composition1, composition2...''. However, the first 3 (2D)
     * or 6 (3D) composition fields correspond to components of the elastic
     * stress tensor and their material values will not contribute to the volume
     * fractions. If a single value is given, then all the compositional fields
     * are given that value. Other lengths of lists are not allowed. For a given
     * compositional field the material parameters are treated as constant,
     * except density, which varies linearly with temperature according to the
     * thermal expansivity.
     *
     * When more than one compositional field is present at a point, they are
     * averaged arithmetically. An exception is viscosity, which may be averaged
     * arithmetically, harmonically, geometrically, or by selecting the
     * viscosity of the composition with the greatest volume fraction.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class ViscoelasticPlastic : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * This model is not compressible, so this returns false.
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;
        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;


      private:

        /**
         * Reference temperature for thermal expansion. All components use
         * the same value.
         */
        double reference_temperature;

        /**
         * Minimum strain-rate (second invariant) to stabilize strain-rate dependent viscosity
         */
        double minimum_strain_rate;

        /**
         * Strain-rate invariant used on the first iteration of the first time step
         * before the velocity is calculated.
         */
        double reference_strain_rate;

        /**
         * Lower cutoff for effective viscosity
        */
        double minimum_viscosity;

        /**
         * Upper cutoff for effective viscosity
         */
        double maximum_viscosity;

        /**
         * Reference viscosity for scaling the pressure in the Stokes equation
         */
        double input_reference_viscosity;

        /**
         * Enumeration for selecting which averaging scheme to use. Select
         * between harmonic, arithmetic, geometric, and maximum_composition.
         * The max composition scheme simply uses the parameter of whichever
         * field has the highest volume fraction.
         */
        enum AveragingScheme
        {
          harmonic,
          arithmetic,
          geometric,
          maximum_composition
        };


        AveragingScheme viscosity_averaging;

        double average_value (const std::vector<double> &volume_fractions,
                              const std::vector<double> &parameter_values,
                              const enum AveragingScheme &average_type) const;

        /**
         * Vector for field densities, read from parameter file.
         */
        std::vector<double> densities;

        /**
         * Vector for field thermal expansivities, read from parameter file.
         */
        std::vector<double> thermal_expansivities;

        /**
         * Vector for field thermal conductivities, read from parameter file.
         */
        std::vector<double> thermal_conductivities;

        /**
         * Vector for field specific heats, read from parameter file.
         */
        std::vector<double> specific_heats;


        /**
         * Vector for field linear (fixed) viscosities, read from parameter file.
         */
        std::vector<double> linear_viscosities;

        /**
         * Vector for field angles of internal friction, read from parameter file.
         */
        std::vector<double> angles_internal_friction;

        /**
         * Vector for cohesions, read from parameter file.
         */
        std::vector<double> cohesions;

        /**
         * Vector for field elastic shear moduli, read from parameter file.
         */
        std::vector<double> elastic_shear_moduli;

        /**
         * Bool indicating whether to use a fixed material time scale in the
         * viscoelastic rheology for all time steps (if true) or to use the
         * actual (variable) advection time step of the model (if false). Read
         * from parameter file.
         */
        bool use_fixed_elastic_time_step;

        /**
         * Bool indicating whether to use a stress averaging scheme to account
         * for differences between the numerical and fixed elastic time step
         * (if true). When set to false, the viscoelastic stresses are not
         * modified to account for differences between the viscoelastic time
         * step and the numerical time step. Read from parameter file.
         */
        bool use_stress_averaging;

        /**
         * Double for fixed elastic time step value, read from parameter file
         */
        double fixed_elastic_time_step;

    };

  }
}

#endif
