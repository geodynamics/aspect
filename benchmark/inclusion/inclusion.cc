#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

      /**
       * A material model that describes the "Pure shear/Inclusion" benchmark
       * of the paper cited in the documentation of the DuretzEtAl namespace.
       *
       * @note This benchmark only talks about the flow field, not about a
       * temperature field. All quantities related to the temperature are
       * therefore set to zero in the implementation of this class.
       *
       * @ingroup MaterialModels
       */
     template <int dim>
      class Inclusion : public MaterialModel::InterfaceCompatibility<dim>
      {
        public:
          /**
           * @name Physical parameters used in the basic equations
           * @{
           */
          virtual double viscosity (const double                  temperature,
                                    const double                  pressure,
                                    const std::vector<double>    &compositional_fields,
                                    const SymmetricTensor<2,dim> &strain_rate,
                                    const Point<dim>             &position) const;

          virtual double density (const double temperature,
                                  const double pressure,
                                  const std::vector<double> &compositional_fields,
                                  const Point<dim> &position) const;

          virtual double compressibility (const double temperature,
                                          const double pressure,
                                          const std::vector<double> &compositional_fields,
                                          const Point<dim> &position) const;

          virtual double specific_heat (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

          virtual double thermal_expansion_coefficient (const double      temperature,
                                                        const double      pressure,
                                                        const std::vector<double> &compositional_fields,
                                                        const Point<dim> &position) const;

          virtual double thermal_conductivity (const double temperature,
                                               const double pressure,
                                               const std::vector<double> &compositional_fields,
                                               const Point<dim> &position) const;
          /**
           * @}
           */

          /**
           * @name Qualitative properties one can ask a material model
           * @{
           */

          /**
           * Return true if the viscosity() function returns something that
           * may depend on the variable identifies by the argument.
           */
          virtual bool
          viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the density() function returns something that may
           * depend on the variable identifies by the argument.
           */
          virtual bool
          density_depends_on (const NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the compressibility() function returns something
           * that may depend on the variable identifies by the argument.
           *
           * This function must return false for all possible arguments if the
           * is_compressible() function returns false.
           */
          virtual bool
          compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the specific_heat() function returns something
           * that may depend on the variable identifies by the argument.
           */
          virtual bool
          specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

          /**
           * Return true if the thermal_conductivity() function returns
           * something that may depend on the variable identifies by the
           * argument.
           */
          virtual bool
          thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

          /**
           * Return whether the model is compressible or not.
           * Incompressibility does not necessarily imply that the density is
           * constant; rather, it may still depend on temperature or pressure.
           * In the current context, compressibility means whether we should
           * solve the contuity equation as $\nabla \cdot (\rho \mathbf u)=0$
           * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
           * (incompressible Stokes).
           */
          virtual bool is_compressible () const;
          /**
           * @}
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
           * @name Reference quantities
           * @{
           */
          virtual double reference_viscosity () const;

          virtual double reference_density () const;

          virtual double reference_thermal_expansion_coefficient () const;

//TODO: should we make this a virtual function as well? where is it used?
          double reference_thermal_diffusivity () const;

          double reference_cp () const;
          /**
           * @}
           */
          /**
           * Returns the viscosity value in the inclusion
           */
          double get_eta_B() const;

        private:
          /**
           * viscosity value in the inclusion
           */
          double eta_B;
      };
    
      template <int dim>
      double
      Inclusion<dim>::
      viscosity (const double,
                 const double,
                 const std::vector<double> &,       /*composition*/
                 const SymmetricTensor<2,dim> &,
                 const Point<dim> &p) const
      {
        const double r2 = (p(0)-1.0)*(p(0)-1.0) + (p(1)-1.0)*(p(1)-1.0);
        return (r2<0.2*0.2)? eta_B : 1.0;
      }


      template <int dim>
      double
      Inclusion<dim>::
      reference_viscosity () const
      {
        return 1;
      }

      template <int dim>
      double
      Inclusion<dim>::
      reference_density () const
      {
        return 0;
      }

      template <int dim>
      double
      Inclusion<dim>::
      reference_thermal_expansion_coefficient () const
      {
        return 0;
      }

      template <int dim>
      double
      Inclusion<dim>::
      specific_heat (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
      {
        return 0;
      }

      template <int dim>
      double
      Inclusion<dim>::
      reference_cp () const
      {
        return 0;
      }

      template <int dim>
      double
      Inclusion<dim>::
      thermal_conductivity (const double,
                            const double,
                            const std::vector<double> &, /*composition*/
                            const Point<dim> &) const
      {
        return 0;
      }

      template <int dim>
      double
      Inclusion<dim>::
      reference_thermal_diffusivity () const
      {
        return 0;
      }

      template <int dim>
      double
      Inclusion<dim>::
      density (const double,
               const double,
               const std::vector<double> &, /*composition*/
               const Point<dim> &p) const
      {
        return 0;
      }


      template <int dim>
      double
      Inclusion<dim>::
      thermal_expansion_coefficient (const double temperature,
                                     const double,
                                     const std::vector<double> &, /*composition*/
                                     const Point<dim> &) const
      {
        return 0;
      }


      template <int dim>
      double
      Inclusion<dim>::
      compressibility (const double,
                       const double,
                       const std::vector<double> &, /*composition*/
                       const Point<dim> &) const
      {
        return 0.0;
      }



      template <int dim>
      bool
      Inclusion<dim>::
      viscosity_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }


      template <int dim>
      bool
      Inclusion<dim>::
      density_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      Inclusion<dim>::
      compressibility_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      Inclusion<dim>::
      specific_heat_depends_on (const NonlinearDependence::Dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      Inclusion<dim>::
      thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
      {
        return false;
      }

      template <int dim>
      bool
      Inclusion<dim>::
      is_compressible () const
      {
        return false;
      }

      template <int dim>
      void
      Inclusion<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Inclusion");
          {
            prm.declare_entry ("Viscosity jump", "1e3",
                               Patterns::Double (0),
                               "Viscosity in the Inclusion.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      Inclusion<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Material model");
        {
          prm.enter_subsection("Inclusion");
          {
            eta_B = prm.get_double ("Viscosity jump");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      double
      Inclusion<dim>::get_eta_B() const
      {
        return eta_B;
      }

    

  }
}








// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
 ASPECT_REGISTER_MATERIAL_MODEL(Inclusion,
                                     "Inclusion",
                                     "A material model that corresponds to the 'Inclusion' benchmark "
                                     "defined in Duretz et al., G-Cubed, 2011.")
 
  }
}
