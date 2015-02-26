

#include <aspect/postprocess/interface.h>
#include <aspect/material_model/simple.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  using namespace dealii;



  /**
   * This benchmark is from the article
   * @code
   *  @article{tan2007compressible,
   *    title={Compressible thermochemical convection and application to lower mantle structures},
   *    author={Tan, E. and Gurnis, M.},
   *    journal={JOURNAL OF GEOPHYSICAL RESEARCH-ALL SERIES-},
   *    volume={112},
   *    number={B6},
   *    pages={6304},
   *    year={2007},
   *    publisher={AGU AMERICAN GEOPHYSICAL UNION}
   *    }
   * @endcode
   *
   * @ingroup Postprocessing
   */



  namespace MaterialModel
  {

    template <int dim>
    class TanGurnis : public MaterialModel::InterfaceCompatibility<dim>
    {
      public:

        TanGurnis();

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
         * Return true if the viscosity() function returns something that may
         * depend on the variable identifies by the argument.
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
         * Return true if the specific_heat() function returns something that
         * may depend on the variable identifies by the argument.
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
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
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

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

//TODO: should we make this a virtual function as well? where is it used?
        double reference_thermal_diffusivity () const;

        double reference_cp () const;

        double parameter_a() const;
        double parameter_wavenumber() const;
        double parameter_Di() const;
        double parameter_gamma() const;

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

      private:

        double a;
        double wavenumber;
        double Di;
        double gamma;

        double reference_rho;
        double reference_T;
        double eta;
        double thermal_alpha;
        double reference_specific_heat;

        /**
         * The thermal conductivity.
         */
        double k_value;
    };

    template <int dim>
    TanGurnis<dim>::TanGurnis()
    {
      a=0; // 0 or 2

      //BA:
      //Di=0;gamma=10000; //=inf
      //EBA:
      //Di=0.5;gamma=inf;
      //TALA:
      Di=0.5;
      gamma=1.0;

      wavenumber=1;
    }


    template <int dim>
    double
    TanGurnis<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &pos) const
    {
      const double depth = 1.0-pos(dim-1);
      return exp(a*depth);
    }


    template <int dim>
    double
    TanGurnis<dim>::
    reference_viscosity () const
    {
      return 1.0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_density () const
    {
      return 1.0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 1.0;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 1250;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_cp () const
    {
      return 1250;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 2e-5;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }



    template <int dim>
    double
    TanGurnis<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &pos) const
    {
      const double depth = 1.0-pos(dim-1);
      const double temperature = sin(numbers::PI*pos(dim-1))*cos(numbers::PI*wavenumber*pos(0));
      return -1.0*temperature*exp(Di/gamma*(depth));
    }



    template <int dim>
    double
    TanGurnis<dim>::
    thermal_expansion_coefficient (const double,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return thermal_alpha;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &pos) const
    {
      return Di/gamma / density(temperature, pressure, compositional_fields, pos);
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    density_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }



    template <int dim>
    bool
    TanGurnis<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_a() const
    {
      return a;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_wavenumber() const
    {
      return wavenumber;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_Di() const
    {
      return Di;
    }



    template <int dim>
    double
    TanGurnis<dim>::
    parameter_gamma() const
    {
      return gamma;
    }



    template <int dim>
    void
    TanGurnis<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tan Gurnis model");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity. Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.declare_entry ("a", "0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Di", "0.5",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("gamma", "1",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("wavenumber", "1",
                             Patterns::Double (0),
                             "");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    TanGurnis<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tan Gurnis model");
        {
          reference_rho     = prm.get_double ("Reference density");
          reference_T = prm.get_double ("Reference temperature");
          eta                   = prm.get_double ("Viscosity");
          k_value               = prm.get_double ("Thermal conductivity");
          reference_specific_heat = prm.get_double ("Reference specific heat");
          thermal_alpha = prm.get_double ("Thermal expansion coefficient");
          a = prm.get_double("a");
          Di = prm.get_double("Di");
          gamma = prm.get_double("gamma");
          wavenumber = prm.get_double("wavenumber");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }


  /**
   * A class that implements a temperature boundary condition for the
   * tan/gurnis benchmark in a box geometry.
   *
   * @ingroup BoundaryTemperatures
   */
  template <int dim>
  class TanGurnisBoundary : public BoundaryTemperature::Interface<dim>
  {
    public:
      /**
       * Return the temperature that is to hold at a particular location on
       * the boundary of the domain. This function returns constant
       * temperatures at the left and right boundaries.
       *
       * @param geometry_model The geometry model that describes the domain.
       * This may be used to determine whether the boundary temperature
       * model is implemented for this geometry.
       * @param boundary_indicator The boundary indicator of the part of the
       * boundary of the domain on which the point is located at which we
       * are requesting the temperature.
       * @param location The location of the point at which we ask for the
       * temperature.
       */
      virtual
      double temperature (const GeometryModel::Interface<dim> &geometry_model,
                          const unsigned int                   boundary_indicator,
                          const Point<dim>                    &location) const;

      /**
       * Return the minimal the temperature on that part of the boundary on
       * which Dirichlet conditions are posed.
       *
       * This value is used in computing dimensionless numbers such as the
       * Nusselt number indicating heat flux.
       */
      virtual
      double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;

      /**
       * Return the maximal the temperature on that part of the boundary on
       * which Dirichlet conditions are posed.
       *
       * This value is used in computing dimensionless numbers such as the
       * Nusselt number indicating heat flux.
       */
      virtual
      double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;
  };

  template <int dim>
  double
  TanGurnisBoundary<dim>::
  temperature (const GeometryModel::Interface<dim> &geometry_model,
               const unsigned int                   boundary_indicator,
               const Point<dim>                    &location) const
  {
    // verify that the geometry is in fact a box since only
    // for this geometry do we know for sure what boundary indicators it
    // uses and what they mean
    Assert (dynamic_cast<const GeometryModel::Box<dim>*>(&geometry_model)
            != 0,
            ExcMessage ("This boundary model is only implemented if the geometry is "
                        "in fact a box."));

    double wavenumber=1;
    return sin(numbers::PI*location(dim-1))*cos(numbers::PI*wavenumber*location(0));
  }


  template <int dim>
  double
  TanGurnisBoundary<dim>::
  minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
  {
    return 0;
  }



  template <int dim>
  double
  TanGurnisBoundary<dim>::
  maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
  {
    return 1;
  }



  /*
   * A postprocessor that evaluates the accuracy of the solution of the
   * aspect::MaterialModel::TanGurnis material model.
   *
   * The implementation writes out the solution to be read in by a matlab
   * script.
   */
  template <int dim>
  class TanGurnisPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      virtual
      std::pair<std::string,std::string>
      execute (TableHandler &statistics);
  };

  template <int dim>
  std::pair<std::string,std::string>
  TanGurnisPostprocessor<dim>::execute (TableHandler &statistics)
  {
    AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                ExcNotImplemented());

    const MaterialModel::TanGurnis<dim> *
    material_model = dynamic_cast<const MaterialModel::TanGurnis<dim> *>(&this->get_material_model());

    AssertThrow(material_model!=NULL, ExcMessage("tan gurnis postprocessor only works with tan gurnis material model"));

    double ref=1.0/this->get_triangulation().begin_active()->minimum_vertex_distance();
    std::ofstream f ((this->get_output_directory() + "vel_" +
                      Utilities::int_to_string(static_cast<unsigned int>(ref)) +
                      ".csv").c_str());
    f.precision (16);
    f << material_model->parameter_Di() << ' '
      << material_model->parameter_gamma() << ' '
      << material_model->parameter_wavenumber() << ' '
      << material_model->parameter_a() << ' '
      << " -1 -1 -1" << std::endl; //pad to 7 values, so matlab is happy

    f << std::scientific;


    const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

    const unsigned int n_q_points =  quadrature_formula.size();
    FEValues<dim> fe_values (this->get_mapping(), this->get_fe(),  quadrature_formula,
                             update_JxW_values | update_values | update_quadrature_points);

    std::vector<Tensor<1, dim> > velocity_values (quadrature_formula.size());
    std::vector<double>         temperature_values (quadrature_formula.size());
    std::vector<double>         pressure_values (quadrature_formula.size());

    typename DoFHandler<dim>::active_cell_iterator
    cell = this->get_dof_handler().begin_active(),
    endc = this->get_dof_handler().end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit (cell);
        fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(), velocity_values);
        fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(), pressure_values);
        fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(), temperature_values);

        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            f
                <<  fe_values.quadrature_point (q) (0)
                << ' ' << fe_values.quadrature_point (q) (1)
                << ' ' << velocity_values[q][0]
                << ' ' << velocity_values[q][1]
                << ' ' << fe_values.JxW (q)
                << ' ' << pressure_values[q]
                << ' ' << temperature_values[q]
                << std::endl;
          }
      }

    return std::make_pair("writing:", "output.csv");
  }

}



// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_POSTPROCESSOR(TanGurnisPostprocessor,
                                "Tan Gurnis error",
                                "A postprocessor that compares the solution of the benchmarks from "
                                "the Tan/Gurnis (2007) paper with the one computed by ASPECT "
                                "by outputing data that is compared using a matlab script.")

  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(TanGurnis,
                                   "Tan Gurnis",
                                   "A simple compressible material model based on a benchmark"
                                   " from the paper of Tan/Gurnis (2007). This does not use the"
                                   " temperature equation, but has a hardcoded temperature.")
  }

  ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(TanGurnisBoundary,
                                             "Tan Gurnis",
                                             "A model for the Tan/Gurnis benchmark.")
}
