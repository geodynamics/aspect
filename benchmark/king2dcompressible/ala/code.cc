

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
   * @endcode
   *
   * @ingroup Postprocessing
   */


  namespace MaterialModel
  {

    template <int dim>
    class Material : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:



        /**
        * Evaluate material properties.
        */
        virtual void evaluate(const MaterialModelInputs<dim> &in,
                              MaterialModelOutputs<dim> &out) const
        {

          for (unsigned int i=0; i < in.temperature.size(); ++i)
            {
              const Point<dim> position = in.position[i];
              const double temperature = in.temperature[i];
              const double pressure = in.pressure[i];

              const double depth = 1.0-position(dim-1);

              out.viscosities[i] = ((Di==0.0)?1.0:Di)/Ra*exp(-b*(temperature- this->get_adiabatic_conditions().temperature(position))+c*depth);

              out.specific_heat[i] = reference_specific_heat;
              out.thermal_conductivities[i] = 1.0;
              out.thermal_expansion_coefficients[i] = (Di==0.0)?1.0:Di;

              double rho = reference_rho;
              rho *= 1.0 - out.thermal_expansion_coefficients[i] * (temperature - this->get_adiabatic_conditions().temperature(position));
	      //rho*= exp(depth);

              out.densities[i] = rho;

              out.compressibilities[i] = 0.0;//(d==0) ? 1.0 : (Di/gamma / d);//reference_compressibility; // 1/rho drho/dp
              out.entropy_derivative_pressure[i] = 0.0;
              out.entropy_derivative_temperature[i] = 0.0;
              // Change in composition due to chemical reactions at the
              // given positions. The term reaction_terms[i][c] is the
              // change in compositional field c at point i.
              for (unsigned int c=0; c<in.composition[i].size(); ++c)
                out.reaction_terms[i][c] = 0.0;
            }
        }

        virtual bool is_compressible () const
        {
          return false;
        }



        virtual double reference_viscosity () const
        {
          return (Di==0.0?1.0:Di)/Ra;
        }


        virtual double reference_density () const
        {

          return 1.0;
        }


        virtual double reference_thermal_expansion_coefficient () const
        {
          return 1.0;
        }


        double reference_thermal_diffusivity () const
        {
          return 1/(reference_rho*reference_specific_heat);
        }


        double reference_cp () const
        {

          return 1;
        }

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

        /**
         * blankenbach parameters
         */
        double b, c;


        /**
         * The reference surface temperature
         */
        double reference_rho;

        /**
         * The constant viscosity
         */
        double Di, Ra;

        /**
         * The constant specific heat
         */
        double reference_specific_heat;

        /**
         * The constant compressibility.
         */
        double reference_compressibility;


    };



    template <int dim>
    void
    Material<dim>::declare_parameters (ParameterHandler &prm)
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
          prm.declare_entry ("Ra", "1e4",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Di", "0.0",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("b", "6.907755279",
                             Patterns::Double (0),
                             "");
          prm.declare_entry ("c", "0",
                             Patterns::Double (0),
                             "");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Material<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Tan Gurnis model");
        {
          reference_rho     = prm.get_double ("Reference density");
          //          reference_T = prm.get_double ("Reference temperature");
          //eta                   = prm.get_double ("Viscosity");
          b               = prm.get_double ("b");
          c               = prm.get_double ("c");
          Di               = prm.get_double ("Di");
          Ra               = prm.get_double ("Ra");

          reference_specific_heat = prm.get_double ("Reference specific heat");
          //thermal_alpha = prm.get_double ("Thermal expansion coefficient");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature;
      this->model_dependence.density = NonlinearDependence::none;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
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
                          const types::boundary_id                   boundary_indicator,
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
               const types::boundary_id                   boundary_indicator,
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
    /*
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
    */
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
    ASPECT_REGISTER_MATERIAL_MODEL(Material,
                                   "Material",
                                   "A simple compressible material model based on a benchmark"
                                   " from the paper of Tan/Gurnis (2007). This does not use the"
                                   " temperature equation, but has a hardcoded temperature.")
  }

  ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(TanGurnisBoundary,
                                             "Tan Gurnis",
                                             "A model for the Tan/Gurnis benchmark.")
}
