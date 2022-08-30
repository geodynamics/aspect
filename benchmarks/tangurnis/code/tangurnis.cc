/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/postprocess/interface.h>
#include <aspect/material_model/simple.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/heating_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/signaling_nan.h>



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
    class TanGurnis : public MaterialModel::Interface<dim>
    {
      public:

        TanGurnis();

        /**
         * @name Physical parameters used in the basic equations
         * @{
         */

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              /**
               * @name Physical parameters used in the basic equations
               * @{
               */

              const Point<dim> &pos = in.position[i];
              const double depth = 1.0 - pos[dim-1];
              const double temperature = sin(numbers::PI*pos(dim-1))*cos(numbers::PI*wavenumber*pos(0));

              out.viscosities[i] = ( Di==0.0 ? 1.0 : Di ) * exp( a * depth );
              out.densities[i] = ( Di==0.0 ? 1.0 : Di ) * (-1.0 * temperature ) * exp( Di/gamma * (depth) );
              out.specific_heat[i] = 1.0;
              out.thermal_conductivities[i] = 1.0;
              out.thermal_expansion_coefficients[i] = ( Di==0.0 ) ? 1.0 : Di;
              out.compressibilities[i] = numbers::signaling_nan<double>();;
            }
        }

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */


        double parameter_a() const;
        double parameter_wavenumber() const;
        double parameter_Di() const;
        double parameter_gamma() const;


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

    };

    template <int dim>
    TanGurnis<dim>::TanGurnis()
    {
      // a=0 for a constant viscosity or
      // a=2 for a depth-dependent viscosity (case tala_c)
      a=0;

      // for BA: Di=0; gamma=inf
      // however, we do not use gamma if Di=0, so the actual value does not matter

      //TALA:
      Di=0.5;
      gamma=1.0;

      wavenumber=1;
    }


    template <int dim>
    bool
    TanGurnis<dim>::
    is_compressible () const
    {
      return Di != 0.0;
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
          a                       = prm.get_double("a");
          Di                      = prm.get_double("Di");
          gamma                   = prm.get_double("gamma");
          wavenumber              = prm.get_double("wavenumber");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
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
      virtual
      double boundary_temperature (const types::boundary_id /*boundary_indicator*/,
                                   const Point<dim> &position) const
      {
        double wavenumber=1;
        return sin(numbers::PI*position(dim-1))*cos(numbers::PI*wavenumber*position(0));
      }

      virtual
      double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;

      virtual
      double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;
  };

  template <int dim>
  double
  TanGurnisBoundary<dim>::
  minimal_temperature (const std::set<types::boundary_id> &/*fixed_boundary_ids*/) const
  {
    return 0;
  }



  template <int dim>
  double
  TanGurnisBoundary<dim>::
  maximal_temperature (const std::set<types::boundary_id> &/*fixed_boundary_ids*/) const
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
  TanGurnisPostprocessor<dim>::execute (TableHandler &/*statistics*/)
  {
    AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                ExcNotImplemented());

    AssertThrow(Plugins::plugin_type_matches<const MaterialModel::TanGurnis<dim>>(this->get_material_model()),
                ExcMessage("tan gurnis postprocessor only works with tan gurnis material model"));

    const MaterialModel::TanGurnis<dim> &
    material_model = Plugins::get_plugin_as_type<const MaterialModel::TanGurnis<dim>>(this->get_material_model());

    double ref=1.0/this->get_triangulation().begin_active()->minimum_vertex_distance();
    std::ofstream f ((this->get_output_directory() + "vel_" +
                      Utilities::int_to_string(static_cast<unsigned int>(ref)) +
                      ".csv").c_str());
    f.precision (16);
    f << material_model.parameter_Di() << ' '
      << material_model.parameter_gamma() << ' '
      << material_model.parameter_wavenumber() << ' '
      << material_model.parameter_a();

    // pad the first line to the same number of columns as the data below to make MATLAB happy
    for (unsigned int i=4; i<7+this->get_heating_model_manager().get_active_heating_models().size(); ++i)
      f << " -1";

    f << std::endl;
    f << std::scientific;


    const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities+2);

    const unsigned int n_q_points =  quadrature_formula.size();
    FEValues<dim> fe_values (this->get_mapping(), this->get_fe(),  quadrature_formula,
                             update_JxW_values | update_values    |
                             update_gradients  | update_quadrature_points);

    MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
    MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());

    const std::list<std::unique_ptr<HeatingModel::Interface<dim>>> &heating_model_objects = this->get_heating_model_manager().get_active_heating_models();

    std::vector<HeatingModel::HeatingModelOutputs> heating_model_outputs (heating_model_objects.size(),
                                                                          HeatingModel::HeatingModelOutputs (n_q_points, this->n_compositional_fields()));

    typename DoFHandler<dim>::active_cell_iterator
    cell = this->get_dof_handler().begin_active(),
    endc = this->get_dof_handler().end();
    for (; cell != endc; ++cell)
      {
        fe_values.reinit (cell);
        in.reinit(fe_values,cell,this->introspection(),this->get_solution(),true);

        this->get_material_model().evaluate(in, out);

        if (this->get_parameters().formulation_temperature_equation ==
            Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
          {
            for (unsigned int q=0; q<n_q_points; ++q)
              {
                out.densities[q] = this->get_adiabatic_conditions().density(in.position[q]);
              }
          }

        unsigned int index = 0;
        for (typename std::list<std::unique_ptr<HeatingModel::Interface<dim>>>::const_iterator
             heating_model = heating_model_objects.begin();
             heating_model != heating_model_objects.end(); ++heating_model, ++index)
          {
            (*heating_model)->evaluate(in, out, heating_model_outputs[index]);
          }


        for (unsigned int q = 0; q < n_q_points; ++q)
          {
            f
                <<  fe_values.quadrature_point (q) (0)
                << ' ' << fe_values.quadrature_point (q) (1)
                << ' ' << in.velocity[q][0]
                << ' ' << in.velocity[q][1]
                << ' ' << fe_values.JxW (q)
                << ' ' << in.pressure[q]
                << ' ' << in.temperature[q];

            for (unsigned int i = 0; i < heating_model_objects.size(); ++i)
              f << ' ' << heating_model_outputs[i].heating_source_terms[q];

            f  << std::endl;
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
