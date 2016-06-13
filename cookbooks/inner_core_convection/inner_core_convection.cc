#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>
#include <aspect/traction_boundary_conditions/interface.h>
#include <aspect/global.h>
#include <aspect/heating_model/interface.h>
#include <aspect/simulator.h>
#include <aspect/assembly.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/parsed_function.h>

#include <typeinfo>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class InnerCore : public MaterialModel::Simple<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    InnerCore<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // First, we use the material descriptions of the 'simple' material model to fill all of the material
      // model outputs. Below, we will then overwrite selected properties (the specific heat) to make the
      // product of density and specifi heat a constant.
      Simple<dim>::evaluate(in, out);

      // we want density * cp to be 1
      for (unsigned int q=0; q < in.position.size(); ++q)
        out.specific_heat[q] /= out.densities[q];
    }

  }
}

namespace aspect
{
  namespace TractionBoundaryConditions
  {
    using namespace dealii;

    /**
     * A class that implements traction boundary conditions based on a
     * functional description provided in the input file.
     *
     * @ingroup TractionBoundaryConditionsModels
     */
    template <int dim>
    class NormalFunction : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        NormalFunction ();

        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        virtual
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        virtual
        void
        update ();

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * A function object representing the components of the traction.
         */
        Functions::ParsedFunction<dim> boundary_traction_function;
    };
  }
}

namespace aspect
{
  namespace TractionBoundaryConditions
  {
    template <int dim>
    NormalFunction<dim>::NormalFunction ()
      :
      boundary_traction_function (1)
    {}



    template <int dim>
    Tensor<1,dim>
    NormalFunction<dim>::
    boundary_traction (const types::boundary_id,
                       const Point<dim> &p,
                       const Tensor<1,dim> &normal_vector) const
    {
      return boundary_traction_function.value(p) * normal_vector;;
    }


    template <int dim>
    void
    NormalFunction<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_traction_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_traction_function.set_time (this->get_time());
    }


    template <int dim>
    void
    NormalFunction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Normal function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    NormalFunction<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Normal function");
        try
          {
            boundary_traction_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary traction model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a constant radiogenic heating rate.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class ConstantCoreHeating : public Interface<dim>
    {
      public:
        /**
         * Return the heating terms. For the current class, this
         * function obviously simply returns a constant value.
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

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
        double radiogenic_heating_rate;
    };
  }
}

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ConstantCoreHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &/*material_model_inputs*/,
              const MaterialModel::MaterialModelOutputs<dim> &/*material_model_outputs*/,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // return a constant value
          heating_model_outputs.heating_source_terms[q] = radiogenic_heating_rate;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ConstantCoreHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Constant core heating");
        {
          prm.declare_entry ("Radiogenic heating rate", "0e0",
                             Patterns::Double (0),
                             "The specific rate of heating due to radioactive decay (or other bulk sources "
                             "you may want to describe). This parameter corresponds to the variable "
                             "$H$ in the temperature equation stated in the manual, and the heating "
                             "term is $\rho H$. Units: W/kg.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ConstantCoreHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Constant core heating");
        {
          radiogenic_heating_rate    = prm.get_double ("Radiogenic heating rate");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

namespace aspect
{
    template <int dim>
    class PhaseBoundaryAssembler :
      public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>,
            public SimulatorAccess<dim>
    {
      public:

        virtual
        void
        phase_change_boundary_conditions (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          const unsigned int                                    face_no,
                                          internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                          internal::Assembly::CopyData::StokesSystem<dim>      &data) const
        {
          const Introspection<dim> &introspection = this->introspection();

          // see if any of the faces are traction boundaries for which
          // we need to assemble force terms for the right hand side
          const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
          if (this->get_traction_boundary_conditions()
              .find (cell->face(face_no)->boundary_id())
              !=
              this->get_traction_boundary_conditions().end())
            {
              scratch.face_finite_element_values.reinit (cell, face_no);

              for (unsigned int q=0; q<scratch.face_finite_element_values.n_quadrature_points; ++q)
                {
                  const Tensor<1,dim> traction
                    = this->get_traction_boundary_conditions().find(
                        cell->face(face_no)->boundary_id()
                      )->second
                      ->boundary_traction (cell->face(face_no)->boundary_id(),
                                           scratch.face_finite_element_values.quadrature_point(q),
                                           scratch.face_finite_element_values.normal_vector(q));

                  // boundary term: P*u*n*v*n*JxW(q);
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                  data.local_matrix(i,j) += scratch.face_finite_element_values[introspection.extractors.velocities].value(i,q) *
                                            traction *
                                            scratch.face_finite_element_values[introspection.extractors.velocities].value(j,q) *
                                            scratch.face_finite_element_values.normal_vector(q) *
                                            scratch.face_finite_element_values.JxW(q);
                }
            }
        }

    };

    template <int dim>
    void set_assemblers_phase_boundary(const SimulatorAccess<dim> &,
                                       internal::Assembly::AssemblerLists<dim> &assemblers,
                                       std::vector<dealii::std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> > > &assembler_objects)
    {
      std::cout << "* set_assemblers()" << std::endl;

      PhaseBoundaryAssembler<dim> *phase_boundary_assembler = new PhaseBoundaryAssembler<dim>();
      assembler_objects.push_back(std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> >(phase_boundary_assembler));

      // add the terms for traction boundary conditions
      assemblers.local_assemble_stokes_system_on_boundary_face
      .connect (std_cxx11::bind(&PhaseBoundaryAssembler<dim>::phase_change_boundary_conditions,
                                std_cxx11::cref (*phase_boundary_assembler),
                                std_cxx11::_1,
                                std_cxx11::_2,
                                // discard pressure_scaling,
                                // discard rebuild_stokes_matrix,
                                std_cxx11::_5,
                                std_cxx11::_6));


    }


}

template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
  std::cout << "* Connecting signals" << std::endl;
  signals.set_assemblers.connect (&aspect::set_assemblers_phase_boundary<dim>);

}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(InnerCore,
                                   "inner core material",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a constant $\rho c_p$.")
  }

  namespace TractionBoundaryConditions
  {
    ASPECT_REGISTER_TRACTION_BOUNDARY_CONDITIONS(NormalFunction,
                                                 "normal function",
                                                 "Implementation of a model in which the boundary "
                                                 "traction is given in terms of an explicit formula "
                                                 "that is elaborated in the parameters in section "
                                                 "``Boundary traction model|Function''. "
                                                 "\n\n"
                                                 "The formula you describe in the mentioned "
                                                 "section constains only a normal component of the traction."
                                                 "for each of the $d$ components of the traction vector. "
                                                 "The formula is interpreted as having units Pa.")
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ConstantCoreHeating,
                                  "constant core heating",
                                  "Implementation of a model in which the heating "
                                  "rate is constant.")
  }
}
