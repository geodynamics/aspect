#include <aspect/material_model/simple.h>
#include <aspect/heating_model/interface.h>
#include <deal.II/base/parsed_function.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/simulator.h>
#include <aspect/assembly.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

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

        /**
         * A function object representing resistance to phase change at the
         * inner core boundary as a function the position (and, optionally,
         * the model time).
         */
        Functions::ParsedFunction<dim> resistance_to_phase_change;
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
      // product of density and specific heat a constant.
      Simple<dim>::evaluate(in, out);

      // We want the right-hand side of the momentum equation to be (- Ra T gravity) and
      // density * cp to be 1
      for (unsigned int q=0; q < in.position.size(); ++q)
        {
          out.densities[q] = - out.thermal_expansion_coefficients[q] * in.temperature[q];
          if (std::abs(out.densities[q]) > 0.0)
            out.specific_heat[q] /= out.densities[q];
        }
    }


    template <int dim>
    void
    InnerCore<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        resistance_to_phase_change.set_time (this->get_time() / year_in_seconds);
      else
        resistance_to_phase_change.set_time (this->get_time());
    }


    template <int dim>
    void
    InnerCore<dim>::declare_parameters (ParameterHandler &prm)
    {
      Simple<dim>::declare_parameters (prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Inner core");
        {
          prm.enter_subsection("Phase change resistance function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    InnerCore<dim>::parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Inner core");
        {
          prm.enter_subsection("Phase change resistance function");
          try
            {
              resistance_to_phase_change.parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Phase boundary model.Function'\n"
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
  /**
   * A new assembler class that implements boundary conditions for the
   * normal stress and the normal velocity that take into account the
   * rate of phase change (melting/freezing) at the inner-outer core
   * boundary. The model is based on Deguen, Alboussiere, and Cardin
   * (2013), Thermal convection in Earth’s inner core with phase change
   * at its boundary. GJI, 194, 1310-133.
   *
   * The mechanical boundary conditions for the inner core are
   * tangential stress-free and continuity of the normal stress at the
   * inner-outer core boundary. For the non-dimensional equations, that
   * means that we can define a 'phase change number' $\mathcal{P}$ so
   * that the normal stress at the boundary is $-\mathcal{P} u_r$ with
   * the radial velocity $u_r$. This number characterizes the resistance
   * to phase change at the boundary, with $\mathcal{P}\rightarrow\infty$
   * corresponding to infinitely slow melting/freezing (free slip
   * boundary), and $\mathcal{P}\rightarrow0$ corresponding to
   * instantaneous melting/freezing (zero normal stress, open boundary).
   *
   * In the weak form, this results in boundary conditions of the form
   * of a surface integral:
   * $$\int_S \mathcal{P} (\mathbf u \cdot \mathbf n) (\mathbf v \cdot \mathbf n) dS,$$,
   * with the normal vector $\mathbf n$.
   *
   * The function value of $\mathcal{P}$ is taken from the inner core
   * material model.
   */
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
        const FiniteElement<dim> &fe = this->get_fe();
        const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
        const unsigned int n_q_points           = scratch.face_finite_element_values.n_quadrature_points;

        // assemble force terms for the matrix for all boundary faces
        if (cell->face(face_no)->at_boundary())
          {
            scratch.face_finite_element_values.reinit (cell, face_no);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double P = dynamic_cast<const MaterialModel::InnerCore<dim>&>
                                 (this->get_material_model()).resistance_to_phase_change
                                 .value(scratch.material_model_inputs.position[q]);

                for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
                  {
                    if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                      {
                        scratch.phi_u[i_stokes] = scratch.face_finite_element_values[introspection
                                                                                     .extractors.velocities].value(i, q);
                        ++i_stokes;
                      }
                    ++i;
                  }

                const Tensor<1,dim> normal_vector = scratch.face_finite_element_values.normal_vector(q);
                const double JxW = scratch.face_finite_element_values.JxW(q);

                // boundary term: P*u*n*v*n*JxW(q)
                for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                  for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                    data.local_matrix(i,j) += P *
                                              scratch.phi_u[i] *
                                              normal_vector *
                                              scratch.phi_u[j] *
                                              normal_vector *
                                              JxW;
              }
          }
      }
  };

  template <int dim>
  void set_assemblers_phase_boundary(const SimulatorAccess<dim> &simulator_access,
                                     internal::Assembly::AssemblerLists<dim> &assemblers,
                                     std::vector<dealii::std_cxx11::shared_ptr<
                                     internal::Assembly::Assemblers::AssemblerBase<dim> > > &assembler_objects)
  {
    AssertThrow (dynamic_cast<const MaterialModel::InnerCore<dim>*>
                 (&simulator_access.get_material_model()) != 0,
                 ExcMessage ("The phase boundary assembler can only be used with the "
                             "material model 'inner core material'!"));

    std_cxx11::shared_ptr<PhaseBoundaryAssembler<dim> > phase_boundary_assembler
    (new PhaseBoundaryAssembler<dim>());
    assembler_objects.push_back (phase_boundary_assembler);

    // add the terms for phase change boundary conditions
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
                                   "'Simple' model, but has a constant $\rho c_p$, "
                                   "and implements a function that characterizes the "
                                   "resistance to melting/freezing at the inner core "
                                   "boundary.")
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ConstantCoreHeating,
                                  "constant core heating",
                                  "Implementation of a model in which the heating "
                                  "rate is constant.")
  }
}
