#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/simulator.h>
#include <aspect/assembly.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


#include "inner_core_convection.cc"

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
        const FiniteElement<dim> &fe            = this->get_fe();
        const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
        const unsigned int n_q_points           = scratch.face_finite_element_values.n_quadrature_points;

        //assemble force terms for the matrix for all boundary faces
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
