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

      //assemble force terms for the matrix for all boundary faces
      const unsigned int dofs_per_cell = scratch.finite_element_values.get_fe().dofs_per_cell;
      if (cell->face(face_no)->at_boundary())
        {
          scratch.face_finite_element_values.reinit (cell, face_no);

          for (unsigned int q=0; q<scratch.face_finite_element_values.n_quadrature_points; ++q)
            {
              const double P = dynamic_cast<const MaterialModel::InnerCore<dim>&>
                               (this->get_material_model()).resistance_to_phase_change.value
                               (scratch.material_model_inputs.position[q]);

              // boundary term: P*u*n*v*n*JxW(q)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  data.local_matrix(i,j) += P *
                                            scratch.face_finite_element_values
                                              [introspection.extractors.velocities].value(i,q) *
                                            scratch.face_finite_element_values.normal_vector(q) *
                                            scratch.face_finite_element_values
                                              [introspection.extractors.velocities].value(j,q) *
                                            scratch.face_finite_element_values.normal_vector(q) *
                                            scratch.face_finite_element_values.JxW(q);
            }
        }
    }

};

template <int dim>
void set_assemblers_phase_boundary(const SimulatorAccess<dim> &simulator_access,
                                   internal::Assembly::AssemblerLists<dim> &assemblers,
                                   std::vector<dealii::std_cxx11::shared_ptr
                                   <internal::Assembly::Assemblers::AssemblerBase<dim> > > &assembler_objects)
{
  std::cout << "* set_assemblers()" << std::endl;

  AssertThrow (dynamic_cast<const MaterialModel::InnerCore<dim>*>(&simulator_access.get_material_model()) != 0,
               ExcMessage ("The phase boundary assembler can only be used with "
                           "the material model 'inner core'!"));

  PhaseBoundaryAssembler<dim> *phase_boundary_assembler = new PhaseBoundaryAssembler<dim>();
  assembler_objects.push_back(std_cxx11::shared_ptr<internal::Assembly::Assemblers::AssemblerBase<dim> >
                              (phase_boundary_assembler));

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

template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &signals)
{
std::cout << "* Connecting signals" << std::endl;
signals.set_assemblers.connect (&aspect::set_assemblers_phase_boundary<dim>);

}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                signal_connector<3>)
