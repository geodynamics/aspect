#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class Compressibility : public MaterialModel::Simple<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
          * Return true if the compressibility() function returns something that
          * is not zero.
          */
        virtual bool
        is_compressible () const;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    Compressibility<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);
      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          const double pressure = in.pressure[i];
          out.densities[i] = 1.0 + pressure;
          out.compressibilities[i] = 1.0 / (1. + pressure);
        }
    }

    template <int dim>
    bool
    Compressibility<dim>::
    is_compressible () const
    {
      return true;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Compressibility,
                                   "compressibility",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a non-zero compressibility.")
  }
}



#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class Compressibility : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    Compressibility<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      // be defensive about determining that what we think is the temperature
      // element, is it in fact
      Assert (this->get_fe().n_base_elements() == 3+(this->n_compositional_fields()>0 ? 1 : 0),
              ExcNotImplemented());
      const QGauss<dim-1> quadrature_formula (this->get_fe().base_element(2).degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_gradients      | update_values |
                                        update_q_points       | update_JxW_values);

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));
      std::vector<Tensor<1,dim> > velocity_values(quadrature_formula.size());

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      MaterialModel::MaterialModelInputs<dim> in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

      // compute the integral of the viscosity. since we're on a unit box,
      // this also is the average value
      double bottom_flux_integral = 0;
      double top_flux_integral = 0;
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<2*dim; ++f)
            if (cell->at_boundary(f)
                &&
                ((cell->face(f)->boundary_indicator() == 2)
                 ||
                 (cell->face(f)->boundary_indicator() == 3)))
              {
                fe_face_values.reinit (cell,f);
                fe_face_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                    in.temperature);
                fe_face_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                               in.pressure);
                fe_face_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                    in.strain_rate);
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  fe_face_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                      composition_values[c]);
                fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                    velocity_values);

                in.position = fe_face_values.get_quadrature_points();
                for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
                  {
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      in.composition[i][c] = composition_values[c][i];
                  }

                this->get_material_model().evaluate(in, out);

                if (cell->face(f)->boundary_indicator() == 2)
                  for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                    bottom_flux_integral += out.densities[q] * velocity_values[q][1] * fe_face_values.JxW(q);
                if (cell->face(f)->boundary_indicator() == 3)
                  for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                    top_flux_integral += out.densities[q] * velocity_values[q][1] * fe_face_values.JxW(q);
              }

      std::ostringstream screen_text1;
      std::ostringstream screen_text2;
      screen_text1.precision(4);
      screen_text2.precision(4);
      screen_text1 << Utilities::MPI::sum(bottom_flux_integral, this->get_mpi_communicator());
      screen_text2 << Utilities::MPI::sum(top_flux_integral, this->get_mpi_communicator());

      return std::pair<std::string, std::string> ("Top/bottom flux:",
                                                  screen_text2.str() + "/" + screen_text1.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Compressibility,
                                  "compressibility",
                                  "A postprocessor that computes some statistics about "
                                  "the mass fluxes.")
  }
}
