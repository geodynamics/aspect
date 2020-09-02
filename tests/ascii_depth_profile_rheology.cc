#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class AsciiDepthProfileRheology : public MaterialModel::Simple<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        static
        void declare_parameters (ParameterHandler &prm);

        void parse_parameters (ParameterHandler &prm);


      private:
        std::unique_ptr<Rheology::AsciiDepthProfile<dim> > depth_dependent_rheology;
        bool use_depth_dependent_viscosity;
    };
  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    AsciiDepthProfileRheology<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          if (use_depth_dependent_viscosity)
            out.viscosities[i] = depth_dependent_rheology->get_viscosity(in.position[i]);
        }
    }



    template <int dim>
    void
    AsciiDepthProfileRheology<dim>::declare_parameters (ParameterHandler &prm)
    {
      Simple<dim>::declare_parameters (prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          prm.declare_entry ("Use depth dependent viscosity", "false",
                             Patterns::Bool (),
                             "This parameter value determines if we want to use the layered depth dependent "
                             "rheology, which is input as an ascii data file.");

          // Depth-dependent viscosity parameters
          Rheology::AsciiDepthProfile<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Simple<dim>::parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {

          // Parse depth-dependent viscosity parameters
          if (use_depth_dependent_viscosity)
            {
              depth_dependent_rheology = std_cxx14::make_unique<Rheology::AsciiDepthProfile<dim>>();
              depth_dependent_rheology->initialize();
              depth_dependent_rheology->initialize_simulator (this->get_simulator());
              depth_dependent_rheology->parse_parameters(prm);
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(AsciiDepthProfileRheology,
                                   "ascii depth profile rheology",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a depth-dependent "
                                   "viscosity input from an ascii data file.")
  }
}
