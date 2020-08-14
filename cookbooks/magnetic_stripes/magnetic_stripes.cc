#include <aspect/material_model/composition_reaction.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class MagneticStripes : public MaterialModel::CompositionReaction<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        static void declare_parameters (ParameterHandler &prm);
        virtual void parse_parameters (ParameterHandler &prm);

      private:
       std::vector<double> reversal_times;
    };


    template <int dim>
    void
    MagneticStripes<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      CompositionReaction<dim>::evaluate(in, out);

      int magnetic_orientation = 1;
      for (unsigned int i=0; i < reversal_times.size(); ++i)
        if (this->get_time() < reversal_times[i])
          {
            magnetic_orientation = -2 * (i % 2) + 1;
            break;
          }

      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          const double depth = this->get_geometry_model().depth(in.position[i]);
          const double reaction_depth = 7000.0;

          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              if (depth < reaction_depth && in.velocity[i][1] > std::abs(in.velocity[i][0])) 
                out.reaction_terms[i][c] = -in.composition[i][0] + magnetic_orientation;
              else
                out.reaction_terms[i][c] = 0.0;
            }
        }
    }

    template <int dim>
    void
    MagneticStripes<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Magnetic stripes model");
        {
          CompositionReaction<dim>::declare_parameters (prm);

          prm.declare_entry ("Reversal times", "5.0, 5.0, 2.0, 2.0, 2.092, 2.419, 2.419",
                             Patterns::List(Patterns::Double(0)),
                             "Reversal times of the magnatic field."
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MagneticStripes<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Magnetic stripes model");
        {
          CompositionReaction<dim>::parse_parameters (prm);

          reversal_times = Utilities::string_to_double (Utilities::split_string_list(prm.get ("Reversal times")));
          std::sort (reversal_times.begin(), reversal_times.end());

          if (this->get_parameters().convert_to_years == true)
            for (unsigned int i=0; i<reversal_times.size(); ++i)
              reversal_times[i] *= constants::year_in_seconds;
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
    ASPECT_REGISTER_MATERIAL_MODEL(MagneticStripes,
                                   "magnetic stripes",
                                   "A simple material model that is like the "
                                   "'melt simple' model, but has a constant reaction term.")
  }
}

