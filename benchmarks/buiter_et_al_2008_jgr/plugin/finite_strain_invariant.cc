#include <aspect/material_model/visco_plastic.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class FiniteStrainInvariant : public MaterialModel::ViscoPlastic<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
        virtual void parse_parameters(ParameterHandler &prm);
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    FiniteStrainInvariant<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // First, we use the material descriptions of the 'visco_plastic' material model to fill all of the material
      // model outputs. Below, we will then overwrite selected properties (the reaction terms), which are
      // needed to track the finite strain invariant.
      ViscoPlastic<dim>::evaluate(in, out);

      if (in.cell && this->get_timestep_number() > 0)
        {

          // Loop through quadrature points
          for (unsigned int q=0; q < in.position.size(); ++q)
            {
              const double edot_ii = std::max(sqrt(std::fabs(second_invariant(deviator(in.strain_rate[q])))),this->min_strain_rate);

              // New strain invariant is old strain invariant plut the strain invariant of the current time step
              const double e_ii = in.composition[q][0] + edot_ii*this->get_timestep();

              // Update reaction term
              out.reaction_terms[q][0] = - in.composition[q][0] + e_ii;

            }
        }
    }


    template <int dim>
    void
    FiniteStrainInvariant<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      ViscoPlastic<dim>::parse_parameters (prm);

      AssertThrow(this->n_compositional_fields() >= 1,
                  ExcMessage("There must be at least as one compositional field. "));
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(FiniteStrainInvariant,
                                   "finite strain invariant",
                                   "A simple material model that is like the "
                                   "'ViscoPlastic' model, but tracks the finite strain invariant as a "
                                   "a compositional field. More precisely, the model assumes that the first "
                                   "compositional field contains the finite strain invariant, which "
                                   "is calculated using the second invariant of the deviatoric strain rate "
                                   "tensor ($\\dot{\\varepsilon}_{ii}$) and the current time step ($dt$). "
                                   "At time step $n$, the finite strain invariant is  "
                                   "$\\varepsilon_{ii}^{n} = \\varepsilon_{ii}^{n-1} + $ "
                                   "$dt*\\dot{\\varepsilon}_{ii}$. The latter terms are stored in the "
                                   "reaction term, which then updates the compositional field value.")

  }
}
