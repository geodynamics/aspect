#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class Anisotropic : public MaterialModel::Simple<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
        static void declare_parameters(ParameterHandler &prm);
        virtual void parse_parameters(ParameterHandler &prm);

        /**
          * Return true if the compressibility() function returns something that
          * is not zero.
          */
        virtual bool
        is_compressible () const;

      private:
        SymmetricTensor<4,dim> C; // Constitutive tensor
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    Anisotropic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      Simple<dim>::evaluate(in, out);
      Point<dim> center;
      center[0] = 0.5;
      center[1] = 0.5;
      if (dim == 3)
        center[2] = 0.5;
      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          const double pressure = in.pressure[i];
          out.densities[i] = 1.0 + pressure;
          out.compressibilities[i] = 1.0 / (1. + pressure);
          if ((in.position[i]-center).norm() < 0.25)
            {
              out.stress_strain_directors[i] = C;
            }
        }
    }

    template <int dim>
    bool
    Anisotropic<dim>::
    is_compressible () const
    {
      return true;
    }

    template <int dim>
    void
    Anisotropic<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      Simple<dim>::declare_parameters (prm);
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("Anisotropic");
        {
          if (dim == 2)
            prm.declare_entry ("Viscosity tensor",
                               "1, 0, 0,"
                               "0, 1, 0,"
                               "0, 0,.5",
                               Patterns::List(Patterns::Double()),
                               "Viscosity-scaling tensor in Voigt notation.");
          else
            prm.declare_entry ("Viscosity tensor",
                               "1, 0, 0, 0, 0, 0,"
                               "0, 1, 0, 0, 0, 0,"
                               "0, 0, 1, 0, 0, 0,"
                               "0, 0, 0,.5, 0, 0,"
                               "0, 0, 0, 0,.5, 0,"
                               "0, 0, 0, 0, 0,.5",
                               Patterns::List(Patterns::Double()),
                               "Viscosity-scaling tensor in Voigt notation.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Anisotropic<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("Anisotropic");
        {
          const int size_voigt = (dim == 3 ? 6 : 3);
          const std::vector<double> tmp_tensor =
            Utilities::string_to_double(Utilities::split_string_list(prm.get ("Viscosity tensor")));
          Assert(tmp_tensor.size() == size_voigt*size_voigt,
                 ExcMessage("Constitutive voigt matrix must have 9 components in 2D, or 36 components in 3d"));

          std::vector<std::vector<double> > voigt_visc_tensor (size_voigt);
          for (unsigned int i=0; i<size_voigt; ++i)
            {
              voigt_visc_tensor[i].resize(size_voigt);
              for (unsigned int j=0; j<size_voigt; ++j)
                voigt_visc_tensor[i][j] = tmp_tensor[i*size_voigt+j];
            }

          // Voigt indices (For mapping back to real tensor)
          const unsigned int vi3d0[] = {0, 1, 2, 1, 0, 0};
          const unsigned int vi3d1[] = {0, 1, 2, 2, 2, 1};
          const unsigned int vi2d0[] = {0, 1, 0};
          const unsigned int vi2d1[] = {0, 1, 1};

          // Fill the constitutive tensor with values from the Voigt tensor
          for (unsigned int i=0; i<size_voigt; ++i)
            for (unsigned int j=0; j<size_voigt; ++j)
              if (dim == 2)
                C[vi2d0[i]][vi2d1[i]][vi2d0[j]][vi2d1[j]] = voigt_visc_tensor[i][j];
              else
                C[vi3d0[i]][vi3d1[i]][vi3d0[j]][vi3d1[j]] = voigt_visc_tensor[i][j];
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Anisotropic,
                                   "anisotropic",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a non-zero compressibility "
                                   "and always has a blob of anisotropic material in "
                                   "the center.")
  }
}
