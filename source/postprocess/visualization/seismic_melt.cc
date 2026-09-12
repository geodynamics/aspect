/*
  Copyright (C) 2015 - 2018 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/seismic_melt.h>
#include <aspect/melt.h>
#include <aspect/utilities.h>
#include <aspect/simulator.h>
#include <aspect/material_model/interface.h>

#include <deal.II/numerics/data_out.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

      template <int dim>
      MeltSeismicProperties<dim>::
      MeltSeismicProperties ()
        :
        DataPostprocessor<dim> ()
      {}

      template <int dim>
      std::vector<std::string>
      MeltSeismicProperties<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;

        for (unsigned int i=0; i<property_names.size(); ++i)
          {
            solution_names.push_back(property_names[i]);
            std::replace(solution_names.back().begin(),solution_names.back().end(),' ', '_');
          }
        return solution_names;
      }

      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      MeltSeismicProperties<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
        for (unsigned int i=0; i<property_names.size(); ++i)
          {
            interpretation.push_back (DataComponentInterpretation::component_is_scalar);
          }
        return interpretation;
      }

      template <int dim>
      UpdateFlags
      MeltSeismicProperties<dim>::
      get_needed_update_flags () const
      {
        return  update_values  | update_q_points;
      }

      template <int dim>
      void
      MeltSeismicProperties<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true,
                    ExcMessage("'Include melt transport' has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());

        // Set use_strain_rates to true since the compaction viscosity might also depend on the strain rate.
        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection(), true);
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        MeltHandler<dim>::create_material_model_outputs(out);

        this->get_material_model().evaluate(in, out);
        MaterialModel::MeltOutputs<dim> *melt_outputs = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
        AssertThrow(melt_outputs != NULL,
                    ExcMessage("Need MeltOutputs from the material model for computing the melt properties."));

        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            unsigned output_index = 0;
            for (unsigned int i=0; i<property_names.size(); ++i, ++output_index)
              {
                double porosity = std::max(in.composition[q][porosity_idx],0.0);

                if (porosity < 0)
                  porosity = 0.0;

                double Vp = 0.0;
                double Vs = 0.0;
                double contiguity = 0.0;

                const double solid_bulk_modulus = solid_bulk_modulus_GPa * std::pow(10,9);
                const double liquid_bulk_modulus = liquid_bulk_modulus_GPa * std::pow(10,9);
                const double shear_modulus = shear_modulus_GPa * std::pow(10,9);

                double n_k = 0.0;
                double n_mu = 0.0;

                double solid_bulk_modulus_sk = 0.0;
                double shear_modulus_sk = 0.0;

                double skeleton_solid_bulk_modulus = 0.0;
                double skeleton_shear_modulus = 0.0;

                double a_1 = 0.0;
                double a_2 = 0.0;
                double a_3 = 0.0;
                double b_1 = 0.0;
                double b_2 = 0.0;

                if (porosity > 0)
                  {
                    //contiguity relation from Wimert and Hier-Majumder 2011
                    contiguity = p1 * std::pow(porosity,5) + p2 * std::pow(porosity,4) + p3 * std::pow(porosity,3) + p4 * std::pow(porosity,2) + p5 * porosity + p6 ;
                    if (contiguity < 0)
                      contiguity = 0.0;

                    // following from Takei (1998), Takei (2000), and Takei (2001)
                    const double a_hat [3][4] =
                    {
                      { 0.318, 6.780, 57.560, 0.182 },
                      { 0.164, 4.290, 26.658, 0.464 },
                      { 1.549, 4.814, 8.777, -0.290 }
                    };

                    const double b_hat [2][2] =
                    {
                      { -0.3238, 0.2341},
                      { -0.1819, 0.5103},
                    };

                    a_1 = a_hat[0][0] * std::exp(((a_hat[0][1] * (poisson - 0.25)) + a_hat[0][2] * std::pow((poisson - 0.25),3))) + a_hat[0][3] ;
                    a_2 = a_hat[1][0] * std::exp(((a_hat[1][1] * (poisson - 0.25)) + a_hat[1][2] * std::pow((poisson - 0.25),3))) + a_hat[1][3] ;
                    a_3 = a_hat[2][0] * std::exp(((a_hat[2][1] * (poisson - 0.25)) + a_hat[2][2] * std::pow((poisson - 0.25),3))) + a_hat[2][3] ;

                    b_1 = (b_hat[0][0] * poisson) + b_hat[0][1];
                    b_2 = (b_hat[1][0] * poisson) + b_hat[1][1];

                    n_k = (a_1 * contiguity) + (a_2 * (1 - contiguity)) + (a_3 * contiguity * (1 - contiguity) * (0.5 - contiguity)) ;
                    n_mu = (b_1 * contiguity) + (b_2 * (1 - contiguity));

                    solid_bulk_modulus_sk = solid_bulk_modulus * std::pow(contiguity,n_k);
                    shear_modulus_sk = shear_modulus * std::pow(contiguity,n_mu);

                    skeleton_solid_bulk_modulus = (1 - porosity) * solid_bulk_modulus_sk;
                    skeleton_shear_modulus = (1 - porosity) * shear_modulus_sk;

                    double H = skeleton_solid_bulk_modulus + (4 * skeleton_shear_modulus / 3) + ((solid_bulk_modulus * (std::pow((1 - (skeleton_solid_bulk_modulus / solid_bulk_modulus)),2))) / (1 - porosity - (skeleton_solid_bulk_modulus / solid_bulk_modulus) + (porosity * solid_bulk_modulus / liquid_bulk_modulus)));

                    double total_density = (1 - porosity) * out.densities[q] + porosity * melt_outputs->fluid_densities[q] ;

                    Vp = std::sqrt(H / total_density);
                    Vs = std::sqrt(skeleton_shear_modulus / total_density);
                  }

                else
                  {
                    contiguity = 1.0;
                    const double H_no_melt = solid_bulk_modulus + (4 * shear_modulus / 3);
                    Vp = std::sqrt( H_no_melt / out.densities[q]);
                    Vs = std::sqrt(shear_modulus / out.densities[q]);
                  }

                if (property_names[i] == "Vp")
                  computed_quantities[q][output_index] = Vp;
                else if (property_names[i] == "Vs")
                  computed_quantities[q][output_index] = Vs;
                else if (property_names[i] == "Vp_anomaly")
                  {
                    const double Vp_melt = Vp;
                    const double H_no_melt = solid_bulk_modulus + (4 * shear_modulus / 3);
                    const double Vp_no_melt = std::sqrt( H_no_melt / out.densities[i]);
                    const double Vp_anomaly = (Vp_melt - Vp_no_melt)/Vp_no_melt*1e2;
                    computed_quantities[q][output_index] = Vp_anomaly;
                  }
                else if (property_names[i] == "Vs_anomaly")
                  {
                    const double Vs_melt = Vs;
                    const double Vs_no_melt = std::sqrt(shear_modulus / out.densities[i]);
                    const double Vs_anomaly = (Vs_melt - Vs_no_melt)/Vs_no_melt*1e2;
                    computed_quantities[q][output_index] = Vs_anomaly;
                  }
                else if (property_names[i] == "contiguity")
                  {
                    computed_quantities[q][output_index] = contiguity;
                  }
                else
                  AssertThrow(false, ExcNotImplemented());
              }
          }
      }

      template <int dim>
      void
      MeltSeismicProperties<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt seismic properties");
            {
              const std::string pattern_of_names
                = "Vp|Vs|Vp_anomaly|Vs_anomaly|contiguity";

              prm.declare_entry("List of properties",
                                "Vp, Vs, Vp_anomaly, Vs_anomaly, contiguity",
                                Patterns::MultipleSelection(pattern_of_names),
                                "A comma separated list of melt properties that should be "
                                "written whenever writing graphical output. "
                                "The following seismic properties are available:\n\n"
                                +
                                pattern_of_names);

              prm.declare_entry ("Solid bulk modulus", "124",
                                 Patterns::Double (),
                                 "The bulk modulus of the solid material "
                                 "assuming no porosity. "
                                 "Units: $GPa$.");

              prm.declare_entry ("Liquid bulk modulus", "40",
                                 Patterns::Double (),
                                 "The bulk modulus of the liquid. "
                                 "Units: $GPa$.");

              prm.declare_entry ("Shear modulus", "64",
                                 Patterns::Double (),
                                 "The shear modulus of the solid material "
                                 "assuming no porosity. "
                                 "Units: $GPa$.");

              prm.declare_entry ("Poisson ratio", "0.25",
                                 Patterns::Double (),
                                 "The Poisson ratio of the solid material. "
                                 "Units: non-dimensional.");

              prm.declare_entry ("p1", "-8065.00",
                                 Patterns::Double (),
                                 "The first coefficient, associated with "
                                 "porosity to the fifth power, in the fifth "
                                 "order polynomial contiguity equation "
                                 "of Wimert and Hier-Majumder (2012).");

              prm.declare_entry ("p2", "6149.00",
                                 Patterns::Double (),
                                 "The second coefficient, associated with "
                                 "porosity to the fourth power, in the fifth "
                                 "order polynomial contiguity equation "
                                 "of Wimert and Hier-Majumder (2012).");

              prm.declare_entry ("p3", "-1778.00",
                                 Patterns::Double (),
                                 "The third coefficient, associated with "
                                 "porosity to the third power, in the fifth "
                                 "order polynomial contiguity equation "
                                 "of Wimert and Hier-Majumder (2012).");

              prm.declare_entry ("p4", "249.00",
                                 Patterns::Double (),
                                 "The fourth coefficient, associated with "
                                 "porosity to the second power, in the fifth "
                                 "order polynomial contiguity equation "
                                 "of Wimert and Hier-Majumder (2012).");

              prm.declare_entry ("p5", "-19.77",
                                 Patterns::Double (),
                                 "The fifth coefficient, associated with "
                                 "porosity to the first power, in the fifth "
                                 "order polynomial contiguity equation "
                                 "of Wimert and Hier-Majumder (2012).");

              prm.declare_entry ("p6", "1.00",
                                 Patterns::Double (),
                                 "The sixth coefficient in the fifth "
                                 "order polynomial contiguity equation "
                                 "of Wimert and Hier-Majumder (2012).");

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      MeltSeismicProperties<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt seismic properties");
            {
              property_names = Utilities::split_string_list(prm.get ("List of properties"));
              AssertThrow(Utilities::has_unique_entries(property_names),
                          ExcMessage("The list of strings for the parameter "
                                     "'Postprocess/Visualization/Melt seismic properties/List of properties' contains entries more than once. "
                                     "This is not allowed. Please check your parameter file."));

              solid_bulk_modulus_GPa    = prm.get_double ("Solid bulk modulus");
              liquid_bulk_modulus_GPa   = prm.get_double ("Liquid bulk modulus");
              shear_modulus_GPa         = prm.get_double ("Shear modulus");
              poisson                   = prm.get_double ("Poisson ratio");
              p1                        = prm.get_double ("p1");
              p2                        = prm.get_double ("p2");
              p3                        = prm.get_double ("p3");
              p4                        = prm.get_double ("p4");
              p5                        = prm.get_double ("p5");
              p6                        = prm.get_double ("p6");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltSeismicProperties,
                                                  "melt seismic properties",
                                                  "A visualization output object that generates output "
                                                  "for seismic velocities of material models "
                                                  "involving melt. Seismic velocities are dependent on "
                                                  "porosity and calculated from the equations of "
                                                  "Takei [1998, 2000, 2001]. ")

    }
  }
}
