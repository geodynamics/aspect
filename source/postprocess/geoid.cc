/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/geoid.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      HarmonicCoefficients::HarmonicCoefficients(const unsigned int max_degree)
      {
        unsigned int k= max_degree*(max_degree+1)/2;
        sine_coefficients.resize(k);
        cosine_coefficients.resize(k);
      }

      template <int dim>
      SphericalHarmonicsExpansion<dim>::SphericalHarmonicsExpansion(const unsigned int max_degree)
      :
      max_degree(max_degree),
      coefficients(max_degree)
      {}

      template <int dim>
      void
      SphericalHarmonicsExpansion<dim>::add_data_point (const Point<dim> &position,
                                                        const double value)
      {
        const std_cxx11::array<double,dim> spherical_position =
            Utilities::spherical_coordinates(position);

        for (unsigned int i = 0, k = 0; i < max_degree; ++i)
          for (unsigned int j = 0; j <= i; ++j, ++k)
          {
              coefficients.sine_coefficients[k] += value
                                                  * boost::math::spherical_harmonic_r(i,j,spherical_position[2],spherical_position[1]);
                                                  //* boost::math::legendre_p(i,j,cos(spherical_position[2]))
                                                  // * cos(j*spherical_position[1]);
              coefficients.cosine_coefficients[k] += value
                                                  * boost::math::spherical_harmonic_i(i,j,spherical_position[2],spherical_position[1]);
                                                  //* boost::math::legendre_p(i,j,cos(spherical_position[2]))
                                                  // * sin(j*spherical_position[1]);
          }
      }

      template <int dim>
      HarmonicCoefficients
      SphericalHarmonicsExpansion<dim>::get_coefficients () const
      {
        return coefficients;
      }

      template <int dim>
      void
      SphericalHarmonicsExpansion<dim>::mpi_sum_coefficients (MPI_Comm mpi_communicator)
      {
        dealii::Utilities::MPI::sum(coefficients.sine_coefficients,mpi_communicator,coefficients.sine_coefficients);
        dealii::Utilities::MPI::sum(coefficients.cosine_coefficients,mpi_communicator,coefficients.cosine_coefficients);
      }
    }

    template <int dim>
    std::pair<std::string,std::string>
    Geoid<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      const QMidpoint<dim> quadrature_formula;
      const unsigned int n_q_points = quadrature_formula.size();

      Assert(quadrature_formula.size()==1, ExcInternalError());

      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                    (&this->get_geometry_model());

      AssertThrow (geometry_model != 0,
                   ExcMessage("The geoid postprocessor is currently only implemented for "
                       "the spherical shell geometry model."));

      // TODO AssertThrow (no_free_surface);

      buoyancy_expansions.resize(number_of_layers);
      for (unsigned int i = 0; i < number_of_layers; ++i)
        buoyancy_expansions[i].reset(new internal::SphericalHarmonicsExpansion<dim>(max_degree));

      surface_topography_expansion.reset(new internal::SphericalHarmonicsExpansion<dim>(max_degree));
      bottom_topography_expansion.reset(new internal::SphericalHarmonicsExpansion<dim>(max_degree));

      //TODO: remove this when finished, it is only there for benchmark purposes
      internal::SphericalHarmonicsExpansion<dim> example_expansion(max_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_q_points);

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      std::vector<double> average_densities(number_of_layers);
      this->get_depth_average_density(average_densities);

      // Some constant that are used several times
      const double inner_radius = geometry_model->inner_radius();
      const double outer_radius = geometry_model->outer_radius();
      const double layer_thickness = geometry_model->maximal_depth() / number_of_layers;
      const double gravitational_constant = 6.67384e-11;

      // loop over all of the surface cells and if one less than h/3 away from
      // the top surface, evaluate the stress at its center
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);

              // get the various components of the solution, then
              // evaluate the material properties there
              fe_values[this->introspection().extractors.temperature]
              .get_function_values (this->get_solution(), in.temperature);
              fe_values[this->introspection().extractors.pressure]
              .get_function_values (this->get_solution(), in.pressure);
              fe_values[this->introspection().extractors.velocities]
              .get_function_symmetric_gradients (this->get_solution(), in.strain_rate);

              in.position = fe_values.get_quadrature_points();

              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]]
                .get_function_values(this->get_solution(),
                                     composition_values[c]);
              for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
                {
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    in.composition[i][c] = composition_values[c][i];
                }

              this->get_material_model().evaluate(in, out);

              // see if the cell is at the *top* or *bottom* boundary
              bool surface_cell = false;
              bool bottom_cell = false;

              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                  {
                    surface_cell = true;
                    break;
                  }
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) > (outer_radius - cell->face(f)->minimum_vertex_distance()/3))
                  {
                    bottom_cell = true;
                    break;
                  }
                }
              // for each of the quadrature points, evaluate the
              // density and add its contribution to the spherical harmonics

              for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                {
                  const Point<dim> location = fe_values.quadrature_point(q);

                  const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                  const Tensor<1,dim> gravity_direction = gravity/gravity.norm();
                  const unsigned int layer_id = static_cast<unsigned int> (geometry_model->depth(location) / geometry_model->maximal_depth() * (number_of_layers-1));

                  const double density   = out.densities[q];
                  const double buoyancy = -1.0 * (density - average_densities[layer_id]) * gravity.norm();

                  buoyancy_expansions[layer_id]->add_data_point(location,buoyancy);

                  //TODO: remove this when finished, it is only there for benchmark purposes
                  const std_cxx1x::array<double,dim> spherical_position =
                      Utilities::spherical_coordinates(location);
                  example_expansion.add_data_point(location,boost::math::spherical_harmonic_r(3,2,spherical_position[2],spherical_position[1]));

                  // if this is a cell at the surface, add the topography to
                  // the topography expansion
                  if (surface_cell)
                    {
                      const double viscosity = out.viscosities[q];

                      const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                      const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                      // Subtract the adiabatic pressure
                      const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(location);
                      const double sigma_rr           = gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure;
                      const double dynamic_topography = - sigma_rr / gravity.norm() / (density - density_above);

                      // Add topography contribution
                      surface_topography_expansion->add_data_point(location,dynamic_topography);
                    }

                  // if this is a cell at the bottom, add the topography to
                  // the bottom expansion
                  if (bottom_cell)
                    {
                      const double viscosity = out.viscosities[q];

                      const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                      const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                      // Subtract the adiabatic pressure
                      const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(location);
                      const double sigma_rr           = gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure;
                      const double dynamic_topography = - sigma_rr / gravity.norm() / (density_below - density);

                      // Add topography contribution
                      bottom_topography_expansion->add_data_point(location,dynamic_topography);
                    }
                }
            }

      // From here on we consider the gravity constant at the value of the surface
      const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(geometry_model->representative_point(0.0));

      const double scaling = 4.0 * numbers::PI * outer_radius * gravitational_constant / gravity.norm();

      internal::HarmonicCoefficients buoyancy_contribution(max_degree);

      internal::HarmonicCoefficients surface_geoid_expansion(max_degree);
      internal::HarmonicCoefficients bottom_geoid_expansion(max_degree);

      for (unsigned int layer_id = 0; layer_id < number_of_layers; ++layer_id)
        {
          buoyancy_expansions[layer_id]->mpi_sum_coefficients(this->get_mpi_communicator());
          const internal::HarmonicCoefficients layer_coefficients = buoyancy_expansions[layer_id]->get_coefficients();

          const double layer_radius = inner_radius + (layer_id + 0.5) * layer_thickness;

          for (unsigned int i = 0, k = 0; i < max_degree; ++i)
            {
              const double con1 = scaling * layer_thickness / (2 * i + 1.0);
              const double cont = pow (layer_radius/outer_radius,i+2);
              const double conb = layer_radius / outer_radius * pow(inner_radius / layer_radius, i);

              for (unsigned int j = 0; j <= i; ++j, ++k)
                {
                  surface_geoid_expansion.sine_coefficients[k] += con1 * cont * layer_coefficients.sine_coefficients[k];
                  surface_geoid_expansion.cosine_coefficients[k] += con1 * cont * layer_coefficients.cosine_coefficients[k];

                  buoyancy_contribution.sine_coefficients[k] += con1 * cont * layer_coefficients.sine_coefficients[k];
                  buoyancy_contribution.cosine_coefficients[k] += con1 * cont * layer_coefficients.cosine_coefficients[k];

                  //TODO: bottom_geoid
                  //bottom_geoid_expansion.sine_coefficients[index] += con1 * conb * layer_coefficients.sine_coefficients[index];
                  //bottom_geoid_expansion.cosine_coefficients[index] += con1 * conb * layer_coefficients.cosine_coefficients[index];

                }
            }
        }


      internal::HarmonicCoefficients topography_contribution(max_degree);
      if (include_topography_contribution)
        {
          surface_topography_expansion->mpi_sum_coefficients(this->get_mpi_communicator());
          const internal::HarmonicCoefficients topography_coefficients = surface_topography_expansion->get_coefficients();

          // Use the average density at the surface to compute density jump across the surface
          const double density_contrast = average_densities.front() - density_above;

          for (unsigned int i = 0, k = 0; i < max_degree; i++)
            {
              const double con1 = scaling * density_contrast / (2.0 * i + 1.0);

              for (unsigned int j = 0; j <= i; ++j, ++k)
                {
                  surface_geoid_expansion.sine_coefficients[k] += con1 * topography_coefficients.sine_coefficients[k];
                  surface_geoid_expansion.cosine_coefficients[k] += con1 * topography_coefficients.cosine_coefficients[k];

                  topography_contribution.sine_coefficients[k] += con1 * topography_coefficients.sine_coefficients[k];
                  topography_contribution.cosine_coefficients[k] += con1 * topography_coefficients.cosine_coefficients[k];
                }
            }
        }

      if (include_topography_contribution)
        {
          bottom_topography_expansion->mpi_sum_coefficients(this->get_mpi_communicator());
          const internal::HarmonicCoefficients topography_coefficients = bottom_topography_expansion->get_coefficients();

          // Use the average density at the bottom to compute density jump across the bottom surface
          const double density_contrast = density_below - average_densities.back();

          for (unsigned int i = 0, k = 0; i < max_degree; i++)
            {
              const double con1 = density_contrast * scaling / (2.0 * i + 1.0);
              const double con2 = con1 * pow(inner_radius/outer_radius, i+2);

              for (unsigned int j = 0; j <= i; ++j, ++k)
                {
                  surface_geoid_expansion.sine_coefficients[k] += con2 * topography_coefficients.sine_coefficients[k];
                  surface_geoid_expansion.cosine_coefficients[k] += con2 * topography_coefficients.cosine_coefficients[k];
                }
            }
        }


      const std::string filename = this->get_output_directory() +
                                   "surface_geoid." +
                                   dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

      // On process 0 write output file
      if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                                   this->get_time() / year_in_seconds :
                                                   this->get_time());
          std::ofstream file (filename.c_str());

          file << "# Timestep Maximum_degree Time" << std::endl;
          file << this->get_timestep_number() << " " << max_degree << " " << time_in_years_or_seconds << std::endl;
          file << "# degree order geoid_sine_coefficient geoid_cosine_coefficient buoyancy_sine buoyancy_cosine topography_sine topography_cosine" << std::endl;
          // Write the solution to an output file
          for (unsigned int i=0, k=0; i < max_degree; ++i)
            {
              for (unsigned int j = 0; j <= i; ++j, ++k)
                {
                  file << i << " " << j << " "
                       << surface_geoid_expansion.sine_coefficients[k] << " "
                       << surface_geoid_expansion.cosine_coefficients[k] << " "
                       << buoyancy_contribution.sine_coefficients[k] << " "
                       << buoyancy_contribution.cosine_coefficients[k] << " "
                       << topography_contribution.sine_coefficients[k] << " "
                       << topography_contribution.cosine_coefficients[k] << std::endl;
                }
            }
        }

      return std::pair<std::string,std::string>("Writing geoid:",
                                                filename);
    }

    template <int dim>
    void
    Geoid<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          prm.declare_entry ("Include topography contribution", "false",
                             Patterns::Bool (),
                             "Option to include the contribution of dynamic "
                             "topography to the geoid.");
          prm.declare_entry ("Number of layers", "20",
                             Patterns::Integer (1),
                             "The geoid contribution is added on a per-layer basis. This parameter "
                             "sets the number of layers. Similar to the depth-average "
                             "postprocessor, the number of layers should correspond roughly to "
                             "the available model resolution.");
          prm.declare_entry ("Density below", "8000",
                             Patterns::Double(0),
                             "");
          prm.declare_entry ("Density above", "0",
                             Patterns::Double(),
                             "");
          prm.declare_entry ("Maximum degree of expansion", "7",
                             Patterns::Integer (1),
                             "");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Geoid<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          include_topography_contribution   = prm.get_bool("Include topography contribution");
          number_of_layers                  = prm.get_integer("Number of layers");
          density_below                     = prm.get_double("Density below");
          density_above                     = prm.get_double("Density above");
          max_degree                        = prm.get_integer("Maximum degree of expansion");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Geoid,
                                  "geoid",
                                  "A postprocessor that computes a measure of geoid height "
                                  "based on the internal buoyancy and top and bottom topography")
  }
}
