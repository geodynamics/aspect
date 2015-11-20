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
#include <aspect/lateral_averaging.h>
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
        unsigned int k=0;
        for (unsigned int i=0; i<=max_degree; i++)
          {
            for (unsigned int j=0; j<=i; ++j)
              {
                ++k;
              }
          }

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

        for (unsigned int i = 0, k = 0; i <= max_degree; ++i)
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
    Geoid<dim>::execute (TableHandler &/*statistics*/)
    {
      // create a quadrature formula based on the velocity element alone.
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);

      //Assert(quadrature_formula.size()==1, ExcInternalError());

      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                                 (&this->get_geometry_model());

      AssertThrow (geometry_model != 0,
                   ExcMessage("The geoid postprocessor is currently only implemented for "
                              "the spherical shell geometry model."));

      // TODO AssertThrow (no_free_surface);

      density_expansions.resize(number_of_layers);
      for (unsigned int i = 0; i < number_of_layers; ++i)
        density_expansions[i].reset(new internal::SphericalHarmonicsExpansion<dim>(max_degree));

      surface_topography_expansion.reset(new internal::SphericalHarmonicsExpansion<dim>(max_degree));
      bottom_topography_expansion.reset(new internal::SphericalHarmonicsExpansion<dim>(max_degree));

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_JxW_values |
                               update_values |
                               update_gradients |
                               update_q_points);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_JxW_values |
                                        update_values |
                                        update_gradients |
                                        update_q_points);

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      typename MaterialModel::Interface<dim>::MaterialModelInputs face_in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs face_out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > face_composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula_face.size()));


      std::vector<double> average_densities(number_of_layers);
      this->get_lateral_averaging().get_density_averages(average_densities);

      // Some constant that are used several times
      const double inner_radius = geometry_model->inner_radius();
      const double outer_radius = geometry_model->outer_radius();
      const double layer_thickness = geometry_model->maximal_depth() / number_of_layers;

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
            unsigned int top_face_idx = numbers::invalid_unsigned_int;
            unsigned int bottom_face_idx = numbers::invalid_unsigned_int;

            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                  {
                    top_face_idx = f;
                    break;
                  }
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) > (outer_radius - inner_radius - cell->face(f)->minimum_vertex_distance()/3))
                  {
                    bottom_face_idx = f;
                    break;
                  }
              }

            // in case we need to calculate the dynamic topography, store the
            // radial stress in here
            double sigma_rr(0.0);
            double volume(0.0);

            // for each of the quadrature points, evaluate the
            // density and add its contribution to the spherical harmonics
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              {
                const Point<dim> position = fe_values.quadrature_point(q);
                const double radius = position.norm();

                const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(position);
                const Tensor<1,dim> gravity_direction = gravity/gravity.norm();
                const unsigned int layer_id = rint (geometry_model->depth(position) / geometry_model->maximal_depth() * (number_of_layers-1));

                const double layer_volume = 4.0 / 3.0 * numbers::PI * (pow(radius+layer_thickness/2.0,3) - pow(radius-layer_thickness/2.0,3));

                const double density  = out.densities[q];
                const double buoyancy = (density - average_densities[layer_id]);

                density_expansions[layer_id]->add_data_point(position,buoyancy * fe_values.JxW(q) / layer_volume);

                // If this cell is at one boundary calculate cell integrated
                // radial stress deviation
                if ((top_face_idx != numbers::invalid_unsigned_int)
                    || (bottom_face_idx != numbers::invalid_unsigned_int))
                  {
                    const double viscosity = out.viscosities[q];

                    const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                    const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                    // Subtract the adiabatic pressure
                    const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(position);

                    sigma_rr += (gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure) * fe_values.JxW(q);
                    volume += fe_values.JxW(q);
                  }
              }

            // If this cell is at one boundary calculate boundary contribution
            // from cell averaged radial stress deviation
            if ((top_face_idx != numbers::invalid_unsigned_int)
                || (bottom_face_idx != numbers::invalid_unsigned_int))
              {
                const unsigned int face_index = (top_face_idx != numbers::invalid_unsigned_int)
                                                ?
                                                top_face_idx
                                                :
                                                bottom_face_idx;

                fe_face_values.reinit (cell,face_index);

                // get the various components of the solution, then
                // evaluate the material properties there
                fe_face_values[this->introspection().extractors.temperature]
                .get_function_values (this->get_solution(), face_in.temperature);
                fe_face_values[this->introspection().extractors.pressure]
                .get_function_values (this->get_solution(), face_in.pressure);
                fe_face_values[this->introspection().extractors.velocities]
                .get_function_symmetric_gradients (this->get_solution(), face_in.strain_rate);

                face_in.position = fe_face_values.get_quadrature_points();

                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  fe_face_values[this->introspection().extractors.compositional_fields[c]]
                  .get_function_values(this->get_solution(),
                                       face_composition_values[c]);
                for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
                  {
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      face_in.composition[i][c] = face_composition_values[c][i];
                  }

                this->get_material_model().evaluate(face_in, face_out);

                // average radial stress by dividing integrated stress by volume
                sigma_rr /= volume;

                for (unsigned int q=0; q<quadrature_formula_face.size(); ++q)
                  {
                    const Point<dim> position = fe_face_values.quadrature_point(q);
                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(position);
                    const double surface = fe_face_values.JxW(q);

                    // if this is a cell at the top, add to surface expansion
                    if (top_face_idx != numbers::invalid_unsigned_int)
                      {
                        const double dynamic_topography = - sigma_rr / gravity.norm() / (face_out.densities[q] - density_above);
                        const double total_surface = 4 * numbers::PI * outer_radius * outer_radius;

                        // Add topography contribution
                        surface_topography_expansion->add_data_point(position,dynamic_topography * surface / total_surface);
                      }
                    // if this is a cell at the bottom, add to bottom expansion
                    else
                      {
                        const double dynamic_topography = - sigma_rr / gravity.norm() / (face_out.densities[q] - density_below);
                        const double total_surface = 4 * numbers::PI * inner_radius * inner_radius;

                        // Add topography contribution
                        bottom_topography_expansion->add_data_point(position,dynamic_topography * surface / total_surface);
                      }
                  }
              }
          }

      // From here on we consider the gravity constant at the value of the surface
      const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(geometry_model->representative_point(0.0));

      const double scaling = 4.0 * numbers::PI * outer_radius * constants::big_g / gravity.norm();

      internal::HarmonicCoefficients buoyancy_contribution(max_degree);
      internal::HarmonicCoefficients bottom_geoid_expansion(max_degree);
      surface_geoid_expansion.reset(new internal::HarmonicCoefficients(max_degree));

      for (unsigned int layer_id = 0; layer_id < number_of_layers; ++layer_id)
        {
          density_expansions[layer_id]->mpi_sum_coefficients(this->get_mpi_communicator());
          const internal::HarmonicCoefficients layer_coefficients = density_expansions[layer_id]->get_coefficients();

          const double layer_radius = inner_radius + (layer_id + 0.5) * layer_thickness;

          if (outer_radius - layer_radius > density_exclusion_depth)
            for (unsigned int i = 0, k = 0; i <= max_degree; ++i)
              {
                const double con1 = scaling * layer_thickness / (2 * i + 1.0);
                const double cont = pow (layer_radius/outer_radius,i+2);
                //const double conb = layer_radius / outer_radius * pow(inner_radius / layer_radius, i);

                for (unsigned int j = 0; j <= i; ++j, ++k)
                  {
                    surface_geoid_expansion->sine_coefficients[k] += con1 * cont * layer_coefficients.sine_coefficients[k];
                    surface_geoid_expansion->cosine_coefficients[k] += con1 * cont * layer_coefficients.cosine_coefficients[k];

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

          for (unsigned int i = 0, k = 0; i <= max_degree; i++)
            {
              const double con1 = scaling * density_contrast / (2.0 * i + 1.0);

              for (unsigned int j = 0; j <= i; ++j, ++k)
                {
                  surface_geoid_expansion->sine_coefficients[k] += con1 * topography_coefficients.sine_coefficients[k];
                  surface_geoid_expansion->cosine_coefficients[k] += con1 * topography_coefficients.cosine_coefficients[k];

                  topography_contribution.sine_coefficients[k] += con1 * topography_coefficients.sine_coefficients[k];
                  topography_contribution.cosine_coefficients[k] += con1 * topography_coefficients.cosine_coefficients[k];
                }
            }
        }

      internal::HarmonicCoefficients bottom_topography_contribution(max_degree);
      if (include_topography_contribution)
        {
          bottom_topography_expansion->mpi_sum_coefficients(this->get_mpi_communicator());
          const internal::HarmonicCoefficients topography_coefficients = bottom_topography_expansion->get_coefficients();

          // Use the average density at the bottom to compute density jump across the bottom surface
          const double density_contrast = density_below - average_densities.back();

          for (unsigned int i = 0, k = 0; i <= max_degree; i++)
            {
              const double con1 = density_contrast * scaling / (2.0 * i + 1.0);
              const double con2 = con1 * pow(inner_radius/outer_radius, i+2);

              for (unsigned int j = 0; j <= i; ++j, ++k)
                {
                  surface_geoid_expansion->sine_coefficients[k] += con2 * topography_coefficients.sine_coefficients[k];
                  surface_geoid_expansion->cosine_coefficients[k] += con2 * topography_coefficients.cosine_coefficients[k];

                  bottom_topography_contribution.sine_coefficients[k] += con2 * topography_coefficients.sine_coefficients[k];
                  bottom_topography_contribution.cosine_coefficients[k] += con2 * topography_coefficients.cosine_coefficients[k];
                  surface_geoid_expansion->sine_coefficients[k] += con2 * topography_coefficients.sine_coefficients[k];
                  surface_geoid_expansion->cosine_coefficients[k] += con2 * topography_coefficients.cosine_coefficients[k];
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
          file << "# degree order geoid_sine_coefficient geoid_cosine_coefficient buoyancy_sine buoyancy_cosine topography_sine topography_cosine bottom_topography_sine bottom_topography_cosine" << std::endl;
          // Write the solution to an output file
          for (unsigned int i=0, k=0; i <= max_degree; ++i)
            {
              for (unsigned int j = 0; j <= i; ++j, ++k)
                {
                  file << i << " " << j << " "
                       << surface_geoid_expansion->sine_coefficients[k] << " "
                       << surface_geoid_expansion->cosine_coefficients[k] << " "
                       << buoyancy_contribution.sine_coefficients[k] << " "
                       << buoyancy_contribution.cosine_coefficients[k] << " "
                       << topography_contribution.sine_coefficients[k] << " "
                       << topography_contribution.cosine_coefficients[k] << " "
                       << bottom_topography_contribution.sine_coefficients[k] << " "
                       << bottom_topography_contribution.cosine_coefficients[k] <<std::endl;
                }
            }
        }

      return std::pair<std::string,std::string>("Writing geoid:",
                                                filename);
    }

    template <int dim>
    const internal::HarmonicCoefficients &
    Geoid<dim>::get_geoid_coefficients () const
    {
      return *surface_geoid_expansion;
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
          prm.declare_entry ("Density below", "9900",
                             Patterns::Double(0),
                             "This constant density is assumed for the material underlying the "
                             "model to calculate its topography contribution. The default value "
                             "is chosen to match the outer core density.");
          prm.declare_entry ("Density above", "1000",
                             Patterns::Double(0),
                             "This constant density is assumed for the material overlying the "
                             "model to calculate its topography contribution. The default value "
                             "is chosen to match the density of water.");
          prm.declare_entry ("Maximum degree of expansion", "7",
                             Patterns::Integer (1),
                             "The geoid will be expanded in spherical harmonics up to this degree. "
                             "If this degree of expansion is set too high compared to the model "
                             "resolution the output will be unreliable.");
          prm.declare_entry ("Remove density heterogeneity down to specified depth", "0",
                             Patterns::Double (0),
                             "The density variation of the uppermost layer is included in the "
                             "topography contribution. Therefore it can make sense to remove "
                             "it from the buoyancy expansion up to a certain depth.");
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
          density_exclusion_depth           = prm.get_double("Remove density heterogeneity down to specified depth");
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
