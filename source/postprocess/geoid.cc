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
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <aspect/utilities.h>
#include <aspect/postprocess/geoid.h>
#include <aspect/geometry_model/spherical_shell.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      template <int dim>
      HarmonicCoefficients<dim>::HarmonicCoefficients(const unsigned int max_degree)
      {
        // Allocate vectors of the appropriate size for 2D and 3D.  For 2D we
        // need cylindrical harmonics, for 3D we need spherical harmonics.

        unsigned int k= (dim == 2 ?
                         2 * (max_degree + 1 ) :             //cylindrical harmonics
                         (max_degree+1)*(max_degree+2)/2 ) ; //spherical harmonics
        sine_coefficients.resize(k);
        cosine_coefficients.resize(k);
      }

      // Constructor.  Just need the initialize the coefficients struct.
      template <int dim>
      MultipoleExpansion<dim>::MultipoleExpansion(const unsigned int max_degree)
        :
        max_degree(max_degree),
        coefficients(max_degree)
      {}

      // Add a quadrature evaluation to the multipole expansion in 3D using spherical
      // multipole moments. It will also expand in internal or external harmonics,
      //depending upon whether the position vector is at a larger radius than evaluation_radius.
      template <>
      void
      MultipoleExpansion<2>::add_quadrature_point (const Point<2> &position,
                                                   const double value,
                                                   const double JxW,
                                                   const double evaluation_radius)
      {
        const double r = position.norm();
        const double theta = std::atan2(position[1],position[0]);
        const bool is_external = ( r<= evaluation_radius);  //determine whether external or internal

        if ( is_external && evaluation_radius > 0.)
          {
            for (unsigned int n = 2; n <= max_degree; ++n)
              {
                const double factor = value * std::pow( r / evaluation_radius, static_cast<double>(n) )
                                      / static_cast<double>(n) * JxW;
                coefficients.cosine_coefficients[n] += factor * std::cos( static_cast<double>(n) * theta);
                coefficients.sine_coefficients[n] +=   factor * std::sin( static_cast<double>(n) * theta);
              }
          }
        else if ( !is_external && r > 0. )
          {
            for (unsigned int n = 2; n <= max_degree; ++n)
              {
                const double factor = value * std::pow( evaluation_radius/r, static_cast<double>(n) )
                                      / static_cast<double>(n) * JxW;
                coefficients.cosine_coefficients[n] += factor * std::cos( static_cast<double>(n) * theta);
                coefficients.sine_coefficients[n] +=   factor * std::sin( static_cast<double>(n) * theta);
              }
          }
      }

      // Add a quadrature evaluation to the multipole expansion in 3D using spherical
      // multipole moments. It will also expand in internal or external harmonics,
      //depending upon whether the position vector is at a larger radius than evaluation_radius.
      template <>
      void
      MultipoleExpansion<3>::add_quadrature_point (const Point<3> &position,
                                                   const double value,
                                                   const double JxW,
                                                   const double evaluation_radius)
      {
        const double r = position.norm();
        const double phi = std::atan2(position[1],position[0]);
        const double theta = std::acos(position[2]/r);
        const bool is_external = ( r<= evaluation_radius);  //determine whether external or internal

        if ( is_external && evaluation_radius > 0.)
          {
            for (unsigned int l = 2, k = 0; l <= max_degree; ++l)
              for (unsigned int m = 0; m <= l; ++m, ++k)
                {
                  std::pair<double,double> coeffs = aspect::Utilities::real_spherical_harmonic( l, m, theta, phi );

                  coefficients.cosine_coefficients[k] += value*std::pow(r/evaluation_radius,static_cast<double>(l))/evaluation_radius
                                                         * coeffs.first * JxW;
                  coefficients.sine_coefficients[k] += value*std::pow(r/evaluation_radius,static_cast<double>(l))/evaluation_radius
                                                       * coeffs.second * JxW;
                }
          }
        else if ( !is_external && r > 0. )
          {
            for (unsigned int l = 2, k = 0; l <= max_degree; ++l)
              for (unsigned int m = 0; m <= l; ++m, ++k)
                {
                  std::pair<double,double> coeffs = aspect::Utilities::real_spherical_harmonic( l, m, theta, phi );

                  coefficients.cosine_coefficients[k] += value*std::pow(evaluation_radius/r,static_cast<double>(l)) / r
                                                         * coeffs.first * JxW;
                  coefficients.sine_coefficients[k] += value*std::pow(evaluation_radius/r,static_cast<double>(l)) / r
                                                       * coeffs.second * JxW;
                }
          }

      }

      //Return a reference to the internal representation of the multipole expansion
      template <int dim>
      const HarmonicCoefficients<dim> &
      MultipoleExpansion<dim>::get_coefficients () const
      {
        return coefficients;
      }

      //Zero out the expansion in 2D
      template <>
      void MultipoleExpansion<2>::clear()
      {
        for (unsigned int n = 2; n <= max_degree; ++n)
          {
            coefficients.sine_coefficients[n] = 0.0;
            coefficients.cosine_coefficients[n] = 0.0;
          }
      }

      //Zero out the expansion in 3D
      template <>
      void MultipoleExpansion<3>::clear()
      {
        for (unsigned int l = 2, k = 0; l <= max_degree; ++l)
          for (unsigned int m = 0; m <= l; ++m, ++k)
            {
              coefficients.sine_coefficients[k] = 0.0;
              coefficients.cosine_coefficients[k] = 0.0;
            }
      }

      template <>
      void MultipoleExpansion<2>::sadd( double s, double a, const MultipoleExpansion &M )
      {
        AssertThrow( coefficients.sine_coefficients.size() == M.get_coefficients().sine_coefficients.size() ,
                     ExcInternalError() );

        for (unsigned int n = 2; n <= max_degree; ++n)
          {
            coefficients.sine_coefficients[n] = s * coefficients.sine_coefficients[n] +
                                                a * M.get_coefficients().sine_coefficients[n];
            coefficients.cosine_coefficients[n] = s * coefficients.cosine_coefficients[n] +
                                                  a * M.get_coefficients().cosine_coefficients[n];
          }
      }

      template <>
      void MultipoleExpansion<3>::sadd( double s, double a, const MultipoleExpansion &M )
      {
        AssertThrow( coefficients.sine_coefficients.size() == M.get_coefficients().sine_coefficients.size() ,
                     ExcInternalError() );

        for (unsigned int l = 2, k = 0; l <= max_degree; ++l)
          for (unsigned int m = 0; m <= l; ++m, ++k)
            {
              coefficients.sine_coefficients[k] = s * coefficients.sine_coefficients[k] +
                                                  a * M.get_coefficients().sine_coefficients[k];
              coefficients.cosine_coefficients[k] = s * coefficients.cosine_coefficients[k] +
                                                    a * M.get_coefficients().cosine_coefficients[k];
            }
      }

      //Version of sadd which has a different factor for each degree.  As such, s and a should have sizes of (max_degree+1)
      template <>
      void MultipoleExpansion<2>::sadd( const std::vector<double> &s, const std::vector<double> &a, const MultipoleExpansion &M )
      {
        AssertThrow( coefficients.sine_coefficients.size() == M.get_coefficients().sine_coefficients.size() ,
                     ExcInternalError() );

        AssertThrow( s.size() == max_degree+1 , ExcInternalError() );
        AssertThrow( a.size() == max_degree+1 , ExcInternalError() );

        for (unsigned int n = 2; n <= max_degree; ++n)
          {
            coefficients.sine_coefficients[n] = s[n] * coefficients.sine_coefficients[n] +
                                                a[n] * M.get_coefficients().sine_coefficients[n];
            coefficients.cosine_coefficients[n] = s[n] * coefficients.cosine_coefficients[n] +
                                                  a[n] * M.get_coefficients().cosine_coefficients[n];
          }
      }

      //Version of sadd which has a different factor for each degree.  As such, s and a should have sizes of (max_degree+1)
      template <>
      void MultipoleExpansion<3>::sadd( const std::vector<double> &s, const std::vector<double> &a, const MultipoleExpansion &M )
      {
        AssertThrow( coefficients.sine_coefficients.size() == M.get_coefficients().sine_coefficients.size() ,
                     ExcInternalError() );

        AssertThrow( s.size() == max_degree+1 , ExcInternalError() );
        AssertThrow( a.size() == max_degree+1 , ExcInternalError() );

        for (unsigned int l = 2, k = 0; l <= max_degree; ++l)
          for (unsigned int m = 0; m <= l; ++m, ++k)
            {
              coefficients.sine_coefficients[k] = s[l] * coefficients.sine_coefficients[k] +
                                                  a[l] * M.get_coefficients().sine_coefficients[k];
              coefficients.cosine_coefficients[k] = s[l] * coefficients.cosine_coefficients[k] +
                                                    a[l] * M.get_coefficients().cosine_coefficients[k];
            }
      }

      template <int dim>
      void
      MultipoleExpansion<dim>::mpi_sum_coefficients (MPI_Comm mpi_communicator)
      {
        dealii::Utilities::MPI::sum(coefficients.sine_coefficients,mpi_communicator,coefficients.sine_coefficients);
        dealii::Utilities::MPI::sum(coefficients.cosine_coefficients,mpi_communicator,coefficients.cosine_coefficients);
      }
    }

    template <int dim>
    std::pair<std::string,std::string>
    Geoid<dim>::execute (TableHandler &/*statistics*/)
    {
      internal_density_expansion_surface->clear();
      internal_density_expansion_bottom->clear();
      surface_topography_expansion->clear();
      bottom_topography_expansion->clear();
      surface_potential_expansion->clear();
      bottom_potential_expansion->clear();

      boundary_pressure_postprocessor = this->template find_postprocessor<Postprocess::BoundaryPressures<dim> >();
      boundary_density_postprocessor = this->template find_postprocessor<Postprocess::BoundaryDensities<dim> >();

      AssertThrow(boundary_pressure_postprocessor != NULL,
                  ExcMessage("Could not find BoundaryPressures postprocessor"));
      AssertThrow(boundary_density_postprocessor != NULL,
                  ExcMessage("Could not find BoundaryDensities postprocessor"));

      compute_internal_density_expansions();
      compute_topography_expansions();
      compute_geoid_expansions();

      output_geoid_information();
      return std::pair<std::string,std::string>("Writing geoid file", "");
    }

    template <int dim>
    void
    Geoid<dim>::compute_internal_density_expansions()
    {
      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim> cell_quadrature_formula (this->get_fe()
                                                 .base_element(this->introspection().base_elements.temperature)
                                                 .degree+1); //Need to do the volume integration with this quadrature

      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                                 (&this->get_geometry_model());
      AssertThrow (geometry_model != 0,
                   ExcMessage("The geoid postprocessor is currently only implemented for "
                              "the spherical shell geometry model."));

      // TODO AssertThrow (no_free_surface);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               cell_quadrature_formula,
                               update_values |
                               update_gradients |
                               update_q_points |
                               update_JxW_values);


      //Material model in/out for the gauss quadrature rule evaluations
      typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (cell_quadrature_formula.size()));

      // Some constant that are used several times
      const double inner_radius = geometry_model->inner_radius();
      const double outer_radius = geometry_model->outer_radius();

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

            // for each of the quadrature points, evaluate the
            // density and add its contribution to the spherical harmonics

            for (unsigned int q=0; q<cell_quadrature_formula.size(); ++q)
              {
                const Point<dim> location = fe_values.quadrature_point(q);
                const double density   = out.densities[q];

                internal_density_expansion_surface->add_quadrature_point( location, density, fe_values.JxW(q), outer_radius );
                internal_density_expansion_bottom->add_quadrature_point( location, density, fe_values.JxW(q), inner_radius );
              }

          }

      internal_density_expansion_surface->mpi_sum_coefficients( this->get_mpi_communicator() );
      internal_density_expansion_bottom->mpi_sum_coefficients( this->get_mpi_communicator() );
    }

    template <int dim>
    void
    Geoid<dim>::compute_topography_expansions()
    {
      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                                 (&this->get_geometry_model());
      AssertThrow (geometry_model != 0,
                   ExcMessage("The geoid postprocessor is currently only implemented for "
                              "the spherical shell geometry model."));

      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      const QGauss<dim> center_quadrature_formula(quadrature_degree);  //need to retrieve viscous stress here
      const QGauss<dim-1> face_quadrature_formula(quadrature_degree);  //need to retrieve pressure here

      FEValues<dim> fe_center_values (this->get_mapping(),
                                      this->get_fe(),
                                      center_quadrature_formula,
                                      update_values |
                                      update_gradients |
                                      update_q_points |
                                      update_JxW_values);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        face_quadrature_formula,
                                        update_values |
                                        update_gradients |
                                        update_q_points |
                                        update_JxW_values);

      //AssertThrow( fe_center_values.n_quadrature_points == fe_face_values.n_quadrature_points == 1,
      //             ExcInternalError() );

      //Material model in/out for the gauss quadrature rule evaluations
      typename MaterialModel::Interface<dim>::MaterialModelInputs in_face(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelInputs in_center(fe_center_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out_face(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out_center(fe_center_values.n_quadrature_points, this->n_compositional_fields());
      std::vector<std::vector<double> > face_composition_values (this->n_compositional_fields(),std::vector<double> (fe_face_values.n_quadrature_points /*q_point*/));
      std::vector<std::vector<double> > center_composition_values (this->n_compositional_fields(),std::vector<double> (fe_center_values.n_quadrature_points /*q_point*/));

      const double outer_radius = geometry_model->outer_radius();
      const double inner_radius = geometry_model->inner_radius();

      const double scale_surface_area = (dim == 2 ? outer_radius : outer_radius * outer_radius );
      const double scale_bottom_area = (dim == 2 ? inner_radius : inner_radius * inner_radius );

      const double surface_pressure = boundary_pressure_postprocessor->pressure_at_top();
      const double bottom_pressure = boundary_pressure_postprocessor->pressure_at_bottom();

      Tensor<1,dim> gravity;
      Tensor<1,dim> gravity_direction;

      // loop over all of the surface cells and if one less than h/3 away from
      // the top surface, evaluate the stress at its center
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {

            // see if the cell is at the *top* or *bottom* boundary
            bool surface_cell = false;
            bool bottom_cell = false;

            unsigned int f = 0;
            for (; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                  {
                    surface_cell = true;
                    break;
                  }
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) > (outer_radius - inner_radius - cell->face(f)->minimum_vertex_distance()/3))
                  {
                    bottom_cell = true;
                    break;
                  }
              }

            if (surface_cell || bottom_cell)
              {
                fe_face_values.reinit (cell, f);

                // get the various components of the solution, then
                // evaluate the material properties there
                fe_face_values[this->introspection().extractors.temperature]
                .get_function_values (this->get_solution(), in_face.temperature);
                fe_face_values[this->introspection().extractors.pressure]
                .get_function_values (this->get_solution(), in_face.pressure);
                fe_face_values[this->introspection().extractors.velocities]
                .get_function_symmetric_gradients (this->get_solution(), in_face.strain_rate);

                in_face.position = fe_face_values.get_quadrature_points();

                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  fe_face_values[this->introspection().extractors.compositional_fields[c]]
                  .get_function_values(this->get_solution(),
                                       face_composition_values[c]);
                for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
                  {
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      in_face.composition[i][c] = face_composition_values[c][i];
                  }

                this->get_material_model().evaluate(in_face, out_face);

                fe_center_values.reinit(cell);

                // get the various components of the solution, then
                // evaluate the material properties there
                fe_center_values[this->introspection().extractors.temperature]
                .get_function_values (this->get_solution(), in_center.temperature);
                fe_center_values[this->introspection().extractors.pressure]
                .get_function_values (this->get_solution(), in_center.pressure);
                fe_center_values[this->introspection().extractors.velocities]
                .get_function_symmetric_gradients (this->get_solution(), in_center.strain_rate);

                in_center.position = fe_center_values.get_quadrature_points();

                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  fe_center_values[this->introspection().extractors.compositional_fields[c]]
                  .get_function_values(this->get_solution(),
                                       center_composition_values[c]);
                for (unsigned int i=0; i<fe_center_values.n_quadrature_points; ++i)
                  {
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      in_face.composition[i][c] = center_composition_values[c][i];
                  }

                this->get_material_model().evaluate(in_center, out_center);

                //Calculate the average pressure, density, gravity over the relevant face of the cell
                double face_area = 0.;
                double face_pressure = 0.;
                double face_density = 0.;
                double face_gravity = 0.;
                for ( unsigned int q = 0; q < face_quadrature_formula.size(); ++q)
                  {
                    face_area += fe_face_values.JxW(q);
                    face_pressure += in_face.pressure[q] * fe_face_values.JxW(q);
                    face_density += out_face.densities[q] * fe_face_values.JxW(q);

                    gravity = this->get_gravity_model().gravity_vector(fe_face_values.quadrature_point(q));
                    face_gravity += gravity.norm() * fe_face_values.JxW(q);
                  }
                face_pressure /= face_area;
                face_density /= face_area;
                face_gravity /= face_area;

                //Calculate the average viscous stress over the volume of the surface cell.
                //The volume avoids certain polynomial overshoot problems with evaluating
                //second order shape functions at the surface
                double cell_volume = 0.;
                double cell_surface_normal_viscous_stress = 0.;
                for ( unsigned int q = 0; q < center_quadrature_formula.size(); ++q)
                  {
                    const double viscosity = out_center.viscosities[q];

                    const SymmetricTensor<2,dim> strain_rate = in_center.strain_rate[q] - 1./3 * trace(in_center.strain_rate[q]) * unit_symmetric_tensor<dim>();
                    const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                    gravity = this->get_gravity_model().gravity_vector(fe_center_values.quadrature_point(q));
                    gravity_direction = (gravity.norm() == 0.0 ? Tensor<1,dim>() : gravity/gravity.norm());

                    cell_volume += fe_center_values.JxW(q);
                    cell_surface_normal_viscous_stress += gravity_direction * (shear_stress * gravity_direction) * fe_center_values.JxW(q);
                  }
                cell_surface_normal_viscous_stress /= cell_volume;

                // if this is a cell at the surface, add the topography to
                // the topography expansion
                if (surface_cell)
                  {
                    const double dynamic_pressure = face_pressure - surface_pressure;
                    const double sigma_rr           = cell_surface_normal_viscous_stress - dynamic_pressure;
                    const double dynamic_topography = - sigma_rr / face_gravity / (face_density - density_above);

                    // Add topography contribution
                    const Point<dim> location = cell->face(f)->center();
                    surface_topography_expansion->add_quadrature_point(location/location.norm(), dynamic_topography, face_area/scale_surface_area, 1.0);
                  }

                // if this is a cell at the bottom, add the topography to
                // the bottom expansion
                if (bottom_cell)
                  {
                    const double dynamic_pressure   = face_pressure - bottom_pressure;
                    const double sigma_rr           = cell_surface_normal_viscous_stress - dynamic_pressure;
                    const double dynamic_topography = - sigma_rr / face_gravity / (face_density - density_below);

                    // Add topography contribution
                    const Point<dim> location = cell->face(f)->center();
                    bottom_topography_expansion->add_quadrature_point(location/location.norm(), dynamic_topography, face_area/scale_bottom_area, 1.0);
                  }
              }
          }
      surface_topography_expansion->mpi_sum_coefficients( this->get_mpi_communicator() );
      bottom_topography_expansion->mpi_sum_coefficients( this->get_mpi_communicator() );
    }


    template <int dim>
    void
    Geoid<dim>::compute_geoid_expansions()
    {
      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                                 (&this->get_geometry_model());
      AssertThrow (geometry_model != 0,
                   ExcMessage("The geoid postprocessor is currently only implemented for "
                              "the spherical shell geometry model."));

      const double outer_radius = geometry_model->outer_radius();
      const double inner_radius = geometry_model->inner_radius();
      const double G = constants::big_g;
      const double surface_density = boundary_density_postprocessor->density_at_top();
      const double bottom_density = boundary_density_postprocessor->density_at_bottom();

      const double delta_rho_top = surface_density-density_above;
      const double delta_rho_bottom = density_below-bottom_density;

      std::vector<double> s(max_degree+1);

      std::vector<double> surface_potential_at_surface(max_degree+1);
      std::vector<double> surface_potential_at_bottom(max_degree+1);
      std::vector<double> bottom_potential_at_bottom(max_degree+1);
      std::vector<double> bottom_potential_at_surface(max_degree+1);

      if ( dim == 3)
        {
          for ( unsigned int l = 2; l <= max_degree; ++l)
            {
              s[l] = 1.0;
              bottom_potential_at_surface[l] = -G * std::pow(inner_radius/outer_radius, static_cast<double>(l+1) ) * inner_radius * delta_rho_bottom;
              surface_potential_at_bottom[l] = -G * std::pow(inner_radius/outer_radius, static_cast<double>(l) ) * outer_radius * delta_rho_top;
            }

          surface_potential_expansion->clear();
          surface_potential_expansion->sadd( 1.0, -G , *internal_density_expansion_surface );
          surface_potential_expansion->sadd( 1.0, -G * outer_radius * delta_rho_top, *surface_topography_expansion );
          surface_potential_expansion->sadd( s, bottom_potential_at_surface, *bottom_topography_expansion );

          bottom_potential_expansion->clear();
          bottom_potential_expansion->sadd( 1.0, -G , *internal_density_expansion_bottom );
          bottom_potential_expansion->sadd( 1.0, -G * inner_radius * delta_rho_bottom, *bottom_topography_expansion );
          bottom_potential_expansion->sadd( s, surface_potential_at_bottom, *surface_topography_expansion );
        }
      else
        {
          const double gravity_constant = 4./3. * G;
          for ( unsigned int n = 2; n <= max_degree; ++n)
            {

              s[n] = 1.0;
              bottom_potential_at_surface[n] = -gravity_constant * std::pow(inner_radius/outer_radius, static_cast<double>(n)) *
                                               delta_rho_bottom * inner_radius / static_cast<double>(n);
              surface_potential_at_bottom[n] = -gravity_constant * std::pow(inner_radius/outer_radius, static_cast<double>(n)) *
                                               delta_rho_top * outer_radius / static_cast<double>(n);
              bottom_potential_at_bottom[n] = -gravity_constant * inner_radius * delta_rho_bottom / static_cast<double>(n);
              surface_potential_at_surface[n] = -gravity_constant * outer_radius * delta_rho_top / static_cast<double>(n);
            }

          surface_potential_expansion->clear();
          surface_potential_expansion->sadd( 1.0, -gravity_constant , *internal_density_expansion_surface );
          surface_potential_expansion->sadd( s, surface_potential_at_surface, *surface_topography_expansion );
          surface_potential_expansion->sadd( s, bottom_potential_at_surface, *bottom_topography_expansion );

          bottom_potential_expansion->clear();
          bottom_potential_expansion->sadd( 1.0, -gravity_constant , *internal_density_expansion_bottom );
          bottom_potential_expansion->sadd( s, bottom_potential_at_bottom, *bottom_topography_expansion );
          bottom_potential_expansion->sadd( s, surface_potential_at_bottom, *surface_topography_expansion );
        }
    }

    template <int dim>
    void
    Geoid<dim>::output_geoid_information()
    {
      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                                 (&this->get_geometry_model());
      AssertThrow (geometry_model != 0,
                   ExcMessage("The geoid postprocessor is currently only implemented for "
                              "the spherical shell geometry model."));

      const std::string filename = this->get_output_directory() +
                                   "geoid." +
                                   dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

      const double gravity_at_surface = this->get_gravity_model().gravity_vector(
                                          ( this->get_geometry_model().representative_point( 0.0 ) ) ).norm();
      const double gravity_at_bottom = this->get_gravity_model().gravity_vector(
                                         this->get_geometry_model().representative_point(
                                           this->get_geometry_model().maximal_depth() )  ).norm();

      // On process 0 write output file
      if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                                   this->get_time() / year_in_seconds :
                                                   this->get_time());
          std::ofstream file (filename.c_str());

          file << "# Timestep Maximum_degree Time" << std::endl;
          file << this->get_timestep_number() << " " << max_degree << " " << time_in_years_or_seconds << std::endl;
          if ( dim == 3 )
            {
              file << "# degree   order geoid_surf_s geoid_surf_c  geoid_cmb_s  cmb_geoid_c  dens_surf_s  dens_surf_c   dens_cmb_s   dens_cmb_c  topo_surf_s  topo_surf_c   topo_cmb_s   topo_cmb_c" << std::endl;
              // Write the solution to an output file
              for (unsigned int l=2, k=0; l <= max_degree; ++l)
                {
                  for (unsigned int m = 0; m <= l; ++m, ++k)
                    {
                      file << std::setw(8) << l << std::setw(8) << m
                           << std::setw(13) << surface_potential_expansion->get_coefficients().sine_coefficients[k]/gravity_at_surface
                           << std::setw(13) << surface_potential_expansion->get_coefficients().cosine_coefficients[k]/gravity_at_surface
                           << std::setw(13) << bottom_potential_expansion->get_coefficients().sine_coefficients[k]/gravity_at_bottom
                           << std::setw(13) << bottom_potential_expansion->get_coefficients().cosine_coefficients[k]/gravity_at_bottom
                           << std::setw(13) << internal_density_expansion_surface->get_coefficients().sine_coefficients[k]
                           << std::setw(13) << internal_density_expansion_surface->get_coefficients().cosine_coefficients[k]
                           << std::setw(13) << internal_density_expansion_bottom->get_coefficients().sine_coefficients[k]
                           << std::setw(13) << internal_density_expansion_bottom->get_coefficients().cosine_coefficients[k]
                           << std::setw(13) << surface_topography_expansion->get_coefficients().sine_coefficients[k]
                           << std::setw(13) << surface_topography_expansion->get_coefficients().cosine_coefficients[k]
                           << std::setw(13) << bottom_topography_expansion->get_coefficients().sine_coefficients[k]
                           << std::setw(13) << bottom_topography_expansion->get_coefficients().cosine_coefficients[k] << std::endl;
                    }
                }
            }
          else
            {
              file << "# degree geoid_surf_s geoid_surf_c  geoid_cmb_s  cmb_geoid_c  dens_surf_s  dens_surf_c   dens_cmb_s   dens_cmb_c  topo_surf_s  topo_surf_c   topo_cmb_s   topo_cmb_c" << std::endl;
              for (unsigned int n=2; n <= max_degree; ++n)
                {
                  file << std::setw(8) << n
                       << std::setw(13) << surface_potential_expansion->get_coefficients().sine_coefficients[n]/gravity_at_surface
                       << std::setw(13) << surface_potential_expansion->get_coefficients().cosine_coefficients[n]/gravity_at_surface
                       << std::setw(13) << bottom_potential_expansion->get_coefficients().sine_coefficients[n]/gravity_at_bottom
                       << std::setw(13) << bottom_potential_expansion->get_coefficients().cosine_coefficients[n]/gravity_at_bottom
                       << std::setw(13) << internal_density_expansion_surface->get_coefficients().sine_coefficients[n]
                       << std::setw(13) << internal_density_expansion_surface->get_coefficients().cosine_coefficients[n]
                       << std::setw(13) << internal_density_expansion_bottom->get_coefficients().sine_coefficients[n]
                       << std::setw(13) << internal_density_expansion_bottom->get_coefficients().cosine_coefficients[n]
                       << std::setw(13) << surface_topography_expansion->get_coefficients().sine_coefficients[n]
                       << std::setw(13) << surface_topography_expansion->get_coefficients().cosine_coefficients[n]
                       << std::setw(13) << bottom_topography_expansion->get_coefficients().sine_coefficients[n]
                       << std::setw(13) << bottom_topography_expansion->get_coefficients().cosine_coefficients[n] << std::endl;
                }
            }
        }
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
          prm.declare_entry ("Density below", "8000",
                             Patterns::Double(0),
                             "Density of the fluid beneath the domain, "
                             "most likely that of the iron-alloy core, ");
          prm.declare_entry ("Density above", "0",
                             Patterns::Double(0),
                             "Density of the fluid above the domain, "
                             "most likely either air or water.");
          prm.declare_entry ("Expansion degree", "7",
                             Patterns::Integer (0),
                             "Degree of the spherical harmonic expansion. "
                             "The expansion into spherical harmonics can be "
                             "expensive, especially for high degrees.");
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
          density_below                     = prm.get_double("Density below");
          density_above                     = prm.get_double("Density above");
          max_degree                        = prm.get_integer("Expansion degree");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      internal_density_expansion_surface.reset(new internal::MultipoleExpansion<dim>(max_degree) );
      internal_density_expansion_bottom.reset(new internal::MultipoleExpansion<dim>(max_degree) );
      surface_topography_expansion.reset(new internal::MultipoleExpansion<dim>(max_degree) );
      bottom_topography_expansion.reset(new internal::MultipoleExpansion<dim>(max_degree) );

      surface_potential_expansion.reset(new internal::MultipoleExpansion<dim>(max_degree) );
      bottom_potential_expansion.reset(new internal::MultipoleExpansion<dim>(max_degree) );
    }

    template <int dim>
    std::list<std::string>
    Geoid<dim>::required_other_postprocessors() const
    {
      std::list<std::string> deps;
      deps.push_back("boundary pressures");
      deps.push_back("boundary densities");
      return deps;
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
