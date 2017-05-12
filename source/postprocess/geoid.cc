/*
 Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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


#include <aspect/utilities.h>
#include <aspect/postprocess/geoid.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <aspect/geometry_model/spherical_shell.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::vector<double>,std::vector<double> >
    Geoid<dim>::to_spherical_harmonic_coefficients(const std::vector<std::vector<double> > &spherical_function) const
    {
      std::vector<double> cosi(spherical_function.size(),0);
      std::vector<double> sini(spherical_function.size(),0);
      std::vector<double> coecos;
      std::vector<double> coesin;

      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
        {
          for (unsigned int iord = 0; iord < ideg+1; iord++)
            {
              // do the sphercial harmonic expansion
              for (unsigned int ds_num = 0; ds_num < spherical_function.size(); ds_num++)
                {
                  // normalization after Dahlen and Tromp, 1986, Appendix B.6
                  const std::pair<double,double> sph_harm_vals = aspect::Utilities::real_spherical_harmonic(ideg,iord,spherical_function.at(ds_num).at(0),spherical_function.at(ds_num).at(1));
                  const double cos_component = sph_harm_vals.first; //real / cos part
                  const double sin_component = sph_harm_vals.second; //imaginary / sine part

                  cosi.at(ds_num) = (spherical_function.at(ds_num).at(3) * cos_component);
                  sini.at(ds_num) = (spherical_function.at(ds_num).at(3) * sin_component);
                }
              // integrate the contribution of each spherical infinitesimal
              double cosii = 0;
              double sinii = 0;
              for (unsigned int ds_num = 0; ds_num < spherical_function.size(); ds_num++)
                {
                  cosii += cosi.at(ds_num) * spherical_function.at(ds_num).at(2);
                  sinii += sini.at(ds_num) * spherical_function.at(ds_num).at(2);
                }
              coecos.push_back(cosii);
              coesin.push_back(sinii);
            }
        }
      // sum over each processor
      dealii::Utilities::MPI::sum (coecos,this->get_mpi_communicator(),coecos);
      dealii::Utilities::MPI::sum (coesin,this->get_mpi_communicator(),coesin);

      return std::make_pair(coecos,coesin);
    }

    template <int dim>
    std::pair<std::vector<double>,std::vector<double> >
    Geoid<dim>::density_contribution (const double &outer_radius) const
    {
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      // need to evaluate density contribution of each volume quadrature point
      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_q_points |
                               update_JxW_values |
                               update_gradients);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      // Directly do the global 3D intergral over each quadrature point of every cell (different from traditional way to do layer integral).
      // This work around ASPECT's adpative mesh refinment feature.
      std::vector<double> SH_density_coecos;
      std::vector<double> SH_density_coesin;
      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
        {
          for (unsigned int iord = 0; iord < ideg+1; iord++)
            {
              // initialization of the density contribution integral per degree, order.
              double integrated_density_cos_component = 0;
              double integrated_density_sin_component = 0;

              // loop over all of the cells
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
                    .get_function_values (this->get_solution(), in.velocity);
                    fe_values[this->introspection().extractors.velocities]
                    .get_function_symmetric_gradients (this->get_solution(), in.strain_rate);
                    fe_values[this->introspection().extractors.pressure]
                    .get_function_gradients (this->get_solution(), in.pressure_gradient);

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
                    in.cell = &cell;

                    this->get_material_model().evaluate(in, out);

                    // Compute the integral of the density function
                    // over the cell, by looping over all quadrature points
                    for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                      {
                        // convert coordinates from [x,y,z] to [r, phi, theta]
                        const std_cxx11::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(in.position[q]);

                        // normalization after Dahlen and Tromp, 1986, Appendix B.6
                        const std::pair<double,double> sph_harm_vals = aspect::Utilities::real_spherical_harmonic(ideg,iord,scoord[2],scoord[1]);
                        const double cos_component = sph_harm_vals.first; //real / cos part
                        const double sin_component = sph_harm_vals.second; //imaginary / sine part

                        const double density = out.densities[q];
                        const double r_q = in.position[q].norm();

                        integrated_density_cos_component += density * (1./r_q) * std::pow(r_q/outer_radius,ideg+1) * cos_component * fe_values.JxW(q);
                        integrated_density_sin_component += density * (1./r_q) * std::pow(r_q/outer_radius,ideg+1) * sin_component * fe_values.JxW(q);
                      }
                  }
              SH_density_coecos.push_back(integrated_density_cos_component);
              SH_density_coesin.push_back(integrated_density_sin_component);
            }
        }
      // sum over each processor
      dealii::Utilities::MPI::sum (SH_density_coecos,this->get_mpi_communicator(),SH_density_coecos);
      dealii::Utilities::MPI::sum (SH_density_coesin,this->get_mpi_communicator(),SH_density_coesin);

      return std::make_pair(SH_density_coecos,SH_density_coesin);
    }

    template <int dim>
    std::pair<std::pair<double, std::pair<std::vector<double>,std::vector<double> > >, std::pair<double, std::pair<std::vector<double>,std::vector<double> > > >
    Geoid<dim>::dynamic_topography_contribution(const double &outer_radius,
                                                const double &inner_radius) const
    {
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      const QGauss<dim> quadrature_formula(quadrature_degree); // need to retrieve normal shear stess, pressure, density to calculate dynamic topography here
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree); // need to grab the infinitesimal area of each quadrature points on every boundary face here

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_q_points |
                               update_JxW_values);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_gradients |
                                        update_q_points |
                                        update_JxW_values);

      // Material model in/out for the gauss quadrature rule evaluations
      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelInputs<dim> in_face(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out_face(fe_face_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));
      std::vector<std::vector<double> > face_composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula_face.size()));

      // initialization of the surface and CMB integrated dynamic topography, density, area, and volume, later used to calculate average dynamic topography and boundary layer density
      double integrated_surface_topography = 0;
      double integrated_surface_area = 0;
      double integrated_CMB_topography = 0;
      double integrated_CMB_area = 0;
      double integrated_top_layer_density = 0;
      double integrated_top_layer_volume = 0;
      double integrated_bottom_layer_density = 0;
      double integrated_bottom_layer_volume = 0;

      // vectors to store the location, infitesimal area, and dynamic topography associated with each quadrature point of each surface and bottom cell respectively.
      std::vector<std::pair<Point<dim>,std::pair<double,double> > > surface_stored_values;
      std::vector<std::pair<Point<dim>,std::pair<double,double> > > CMB_stored_values;

      // loop over all of the boundary cells and if one is at
      // surface or CMB, evaluate the stress at its center
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          if (cell->at_boundary())
            {
              // see if the cell is at the *top* boundary or CMB, not just any boundary
              unsigned int top_face_idx = numbers::invalid_unsigned_int;
              unsigned int CMB_face_idx = numbers::invalid_unsigned_int;
              {
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                  {
                    if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                      {
                        // if the cell is at top boundary, assign the correspongding upper face index to top_face_idx but keep the CMB_face_idx still invalid_unsigned_int
                        top_face_idx = f;
                        break;
                      }
                    else if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center())
                             > (this->get_geometry_model().maximal_depth() - cell->face(f)->minimum_vertex_distance()/3.))
                      {
                        // if the cell is at bottom boundary, assign the correspongding lower face index to CMB_face_idx but keep the top_face_idx still invalid_unsigned_int
                        CMB_face_idx = f;
                        break;
                      }
                    else
                      continue;
                  }

                // if the cell is not at the boundary, keep the search loop
                if (top_face_idx == numbers::invalid_unsigned_int && CMB_face_idx == numbers::invalid_unsigned_int)
                  continue;
              }

              // focus on the boundary cell
              fe_values.reinit (cell);
              // get the various components of the solution, then
              // evaluate the material properties there
              fe_values[this->introspection().extractors.temperature]
              .get_function_values (this->get_solution(), in.temperature);
              fe_values[this->introspection().extractors.pressure]
              .get_function_values (this->get_solution(), in.pressure);
              fe_values[this->introspection().extractors.velocities]
              .get_function_values (this->get_solution(), in.velocity);
              fe_values[this->introspection().extractors.velocities]
              .get_function_symmetric_gradients (this->get_solution(), in.strain_rate);
              fe_values[this->introspection().extractors.pressure]
              .get_function_gradients (this->get_solution(), in.pressure_gradient);

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
              in.cell = &cell;

              this->get_material_model().evaluate(in, out);


              // Calculate the average dynamic topography of the cell.
              // For each of the quadrature points, evaluate the density, dynamic pressure, and shear stress in direction of the gravity vector.
              // Compute the integral of the dynamic topograhy function over the entire cell, by looping over all quadrature points.
              double dynamic_topography_x_volume = 0.;
              double cell_volume = 0.;
              double density_x_volume = 0.;

              for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                {
                  const Point<dim> location = fe_values.quadrature_point(q);
                  const double viscosity = out.viscosities[q];
                  const double density   = out.densities[q];

                  const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                  const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                  const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                  const Tensor<1,dim> gravity_direction = gravity/gravity.norm();

                  // Subtract the dynamic pressure to get the normal stress
                  const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(location);
                  const double sigma_rr           = gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure;

                  // Compute dynamic topography either at top or bottom cell
                  double dynamic_topography;
                  if (top_face_idx != numbers::invalid_unsigned_int)
                    dynamic_topography = - sigma_rr / gravity.norm() / (density - density_above);
                  else if (CMB_face_idx != numbers::invalid_unsigned_int)
                    dynamic_topography = sigma_rr / gravity.norm() / (density_below - density);
                  else
                    dynamic_topography = std::numeric_limits<double>::signaling_NaN();

                  Assert(std::isnan(dynamic_topography) == false,
                         ExcMessage("The cells not at the top and bottom boundaries are mistakenly used to calculate dynamic topography contribution to the geoid!"));

                  density_x_volume += density * fe_values.JxW(q);
                  dynamic_topography_x_volume += dynamic_topography * fe_values.JxW(q);
                  cell_volume += fe_values.JxW(q);
                }

              // get the average dynamic topograhy of the cell.
              const double dynamic_topography_cell_average = dynamic_topography_x_volume / cell_volume;

              // Get the boundary face index and check again to make sure the cell is at the boundary.
              unsigned int face_idx;
              if (top_face_idx != numbers::invalid_unsigned_int)
                face_idx = top_face_idx;
              else if (CMB_face_idx != numbers::invalid_unsigned_int)
                face_idx = CMB_face_idx;
              else
                face_idx = numbers::invalid_unsigned_int;

              Assert(face_idx != numbers::invalid_unsigned_int,
                     ExcMessage("The cells not at the top and bottom boundaries are mistakenly used to calculate dynamic topography contribution to the geoid!"));

              // focus on the boundary cell's upper face if on the top boundary and lower face if on the bottom boundary
              fe_face_values.reinit(cell, face_idx);
              // get the various components of the solution, then
              // evaluate the material properties there
              fe_face_values[this->introspection().extractors.temperature]
              .get_function_values (this->get_solution(), in_face.temperature);
              fe_face_values[this->introspection().extractors.pressure]
              .get_function_values (this->get_solution(), in_face.pressure);
              fe_face_values[this->introspection().extractors.velocities]
              .get_function_values (this->get_solution(), in_face.velocity);
              fe_face_values[this->introspection().extractors.velocities]
              .get_function_symmetric_gradients (this->get_solution(), in_face.strain_rate);
              fe_face_values[this->introspection().extractors.pressure]
              .get_function_gradients (this->get_solution(), in_face.pressure_gradient);

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

              // if the cell at top boundary, calculate dynamic topography and add the results into surface dynamic topogaphy storage vector
              if (top_face_idx != numbers::invalid_unsigned_int)
                {
                  // Although we calculate the average dynamic topograhy of the cell's boundary face, we store this same dynamic topography with every
                  // surface quadrature point's associated area into the storage vector. The reason to do this is that later in the spherical
                  // harmonic expansion, we will calculate sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2. The accuracy of this relation
                  // is improved as infinitesimal_area is closer to zero, so using every surface quadrature point's associated area of each
                  // surface cell will lead to better accuracy in spherical harmonic expansion, especially in the coarse mesh.
                  double face_area = 0.;
                  for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                    {
                      surface_stored_values.push_back (std::make_pair(in_face.position[q], std::make_pair(fe_face_values.JxW(q),dynamic_topography_cell_average)));
                      face_area += fe_face_values.JxW(q);
                    }
                  // Contribute the density, volume, dynamic topography, and surface area of the cell to the integration of the current processor.
                  integrated_top_layer_density += density_x_volume;
                  integrated_top_layer_volume += cell_volume;
                  integrated_surface_topography += dynamic_topography_cell_average * face_area;
                  integrated_surface_area += face_area;
                }

              // if the cell at bottom boundary, calculate dynamic topography and add the results into bottom dynamic topogaphy storage vector
              if (CMB_face_idx != numbers::invalid_unsigned_int)
                {
                  // Although we calculate the average dynamic topograhy of the cell's boundary face, we store this same dynamic topography with every
                  // bottom quadrature point's associated area into the storage vector. The reason to do this is that later in the spherical
                  // harmonic expansion, we will calculate sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2. The accuracy of this relation
                  // is improved as infinitesimal_area is closer to zero, so using every bottom quadrature point's associated area of each
                  // bottom cell will lead to better accuracy in spherical harmonic expansion, especially in the coarse mesh.
                  double face_area = 0.;
                  for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                    {
                      CMB_stored_values.push_back (std::make_pair(in_face.position[q], std::make_pair(fe_face_values.JxW(q),dynamic_topography_cell_average)));
                      face_area += fe_face_values.JxW(q);
                    }
                  // Contribute the density, volume, dynamic topography, and CMB area of the cell to the integration of the current processor.
                  integrated_bottom_layer_density += density_x_volume;
                  integrated_bottom_layer_volume += cell_volume;
                  integrated_CMB_topography += dynamic_topography_cell_average * face_area;
                  integrated_CMB_area += face_area;
                }
            }

      // Calculate the surface and CMB weighted average dynamic topography
      const double surface_average_topography = dealii::Utilities::MPI::sum (integrated_surface_topography,this->get_mpi_communicator()) / dealii::Utilities::MPI::sum (integrated_surface_area,this->get_mpi_communicator());
      const double CMB_average_topography = dealii::Utilities::MPI::sum (integrated_CMB_topography,this->get_mpi_communicator()) / dealii::Utilities::MPI::sum (integrated_CMB_area,this->get_mpi_communicator());

      // Calculate the surface and CMB layer weighted average density.
      const double top_layer_average_density = dealii::Utilities::MPI::sum (integrated_top_layer_density,this->get_mpi_communicator()) / dealii::Utilities::MPI::sum (integrated_top_layer_volume,this->get_mpi_communicator());
      const double bottom_layer_average_density = dealii::Utilities::MPI::sum (integrated_bottom_layer_density,this->get_mpi_communicator()) / dealii::Utilities::MPI::sum (integrated_bottom_layer_volume,this->get_mpi_communicator());

      // Subtract the average dynamic topography.
      // Transfer the geocentric coordinates to the spherical coordinates.
      // Prepare the value of spherical infinitesimal, i.e., sin(theta)*d_theta*d_phi, for the later spherical harmonic expansion.
      std::vector<std::vector<double> > surface_topo_spherical_function; // store theta, phi, spherical infinitesimal, and surface dynamic topography
      std::vector<std::vector<double> > CMB_topo_spherical_function; // store theta, phi, spherical infinitesimal, and CMB dynamic topography
      for (unsigned int i=0; i<surface_stored_values.size(); ++i)
        {
          surface_stored_values.at(i).second.second -= surface_average_topography;
          const std_cxx11::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(surface_stored_values.at(i).first);
          const double theta = scoord[2];
          const double phi = scoord[1];
          // calculate spherical infinitesimal sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2
          const double infinitesimal = surface_stored_values.at(i).second.first/(outer_radius*outer_radius);
          // assign the spherical function containing theta, phi, spherical infinitesimal, and surface dynamic topography
          std::vector<double> tmp;
          tmp.push_back(theta);
          tmp.push_back(phi);
          tmp.push_back(infinitesimal);
          tmp.push_back(surface_stored_values.at(i).second.second);
          surface_topo_spherical_function.push_back(tmp);
        }
      for (unsigned int i=0; i<CMB_stored_values.size(); ++i)
        {
          CMB_stored_values.at(i).second.second -= CMB_average_topography;
          const std_cxx11::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(CMB_stored_values.at(i).first);
          const double theta = scoord[2];
          const double phi = scoord[1];
          // calculate spherical infinitesimal sin(theta)*d_theta*d_phi by infinitesimal_area/radius^2
          const double infinitesimal = CMB_stored_values.at(i).second.first/(inner_radius*inner_radius);
          // assign the spherical function containing theta, phi, spherical infinitesimal, and CMB dynamic topography
          std::vector<double> tmp;
          tmp.push_back(theta);
          tmp.push_back(phi);
          tmp.push_back(infinitesimal);
          tmp.push_back(CMB_stored_values.at(i).second.second);
          CMB_topo_spherical_function.push_back(tmp);
        }

      std::pair<double, std::pair<std::vector<double>,std::vector<double> > > SH_surface_dyna_topo_coes;
      SH_surface_dyna_topo_coes = std::make_pair(top_layer_average_density,to_spherical_harmonic_coefficients(surface_topo_spherical_function));
      std::pair<double, std::pair<std::vector<double>,std::vector<double> > > SH_CMB_dyna_topo_coes;
      SH_CMB_dyna_topo_coes = std::make_pair(bottom_layer_average_density,to_spherical_harmonic_coefficients(CMB_topo_spherical_function));
      return std::make_pair(SH_surface_dyna_topo_coes,SH_CMB_dyna_topo_coes);
    }

    template <int dim>
    std::pair<std::string,std::string>
    Geoid<dim>::execute (TableHandler &)
    {
      // Current geoid code only works for spherical shell geometry
      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                                 (&this->get_geometry_model());
      AssertThrow (geometry_model != 0 && dim == 3,
                   ExcMessage("The geoid postprocessor is currently only implemented for the 3D spherical shell geometry model."));

      // Get the value of the outer radius and inner radius
      const double outer_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                                  (this->get_geometry_model()).outer_radius();
      const double inner_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                                  (this->get_geometry_model()).inner_radius();

      // Get the value of the surface gravity acceleration from the gravity model
      Point<dim> surface_point;
      surface_point[0] = outer_radius;
      const double surface_gravity = this->get_gravity_model().gravity_vector(surface_point).norm();

      // Get the value of the universal gravitational constant
      const double G = aspect::constants::big_g;

      // Get the spherical harmonic coefficients of the density contribution.
      std::pair<std::vector<double>,std::vector<double> > SH_density_coes = density_contribution(outer_radius);

      // Get the spherical harmonic coefficients of the surface and CMB dynamic topography
      std::pair<std::pair<double, std::pair<std::vector<double>,std::vector<double> > >, std::pair<double, std::pair<std::vector<double>,std::vector<double> > > > SH_dyna_topo_coes;
      SH_dyna_topo_coes = dynamic_topography_contribution(outer_radius,inner_radius);
      std::pair<double, std::pair<std::vector<double>,std::vector<double> > > SH_surface_dyna_topo_coes = SH_dyna_topo_coes.first;
      std::pair<double, std::pair<std::vector<double>,std::vector<double> > > SH_CMB_dyna_topo_coes = SH_dyna_topo_coes.second;

      // Get the density contrast at the surface and CMB
      const double surface_delta_rho =  SH_surface_dyna_topo_coes.first - density_above;
      const double CMB_delta_rho = density_below - SH_CMB_dyna_topo_coes.first;

      // Compute the spherical harmonic coefficients of geoid anomaly
      std::vector<double> density_anomaly_contribution_coecos; // a vector to store cos terms of density anomaly contribution SH coefficients
      std::vector<double> density_anomaly_contribution_coesin; // a vector to store sin terms of density anomaly contribution SH coefficients
      std::vector<double> surface_dyna_topo_contribution_coecos; // a vector to store cos terms of surface dynamic topography contribution SH coefficients
      std::vector<double> surface_dyna_topo_contribution_coesin; // a vector to store sin terms of surface dynamic topography contribution SH coefficients
      std::vector<double> CMB_dyna_topo_contribution_coecos; // a vector to store cos terms of CMB dynamic topography contribution SH coefficients
      std::vector<double> CMB_dyna_topo_contribution_coesin; // a vector to store sin terms of CMB dynamic topography contribution SH coefficients
      std::vector<double> geoid_coecos;  // a vector to store cos terms of geoid anomaly SH coefficients
      std::vector<double> geoid_coesin;  // a vector to store sin terms of geoid anomaly SH coefficients

      // First compute the spherical harmonic contributions from density anomaly, surface dynamic topography and CMB dynamic topography
      int ind = 0; // coefficients index
      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
        {
          for (unsigned int iord = 0; iord < ideg+1; iord++)
            {
              double coecos_density_anomaly = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1))) * SH_density_coes.first.at(ind);
              double coesin_density_anomaly = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1))) * SH_density_coes.second.at(ind);
              density_anomaly_contribution_coecos.push_back(coecos_density_anomaly);
              density_anomaly_contribution_coesin.push_back(coesin_density_anomaly);

              double coecos_surface_dyna_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                                * surface_delta_rho*SH_surface_dyna_topo_coes.second.first.at(ind)*outer_radius;
              double coesin_surface_dyna_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                                * surface_delta_rho*SH_surface_dyna_topo_coes.second.second.at(ind)*outer_radius;
              surface_dyna_topo_contribution_coecos.push_back(coecos_surface_dyna_topo);
              surface_dyna_topo_contribution_coesin.push_back(coesin_surface_dyna_topo);

              double coecos_CMB_dyna_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                            * CMB_delta_rho*SH_CMB_dyna_topo_coes.second.first.at(ind)*inner_radius*std::pow(inner_radius/outer_radius,ideg+1);
              double coesin_CMB_dyna_topo = (4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1)))
                                            * CMB_delta_rho*SH_CMB_dyna_topo_coes.second.second.at(ind)*inner_radius*std::pow(inner_radius/outer_radius,ideg+1);
              CMB_dyna_topo_contribution_coecos.push_back(coecos_CMB_dyna_topo);
              CMB_dyna_topo_contribution_coesin.push_back(coesin_CMB_dyna_topo);

              ++ind;
            }
        }

      // Then sum the three contributions together to get the spherical harmonic coefficients of geoid anomaly
      ind = 0; // coefficients index
      for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
        {
          for (unsigned int iord = 0; iord < ideg+1; iord++)
            {
              geoid_coecos.push_back(density_anomaly_contribution_coecos.at(ind)+surface_dyna_topo_contribution_coecos.at(ind)+CMB_dyna_topo_contribution_coecos.at(ind));
              geoid_coesin.push_back(density_anomaly_contribution_coesin.at(ind)+surface_dyna_topo_contribution_coesin.at(ind)+CMB_dyna_topo_contribution_coesin.at(ind));
              ind += 1;
            }
        }

      const QMidpoint<dim-1> quadrature_formula_face_center; // need to get each surface cell's upper face midpoint location to calculate geoid in grid from spherical harmonics
      Assert(quadrature_formula_face_center.size() == 1, ExcInternalError());
      FEFaceValues<dim> fe_face_center_values (this->get_mapping(),
                                               this->get_fe(),
                                               quadrature_formula_face_center,
                                               update_values |
                                               update_q_points|
                                               update_JxW_values);

      // define a vector to store the location of the cells along the surface
      std::vector<Point<dim> > surface_cell_locations;

      // loop over all the cells to get the locations of the surface cells to prepare for the geoid computation.
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          if (cell->at_boundary())
            {
              // if the cell is at the *top* boundary, store the cell's upper face midpoint location
              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3.)
                  {
                    fe_face_center_values.reinit(cell,f);
                    const Point<dim> midpoint_at_top_face = fe_face_center_values.get_quadrature_points().at(0);
                    surface_cell_locations.push_back(midpoint_at_top_face);
                    break;
                  }
            }

      // Transfer the geocentric coordinates of the surface cells to the surface spherical coordinates(theta,phi)
      std::vector<std::pair<double,double> > surface_cell_spherical_coordinates;
      for (unsigned int i=0; i<surface_cell_locations.size(); ++i)
        {
          const std_cxx11::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(surface_cell_locations.at(i));
          const double phi = scoord[1];
          const double theta = scoord[2];
          surface_cell_spherical_coordinates.push_back(std::make_pair(theta,phi));
        }

      // Compute the grid geoid anomaly based on spherical harmonics
      std::vector<double> geoid_anomaly;
      for (unsigned int i=0; i<surface_cell_spherical_coordinates.size(); ++i)
        {
          int ind = 0;
          double geoid_value = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
            {
              for (unsigned int iord = 0; iord < ideg+1; iord++)
                {
                  // normalization after Dahlen and Tromp, 1986, Appendix B.6
                  const std::pair<double,double> sph_harm_vals = aspect::Utilities::real_spherical_harmonic(ideg,iord,surface_cell_spherical_coordinates.at(i).first,surface_cell_spherical_coordinates.at(i).second);
                  const double cos_component = sph_harm_vals.first; //real / cos part
                  const double sin_component = sph_harm_vals.second; //imaginary / sine part

                  geoid_value += geoid_coecos.at(ind)*cos_component+geoid_coesin.at(ind)*sin_component;
                  ++ind;
                }
            }
          geoid_anomaly.push_back(geoid_value);
        }

      // The user can get the spherical harmonic coefficients of the density anomaly contribution if needed
      if (also_output_density_anomaly_contribution_SH_coes == true)
        {
          // have a stream into which we write the SH coefficients data from density anomaly contribution.
          // The text stream is then later sent to processor 0
          std::ostringstream output_density_anomaly_contribution_SH_coes;

          // Prepare the output SH coefficients data from density anomaly contribution.
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
            {
              for (unsigned int iord = 0; iord < ideg+1; iord++)
                {
                  output_density_anomaly_contribution_SH_coes << ideg
                                                              << ' '
                                                              << iord
                                                              << ' '
                                                              << density_anomaly_contribution_coecos.at(SH_coes_ind)
                                                              << ' '
                                                              << density_anomaly_contribution_coesin.at(SH_coes_ind)
                                                              << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string density_anomaly_contribution_SH_coes_filename = this->get_output_directory() +
                                                                            "density_anomaly_contribution_SH_coefficients." +
                                                                            dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the SH coefficients from density anomaly contribution, we only need to stop by the processor 0 to get the data.
          // On processor 0, collect all the data and put them into the output density anomaly contribution SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream density_anomaly_contribution_SH_coes_file (density_anomaly_contribution_SH_coes_filename.c_str());
              density_anomaly_contribution_SH_coes_file << "# "
                                                        << "degree order cosine_coefficient sine_coefficient"
                                                        << std::endl;

              // write out the data on processor 0
              density_anomaly_contribution_SH_coes_file << output_density_anomaly_contribution_SH_coes.str();
            }
        }

      // The user can get the spherical harmonic coefficients of the surface dynamic topography contribution if needed
      if (also_output_surface_dynamic_topo_contribution_SH_coes == true)
        {
          // have a stream into which we write the SH coefficients data from surface dynamic topography contribution.
          // The text stream is then later sent to processor 0
          std::ostringstream output_surface_dynamic_topo_contribution_SH_coes;

          // Prepare the output SH coefficients data from surface dynamic topography contribution.
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
            {
              for (unsigned int iord = 0; iord < ideg+1; iord++)
                {
                  output_surface_dynamic_topo_contribution_SH_coes << ideg
                                                                   << ' '
                                                                   << iord
                                                                   << ' '
                                                                   << surface_dyna_topo_contribution_coecos.at(SH_coes_ind)
                                                                   << ' '
                                                                   << surface_dyna_topo_contribution_coesin.at(SH_coes_ind)
                                                                   << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string surface_dynamic_topo_contribution_SH_coes_filename = this->get_output_directory() +
                                                                                 "surface_dynamic_topography_contribution_SH_coefficients." +
                                                                                 dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the SH coefficients from surface dynamic topography contribution, we only need to stop by the processor 0
          // to get the data. On processor 0, collect all the data and put them into the output surface dynamic topography contribution SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream surface_dynamic_topo_contribution_SH_coes_file (surface_dynamic_topo_contribution_SH_coes_filename.c_str());
              surface_dynamic_topo_contribution_SH_coes_file << "# "
                                                             << "degree order cosine_coefficient sine_coefficient"
                                                             << std::endl;
              std::ostringstream output_surface_delta_rho;
              output_surface_delta_rho << surface_delta_rho;
              surface_dynamic_topo_contribution_SH_coes_file << "surface density contrast(kg/m^3): "
                                                             << output_surface_delta_rho.str()
                                                             << std::endl;
              // write out the data on processor 0
              surface_dynamic_topo_contribution_SH_coes_file << output_surface_dynamic_topo_contribution_SH_coes.str();
            }
        }

      // The user can get the spherical harmonic coefficients of the CMB dynamic topography contribution if needed
      if (also_output_CMB_dynamic_topo_contribution_SH_coes == true)
        {
          // have a stream into which we write the SH coefficients data from CMB dynamic topography contribution.
          // The text stream is then later sent to processor 0
          std::ostringstream output_CMB_dynamic_topo_contribution_SH_coes;

          // Prepare the output SH coefficients data from CMB dynamic topography contribution.
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
            {
              for (unsigned int iord = 0; iord < ideg+1; iord++)
                {
                  output_CMB_dynamic_topo_contribution_SH_coes << ideg
                                                               << ' '
                                                               << iord
                                                               << ' '
                                                               << CMB_dyna_topo_contribution_coecos.at(SH_coes_ind)
                                                               << ' '
                                                               << CMB_dyna_topo_contribution_coesin.at(SH_coes_ind)
                                                               << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string CMB_dynamic_topo_contribution_SH_coes_filename = this->get_output_directory() +
                                                                             "CMB_dynamic_topography_contribution_SH_coefficients." +
                                                                             dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the SH coefficients from CMB dynamic topography contribution, we only need to stop by the processor 0
          // to get the data. On processor 0, collect all the data and put them into the output CMB dynamic topography contribution SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream CMB_dynamic_topo_contribution_SH_coes_file (CMB_dynamic_topo_contribution_SH_coes_filename.c_str());
              CMB_dynamic_topo_contribution_SH_coes_file << "# "
                                                         << "degree order cosine_coefficient sine_coefficient"
                                                         << std::endl;
              std::ostringstream output_CMB_delta_rho;
              output_CMB_delta_rho << CMB_delta_rho;
              CMB_dynamic_topo_contribution_SH_coes_file << "CMB density contrast(kg/m^3): "
                                                         << output_CMB_delta_rho.str()
                                                         << std::endl;
              // write out the data on processor 0
              CMB_dynamic_topo_contribution_SH_coes_file << output_CMB_dynamic_topo_contribution_SH_coes.str();
            }
        }

      // The user can get the spherical harmonic coefficients of the geoid anomaly if needed.
      if (also_output_geoid_anomaly_SH_coes == true)
        {
          // have a stream into which we write the geoid anomaly SH coefficients data.
          // The text stream is then later sent to processor 0
          std::ostringstream output_geoid_anomaly_SH_coes;

          // Prepare the output geoid anomaly SH coefficients data
          unsigned int SH_coes_ind = 0;
          for (unsigned int ideg =  min_degree; ideg < max_degree+1; ideg++)
            {
              for (unsigned int iord = 0; iord < ideg+1; iord++)
                {
                  output_geoid_anomaly_SH_coes << ideg
                                               << ' '
                                               << iord
                                               << ' '
                                               << geoid_coecos.at(SH_coes_ind)
                                               << ' '
                                               << geoid_coesin.at(SH_coes_ind)
                                               << std::endl;
                  ++SH_coes_ind;
                }
            }

          const std::string geoid_anomaly_SH_coes_filename = this->get_output_directory() +
                                                             "geoid_anomaly_SH_coefficients." +
                                                             dealii::Utilities::int_to_string(this->get_timestep_number(), 5);

          // Because each processor already held all the geoid anomaly SH coefficients, we only need to stop by the processor 0 to get the data.
          // On processor 0, collect all the data and put them into the output geoid anomaly SH coefficients file.
          if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              std::ofstream geoid_anomaly_SH_coes_file (geoid_anomaly_SH_coes_filename.c_str());
              geoid_anomaly_SH_coes_file << "# "
                                         << "degree order cosine_coefficient sine_coefficient"
                                         << std::endl;

              // write out the data on processor 0
              geoid_anomaly_SH_coes_file << output_geoid_anomaly_SH_coes.str();
            }
        }

      // have a stream into which we write the geoid height data. the text stream is then
      // later sent to processor 0
      std::ostringstream output;

      // Prepare the output data
      if (output_in_lat_lon == true)
        {
          std::vector<std::pair<std::pair<double,double>,double> > stored_values_lon_lat;
          double lon, lat;
          for (unsigned int i=0; i<surface_cell_spherical_coordinates.size(); ++i)
            {
              // Transfer the spherical coordinates to geographical coordinates
              lat = 90. - surface_cell_spherical_coordinates.at(i).first*(180./numbers::PI);
              lon = (surface_cell_spherical_coordinates.at(i).second <= numbers::PI
                     ?
                     surface_cell_spherical_coordinates.at(i).second*(180./numbers::PI)
                     :
                     surface_cell_spherical_coordinates.at(i).second*(180./numbers::PI) - 360.);

              stored_values_lon_lat.push_back(std::make_pair(std::make_pair(lon,lat),geoid_anomaly.at(i)));
            }
          // Write the solution to the stream output
          for (unsigned int i=0; i<stored_values_lon_lat.size(); ++i)
            {
              output << stored_values_lon_lat.at(i).first.first
                     << ' '
                     << stored_values_lon_lat.at(i).first.second
                     << ' '
                     << stored_values_lon_lat.at(i).second
                     << std::endl;
            }
        }
      else
        {
          std::vector<std::pair<Point<dim>,double> > stored_values_xyz;
          for (unsigned int i=0; i<surface_cell_locations.size(); ++i)
            {
              stored_values_xyz.push_back(std::make_pair(surface_cell_locations.at(i),geoid_anomaly.at(i)));
            }
          // Write the solution to the stream output
          for (unsigned int i=0; i<stored_values_xyz.size(); ++i)
            {
              output << stored_values_xyz.at(i).first
                     << ' '
                     << stored_values_xyz.at(i).second
                     << std::endl;
            }
        }

      const std::string filename = this->get_output_directory() +
                                   "geoid_anomaly." +
                                   dealii::Utilities::int_to_string(this->get_timestep_number(), 5);
      const unsigned int max_data_length = dealii::Utilities::MPI::max (output.str().size()+1,
                                                                        this->get_mpi_communicator());
      const unsigned int mpi_tag = 123;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());
          file << "# "
               << ((output_in_lat_lon == true)? "longitude latitude" : "x y z")
               << " geoid_anomaly" << std::endl;

          // first write out the data we have created locally
          file << output.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // then loop through all of the other processors and collect
          // data, then write it to the file
          for (unsigned int p=1; p<dealii::Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // get the data. note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                        this->get_mpi_communicator(), &status);

              // output the string. note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. write only this part by outputting it as a
              // C string object, rather than as a std::string
              file << tmp.c_str();
            }
        }
      else
        // on other processors, send the data to processor zero. include the \0
        // character at the end of the string
        {
          MPI_Send (&output.str()[0], output.str().size()+1, MPI_CHAR, 0, mpi_tag,
                    this->get_mpi_communicator());
        }

      return std::pair<std::string,std::string>("Writing geoid anomaly:",
                                                filename);
    }

    template <int dim>
    void
    Geoid<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          prm.declare_entry("Maximum degree","20",
                            Patterns::Integer (0),
                            "This parameter can be a random positive integral. However, the value normally should not exceed the maximum"
                            "degree of the initial perturbed temperature field. For example, if the initial condition uses S40RTS, the"
                            "maximun degree should not be larger than 40.");
          prm.declare_entry("Minimum degree","2",
                            Patterns::Integer (0),
                            "This parameter normally is set to 2 since the perturbed gravitational potential at degree 1 always vanishes"
                            "in a reference frame with the planetary center of mass as the center of the reference frame.");
          prm.declare_entry("Output data in geographical coordinates", "false",
                            Patterns::Bool(),
                            "Option to output the geoid anomaly in geographical coordinates (latitude and longitude). "
                            "The default is false, so postprocess will output the data in geocentric coordinates (x,y,z) as normally.");
          prm.declare_entry("Density above","0",
                            Patterns::Double (0),
                            "The density value above the surface boundary.");
          prm.declare_entry("Density below","8000",
                            Patterns::Double (0),
                            "The density value below the CMB boundary.");
          prm.declare_entry("Also output the spherical harmonic coefficients of geoid anomaly", "false",
                            Patterns::Bool(),
                            "Option to also output the spherical harmonic coefficients of the geoid anomaly up to the maximum degree. "
                            "The default is false, so postprocess will only output the grid geoid anomaly. ");
          prm.declare_entry("Also output the spherical harmonic coefficients of surface dynamic topography contribution", "false",
                            Patterns::Bool(),
                            "Option to also output the spherical harmonic coefficients of the surface dynamic topography contribution "
                            "to the maximum degree. The default is false. ");
          prm.declare_entry("Also output the spherical harmonic coefficients of CMB dynamic topography contribution", "false",
                            Patterns::Bool(),
                            "Option to also output the spherical harmonic coefficients of the CMB dynamic topography contribution "
                            "to the maximum degree. The default is false. ");
          prm.declare_entry("Also output the spherical harmonic coefficients of density anomaly contribution", "false",
                            Patterns::Bool(),
                            "Option to also output the spherical harmonic coefficients of the density anomaly contribution to the "
                            "maximum degree. The default is false. ");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    Geoid<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          max_degree = prm.get_integer ("Maximum degree");
          min_degree = prm.get_integer ("Minimum degree");
          output_in_lat_lon = prm.get_bool ("Output data in geographical coordinates");
          density_above = prm.get_double ("Density above");
          density_below = prm.get_double ("Density below");
          also_output_geoid_anomaly_SH_coes = prm.get_bool ("Also output the spherical harmonic coefficients of geoid anomaly");
          also_output_surface_dynamic_topo_contribution_SH_coes = prm.get_bool ("Also output the spherical harmonic coefficients of surface dynamic topography contribution");
          also_output_CMB_dynamic_topo_contribution_SH_coes = prm.get_bool ("Also output the spherical harmonic coefficients of CMB dynamic topography contribution");
          also_output_density_anomaly_contribution_SH_coes = prm.get_bool ("Also output the spherical harmonic coefficients of density anomaly contribution");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Geoid,
                                  "geoid",
                                  "A postprocessor that computes a measure of geoid anomaly based on the density anomaly determined "
                                  "from the temperature field in the mantle, and the dynamic topography at the surface and core mantle "
                                  "boundary(CMB). The geoid is computed from the spherical harmonics expansion, so the geometry of the "
                                  "domain needs to be a spherical shell.")

  }
}
