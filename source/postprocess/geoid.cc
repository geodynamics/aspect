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
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <aspect/geometry_model/spherical_shell.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>


namespace aspect
{
  namespace Postprocess
  {
    // Compute the spherical harmonic coefficients
    template <int dim>
    std::pair<std::vector<double>,std::vector<double> >
    Geoid<dim>::sph_fun2coes(const std::vector<std::vector<double> > &spherical_function) const
    {
      std::vector<double> cosi(spherical_function.size(),0);
      std::vector<double> sini(spherical_function.size(),0);
      std::vector<double> coecos;
      std::vector<double> coesin;
      double cosii, sinii;
      double prefact;

      for (int ideg =  min_degree; ideg < max_degree+1; ideg++)
        {
          for (int iord = 0; iord < ideg+1; iord++)
            {
              // normalization after Dahlen and Tromp, 1986, Appendix B.6
              if (iord == 0)
                {
                  prefact = 1.;
                }
              else
                {
                  prefact = sqrt(2.);
                }

              for (unsigned int ds_num = 0; ds_num < spherical_function.size(); ds_num++)
                {
                  const double cos_component = boost::math::spherical_harmonic_r(ideg,iord,spherical_function.at(ds_num).at(0),spherical_function.at(ds_num).at(1)); //real / cos part
                  const double sin_component = boost::math::spherical_harmonic_i(ideg,iord,spherical_function.at(ds_num).at(0),spherical_function.at(ds_num).at(1)); //imaginary / sine part
                       
                  cosi.at(ds_num) = (spherical_function.at(ds_num).at(3) * prefact * cos_component);
                  sini.at(ds_num) = (spherical_function.at(ds_num).at(3) * prefact * sin_component);
                }

              cosii = 0;
              sinii = 0;
              for (unsigned int ds_num = 0; ds_num < spherical_function.size(); ds_num++)
                {
                  cosii += cosi.at(ds_num) * spherical_function.at(ds_num).at(2);
                  sinii += sini.at(ds_num) * spherical_function.at(ds_num).at(2);
                }
              const double all_cosii = Utilities::MPI::sum (cosii,this->get_mpi_communicator());
              const double all_sinii = Utilities::MPI::sum (sinii,this->get_mpi_communicator());
              coecos.push_back(all_cosii);
              coesin.push_back(all_sinii);
            }
        }
      std::pair<std::vector<double>,std::vector<double> > coes;
      coes = std::make_pair(coecos,coesin);
      return coes;
    }

    // Compute the surface dynamic topography in spherical harmonic expansion
    template <int dim>
    std::pair<std::vector<double>,std::vector<double> >
    Geoid<dim>::get_surface_dynamic_topography() const
    {
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);
        
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_q_points);
        
      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_q_points |
                                        update_JxW_values);
        
      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
        
      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      double integrated_topography = 0;
      double integrated_surface_area = 0;
        
      std::vector<std::pair<Point<dim>,std::pair<double,double> > > stored_values;
        
      // loop over all of the boundary cells and if one less than h/3 away from
      // the top surface, evaluate the stress at its center
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
        
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          if (cell->at_boundary())
            {
              // see if the cell is at the *top* boundary, not just any boundary
              unsigned int top_face_idx = numbers::invalid_unsigned_int;
              {
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                  if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                    {
                      top_face_idx = f;
                      break;
                    }
                        
                if (top_face_idx == numbers::invalid_unsigned_int)
                  continue;
              }
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
                    
                    
              // for each of the quadrature points, evaluate the
              // stress and compute the component in direction of the
              // gravity vector
                    
              double dynamic_topography_x_volume = 0;
              double volume = 0;
                    
              // Compute the integral of the dynamic topography function
              // over the entire cell, by looping over all quadrature points
              // (currently, there is only one, but the code is generic).
              for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                {
                  Point<dim> location = fe_values.quadrature_point(q);
                  const double viscosity = out.viscosities[q];
                  const double density   = out.densities[q];
                        
                  const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                  const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;
                        
                  const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                  const Tensor<1,dim> gravity_direction = gravity/gravity.norm();
                        
                  // Subtract the dynamic pressure
                  const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(location);
                  const double sigma_rr           = gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure;
                  const double dynamic_topography = - sigma_rr / gravity.norm() / (density - density_above);
                        
                  // JxW provides the volume quadrature weights. This is a general formulation
                  // necessary for when a quadrature formula is used that has more than one point.
                  dynamic_topography_x_volume += dynamic_topography * fe_values.JxW(q);
                  volume += fe_values.JxW(q);
                }
                    
              const double dynamic_topography_cell_average = dynamic_topography_x_volume / volume;
              // Compute the associated surface area to later compute the surfaces weighted integral at the Earth surface
              fe_face_values.reinit(cell, top_face_idx);
              double surface = 0;
              for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                {
                  surface += 1.0 * fe_face_values.JxW(q);
                }
              const Point<dim> midpoint_at_surface = cell->face(top_face_idx)->center();
                    
              integrated_topography += dynamic_topography_cell_average*surface;
              integrated_surface_area += surface;
                    
              stored_values.push_back (std::make_pair(midpoint_at_surface, std::make_pair(surface,dynamic_topography_cell_average)));
            }
      // Calculate surface weighted average dynamic topography
      const double average_topography = Utilities::MPI::sum (integrated_topography,this->get_mpi_communicator()) / Utilities::MPI::sum (integrated_surface_area,this->get_mpi_communicator());
      // Get the value of outer radius
      const double outer_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
      (this->get_geometry_model()).outer_radius();

      // Subtract the average dynamic topography.
      // Transfer the geocentric coordinates to the spherical coordinates.
      // Prepare the value of infinitesimal, i.e., sin(theta)*d_theta*d_phi, for the later spherical integral to compute SH coefficients.
      double theta, phi;
      std::vector<std::vector<double> > spherical_function;
      for (unsigned int i=0; i<stored_values.size(); ++i)
        {
          stored_values.at(i).second.second -= average_topography;
        
          const double r_point = stored_values.at(i).first.norm();
          theta = std::acos(stored_values.at(i).first[2]/r_point);
          const double phi_pre = std::acos(stored_values.at(i).first[0]/std::sqrt(stored_values.at(i).first[0]*stored_values.at(i).first[0]+stored_values.at(i).first[1]*stored_values.at(i).first[1]));
          if (stored_values.at(i).first[1] >= 0)
            {
              phi = phi_pre;
            }
          else
            {
              phi = 2 * numbers::PI - phi_pre;
            }
          const double infinitesimal = stored_values.at(i).second.first/(outer_radius*outer_radius);
          std::vector<double> tmp = {theta,phi,infinitesimal,stored_values.at(i).second.second};
          spherical_function.push_back(tmp);
        }
      return sph_fun2coes(spherical_function);
    }
                            
                            
    // Compute the CMB dynamic topography in spherical harmonic expansion
    template <int dim>
    std::pair<std::vector<double>,std::vector<double> >
    Geoid<dim>::get_CMB_dynamic_topography() const
    {
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);
        
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_q_points);
        
      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_q_points |
                                        update_JxW_values);
        
      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());
        
      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));
        
      double integrated_topography = 0;
      double integrated_CMB_area = 0;
        
      std::vector<std::pair<Point<dim>,std::pair<double,double> > > stored_values;
        
      // loop over all of the boundary cells and if one is
      // at the CMB, evaluate the stress at its center
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
        
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          if (cell->at_boundary())
            {
              // see if the cell is at the CMB boundary, not just any boundary
              unsigned int CMB_face_idx = numbers::invalid_unsigned_int;
              {
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                  if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center())
                      > (this->get_geometry_model().maximal_depth() - cell->face(f)->minimum_vertex_distance()/3.))
                    {
                      CMB_face_idx = f;
                      break;
                    }
                        
                if (CMB_face_idx == numbers::invalid_unsigned_int)
                  continue;
              }
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
                    
                    
              // for each of the quadrature points, evaluate the
              // stress and compute the component in direction of the
              // gravity vector
                    
              double dynamic_topography_x_volume = 0;
              double volume = 0;
                    
              // Compute the integral of the dynamic topography function
              // over the entire cell, by looping over all quadrature points
              // (currently, there is only one, but the code is generic).
              for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                {
                  Point<dim> location = fe_values.quadrature_point(q);
                  const double viscosity = out.viscosities[q];
                  const double density   = out.densities[q];
                        
                  const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                  const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;
                        
                  const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                  const Tensor<1,dim> gravity_direction = gravity/gravity.norm();
                        
                  // Subtract the dynamic pressure
                  const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(location);
                  const double sigma_rr           = gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure;
                  const double dynamic_topography = sigma_rr / gravity.norm() / (density_below - density);
                        
                  // JxW provides the volume quadrature weights. This is a general formulation
                  // necessary for when a quadrature formula is used that has more than one point.
                  dynamic_topography_x_volume += dynamic_topography * fe_values.JxW(q);
                  volume += fe_values.JxW(q);
                }
                    
              const double dynamic_topography_cell_average = dynamic_topography_x_volume / volume;
              // Compute the associated surface area to later compute the surfaces weighted integral at CMB
              fe_face_values.reinit(cell, CMB_face_idx);
              double surface = 0;
              for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                {
                  surface += 1.0 * fe_face_values.JxW(q);
                }
              const Point<dim> midpoint_at_surface = cell->face(CMB_face_idx)->center();
                    
              integrated_topography += dynamic_topography_cell_average*surface;
              integrated_CMB_area += surface;
                    
              stored_values.push_back (std::make_pair(midpoint_at_surface, std::make_pair(surface,dynamic_topography_cell_average)));
            }
      // Calculate CMB weighted average dynamic topography
      const double average_topography = Utilities::MPI::sum (integrated_topography,this->get_mpi_communicator()) / Utilities::MPI::sum (integrated_CMB_area,this->get_mpi_communicator());

      // Get the value of inner radius
      const double inner_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
      (this->get_geometry_model()).inner_radius();
        
      // Subtract the average dynamic topography.
      // Transfer the geocentric coordinates to the spherical coordinates.
      // Prepare the value of infinitesimal, i.e., sin(theta)*d_theta*d_phi, for the later spherical integral to compute SH coefficients.
      double theta, phi;
      std::vector<std::vector<double> > spherical_function;
      for (unsigned int i=0; i<stored_values.size(); ++i)
        {
          stored_values.at(i).second.second -= average_topography;
            
          double r_point = stored_values.at(i).first.norm();
          theta = std::acos(stored_values.at(i).first[2]/r_point);
          double phi_pre = std::acos(stored_values.at(i).first[0]/std::sqrt(stored_values.at(i).first[0]*stored_values.at(i).first[0]+stored_values.at(i).first[1]*stored_values.at(i).first[1]));
          if (stored_values.at(i).first[1] >= 0)
            {
              phi = phi_pre;
            }
          else
            {
              phi = 2 * numbers::PI - phi_pre;
            }
          const double infinitesimal = stored_values.at(i).second.first/(inner_radius*inner_radius);
          std::vector<double> tmp = {theta,phi,infinitesimal,stored_values.at(i).second.second};
          spherical_function.push_back(tmp);
        }
      return sph_fun2coes(spherical_function);
    }

    // Compute the density anomaly on a global sphere within the mantle in spherical harmonic expansion
    template <int dim>
    std::pair<double,std::pair<std::vector<double>,std::vector<double> > >
    Geoid<dim>::get_sphere_density_anomaly (const double &r_lower_sphere,
                                            const double &r_upper_sphere) const
    {
      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);
      const QMidpoint<dim> quadrature_formula_center;
   
      Assert(quadrature_formula_center.size()==1, ExcInternalError());

        
      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_q_points |
                                        update_gradients |
                                        update_JxW_values);
     
      FEValues<dim> fe_center_values (this->get_mapping(),
                                      this->get_fe(),
                                      quadrature_formula_center,
                                      update_values |
                                      update_q_points |
                                      update_JxW_values |
                                      update_gradients);
        
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

      std::vector<std::pair<Point<dim>,std::pair<double,double> > > stored_values;
        
      double integrated_density = 0;
      double integrated_volume = 0;

      // loop over all of the cells and if the cell's upper face radius is
      // the input sphere radius, evaluate the density anomaly at its upper face
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
        
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_center_values.reinit(cell);
            // Get the geocentric coordinates of the cell midpoint and later transfer them to the spherical coordinates
            const Point<dim> midpoint_of_cell = fe_center_values.get_quadrature_points().at(0);
            const double r_cell_midpoint = midpoint_of_cell.norm();
            bool cell_at_this_layer = false;
            if (r_cell_midpoint > r_lower_sphere && r_cell_midpoint < r_upper_sphere)
              cell_at_this_layer = true;
            if (cell_at_this_layer == false)
              continue;

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

            // for each of the quadrature points, evaluate the density
            // value and then compute the density anomaly
            double density_x_volume = 0;
            double volume = 0;
             
            // Compute the integral of the density anomaly function
            // over the cell's entire upper face, by looping over all quadrature points;
            // Compute the associated surface area
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              {
                double density = out.densities[q];

                // JxW provides the volume quadrature weights. This is a general formulation
                // necessary for when a quadrature formula is used that has more than one point.
                density_x_volume += density * fe_values.JxW(q);
                volume += 1.0 * fe_values.JxW(q);
              }
            // Compute the cell averaged density anomaly
            const double density_cell_average = density_x_volume / volume;
            
            // Contribute the density of the cell to the integrated density of the current processor. So does the volume. 
            integrated_density += density_x_volume;
            integrated_volume += volume;

            // Find the upper and lower face of the cell, set the associated face area to the average of the two faces.
            double min_face_depth = this->get_geometry_model().maximal_depth();
            unsigned int upper_face_idx = numbers::invalid_unsigned_int;
            double max_face_depth = 0.0;
            unsigned int lower_face_idx = numbers::invalid_unsigned_int;
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                if (this->get_geometry_model().depth (cell->face(f)->center()) < min_face_depth)
                  {
                    upper_face_idx = f;
                    min_face_depth = this->get_geometry_model().depth (cell->face(f)->center());
                  }
                if (this->get_geometry_model().depth (cell->face(f)->center()) > max_face_depth)
                  { 
                    lower_face_idx = f;
                    max_face_depth = this->get_geometry_model().depth (cell->face(f)->center());
                  }
              }
            // Compute the upper face area 
            fe_face_values.reinit(cell,upper_face_idx);
            double upper_face_area = 0;
            for (unsigned int q=0; q<quadrature_formula_face.size(); ++q)
              {
                upper_face_area += 1.0 * fe_face_values.JxW(q);
              }
            // Compute the lower face area
            fe_face_values.reinit(cell,lower_face_idx);
            double lower_face_area = 0;             
            for (unsigned int q=0; q<quadrature_formula_face.size(); ++q)
              {
                lower_face_area += 1.0 * fe_face_values.JxW(q);
              }
            // Compute the associated infinitesimal, i.e., sin(theta)*d_theta*d_phi for the later spherical integral to compute SH coefficients.
            const double infinitesimal = 0.5*(upper_face_area/(r_upper_sphere*r_upper_sphere)+lower_face_area/(r_lower_sphere*r_lower_sphere));
            
            stored_values.push_back (std::make_pair(midpoint_of_cell, std::make_pair(infinitesimal,density_cell_average)));
          }
      
      // Compute the value of the average density for this layer.
      const double layer_average_density = Utilities::MPI::sum (integrated_density,this->get_mpi_communicator()) / Utilities::MPI::sum (integrated_volume,this->get_mpi_communicator()); 

      // Subtract the weighted average density of the current layer from the density of each cell to get the density anomaly.
      // Transfer the geocentric coordinates to the spherical coordinates.
      double theta, phi;
      std::vector<std::vector<double> > spherical_function;
      for (unsigned int i=0; i<stored_values.size(); ++i)
        { 
          stored_values.at(i).second.second -= layer_average_density;
          double r_point = stored_values.at(i).first.norm();
          theta = std::acos(stored_values.at(i).first[2]/r_point);
          double phi_pre = std::acos(stored_values.at(i).first[0]/std::sqrt(stored_values.at(i).first[0]*stored_values.at(i).first[0]+stored_values.at(i).first[1]*stored_values.at(i).first[1]));
          if (stored_values.at(i).first[1] >= 0)
            {
              phi = phi_pre;
            }
          else
            {
              phi = 2 * numbers::PI - phi_pre;
            }
          std::vector<double> tmp = {theta,phi,stored_values.at(i).second.first,stored_values.at(i).second.second};
          spherical_function.push_back(tmp);
        }
      return std::make_pair(layer_average_density,sph_fun2coes(spherical_function));
    }

    // Compute the geoid anomaly
    // The number of extracted spherical layers for the integral is determined by the resolution of cells throughout the mantle
    template <int dim>
    std::pair<std::string,std::string>
    Geoid<dim>::execute (TableHandler &)
    {
      // Get the value of the outer radius and inner radius
      const double outer_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
      (this->get_geometry_model()).outer_radius();
      const double inner_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
      (this->get_geometry_model()).inner_radius();
        
      // Get the value of the surface gravity acceleration from the gravity model
      const Point<dim> surface_point = {outer_radius,0,0};
      const double surface_gravity = this->get_gravity_model().gravity_vector(surface_point).norm();

      // Get the value of the universal gravitational constant
      const double G = 6.674e-11;
      
      // Define a vector to store the radius of all the interior spheres
      std::vector<double> R_sphere;
      const double middle_radius = inner_radius*(1+(outer_radius-inner_radius)/(outer_radius+inner_radius));
      R_sphere.push_back(inner_radius);
      R_sphere.push_back(middle_radius);
      R_sphere.push_back(outer_radius);
      const double thickness_lower = middle_radius-inner_radius;
      const double thickness_upper = outer_radius-middle_radius;
      int base = 2;
      for (int id_radius = 0;id_radius < std::pow(base,depth_integral_refinement)-1;id_radius++)
        {
          const double lower_layer_inner_radius = inner_radius + (id_radius+1)*thickness_lower/(std::pow(base,depth_integral_refinement));
          const double upper_layer_inner_radius = middle_radius + (id_radius+1)*thickness_upper/(std::pow(base,depth_integral_refinement));
          R_sphere.push_back(lower_layer_inner_radius);
          R_sphere.push_back(upper_layer_inner_radius);
        }
        
      // reorder R_sphere to make sure the sphere radii are given from the CMB to surface in order
      for (unsigned int i =0; i < R_sphere.size()-1; i++)
        {
          for (unsigned int j = i+1; j < R_sphere.size(); j++)
            {
              if (R_sphere.at(i) > R_sphere.at(j))
               {
                 const double tmp = R_sphere.at(i);
                 R_sphere.at(i) = R_sphere.at(j);
                 R_sphere.at(j) = tmp;
               }
            }
        }

      // Get the spherical harmonic coefficients of the density anomaly for each interior spherical layer.
      std::vector<std::pair<double,std::pair<std::vector<double>,std::vector<double> > > > coes_density_anomaly;
      for (unsigned int i = 0; i < R_sphere.size()-1; i++)
        {
          coes_density_anomaly.push_back(get_sphere_density_anomaly(R_sphere.at(i),R_sphere.at(i+1)));
        }      

      // Get the density contrast at the surface and CMB
      const double surface_delta_rho = coes_density_anomaly.back().first - density_above;
      const double CMB_delta_rho = density_below - coes_density_anomaly.front().first;
        
      // Compute the integral of the density anomaly's spherical harmonic coefficients throughout the mantle
      // The mantle is divided into a series of spherical layers determined by the cell distributions.
      std::pair<std::vector<double>,std::vector<double> > integrated_coes_density_anomaly;
      int ind = 0; // coefficients index
      for (int ideg =  min_degree; ideg < max_degree+1; ideg++)
        {
          for (int iord = 0; iord < ideg+1; iord++)
            {
              double cos_integral = 0;
              double sin_integral = 0;
              for (unsigned int i = 0; i < R_sphere.size()-1; i++)
                {
                  const double r_layer = 0.5*(R_sphere.at(i)+R_sphere.at(i+1));
                  double layer_thickness = R_sphere.at(i+1) - R_sphere.at(i);
                  cos_integral += coes_density_anomaly.at(i).second.first.at(ind)*r_layer*std::pow(r_layer/outer_radius,ideg+1)*layer_thickness;
                  sin_integral += coes_density_anomaly.at(i).second.second.at(ind)*r_layer*std::pow(r_layer/outer_radius,ideg+1)*layer_thickness;
                }
              integrated_coes_density_anomaly.first.push_back(cos_integral);
              integrated_coes_density_anomaly.second.push_back(sin_integral);
              ind += 1;
            }
        }

      // Get the spherical harmonic coefficients of the surface and CMB dynamic topography
      std::pair<std::vector<double>,std::vector<double> > coes_surface_dynamic_topo = get_surface_dynamic_topography();
      std::pair<std::vector<double>,std::vector<double> > coes_CMB_dynamic_topo = get_CMB_dynamic_topography();
        
      // Compute the spherical harmonic coefficients of geoid anomaly
      std::vector<double> geoid_coecos;  // a vector to store cos terms
      std::vector<double> geoid_coesin;  // a vector to store sin terms
      ind = 0; // reinitialize the ind to 0 for the coefficients index
      for (int ideg =  min_degree; ideg < max_degree+1; ideg++)
        {
          for (int iord = 0; iord < ideg+1; iord++)
            {
              double coecos = integrated_coes_density_anomaly.first.at(ind) + surface_delta_rho*coes_surface_dynamic_topo.first.at(ind)*outer_radius + CMB_delta_rho*coes_CMB_dynamic_topo.first.at(ind)*inner_radius*std::pow(inner_radius/outer_radius,ideg+1);
              double coesin = integrated_coes_density_anomaly.second.at(ind) + surface_delta_rho*coes_surface_dynamic_topo.second.at(ind)*outer_radius + CMB_delta_rho*coes_CMB_dynamic_topo.second.at(ind)*inner_radius*std::pow(inner_radius/outer_radius,ideg+1);
              coecos *= 4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1));
              coesin *= 4 * numbers::PI * G / (surface_gravity * (2 * ideg + 1));
              geoid_coecos.push_back(coecos);
              geoid_coesin.push_back(coesin);
              ind += 1;
            }
        }

      const QMidpoint<dim-1> quadrature_formula_face_center;
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
      double theta, phi;
      for (unsigned int i=0; i<surface_cell_locations.size(); ++i)
        {
          double r_cell = surface_cell_locations.at(i).norm();
          theta = std::acos(surface_cell_locations.at(i)[2]/r_cell);
          double phi_pre = std::acos(surface_cell_locations.at(i)[0]/std::sqrt(surface_cell_locations.at(i)[0]*surface_cell_locations.at(i)[0]+surface_cell_locations.at(i)[1]*surface_cell_locations.at(i)[1]));
          if (surface_cell_locations.at(i)[1] >= 0)
            {
              phi = phi_pre;
            }
          else
            {
              phi = 2 * numbers::PI - phi_pre;
            }
          surface_cell_spherical_coordinates.push_back(std::make_pair(theta,phi));
        }
                            
      // Compute the geoid anomaly based on spherical harmonics
      std::vector<double> geoid_anomaly;
      double prefact;
      for (unsigned int i=0; i<surface_cell_spherical_coordinates.size(); ++i)
        {
          int ind = 0;
          double geoid_value = 0;
          for (int ideg =  min_degree; ideg < max_degree+1; ideg++)
            {
              for (int iord = 0; iord < ideg+1; iord++)
                {
                  const double cos_component = boost::math::spherical_harmonic_r(ideg,iord,surface_cell_spherical_coordinates.at(i).first,surface_cell_spherical_coordinates.at(i).second); //real / cos part
                  const double sin_component = boost::math::spherical_harmonic_i(ideg,iord,surface_cell_spherical_coordinates.at(i).first,surface_cell_spherical_coordinates.at(i).second); //imaginary / sine part
                  // normalization after Dahlen and Tromp, 1986, Appendix B.6
                  if (iord == 0)
                    {
                      prefact = 1.;
                    }
                  else
                    {
                      prefact = sqrt(2.);
                    }
                  geoid_value += prefact*(geoid_coecos.at(ind)*cos_component+geoid_coesin.at(ind)*sin_component);
                  ind += 1;
                }
            }
          geoid_anomaly.push_back(geoid_value);
        }

      // have a stream into which we write the data. the text stream is then
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
                                   Utilities::int_to_string(this->get_timestep_number(), 5);
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
               << " geoid anomaly" << std::endl;
          
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
          prm.declare_entry("maximum degree","40",
                            Patterns::Integer (0),
                            "This parameter can be a random positive integral. However, the value normally should not exceed the maximum"
                            "degree of the initial perturbed temperature field. For example, if the initial condition uses S40RTS, the"
                            "maximun degree should not be larger than 40.");
          prm.declare_entry("minimum degree","2",
                            Patterns::Integer (0),
                            "This parameter normally is set to 2 since the perturbed gravitational potential at degree 1 always vanishes"
                            "in a reference frame with the planetary center of mass as the center of the reference frame.");
          prm.declare_entry("Depth integral refinement","0",
                            Patterns::Integer (0),
                            "This parameter controls how many depth layers the code use for the numerical density integral throughout."
                            "the mantle. In fact, it means the global refinement for the desnity integral, i.e., the default 0 will"
                            "divide the whole spherical domain into two layes and 1 for 4 layers.... Each layer's boundary will be consistent"
                            "with the cells' boundary, so this value should not be larger than the initial global refinement performed on the"
                            "mesh (The most efficient way is of course to set this parameter equal to the initial global refinement). This "
                            "approach will align each layer with the cell's boundary and lead to better accuracy of the numerical integral"
                            "with higher inital global refinement. However, this method will break down when the adaptive mesh refinement "
                            "is used because there might be more refined cells only in a regional part of the domain. So for now we should"
                            "always set the initial adaptive refinement to zero when computing geoid postprocess.");
          prm.declare_entry("output data in geographical coordinates", "false",
                            Patterns::Bool(),
                            "Option to output the geoid anomaly in geographical coordinates (latitude and longitude). "
                            "The default is false, so postprocess will output the data in geocentric coordinates (x,y,z) as normally.");
          prm.declare_entry("Density above","0",
                            Patterns::Double (0),
                            "The density value out of the surface boundary.");
          prm.declare_entry("Density below","8000",
                            Patterns::Double (0),
                            "The density value out of the CMB boundary.");
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
          max_degree = prm.get_integer ("maximum degree");
          min_degree = prm.get_integer ("minimum degree");
          depth_integral_refinement = prm.get_integer("Depth integral refinement");
          output_in_lat_lon = prm.get_bool ("output data in geographical coordinates");
          density_above = prm.get_double ("Density above");
          density_below = prm.get_double ("Density below");
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
                                  "boundary(CMB). The geoid is computed form the spherical harmonics expansion, so the geometry of the "
                                  "domain needs to be a spherical shell.")
  }
}
