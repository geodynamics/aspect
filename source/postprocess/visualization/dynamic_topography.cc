/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/dynamic_topography.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      std::pair<std::string, Vector<float> *>
      DynamicTopography<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("dynamic_topography",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree;
        const QGauss<dim> quadrature_formula(quadrature_degree);
        const QGauss<dim-1> quadrature_formula_face(quadrature_degree);

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values   |
                                 update_gradients   |
                                 update_quadrature_points |
                                 update_JxW_values);

        FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                          this->get_fe(),
                                          quadrature_formula_face,
                                          update_JxW_values);

        MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());

        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

        double integrated_uppersurface_topography = 0;
        double integrated_uppersurface_area = 0;
        double integrated_lowersurface_topography = 0;
        double integrated_lowersurface_area = 0;

        // loop over all of the surface cells and if one less than h/3 away from
        // the top surface, evaluate the stress at its center
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            if (cell->at_boundary())
              {
                // see if the cell is at the *top* or *bottom* boundary, not just any boundary
                unsigned int face_idx = numbers::invalid_unsigned_int;

                // if the face is at the upper surface 'at_upper_surface' will be true, if
                // it is at the lower surface 'at_upper_surface' will be false. The default
                // is true and will be changed to false if it's at the lower boundary. If the
                // cell is at neither boundary the loop will continue to the next cell.
                bool at_upper_surface = true;
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                  {
                    const double depth_face_center = this->get_geometry_model().depth (cell->face(f)->center());
                    const double upper_depth_cutoff = cell->face(f)->minimum_vertex_distance()/3.0;
                    const double lower_depth_cutoff = this->get_geometry_model().maximal_depth() - cell->face(f)->minimum_vertex_distance()/3.0;

                    // Check if cell is at upper and lower surface at the same time
                    if (depth_face_center < upper_depth_cutoff && depth_face_center > lower_depth_cutoff)
                      AssertThrow(false, ExcMessage("Your geometry model is so small that the upper and lower boundary of "
                                                    "the domain are bordered by the same cell. "
                                                    "Consider using a higher mesh resolution.") );

                    // Check if the face is at the top boundary
                    if (depth_face_center < upper_depth_cutoff)
                      {
                        face_idx = f;
                        break;
                      }
                    // or at the bottom boundary
                    else if (depth_face_center > lower_depth_cutoff)
                      {
                        face_idx = f;
                        at_upper_surface = false;
                        break;
                      }
                  }

                if (face_idx == numbers::invalid_unsigned_int)
                  {
                    (*return_value.second)(cell_index) = 0;
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

                this->get_material_model().evaluate(in, out);


                // for each of the quadrature points, evaluate the
                // stress and compute the component in direction of the
                // gravity vector
                double dynamic_topography_x_volume = 0;
                double volume = 0;

                // check which way gravity points:
                // gravity_direction_binary is 1 if gravity is pointing in (down) and -1 if it is pointing up (out)
                // This is needed to calculate whether g * n is ||g|| or -||g|| and has been introduced for backward advection
                const Tensor <1,dim> g = this->get_gravity_model().gravity_vector(this->get_geometry_model().representative_point(0));
                const Point<dim> point_surf = this->get_geometry_model().representative_point(0);
                const Point<dim> point_bot = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
                const int gravity_direction_binary =  (g * (point_bot - point_surf) >= 0) ?
                                                      1 :
                                                      -1;

                // Compute the integral of the dynamic topography function
                // over the entire cell, by looping over all quadrature points
                for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                  {
                    const Point<dim> location = fe_values.quadrature_point(q);
                    const double viscosity = out.viscosities[q];
                    const double density   = out.densities[q];

                    const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3.0 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                    const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                    const Tensor<1,dim> gravity_direction = gravity/gravity.norm();

                    // Subtract the dynamic pressure
                    const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(location);
                    const double sigma_rr           = gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure;
                    // Choose the right density contrast depending on whether we're at the top or bottom boundary
                    const double density_outside     = at_upper_surface ?
                                                       density_above :
                                                       density_below;
                    const double dynamic_topography = - sigma_rr / (gravity_direction_binary*gravity.norm()) / (density - density_outside);

                    // JxW provides the volume quadrature weights. This is a general formulation
                    // necessary for when a quadrature formula is used that has more than one point.
                    dynamic_topography_x_volume += dynamic_topography * fe_values.JxW(q);
                    volume += fe_values.JxW(q);
                  }

                const double dynamic_topography_cell_average = dynamic_topography_x_volume / volume;
                // Compute the associated surface area to later compute the surfaces weighted integral
                fe_face_values.reinit(cell, face_idx);
                double surface = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    surface += fe_face_values.JxW(q);
                  }

                if (at_upper_surface == true)
                  {
                    integrated_uppersurface_topography += dynamic_topography_cell_average*surface;
                    integrated_uppersurface_area += surface;
                  }
                else
                  {
                    integrated_lowersurface_topography += dynamic_topography_cell_average*surface;
                    integrated_lowersurface_area += surface;
                  }

                (*return_value.second)(cell_index) = dynamic_topography_cell_average;
              }

        // Calculate surface weighted average dynamic topography
        const double average_upper_topography = Utilities::MPI::sum (integrated_uppersurface_topography,this->get_mpi_communicator())
                                                / Utilities::MPI::sum (integrated_uppersurface_area,this->get_mpi_communicator());
        const double average_lower_topography = Utilities::MPI::sum (integrated_lowersurface_topography,this->get_mpi_communicator())
                                                / Utilities::MPI::sum (integrated_lowersurface_area,this->get_mpi_communicator());

        // subtract the average dynamic topography
        cell_index = 0;
        cell = this->get_dof_handler().begin_active();
        if (subtract_mean_dyn_topography)
          for (; cell!=endc; ++cell,++cell_index)
            if (cell->is_locally_owned()
                && (*return_value.second)(cell_index) != 0
                && cell->at_boundary())
              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                {
                  const double depth_face_center = this->get_geometry_model().depth (cell->face(f)->center());
                  const double upper_depth_cutoff = cell->face(f)->minimum_vertex_distance()/3.0;
                  const double lower_depth_cutoff = this->get_geometry_model().maximal_depth() - cell->face(f)->minimum_vertex_distance()/3.0;
                  if (depth_face_center < upper_depth_cutoff)
                    {
                      (*return_value.second)(cell_index) -= average_upper_topography;
                      break;
                    }
                  else if (depth_face_center > lower_depth_cutoff)
                    {
                      (*return_value.second)(cell_index) -= average_lower_topography;
                      break;
                    }
                }

        return return_value;
      }

      template <int dim>
      void
      DynamicTopography<dim>::
      declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Dynamic Topography");
            {
              prm.declare_entry ("Subtract mean of dynamic topography", "false",
                                 Patterns::Bool (),
                                 "Option to remove the mean dynamic topography "
                                 "in the visualization. The mean dynamic topography is calculated "
                                 "and subtracted separately for the upper and lower boundary. "
                                 "'true' subtracts the mean, 'false' leaves "
                                 "the calculated dynamic topography as is.");
              prm.declare_entry ("Density above","0",
                                 Patterns::Double (0),
                                 "Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. "
                                 "This value depends on the density of material that is moved up or down, i.e. crustal rock, and the "
                                 "density of the material that is displaced (generally water or air). While the density of crustal rock "
                                 "is part of the material model, this parameter `Density above' allows the user to specify the density "
                                 "value of material that is displaced above the solid surface. By default this material is assumed to "
                                 "be air, with a density of 0. "
                                 "Units: $kg/m^3$.");
              prm.declare_entry ("Density below","9900",
                                 Patterns::Double (0),
                                 "Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. "
                                 "This value depends on the density of material that is moved up or down, i.e. mantle above CMB, and the "
                                 "density of the material that is displaced (generally outer core material). While the density of mantle rock "
                                 "is part of the material model, this parameter `Density below' allows the user to specify the density "
                                 "value of material that is displaced below the solid surface. By default this material is assumed to "
                                 "be outer core material with a density of 9900. "
                                 "Units: $kg/m^3$.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      DynamicTopography<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Dynamic Topography");
            {
              subtract_mean_dyn_topography              = prm.get_bool("Subtract mean of dynamic topography");
              density_above = prm.get_double ("Density above");
              density_below = prm.get_double ("Density below");
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DynamicTopography,
                                                  "dynamic topography",
                                                  "A visualization output object that generates output "
                                                  "for the dynamic topography at the top and bottom of the model space. The approach to determine the "
                                                  "dynamic topography requires us to compute the stress tensor and "
                                                  "evaluate the component of it in the direction in which "
                                                  "gravity acts. In other words, we compute "
                                                  "$\\sigma_{rr}={\\hat g}^T(2 \\eta \\varepsilon(\\mathbf u)-\\frac 13 (\\textrm{div}\\;\\mathbf u)I)\\hat g - p_d$ "
                                                  "where $\\hat g = \\mathbf g/\\|\\mathbf g\\|$ is the direction of "
                                                  "the gravity vector $\\mathbf g$ and $p_d=p-p_a$ is the dynamic "
                                                  "pressure computed by subtracting the adiabatic pressure $p_a$ "
                                                  "from the total pressure $p$ computed as part of the Stokes "
                                                  "solve. From this, the dynamic "
                                                  "topography is computed using the formula "
                                                  "$h=\\frac{\\sigma_{rr}}{(\\mathbf g \\cdot \\mathbf n)  \\rho}$ where $\\rho$ "
                                                  "is the density at the cell center. For the bottom surface we chose the convection "
                                                  "that positive values are up (out) and negative values are in (down), analogous to "
                                                  "the deformation of the upper surface. "
                                                  "Note that this implementation takes "
                                                  "the direction of gravity into account, which means that reversing the flow "
                                                  "in backward advection calculations will not reverse the intantaneous topography "
                                                  "because the reverse flow will be divided by the reverse surface gravity."
                                                  "\n\n"
                                                  "Strictly speaking, the dynamic topography is of course a "
                                                  "quantity that is only of interest at the surface. However, "
                                                  "we compute it everywhere to make things fit into the framework "
                                                  "within which we produce data for visualization. You probably "
                                                  "only want to visualize whatever data this postprocessor generates "
                                                  "at the surface of your domain and simply ignore the rest of the "
                                                  "data generated.")
    }
  }
}
