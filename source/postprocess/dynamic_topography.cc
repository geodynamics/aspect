/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/postprocess/dynamic_topography.h>

#include <aspect/postprocess/boundary_pressures.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    DynamicTopography<dim>::execute (TableHandler &)
    {
      const Postprocess::BoundaryPressures<dim> &boundary_pressures =
        this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::BoundaryPressures<dim>>();

      // Get the average pressure at the top and bottom boundaries.
      // This will be used to compute the dynamic pressure at the boundaries.
      const double surface_pressure = boundary_pressures.pressure_at_top();
      const double bottom_pressure = boundary_pressures.pressure_at_bottom();

      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const types::boundary_id bottom_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      // If the gravity vector is pointed *up*, as determined by representative points
      // at the surface and at depth, then we are running backwards advection, and need
      // to reverse the dynamic topography values.
      const bool backward_advection = (this->get_gravity_model().gravity_vector(this->get_geometry_model().representative_point(0.)) *
                                       (this->get_geometry_model().representative_point(0.) -
                                        this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth()))) >= 0.;

      const unsigned int quadrature_degree = this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1;
      // Gauss quadrature in the interior for best accuracy.
      const QGauss<dim> quadrature_formula(quadrature_degree);
      // GLL quadrature on the surface to get a diagonal mass matrix.
      const QGaussLobatto<dim-1> quadrature_formula_face(quadrature_degree);

      const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;
      const unsigned int dofs_per_face = this->get_fe().dofs_per_face;
      const unsigned int n_q_points = quadrature_formula.size();
      const unsigned int n_face_q_points = quadrature_formula_face.size();

      // The CBF method involves both boundary and volume integrals on the
      // cells at the boundary. Construct FEValues objects for each of these integrations.
      FEValues<dim> fe_volume_values (this->get_mapping(),
                                      this->get_fe(),
                                      quadrature_formula,
                                      update_values |
                                      update_gradients |
                                      update_quadrature_points |
                                      update_JxW_values);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_JxW_values |
                                        update_values |
                                        update_gradients |
                                        update_quadrature_points);

      // Storage for shape function values for the current solution.
      // Used for constructing the known side of the CBF system.
      std::vector<Tensor<1,dim>> phi_u (dofs_per_cell);
      std::vector<SymmetricTensor<2,dim>> epsilon_phi_u (dofs_per_cell);
      std::vector<double> div_phi_u (dofs_per_cell);
      std::vector<double> div_solution(n_q_points);

      // Vectors for solving CBF system.
      Vector<double> local_vector(dofs_per_cell);
      Vector<double> local_mass_matrix(dofs_per_cell);

      LinearAlgebra::BlockVector rhs_vector(this->introspection().index_sets.system_partitioning, this->get_mpi_communicator());
      // The mass matrix may be stored in a vector as it is a
      // diagonal matrix.
      LinearAlgebra::BlockVector mass_matrix(this->introspection().index_sets.system_partitioning, this->get_mpi_communicator());

      LinearAlgebra::BlockVector distributed_topo_vector(this->introspection().index_sets.system_partitioning, this->get_mpi_communicator());

      topo_vector.reinit(this->introspection().index_sets.system_partitioning,
                         this->introspection().index_sets.system_relevant_partitioning,
                         this->get_mpi_communicator());
      distributed_topo_vector = 0.;
      topo_vector = 0.;

      // Possibly keep track of the dynamic topography values for
      // later surface output.
      std::vector<std::pair<Point<dim>, double>> stored_values_surface;
      std::vector<std::pair<Point<dim>, double>> stored_values_bottom;
      visualization_values.reinit(this->get_triangulation().n_active_cells());
      visualization_values = 0.;

      // Loop over all of the surface cells and if one less than h/3 away from
      // one of the top or bottom boundaries, assemble CBF system for it.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            // see if the cell is at the *top* or *bottom* boundary, not just any boundary
            unsigned int face_idx = numbers::invalid_unsigned_int;
            for (const unsigned int f : cell->face_indices())
              {
                if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                  {
                    // If the cell is at the top boundary, assign face_idx.
                    face_idx = f;
                    break;
                  }
                else if (cell->at_boundary(f) && cell->face(f)->boundary_id() == bottom_boundary_id)
                  {
                    // If the cell is at the bottom boundary, assign face_idx.
                    face_idx = f;
                    break;
                  }
              }
            // If the cell is not at the boundary, jump to the next cell.
            if (face_idx == numbers::invalid_unsigned_int)
              continue;

            fe_volume_values.reinit (cell);
            fe_face_values.reinit (cell, face_idx);

            local_vector = 0.;
            local_mass_matrix = 0.;

            // Evaluate the material model in the cell volume.
            MaterialModel::MaterialModelInputs<dim> in_volume(fe_volume_values, cell, this->introspection(), this->get_solution());
            MaterialModel::MaterialModelOutputs<dim> out_volume(fe_volume_values.n_quadrature_points, this->n_compositional_fields());
            in_volume.requested_properties = MaterialModel::MaterialProperties::density | MaterialModel::MaterialProperties::viscosity;
            this->get_material_model().evaluate(in_volume, out_volume);

            // Get solution values for the divergence of the velocity, which is not
            // computed by the material model.
            fe_volume_values[this->introspection().extractors.velocities].get_function_divergences (this->get_solution(), div_solution);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                const double eta = out_volume.viscosities[q];
                const double density = out_volume.densities[q];
                const bool is_compressible = this->get_material_model().is_compressible();
                const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(in_volume.position[q]);
                const double JxW = fe_volume_values.JxW(q);

                // Set up shape function values
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                  {
                    phi_u[k] = fe_volume_values[this->introspection().extractors.velocities].value(k,q);
                    epsilon_phi_u[k] = fe_volume_values[this->introspection().extractors.velocities].symmetric_gradient(k,q);
                    div_phi_u[k] = fe_volume_values[this->introspection().extractors.velocities].divergence (k, q);
                  }

                for (unsigned int i = 0; i<dofs_per_cell; ++i)
                  {
                    // Viscous stress part
                    local_vector(i) += 2.0 * eta * ( epsilon_phi_u[i] * in_volume.strain_rate[q]
                                                     - (is_compressible ? 1./3. * div_phi_u[i] * div_solution[q] : 0.0) ) * JxW;
                    // Pressure and compressibility parts
                    local_vector(i) -= div_phi_u[i] * in_volume.pressure[q] * JxW;
                    // Force part
                    local_vector(i) -= density * gravity * phi_u[i] * JxW;
                  }
              }
            // Assemble the mass matrix for cell face. Since we are using GLL
            // quadrature, the mass matrix will be diagonal, and we can just assemble it into a vector.
            for (unsigned int q=0; q < n_face_q_points; ++q)
              {
                const double JxW = fe_face_values.JxW(q);
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    const Tensor<1,dim> phi_u = fe_face_values[this->introspection().extractors.velocities].value(i,q);
                    local_mass_matrix(i) += phi_u *
                                            phi_u *
                                            JxW;
                  }
              }

            cell->distribute_local_to_global(local_vector, rhs_vector);
            cell->distribute_local_to_global(local_mass_matrix, mass_matrix);
          }

      rhs_vector.compress(VectorOperation::add);
      mass_matrix.compress(VectorOperation::add);

      // Since the mass matrix is diagonal, we can just solve for the stress vector by dividing.
      const IndexSet local_elements = mass_matrix.locally_owned_elements();
      for (unsigned int k=0; k<local_elements.n_elements(); ++k)
        {
          const unsigned int global_index = local_elements.nth_index_in_set(k);
          if (mass_matrix[global_index] > 1.e-15)
            distributed_topo_vector[global_index] = rhs_vector[global_index]/mass_matrix[global_index];
        }
      distributed_topo_vector.compress(VectorOperation::insert);
      topo_vector = distributed_topo_vector;

      // Now loop over the cells again and solve for the dynamic topography.
      // We solve for it on the support points of the system, since it can be
      // directly put into a system vector of the right size.
      std::vector<Point<dim-1>> face_support_points = this->get_fe().base_element(this->introspection().base_elements.temperature).get_unit_face_support_points();
      Quadrature<dim-1> support_quadrature(face_support_points);
      FEFaceValues<dim> fe_support_values (this->get_mapping(),
                                           this->get_fe(),
                                           support_quadrature,
                                           update_values | update_normal_vectors
                                           | update_gradients | update_quadrature_points);

      std::vector<Tensor<1,dim>> stress_support_values( support_quadrature.size() );
      std::vector<double> topo_values( support_quadrature.size() );
      std::vector<types::global_dof_index> face_dof_indices (dofs_per_face);

      // Also construct data structures for getting the dynamic topography at the cell face
      // midpoints. This is a more practical thing for text output and visualization.
      const QGauss<dim-1> output_quadrature(quadrature_degree);
      FEFaceValues<dim> fe_output_values (this->get_mapping(),
                                          this->get_fe(),
                                          output_quadrature,
                                          update_values | update_normal_vectors | update_gradients |
                                          update_quadrature_points | update_JxW_values);
      std::vector<Tensor<1,dim>> stress_output_values( output_quadrature.size() );


      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            // see if the cell is at the *top* boundary, not just any boundary
            unsigned int face_idx = numbers::invalid_unsigned_int;
            // if the face is at the upper surface 'at_upper_surface' will be true, if
            // it is at the lower surface 'at_upper_surface' will be false. The default
            // is true and will be changed to false if it's at the lower boundary. If the
            // cell is at neither boundary the loop will continue to the next cell.
            bool at_upper_surface = true;
            for (const unsigned int f : cell->face_indices())
              {
                if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                  {
                    // If the cell is at the top boundary, assign face_idx.
                    face_idx = f;
                    at_upper_surface = true;
                    break;
                  }
                else if (cell->at_boundary(f) && cell->face(f)->boundary_id() == bottom_boundary_id)
                  {
                    // If the cell is at the bottom boundary, assign face_idx.
                    face_idx = f;
                    at_upper_surface = false;
                    break;
                  }
              }
            // If the cell is not at the boundary, jump to the next cell.
            if (face_idx == numbers::invalid_unsigned_int)
              continue;

            fe_support_values.reinit (cell, face_idx);

            // Evaluate the material model on the cell face.
            MaterialModel::MaterialModelInputs<dim> in_support(fe_support_values, cell, this->introspection(), this->get_solution());
            MaterialModel::MaterialModelOutputs<dim> out_support(fe_support_values.n_quadrature_points, this->n_compositional_fields());
            in_support.requested_properties = MaterialModel::MaterialProperties::density;
            this->get_material_model().evaluate(in_support, out_support);

            fe_support_values[this->introspection().extractors.velocities].get_function_values(topo_vector, stress_support_values);
            cell->face(face_idx)->get_dof_indices (face_dof_indices);
            for (unsigned int i = 0; i < face_dof_indices.size(); ++i)
              {
                // Given the face dof, we get the component and overall cell dof index.
                const std::pair<unsigned int, unsigned int> component_index = this->get_fe().face_system_to_component_index(i);
                const unsigned int component = component_index.first;
                const unsigned int support_index = component_index.second;
                // Dynamic topography is stored in the temperature component.
                if (component == this->introspection().component_indices.temperature)
                  {
                    // Compute the local gravity at the support point.
                    const Point<dim> point = fe_support_values.quadrature_point(support_index);
                    const double gravity_norm = this->get_gravity_model().gravity_vector(point).norm();
                    const Tensor<1,dim> normal = fe_support_values.normal_vector(support_index);

                    // Compute the dynamic topography, formulae are slightly different
                    // for upper and lower boundaries.
                    double dynamic_topography;
                    if (at_upper_surface)
                      {
                        const double delta_rho = out_support.densities[support_index] - density_above;
                        AssertThrow(std::abs(delta_rho) > std::numeric_limits<double>::min(),
                                    ExcMessage("delta_rho is close or equal to zero at the surface."));
                        dynamic_topography = (-stress_support_values[support_index]*normal - surface_pressure)
                                             / delta_rho / gravity_norm;
                      }
                    else
                      {
                        const double delta_rho = out_support.densities[support_index] - density_below;
                        AssertThrow(std::abs(delta_rho) > std::numeric_limits<double>::min(),
                                    ExcMessage("delta_rho is close or equal to zero at the bottom."));
                        dynamic_topography = (-stress_support_values[support_index]*normal - bottom_pressure)
                                             / delta_rho / gravity_norm;
                      }
                    distributed_topo_vector[face_dof_indices[i]] = dynamic_topography * (backward_advection ? -1. : 1.);
                  }
              }

            // Also evaluate the dynamic topography on the cell faces. This is more convenient
            // for ASCII output, as well as for use with the visualization postprocessor.
            fe_output_values.reinit(cell, face_idx);

            // Evaluate the material model on the cell face.
            MaterialModel::MaterialModelInputs<dim> in_output(fe_output_values, cell, this->introspection(), this->get_solution());
            MaterialModel::MaterialModelOutputs<dim> out_output(fe_output_values.n_quadrature_points, this->n_compositional_fields());
            in_output.requested_properties = MaterialModel::MaterialProperties::density;
            this->get_material_model().evaluate(in_output, out_output);

            fe_output_values[this->introspection().extractors.velocities].get_function_values(topo_vector, stress_output_values);

            // Compute the average dynamic topography at the cell face.
            double face_area = 0.;
            double dynamic_topography = 0.;
            for (unsigned int q=0; q < output_quadrature.size(); ++q)
              {
                const Point<dim> point = fe_output_values.quadrature_point(q);
                const Tensor<1,dim> normal = fe_output_values.normal_vector(q);
                const double gravity_norm = this->get_gravity_model().gravity_vector(point).norm();
                const double JxW = fe_output_values.JxW(q);

                if (at_upper_surface)
                  {
                    const double delta_rho = out_output.densities[q] - density_above;
                    dynamic_topography += (-stress_output_values[q]*normal - surface_pressure)
                                          / delta_rho / gravity_norm * JxW;
                  }
                else
                  {
                    const double delta_rho = out_output.densities[q] - density_below;

                    dynamic_topography += (-stress_output_values[q]*normal - bottom_pressure)
                                          / delta_rho / gravity_norm * JxW;
                  }
                face_area += JxW;
              }
            // Get the average dynamic topography for the cell
            dynamic_topography = dynamic_topography * (backward_advection ? -1. : 1.) / face_area;

            // Maybe keep track of surface output vector.
            const bool respect_manifold = true;
            if (output_surface && at_upper_surface)
              stored_values_surface.push_back(std::make_pair(cell->face(face_idx)->center(respect_manifold), dynamic_topography));
            // Maybe keep track of bottom output vector.
            if (output_bottom && !at_upper_surface)
              stored_values_bottom.push_back(std::make_pair(cell->face(face_idx)->center(respect_manifold), dynamic_topography));

            // Add the value to the vector for the visualization postprocessor.
            visualization_values(cell->active_cell_index()) = dynamic_topography;
          }
      distributed_topo_vector.compress(VectorOperation::insert);
      topo_vector = distributed_topo_vector;

      // Possibly output the result to file.
      if (output_surface)
        output_to_file(top_boundary_id, stored_values_surface);
      if (output_bottom)
        output_to_file(bottom_boundary_id, stored_values_bottom);

      return std::pair<std::string,std::string>("Computing dynamic topography", "");
    }

    /**
     * Return the topography vector as calculated by CBF formulation
     */
    template <int dim>
    const LinearAlgebra::BlockVector &
    DynamicTopography<dim>::
    topography_vector() const
    {
      return topo_vector;
    }

    /**
     * Return the cellwise topography vector as calculated by CBF formulation
     */
    template <int dim>
    const Vector<float> &
    DynamicTopography<dim>::
    cellwise_topography() const
    {
      return visualization_values;
    }

    /**
     * Register the other postprocessor that we need: BoundaryPressures
     */
    template <int dim>
    std::list<std::string>
    DynamicTopography<dim>::required_other_postprocessors() const
    {
      return {"boundary pressures"};
    }


    /**
     * Output the dynamic topography solution to
     * a file.
     */
    template <int dim>
    void
    DynamicTopography<dim>::output_to_file(const types::boundary_id boundary_id,
                                           const std::vector<std::pair<Point<dim>,double>> &position_and_topography)
    {
      // get boundary name and avoid spaces for file output
      std::string boundary_name = this->get_geometry_model().translate_id_to_symbol_name(boundary_id);
      std::replace(boundary_name.begin(), boundary_name.end(), ' ', '_');

      std::ostringstream output;

      // On processor 0, write header lines
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          output << "# "
                 << ((dim==2)? "x y " : "x y z ")
                 << "dynamic_topography_" << boundary_name << std::endl;
        }

      for (unsigned int i = 0; i < position_and_topography.size(); ++i)
        {
          output << std::setprecision(10)
                 << position_and_topography[i].first
                 << ' '
                 << std::setprecision(10)
                 << position_and_topography[i].second
                 << std::endl;
        }

      std::string filename = this->get_output_directory() +
                             "dynamic_topography_" + boundary_name + "." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      Utilities::collect_and_write_file_content(filename, output.str(), this->get_mpi_communicator());
    }


    /**
     * Declare the parameters for the postprocessor.
     */
    template <int dim>
    void
    DynamicTopography<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Dynamic topography");
        {
          prm.declare_entry ("Density above","0.",
                             Patterns::Double (0.),
                             "Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. "
                             "This value depends on the density of material that is moved up or down, i.e. crustal rock, and the "
                             "density of the material that is displaced (generally water or air). While the density of crustal rock "
                             "is part of the material model, this parameter `Density above' allows the user to specify the density "
                             "value of material that is displaced above the solid surface. By default this material is assumed to "
                             "be air, with a density of 0. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Density below","9900.",
                             Patterns::Double (0.),
                             "Dynamic topography is calculated as the excess or lack of mass that is supported by mantle flow. "
                             "This value depends on the density of material that is moved up or down, i.e. mantle above CMB, and the "
                             "density of the material that is displaced (generally outer core material). While the density of mantle rock "
                             "is part of the material model, this parameter `Density below' allows the user to specify the density "
                             "value of material that is displaced below the solid surface. By default this material is assumed to "
                             "be outer core material with a density of 9900. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Output surface", "true",
                             Patterns::Bool(),
                             "Whether to output a file containing the surface dynamic topography.");
          prm.declare_entry ("Output bottom", "true",
                             Patterns::Bool(),
                             "Whether to output a file containing the bottom (i.e., CMB) dynamic topography.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    /**
     * Declare the parameters for the postprocessor.
     */
    template <int dim>
    void
    DynamicTopography<dim>::parse_parameters (ParameterHandler &prm)
    {
      CitationInfo::add("geoid");

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Dynamic topography");
        {
          density_above = prm.get_double ("Density above");
          density_below = prm.get_double ("Density below");
          output_surface = prm.get_bool ("Output surface");
          output_bottom = prm.get_bool ("Output bottom");
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
    ASPECT_REGISTER_POSTPROCESSOR(DynamicTopography,
                                  "dynamic topography",
                                  "A postprocessor that computes a measure of dynamic topography "
                                  "based on the stress at the boundary. The data is written into text "
                                  "files named `dynamic\\_topography_\\X.NNNNN' in the output directory, "
                                  "where X is the name of the boundary and NNNNN is the number of the time step."
                                  "\n\n"
                                  "The exact approach works as follows: At each selected boundary, we compute "
                                  "the traction that acts normal to the boundary faces using the "
                                  "consistent boundary flux method as described in "
                                  "``Gresho, Lee, Sani, Maslanik, Eaton (1987). "
                                  "The consistent Galerkin FEM for computing derived boundary "
                                  "quantities in thermal and or fluids problems. International "
                                  "Journal for Numerical Methods in Fluids, 7(4), 371-394.'' "
                                  "From this traction, the dynamic topography is computed using the formula "
                                  "$h=\\frac{\\sigma_{n}}{g \\rho}$ where $g$ is the norm of the gravity and $\\rho$ "
                                  "is the density. For the bottom surface we chose the convention "
                                  "that positive values are up and negative values are down, analogous to "
                                  "the deformation of the upper surface. Note that this implementation takes "
                                  "the direction of gravity into account, which means that reversing the flow "
                                  "in backward advection calculations will not reverse the instantaneous topography "
                                  "because the reverse flow will be divided by the reverse surface gravity.  "
                                  "\n"
                                  "The file format then consists of lines with Euclidean coordinates "
                                  "followed by the corresponding topography value.")
  }
}
