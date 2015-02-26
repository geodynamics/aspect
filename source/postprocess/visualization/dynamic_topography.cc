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


#include <aspect/postprocess/visualization/dynamic_topography.h>
#include <aspect/simulator_access.h>

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

        // evaluate a single point per cell
        const QMidpoint<dim> quadrature_formula;
        const QMidpoint<dim-1> quadrature_formula_face;

        Assert(quadrature_formula_face.size()==1, ExcInternalError());

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values   |
                                 update_gradients   |
                                 update_quadrature_points );

        FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                          this->get_fe(),
                                          quadrature_formula_face,
                                          update_JxW_values);

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, this->n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, this->n_compositional_fields());

        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

        double integrated_topography = 0;
        double integrated_surface_area = 0;

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
                // see if the cell is at the *top* boundary, not just any boundary
                unsigned int top_face_idx = numbers::invalid_unsigned_int;
                {
                  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                    if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                      {
                        top_face_idx = f;
                        break;
                      }
                }
                if (top_face_idx == numbers::invalid_unsigned_int)
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
                // stress and compute the component in direction of the
                // gravity vector

                double dynamic_topography_x_volume = 0;
                double volume = 0;

                // Compute the integral of the dynamic topography function
                // over the entire cell, by looping over all quadrature points
                // (currently, there is only one, but the code is generic).
                for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                  {
                    const Point<dim> location = fe_values.quadrature_point(q);
                    const double viscosity = out.viscosities[q];
                    const double density   = out.densities[q];

                    const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q] - 1./3 * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>();
                    const SymmetricTensor<2,dim> shear_stress = 2 * viscosity * strain_rate;

                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                    const Tensor<1,dim> gravity_direction = gravity/gravity.norm();

                    // Subtract the dynamic pressure
                    const double dynamic_pressure   = in.pressure[q] - this->get_adiabatic_conditions().pressure(location);
                    const double sigma_rr           = gravity_direction * (shear_stress * gravity_direction) - dynamic_pressure;
                    const double dynamic_topography = - sigma_rr / gravity.norm() / density;

                    // JxW provides the volume quadrature weights. This is a general formulation
                    // necessary for when a quadrature formula is used that has more than one point.
                    dynamic_topography_x_volume += dynamic_topography * fe_values.JxW(q);
                    volume += fe_values.JxW(q);
                  }

                const double dynamic_topography_cell_average = dynamic_topography_x_volume / volume;
                // Compute the associated surface area to later compute the surfaces weighted integral
                fe_face_values.reinit(cell, top_face_idx);
                const double surface = fe_face_values.JxW(0);

                integrated_topography += dynamic_topography_cell_average*surface;
                integrated_surface_area += surface;

                (*return_value.second)(cell_index) = dynamic_topography_cell_average;
              }

        // Calculate surface weighted average dynamic topography
        const double average_topography = Utilities::MPI::sum (integrated_topography,this->get_mpi_communicator())
                                          / Utilities::MPI::sum (integrated_surface_area,this->get_mpi_communicator());

        // subtract the average dynamic topography
        cell_index = 0;
        cell = this->get_dof_handler().begin_active();
        if (subtract_mean_dyn_topography)
          for (; cell!=endc; ++cell,++cell_index)
            if (cell->is_locally_owned()
                && (*return_value.second)(cell_index) != 0
                && cell->at_boundary())
              {
                for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                  if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                    {
                      (*return_value.second)(cell_index) -= average_topography;
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
                                 "in the outputted data file (not visualization). "
                                 "'true' subtracts the mean, 'false' leaves "
                                 "the calculated dynamic topography as is. ");

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
                                                  "for the dynamic topography. The approach to determine the "
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
                                                  "$h=\\frac{\\sigma_{rr}}{\\|\\mathbf g\\| \\rho}$ where $\\rho$ "
                                                  "is the density at the cell center."
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
