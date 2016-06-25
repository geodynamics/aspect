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


#include <aspect/traction_boundary_conditions/lithostatic_pressure.h>
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <deal.II/base/std_cxx11/array.h>

#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
//#include "two_merged_chunks.h"
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>

namespace aspect
{
  namespace TractionBoundaryConditions
  {

    template <int dim>
    void
    LP<dim>::initialize()
    {
      // Ensure the traction boundary conditions are used
      // There might be other types of traction boundary conditions specified as well,
      // so check only for this type
      const std::map<types::boundary_id,std_cxx11::shared_ptr<TractionBoundaryConditions::Interface<dim> > >
      bvs = this->get_traction_boundary_conditions();
      for (typename std::map<types::boundary_id,std_cxx11::shared_ptr<TractionBoundaryConditions::Interface<dim> > >::const_iterator
           p = bvs.begin();
           p != bvs.end(); ++p)
        {
          if (p->second.get() == this)
            traction_bi.insert(p->first);
        }
      AssertThrow(*(traction_bi.begin()) != numbers::invalid_boundary_id,
                  ExcMessage("Did not find the boundary indicator for the lithostatic pressure plugin."));

      // the below is adapted from adiabatic_conditions/initial_profile.cc
      // but we use the initial temperature and composition and only calculate
      // the pressure profile

      // the spacing of the depth profile
      delta_z = this->get_geometry_model().maximal_depth() / (n_points-1);

      const unsigned int n_compositional_fields = this->n_compositional_fields();

      // the pressure at the surface
      pressure[0]    = this->get_surface_pressure();

      // For spherical(-like) domains, modify the representative point
      Point<dim> spherical_representative_point(representative_point);
      const double degrees_to_radians = dealii::numbers::PI/180.0;

      // check location of representative point and
      // set radius to radius of surface of domain
      // for spherical domains
      if (const GeometryModel::SphericalShell<dim> *gm = dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()))
        {
          spherical_representative_point[1] *= degrees_to_radians;

          AssertThrow(spherical_representative_point[1] >= 0.0 &&
                      spherical_representative_point[1] <= gm->opening_angle(),
                      ExcMessage("Longitude of representative point outside shell domain."));

          // set outer radius
          spherical_representative_point[0] = gm->outer_radius();

          if (dim == 3)
            spherical_representative_point[2] *= degrees_to_radians;
        }
      else if (const GeometryModel::Chunk<dim> *gm = dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()))
        {
          spherical_representative_point[1] *= degrees_to_radians;

          AssertThrow(spherical_representative_point[1] >= gm->west_longitude() &&
                      spherical_representative_point[1] <= gm->east_longitude(),
                      ExcMessage("Longitude of representative point outside chunk domain."));

          if (dim ==3)
            {
              spherical_representative_point[2] *= degrees_to_radians;

              AssertThrow(spherical_representative_point[2] >= gm->south_latitude() &&
                          spherical_representative_point[2] <= gm->north_latitude(),
                          ExcMessage("Latitude of representative point outside chunk domain."));
            }

          // set outer radius
          spherical_representative_point[0] = gm->outer_radius();
        }
//      TODO There will be a LayeredChunk soon :)
//      else if (const GeometryModel::TwoMergedChunks<dim> *gm = dynamic_cast<const GeometryModel::TwoMergedChunks<dim>*> (&this->get_geometry_model()))
//        {
//          spherical_representative_point[1] *= degrees_to_radians;
//
//          AssertThrow(spherical_representative_point[1] >= gm->west_longitude() &&
//                      spherical_representative_point[1] <= gm->east_longitude(),
//                      ExcMessage("Longitude of representative point outside chunk domain."));
//
//          if (dim ==3)
//            {
//              spherical_representative_point[2] *= degrees_to_radians;
//
//              AssertThrow(spherical_representative_point[2] >= gm->south_latitude() &&
//                          spherical_representative_point[2] <= gm->north_latitude(),
//                          ExcMessage("Latitude of representative point outside chunk domain."));
//            }
            // set outer radius
//          spherical_representative_point[0] = gm->outer_radius();
//        }
      else if (const GeometryModel::EllipsoidalChunk<dim> *gm = dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*> (&this->get_geometry_model()))
        {
          // NB: the semi major axis a is only equal to the radius at all surface points if the eccentricity is zero.

          // TODO min/max angles are not available from geometry plugin
//          AssertThrow(//spherical_representative_point[0] >= gm->inner_radius() &&
//                      spherical_representative_point[0] <= gm->get_semi_major_axis_a() &&
//                      spherical_representative_point[1] >= gm->minimum_longitude() &&
//                      spherical_representative_point[1] <= gm->maximum_longitude(),
//                      ExcMessage("Longitude of representative point outside chunk domain."));

          spherical_representative_point[1] *= degrees_to_radians;

          if (dim ==3)
            {

//              AssertThrow(spherical_representative_point[2] >= gm->minimum_latitude() &&
//                          spherical_representative_point[2] <= gm->maximum_latitude(),
//                          ExcMessage("Representative point outside chunk domain."));

              spherical_representative_point[2] *= degrees_to_radians;
            }

          // TODO this is only correct for eccentricity zero
          // or ask for the correct radius from user
          spherical_representative_point[0] = gm->get_semi_major_axis_a();
        }
      else if (const GeometryModel::Sphere<dim> *gm = dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()))
        {
          spherical_representative_point[0] = gm->radius();
          spherical_representative_point[1] *= degrees_to_radians;
          spherical_representative_point[2] *= degrees_to_radians;
        }
      else if (const GeometryModel::Box<dim> *gm = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model()))
        {
          const Point<dim> extents = gm->get_extents();

          for (unsigned int d = 0; d < dim-1; ++d)
            AssertThrow(representative_point[d] <= extents[d], ExcMessage("Coordinate "
                                                                          + dealii::Utilities::int_to_string(d)
                                                                          + " of representative point outside box domain."));

          // set z of surface
          representative_point[dim-1] = extents[dim-1];
        }
      else if (const GeometryModel::TwoMergedBoxes<dim> *gm = dynamic_cast<const GeometryModel::TwoMergedBoxes<dim>*> (&this->get_geometry_model()))
        {
          const Point<dim> extents = gm->get_extents();

          for (unsigned int d = 0; d < dim-1; ++d)
            AssertThrow(representative_point[d] <= extents[d], ExcMessage("Coordinate "
                                                                          + dealii::Utilities::int_to_string(d)
                                                                          +"Representative point outside box domain."));

          // set z of surface
          representative_point[dim-1] = extents[dim-1];
        }
      else
        AssertThrow(false, ExcNotImplemented());

      // set up the input for the density function of the material model
      typename MaterialModel::Interface<dim>::MaterialModelInputs in0(1, n_compositional_fields);
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out0(1, n_compositional_fields);

      // where to calculate the density
      // for spherical domains
      if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0 ||
          dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != 0 ||
          dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != 0 ||
          /*dynamic_cast<const GeometryModel::TwoMergedChunks<dim>*> (&this->get_geometry_model()) != 0*/
          dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*> (&this->get_geometry_model()) != 0 )
        {
          // set spherical surface coordinates
          in0.position[0] = spherical_to_cart(spherical_representative_point);
        }
      // and for cartesian domains
      else
        {
          // set cartesian surface coordinates
          in0.position[0] = representative_point;
        }

      // we need the actual temperature at t0
      in0.temperature[0] = this->get_initial_conditions().initial_temperature(in0.position[0]);

      // the previous pressure
      in0.pressure[0] = pressure[0];

      // we need the actual composition at t0
      for (unsigned int c=0; c<n_compositional_fields; ++c)
        in0.composition[0][c] = this->get_compositional_initial_conditions().initial_composition(in0.position[0], c);

      // we do not need the viscosity
      in0.strain_rate.resize(0);

      // evaluate
      this->get_material_model().evaluate(in0, out0);

      // get the magnitude of gravity. we assume
      // that gravity always points along the depth direction. this
      // may not strictly be true always but is likely a good enough
      // approximation here.
      const double density0 = out0.densities[0];
      const double gravity0 = this->get_gravity_model().gravity_vector(in0.position[0]).norm();

      // now integrate pressure downward using trapezoidal integration
      // p'(z) = rho(p,c,T) * |g| * delta_z
      double sum = delta_z * 0.5 * density0 * gravity0;

      // we assume that the t0 temperature and composition fields at the boundary
      // are representative throughout the model time and model domain
      double z;

      for (unsigned int i=1; i<n_points; ++i)
        {
          AssertThrow (i < pressure.size(), ExcMessage(std::string("The current index ")
                                                       + dealii::Utilities::int_to_string(i)
                                                       + std::string(" is bigger than the size of the pressure vector ")
                                                       + dealii::Utilities::int_to_string(pressure.size())));

          // current depth z
          z = double(i) * delta_z;

          // set up the input for the density function of the material model
          typename MaterialModel::Interface<dim>::MaterialModelInputs in(1, n_compositional_fields);
          typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1, n_compositional_fields);

          // where to calculate the density
          // for spherical domains
          if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0 ||
              dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != 0 ||
              dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != 0 ||
              /*dynamic_cast<const GeometryModel::TwoMergedChunks<dim>*> (&this->get_geometry_model()) != 0*/
              dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*> (&this->get_geometry_model()) != 0)
            {
              // decrease radius with depth increment
              spherical_representative_point[0] -= delta_z;
              in.position[0] = spherical_to_cart(spherical_representative_point);
            }
          // and for cartesian domains
          else
            {
              // decrease z coordinate with depth increment
              representative_point[dim-1] -= delta_z;
              in.position[0] = representative_point;
            }

          // we need the actual temperature at t0
          in.temperature[0] = this->get_initial_conditions().initial_temperature(in.position[0]);
          // the previous pressure
          in.pressure[0] = pressure[i-1];

          // we need the actual composition at t0
          for (unsigned int c=0; c<n_compositional_fields; ++c)
            in.composition[0][c] = this->get_compositional_initial_conditions().initial_composition(in.position[0], c);

          in.strain_rate.resize(0); // we do not need the viscosity
          this->get_material_model().evaluate(in, out);

          // get the magnitude of gravity. we assume
          // that gravity always points along the depth direction. this
          // may not strictly be true always but is likely a good enough
          // approximation here.
          const double density = out.densities[0];
          const double gravity = this->get_gravity_model().gravity_vector(in.position[0]).norm();

          // trapezoid integration
          pressure[i] = sum + delta_z * 0.5 * density * gravity;
          sum += delta_z * density * gravity;

        }

      Assert (*std::min_element (pressure.begin(), pressure.end()) >=
              -std::numeric_limits<double>::epsilon() * pressure.size(),
              ExcInternalError());

    }

    template <int dim>
    Tensor<1,dim>
    LP<dim>::
    traction (const Point<dim> &p,
              const Tensor<1,dim> &normal) const
    {
      // We want to set the normal component to the vertical boundary
      // to the lithostatic pressure, the rest of the traction
      // components are left set to zero
      Tensor<1,dim> traction;

      // assign correct value to traction
      // get the lithostatic pressure from a linear interpolation of
      // the calculated profile
      traction = -get_pressure(p) * normal;
      return traction;
    }

    template <int dim>
    double
    LP<dim>::
    get_pressure (const Point<dim> &p) const
    {
      // depth at which we need the pressure
      const double z = this->get_geometry_model().depth(p);

      // TODO: probably this check is unnecessary as
      // geometry models cut off depth at max_depth
      if (z >= this->get_geometry_model().maximal_depth())
        {
          Assert (z <= this->get_geometry_model().maximal_depth() + delta_z,
                  ExcInternalError());
          // return deepest (last) pressure
          return pressure.back();
        }

      const unsigned int i = static_cast<unsigned int>(z/delta_z);
      Assert ((z/delta_z) >= 0, ExcInternalError());
      Assert (i+1 < pressure.size(), ExcInternalError());

      // now do the linear interpolation
      const double d=1.0+i-z/delta_z;
      Assert ((d>=0) && (d<=1), ExcInternalError());

      return d*pressure[i]+(1-d)*pressure[i+1];
    }

    template <int dim>
    Point<dim>
    LP<dim>::
    spherical_to_cart(const Point<dim >sphere_coord) const
    {
      // Input: radius, longitude, latitude
      Point<dim> cart_coord;

      switch (dim)
        {
          case 2:
          {
            cart_coord[0] = sphere_coord[0] * std::cos(sphere_coord[1]); // X
            cart_coord[1] = sphere_coord[0] * std::sin(sphere_coord[1]); // Y
            break;
          }
          case 3:
          {
            cart_coord[0] = sphere_coord[0] * std::sin(0.5*numbers::PI-sphere_coord[2]) * std::cos(sphere_coord[1]); // X
            cart_coord[1] = sphere_coord[0] * std::sin(0.5*numbers::PI-sphere_coord[2]) * std::sin(sphere_coord[1]); // Y
            cart_coord[2] = sphere_coord[0] * std::cos(0.5*numbers::PI-sphere_coord[2]); // Z
            break;
          }
          default:
            Assert (false, ExcNotImplemented());
            break;
        }

      return cart_coord;
    }

    template <int dim>
    void
    LP<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Lithostatic pressure");
        {
          prm.declare_entry ("Representative point", "",
                             Patterns::List(Patterns::Double()),
                             "The point where the pressure profile will be calculated. "
                             "Cartesian coordinates when geometry is a box, otherwise enter radius, longitude, "
                             "and in 3D latitude. The z or radius coordinate is ignored. "
                             "Units: m or degrees.");
          prm.declare_entry("Number of integration points", "1000",
                            Patterns::Integer(0),
                            "The number of integration points over which we integrate the lithostatic pressure "
                            "downwards. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    LP<dim>::parse_parameters (ParameterHandler &prm)
    {
      unsigned int refinement;
      prm.enter_subsection("Mesh refinement");
      {
        refinement = prm.get_integer("Initial adaptive refinement") + prm.get_integer("Initial global refinement");
      }
      prm.leave_subsection();

      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Lithostatic pressure");
        {
          n_points = prm.get_integer("Number of integration points");
          // The representative point where to calculate the depth profile
          const std::vector<double> rep_point =
            dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Representative point")));
          AssertThrow(rep_point.size() == dim, ExcMessage("Representative point does not have the right dimensions."));
          for (unsigned int d = 0; d<dim; d++)
            representative_point[d] = rep_point[d];
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Check that we have enough integration points for this mesh
      AssertThrow(std::pow(2.0,refinement) <= n_points, ExcMessage("Not enough integration points for this resolution."));

      pressure.resize(n_points,-1);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace TractionBoundaryConditions
  {
    ASPECT_REGISTER_TRACTION_BOUNDARY_CONDITIONS(LP,
                                                 "lithostatic pressure",
                                                 "Implementation of a model in which the boundary "
                                                 "traction is given in terms of a normal traction component "
                                                 "set to the lithostatic pressure "
                                                 "calculated according to the parameters in section "
                                                 "``Boundary traction model|Lithostatic pressure''. "
                                                 "\n\n"
                                                 "The lithostatic pressure is calculated by integrating "
                                                 "the pressure downward based on the initial composition "
                                                 "and temperature. "
                                                 "\n\n"
                                                 "Note that the tangential velocity component(s) should be set "
                                                 "to zero. ")
  }
}
