/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/global.h>
#include <aspect/velocity_boundary_conditions/gplates.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    namespace internal
    {
      GPlatesLookup::GPlatesLookup(const Tensor<1,2> &surface_point_one, const Tensor<1,2> &surface_point_two)
      {
        velocity_vals.reinit(0,0);
        old_velocity_vals.reinit(0,0);

        velocity_values = &velocity_vals;
        old_velocity_values = &old_velocity_vals;

        delta_phi = 0.0;
        delta_theta = 0.0;

        // get the cartesian coordinates of the points the 2D model shall lie in
        // this computation is done also for 3D since it is not expensive and the
        // template dim is currently not used here. Could be changed.
        const Tensor<1,3> point_one = cartesian_surface_coordinates(surface_point_one);
        const Tensor<1,3> point_two = cartesian_surface_coordinates(surface_point_two);

        //Set up the normal vector of an unrotated 2D spherical shell
        // that by default lies in the x-y plane.
        const double normal[3] = {0.0,0.0,1.0};
        const Tensor<1,3> unrotated_normal_vector (normal);

        // Compute the normal vector of the plane that contains
        // the origin and the two user-specified points
        Tensor<1,3> rotated_normal_vector;
        cross_product(rotated_normal_vector,point_one,point_two);
        rotated_normal_vector /= rotated_normal_vector.norm();

        // Calculate the crossing line of the two normals,
        // which will be the rotation axis to transform the
        // into each other
        cross_product(rotation_axis,unrotated_normal_vector,rotated_normal_vector);
        rotation_axis /= rotation_axis.norm();

        // Calculate the rotation angle from the inner product rule
        rotation_angle = std::acos(rotated_normal_vector*unrotated_normal_vector);

      }

      bool
      GPlatesLookup::fexists(const std::string &filename) const
      {
        std::ifstream ifile(filename.c_str());
        return ifile;
      }

      void
      GPlatesLookup::load_file(const std::string &filename)
      {
        using boost::property_tree::ptree;
        ptree pt;

        // Check whether file exists, we do not want to throw
        // an exception in case it does not, because it could be by purpose
        // (i.e. the end of the boundary condition is reached)
        AssertThrow (fexists(filename),
                     ExcMessage (std::string("GPlates file <")
                                 +
                                 filename
                                 +
                                 "> not found!"));

        // populate tree structure pt
        read_xml(filename, pt);

        const int n_points = pt.get_child("gpml:FeatureCollection.gml:featureMember.gpml:VelocityField.gml:domainSet.gml:MultiPoint").size();
        const int n_phi = static_cast<int>(std::sqrt(2.*n_points));
        const int n_theta = n_phi /2;

        if ((delta_theta != 0.0) || (delta_phi != 0.0))
          {
            Assert((delta_theta - numbers::PI / (n_theta-1)) <= 1e-7,  ExcMessage("Resolution changed during boundary conditions. Interpolation is not working correctly."));
            Assert((delta_phi - 2*numbers::PI / n_phi) <= 1e-7,  ExcMessage("Resolution changed during boundary conditions. Interpolation is not working correctly."));
          }

        delta_theta =   numbers::PI / (n_theta-1);
        delta_phi   = 2*numbers::PI / n_phi;

        // swap pointers to old and new values, we overwrite the old ones
        // and the new ones become the old ones
        std::swap( velocity_values, old_velocity_values);

        (*velocity_values).reinit(n_theta,n_phi);

        std::string velos = pt.get<std::string>("gpml:FeatureCollection.gml:featureMember.gpml:VelocityField.gml:rangeSet.gml:DataBlock.gml:tupleList");
        std::stringstream in(velos, std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't find velocities. Is file native gpml format for velocities?")));

        unsigned int i = 0;
        char sep;
        while (!in.eof())
          {
            double vtheta, vphi;
            const double cmyr_si = 3.1557e9;

            in >> vtheta >> sep >> vphi;

            if (in.eof())
              break;

            (*velocity_values)[i%n_theta][i/n_theta][0] = vtheta/cmyr_si;
            (*velocity_values)[i%n_theta][i/n_theta][1] = vphi/cmyr_si;

            i++;
          }
      }

      template <int dim>
      Tensor<1,dim>
      GPlatesLookup::surface_velocity(const Point<dim> &position_, const double time_weight) const
      {

        Tensor<1,dim> tensor_position;
        for (unsigned int i = 0 ; i < dim; i++) tensor_position[i] = position_[i];

        Tensor<1,3> position;
        if (dim == 2)
          {
            position = rotate(convert_tensor<dim,3>(tensor_position),rotation_axis,rotation_angle);
          }
        else
          position = convert_tensor<dim,3>(tensor_position);

        const Tensor<1,3> scoord = spherical_surface_coordinates(position);

        const double idtheta = get_idtheta(scoord[0]);
        const unsigned int idxtheta = static_cast<unsigned int>(idtheta);
        const unsigned int nexttheta = idxtheta + 1;


        const double idphi = get_idphi(scoord[1]);
        const unsigned int idxphi = static_cast<unsigned int>(idphi);
        const unsigned int nextphi = (idxphi+1) % velocity_values->n_cols();


        Assert(idxtheta<velocity_values->n_rows(), ExcMessage("not in range"));
        Assert(idxphi  <velocity_values->n_cols(), ExcMessage("not in range"));

        Assert(nexttheta<velocity_values->n_rows(), ExcMessage("not in range"));
        Assert(nextphi  <velocity_values->n_cols(), ExcMessage("not in range"));

        // compute the coordinates of this point in the
        // reference cell between the data points
        const double xi = idtheta-static_cast<double>(idxtheta);
        const double eta = idphi-static_cast<double>(idxphi);

        Assert ((0 <= xi) && (xi <= 1), ExcInternalError());
        Assert ((0 <= eta) && (eta <= 1), ExcInternalError());
        Assert ((0 <= time_weight) && (time_weight <= 1), ExcInternalError());

        // use xi, eta and time_weight for a trilinear interpolation
        // TODO: interpolation in cartesian probably more accurate


        Tensor<1,2> surf_vel;

        surf_vel[0] = time_weight *
                      ((1-xi)*(1-eta)*(*velocity_values)[idxtheta][idxphi][0] +
                       xi    *(1-eta)*(*velocity_values)[nexttheta][idxphi][0] +
                       (1-xi)*eta    *(*velocity_values)[idxtheta][nextphi][0] +
                       xi    *eta    *(*velocity_values)[nexttheta][nextphi][0]);

        surf_vel[1] = time_weight *
                      ((1-xi)*(1-eta)*(*velocity_values)[idxtheta][idxphi][1] +
                       xi    *(1-eta)*(*velocity_values)[nexttheta][idxphi][1] +
                       (1-xi)*eta    *(*velocity_values)[idxtheta][nextphi][1] +
                       xi    *eta    *(*velocity_values)[nexttheta][nextphi][1]);


        surf_vel[0] += (1-time_weight) *
                       ((1-xi)*(1-eta)*(*old_velocity_values)[idxtheta][idxphi][0] +
                        xi    *(1-eta)*(*old_velocity_values)[nexttheta][idxphi][0] +
                        (1-xi)*eta    *(*old_velocity_values)[idxtheta][nextphi][0] +
                        xi    *eta    *(*old_velocity_values)[nexttheta][nextphi][0]);

        surf_vel[1] += (1-time_weight) *
                       ((1-xi)*(1-eta)*(*old_velocity_values)[idxtheta][idxphi][1] +
                        xi    *(1-eta)*(*old_velocity_values)[nexttheta][idxphi][1] +
                        (1-xi)*eta    *(*old_velocity_values)[idxtheta][nextphi][1] +
                        xi    *eta    *(*old_velocity_values)[nexttheta][nextphi][1]);

        const Tensor<1,3> cart_velo = sphere_to_cart_velocity(surf_vel,scoord);
        Tensor<1,dim> velos;

        if (dim == 2)
          velos = convert_tensor<3,dim>(rotate(cart_velo,rotation_axis,-1.0*rotation_angle));
        else
          velos = convert_tensor<3,dim>(cart_velo);

        return velos;
      }

      Tensor<1,3>
      GPlatesLookup::rotate (const Tensor<1,3> &position, const Tensor<1,3> &rotation_axis, const double angle) const
      {
        Tensor<1,3> cross;
        cross_product(cross,position,rotation_axis);
        Tensor<1,3> newpos = (1-std::cos(angle)) * rotation_axis*(rotation_axis*position) +
                             std::cos(angle) * position + std::sin(angle) * cross;
        return newpos;
      }

      Tensor<1,3>
      GPlatesLookup::spherical_surface_coordinates(const Tensor<1,3> &position) const
      {
        Tensor<1,3> scoord;

        scoord[0] = std::acos(position[2]/std::sqrt(position.norm_square())); // Theta
        scoord[1] = std::atan2(position[1],position[0]); // Phi
        if (scoord[1] < 0.0) scoord[1] = 2*numbers::PI + scoord[1]; // correct phi to [0,2*pi]
        scoord[2] = std::sqrt(position.norm_square()); // R
        return scoord;
      }

      Tensor<1,3>
      GPlatesLookup::cartesian_surface_coordinates(const Tensor<1,2> &sposition) const
      {
        Tensor<1,3> ccoord;

        ccoord[0] = std::sin(sposition[0] * std::cos(sposition[1])); // X
        ccoord[1] = std::sin(sposition[0] * std::sin(sposition[1])); // Y
        ccoord[2] = std::cos(sposition[0]); // Z
        return ccoord;
      }

      Tensor<1,3>
      GPlatesLookup::sphere_to_cart_velocity(const Tensor<1,2> &s_velocities, const Tensor<1,3> &s_position) const
      {
        Tensor<1,3> velocity;

        velocity[0] = -1.0 * s_velocities[1] * std::sin(s_position[1]) + s_velocities[0]*std::cos(s_position[0])*std::cos(s_position[1]);
        velocity[1] = s_velocities[1]*std::cos(s_position[1])+s_velocities[0]*std::cos(s_position[0])*std::sin(s_position[1]);
        velocity[2] = -1.0*s_velocities[0]*std::sin(s_position[0]);
        return velocity;
      }

      template <int in, int out>
      Tensor<1,out>
      GPlatesLookup::convert_tensor (Tensor<1,in> old_tensor) const
      {
        Tensor<1,out> new_tensor;
        for (unsigned int i = 0; i < out; i++)
          if (i < in) new_tensor[i] = old_tensor[i];
          else new_tensor[i] = 0.0;

        return new_tensor;
      }

      double
      GPlatesLookup::get_idphi(const double phi_) const
      {
        const double phi = std::max(std::min(phi_,2*numbers::PI-1e-7),0.0);

        Assert(phi>=0, ExcMessage("not in range"));
        Assert(phi<=2*numbers::PI, ExcMessage("not in range"));
        return phi/delta_phi;
      }

      double
      GPlatesLookup::get_idtheta(const double theta_) const
      {
        double theta = std::max(std::min(theta_,numbers::PI-1e-7),0.0);

        Assert(theta>=0, ExcMessage("not in range"));
        Assert(theta<=numbers::PI, ExcMessage("not in range"));
        return theta/delta_theta;
      }
    }

    template <int dim>
    GPlates<dim>::GPlates ()
      :
      current_time(0.0),
      current_time_step(0),
      time_step(0.0),
      time_weight(0.0),
      time_dependent(true),
      point1("0.0,0.0"),
      point2("0.0,0.0"),
      lookup()
    {}


    template <int dim>
    void
    GPlates<dim>::initialize (const GeometryModel::Interface<dim> &geometry_model)
    {
      int loaded;
      Tensor<1,2> pointone;
      Tensor<1,2> pointtwo;
      char sep;

      std::stringstream streampoint(point1);
      streampoint >> pointone[0] >> sep >> pointone[1];

      std::stringstream streampoint2(point2);
      streampoint2 >> pointtwo[0] >> sep >> pointtwo[1];

      lookup.reset(new internal::GPlatesLookup(pointone,pointtwo));

      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&geometry_model) != 0,
              ExcMessage ("This boundary condition can only be used if the geometry "
                          "is a spherical shell."));
      Assert((dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&geometry_model))->opening_angle()==360, ExcMessage("gplates velocity model just works for opening angle == 360"));

      // TODO: Move the loading of files into set_current_time would remove
      // some double effort if the Start time != 0

      // Load first velocity file
      lookup->load_file(create_filename(current_time_step));

      // Load next velocity file for interpolation
      try
        {
          lookup->load_file(create_filename(current_time_step+1));
          time_dependent = true;
        }
      catch (...)
        // if loading the next time step's file failed, then simply
        // take the last one a second time
        {
          // Next velocity file not found --> Constant velocities
          lookup->load_file(create_filename(current_time_step));
          time_dependent = false;

          // Give warning if first processor
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
            std::cout << std::endl
                      << "Loading second velocity file did not succeed."
                      << "Assuming constant boundary conditions model." << std::endl
                      << std::endl;
        }
    }

    template <int dim>
    std::string
    GPlates<dim>::create_filename (const int timestep) const
    {
      std::string templ = data_directory+velocity_file_name;

      char filename[150];
      snprintf (filename, 150, templ.c_str(), timestep);
      return std::string (filename);
    }

    template <int dim>
    void
    GPlates<dim>::set_current_time (const double time)
    {
      current_time = time-velocity_file_start_time;
      if (time_dependent && (current_time > 0.0))
        {
          const unsigned int old_time_step = current_time_step;

          if (static_cast<unsigned int>(current_time / time_step) > current_time_step)
            {
              current_time_step = static_cast<unsigned int>(current_time / time_step);
              // Load next velocity file for interpolation
              try
                {
                  // If the time step was large enough to move forward more
                  // then one velocity file, we need to load both current files
                  // to stay accurate in interpolation
                  if (current_time_step > old_time_step + 1)
                    lookup->load_file(create_filename(current_time_step));

                  lookup->load_file(create_filename(current_time_step+1));
                }
              catch (...)
                // if loading failed, simply take the previous file a second time
                {
                  // Next velocity file not found --> Constant velocities
                  lookup->load_file(create_filename(old_time_step));
                  // no longer consider the problem time dependent from here on out
                  time_dependent = false;

                  // Give warning if first processor
                  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
                    std::cout << std::endl << "Loading new velocity file did not succeed." << std::endl
                              << "Assuming constant boundary conditions for rest of model." << std::endl
                              << std::endl;
                }
            }
          time_weight = current_time/time_step - current_time_step;
          Assert ((0 <= time_weight) && (time_weight <= 1),
                  ExcMessage("Error in set_current_time. Time_weight has to be in [0,1]"));
        }
    }


    template <int dim>
    Tensor<1,dim>
    GPlates<dim>::
    boundary_velocity (const Point<dim> &position) const
    {
      if (current_time > 0.0)
        return lookup->surface_velocity(position,time_weight);
      else
        return Tensor<1,dim> ();
    }


    template <int dim>
    void
    GPlates<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("GPlates model");
        {
          prm.declare_entry ("Data directory", "data/velocity-boundary-conditions/gplates/",
                             Patterns::DirectoryName (),
                             "The path to the model data.");
          prm.declare_entry ("Velocity file name", "phi_test.%d",
                             Patterns::Anything (),
                             "The file name of the material data. Provide file in format: "
                             "(Velocity file name).%d.gpml where %d is any sprintf integer qualifier, "
                             "specifying the format of the current file number.");
          prm.declare_entry ("Time step", "3.1558e13",
                             Patterns::Anything (),
                             "Time step between following velocity files."
                             " Default is one million years expressed in SI units.");
          prm.declare_entry ("Velocity file start time", "0.0",
                             Patterns::Anything (),
                             "Time at which the velocity file with number 0 shall be loaded."
                             "Previous to this time, a no-slip boundary condition is assumed.");
          prm.declare_entry ("Point one", "1.570796,0.0",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in."
                             " Has to be in the format 'a,b' where a and b are theta (polar angle) "
                             " and phi in radians.");
          prm.declare_entry ("Point two", "1.570796,1.570796",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in."
                             " Has to be in the format 'a,b' where a and b are theta (polar angle) "
                             " and phi in radians.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    GPlates<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("GPlates model");
        {
          data_directory        = prm.get ("Data directory");
          velocity_file_name    = prm.get ("Velocity file name");
          time_step             = prm.get_double ("Time step");
          velocity_file_start_time = prm.get_double ("Velocity file start time");
          point1                = prm.get("Point one");
          point2                = prm.get("Point two");
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
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(GPlates,
                                                 "gplates",
                                                 "Implementation of a model in which the boundary "
                                                 "velocity is derived from GPlates.")
  }
}
