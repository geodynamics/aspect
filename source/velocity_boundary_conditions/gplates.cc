/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/velocity_boundary_conditions/gplates.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    namespace internal
    {
      GPlatesLookup::GPlatesLookup(const Tensor<1,2> &surface_point_one,
                                   const Tensor<1,2> &surface_point_two,
                                   const double interpolation_width_)
        :
        velocity_vals(0,0),
        old_velocity_vals(0,0),
        velocity_positions(0,0),
        velocity_values (&velocity_vals),
        old_velocity_values (&old_velocity_vals),
        delta_phi(0.0),
        delta_theta(0.0),
        interpolation_width(interpolation_width_)
      {
        // get the Cartesian coordinates of the points the 2D model shall lie in
        // this computation is done also for 3D since it is not expensive and the
        // template dim is currently not used here. Could be changed.
        const Tensor<1,3> point_one = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_one));
        const Tensor<1,3> point_two = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_two));

        //Set up the normal vector of an unrotated 2D spherical shell
        // that by default lies in the x-y plane.
        const double normal[3] = {0.0,0.0,1.0};
        const Tensor<1,3> unrotated_normal_vector (normal);

        // Compute the normal vector of the plane that contains
        // the origin and the two user-specified points
        Tensor<1,3> rotated_normal_vector;
        cross_product(rotated_normal_vector,point_one,point_two);
        rotated_normal_vector /= rotated_normal_vector.norm();

        if ((rotated_normal_vector - unrotated_normal_vector).norm() > 1e-3)
          {
            // Calculate the crossing line of the two normals,
            // which will be the rotation axis to transform them
            // into each other
            cross_product(rotation_axis,unrotated_normal_vector,rotated_normal_vector);
            rotation_axis /= rotation_axis.norm();

            // Calculate the rotation angle from the inner product rule
            rotation_angle = std::acos(rotated_normal_vector*unrotated_normal_vector);
          }
        else
          {
            rotation_axis = unrotated_normal_vector;
            rotation_angle = 0.0;
          }
      }

      template <int dim>
      void GPlatesLookup::screen_output(const Tensor<1,2> &surface_point_one,
                                        const Tensor<1,2> &surface_point_two,
                                        const ConditionalOStream &pcout) const
      {
        const Tensor<1,3> point_one = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_one));
        const Tensor<1,3> point_two = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_two));

        std::ostringstream output;

        output << std::setprecision (3) << std::setw(3) << std::fixed << std::endl
               << "   Set up GPlates boundary velocity module."  << std::endl
               << std::endl;
        if (dim == 2)
          {
            output << "   Input point 1 spherical coordinates: " << surface_point_one << std::endl
                   << "   Input point 1 Cartesian output coordinates: " << rotate(point_one,rotation_axis,-rotation_angle)  << std::endl
                   << "   Input point 2 spherical coordinates: " << surface_point_two << std::endl
                   << "   Input point 2 Cartesian output coordinates: " << rotate(point_two,rotation_axis,-rotation_angle)  << std::endl
                   << std::endl <<  std::setprecision(2)
                   << "   Model will be rotated by " << -rotation_angle*180/numbers::PI
                   << " degrees around axis " << rotation_axis << std::endl
                   << std::endl;
          }

        pcout << output.str();
      }

      bool
      GPlatesLookup::fexists(const std::string &filename) const
      {
        std::ifstream ifile(filename.c_str());
        return !(!ifile); // only in c++11 you can convert to bool directly
      }

      void
      GPlatesLookup::load_file(const std::string &filename, const ConditionalOStream &pcout)
      {
        pcout << std::endl << "   Loading GPlates boundary velocity file "
              << filename << "." << std::endl << std::endl;

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

        const unsigned int n_points = pt.get_child("gpml:FeatureCollection.gml:featureMember.gpml:VelocityField.gml:domainSet.gml:MultiPoint").size();

        // These formulas look magic, but they are the proper solution to the equation:
        // n_points = n_theta * n_phi with n_phi = 2 * (n_theta - 1)
        // From the XML information we only know n_points, but need n_theta
        // and n_phi to properly size the arrays and get the grip point positions
        const double dn_theta = 0.5 + std::sqrt(0.25 + n_points/2);
        const unsigned int n_theta = static_cast<unsigned int> (dn_theta);
        const unsigned int n_phi = 2 * (dn_theta - 1);

        AssertThrow(dn_theta - n_theta <= 1e-5,
                    ExcMessage("The velocity file has a grid structure that is not readable. Please refer to the manual for a proper grid structure."));

        if ((delta_theta != 0.0) || (delta_phi != 0.0))
          {
            Assert((delta_theta - numbers::PI / (n_theta-1)) <= 1e-7,  ExcMessage("Resolution changed during boundary conditions. Interpolation is not working correctly."));
            Assert((delta_phi - 2*numbers::PI / n_phi) <= 1e-7,  ExcMessage("Resolution changed during boundary conditions. Interpolation is not working correctly."));
          }

        delta_theta =   numbers::PI / (n_theta-1);
        delta_phi   = 2*numbers::PI / n_phi;

        // swap pointers to old and new values, we overwrite the old ones
        // and the new ones become the old ones
        std::swap (velocity_values, old_velocity_values);

        (*velocity_values).reinit(n_theta,n_phi);
        velocity_positions.reinit(n_theta,n_phi);

        std::string velos = pt.get<std::string>("gpml:FeatureCollection.gml:featureMember.gpml:VelocityField.gml:rangeSet.gml:DataBlock.gml:tupleList");
        std::stringstream in(velos, std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Couldn't find velocities. Is file native gpml format for velocities?")));

        // The lat-lon mesh has changed its starting longitude in gplates1.4
        // correct for this while reading in the velocity data
        unsigned int longitude_correction = 0;
        if (gplates_1_4_or_higher(pt))
          longitude_correction = n_phi/2;

        unsigned int i = 0;
        char sep;
        while (!in.eof())
          {
            Tensor<1,2> spherical_velocities;
            const double cmyr_si = 3.1557e9;

            in >> spherical_velocities[0] >> sep >> spherical_velocities[1];

            if (in.eof())
              break;

            const unsigned int idx_theta = i / n_phi;
            const unsigned int idx_phi = (i + longitude_correction) % n_phi;

            // Currently it would not be necessary to update the grid positions at every timestep
            // since they are not allowed to change. In case we allow this later, do it anyway.
            const Tensor<1,3> spherical_position = get_grid_point_position(idx_theta,idx_phi,false);
            velocity_positions[idx_theta][idx_phi] = cartesian_surface_coordinates(spherical_position);
            (*velocity_values)[idx_theta][idx_phi] = sphere_to_cart_velocity(spherical_velocities,spherical_position)
                                                     / cmyr_si;

            i++;
          }

        AssertThrow(i == n_points,
                    ExcMessage (std::string("Number of read in points does not match number of points in file. File corrupted?")));
      }

      template <int dim>
      Tensor<1,dim>
      GPlatesLookup::surface_velocity(const Point<dim> &position,
                                      const double time_weight) const
      {
        Tensor<1,3> internal_position;
        if (dim == 2)
          internal_position = rotate(convert_tensor<dim,3>(position),rotation_axis,rotation_angle);
        else
          internal_position = convert_tensor<dim,3>(position);

        // Main work, interpolate velocity at this point
        const Tensor<1,3> interpolated_velocity = interpolate(internal_position, time_weight);

        Tensor<1,dim> output_boundary_velocity;

        if (dim == 2)
          // convert_tensor conveniently also handles the projection to the 2D plane by
          // omitting the z-component of velocity (since the 2D model lies in the x-y plane).
          output_boundary_velocity = convert_tensor<3,dim>(rotate(interpolated_velocity,rotation_axis,-1.0*rotation_angle));
        else
          output_boundary_velocity = convert_tensor<3,dim>(interpolated_velocity);

        return output_boundary_velocity;
      }

      Tensor<1,3>
      GPlatesLookup::interpolate (const Tensor<1,3> position, const double time_weight) const
      {
        Tensor<1,3> surf_vel;

        int spatial_index[2] = {0,0};
        calculate_spatial_index(spatial_index,position);

        // always use closest point, ensured by the check_termination = false setting
        double n_interpolation_weight = add_interpolation_point(surf_vel,position,spatial_index,time_weight,false);

        // If interpolation is requested, loop over points in theta, phi direction until we exit a circle of
        // interpolation_width radius around the evaluation point position
        if (interpolation_width > 0)
          {
            unsigned int i = 0;
            unsigned int j = 1;
            do
              {
                do
                  {
                    const int idx[2][2] = { {(int)(spatial_index[0]+i), (int)(spatial_index[1]+j)},
                      {(int)(spatial_index[0]-i), (int)(spatial_index[1]-j)}
                    };

                    const double weight_1 = add_interpolation_point(surf_vel,position,idx[0],time_weight,true);
                    const double weight_2 = add_interpolation_point(surf_vel,position,idx[1],time_weight,true);

                    // termination criterion, our interpolation points are outside of the defined circle
                    if ((std::abs(weight_1 + 1) < 1e-14) || (std::abs(weight_2 + 1) < 1e-14))
                      {
                        if (j == 0)
                          {
                            // If we are outside of the interpolation circle even without extending phi (j), we have
                            // visited all possible interpolation points. Return current surf_vel.

                            // If everything worked as intended the radial component of the velocity should be small.
                            // However we can not do this check for velocity = 0. We assume velocities < 1e-16 m/s
                            // to be essentially zero and pass the check in this case.
                            double residual_normal_velocity(0.0);
                            if (surf_vel.norm() > 1e-16)
                              residual_normal_velocity = (surf_vel * position) / (surf_vel.norm() * position.norm());

                            AssertThrow(residual_normal_velocity < 1e-11,
                                        ExcMessage("Error in velocity boundary module interpolation. "
                                                   "Radial component of velocity should be zero, but is not."));

                            // After checking that the residual normal velocity is small
                            // we might as well remove it to increase compatibility of the
                            // pressure right hand side.
                            const Tensor<1,3> normalized_position = position / position.norm();
                            surf_vel -=  (surf_vel * normalized_position) * normalized_position;

                            return surf_vel / n_interpolation_weight;
                          }
                        else
                          {
                            // We are outside of the interpolation circle, but resetting j to 0 (phi to original index)
                            // and moving forward in theta direction (i) might show us more interpolation points.
                            break;
                          }
                      }

                    n_interpolation_weight += weight_1;
                    n_interpolation_weight += weight_2;

                    j++;
                  }
                while (j < velocity_values->n_cols() / 2);
                j = 0;
                i++;
              }
            while (i < velocity_values->n_rows() / 2);
          }

        // If everything worked as intended the radial component of the velocity should be small.
        // However we can not do this check for velocity close to 0. We assume velocities < 1e-16 m/s
        // to be essentially zero and pass the check in this case.
        double residual_normal_velocity(0.0);
        if (surf_vel.norm() > 1e-16)
          residual_normal_velocity = (surf_vel * position) / (surf_vel.norm() * position.norm());

        AssertThrow( residual_normal_velocity < 1e-11,
                     ExcMessage("Error in velocity boundary module interpolation. "
                                "Radial component of velocity should be zero, but is not."));

        // After checking that the residual normal velocity is small
        // we might as well remove it to increase compatibility of the
        // pressure right hand side.
        const Tensor<1,3> normalized_position = position / position.norm();
        surf_vel -=  (surf_vel * normalized_position) * normalized_position;

        // We have interpolated over the whole dataset of the provided velocities, which is ok.
        // Velocity will be constant for all evaluation points in this case.
        return surf_vel / n_interpolation_weight;
      }

      double
      GPlatesLookup::add_interpolation_point(Tensor<1,3>       &surf_vel,
                                             const Tensor<1,3> &position,
                                             const int          spatial_index[2],
                                             const double       time_weight,
                                             const bool         check_termination) const
      {
        // If the point is extended over the poles, do not use it. It will be found
        // by the check in phi direction.
        if ((spatial_index[0] < 0) || (spatial_index[0] >= static_cast<int>(velocity_values->n_rows())))
          return 0.0;

        int idx[2] = {spatial_index[0],spatial_index[1]};
        reformat_indices(idx);

        const Tensor<1,3> spherical_point = get_grid_point_position(idx[0],idx[1],false);
        const Tensor<1,3> cartesian_point = cartesian_surface_coordinates(spherical_point);

        const double arc_dis = arc_distance(position,position.norm() * cartesian_point);

        // termination criterion, our interpolation points are outside of the defined circle
        if (check_termination && (arc_dis > interpolation_width))
          {
            return -1;
          }

        // sin(theta) accounts for the different area covered by grid points at varying latitudes
        // the minimal value is chosen to weight points at the poles according to their number and area
        const double point_weight = std::max(std::sin(spherical_point[0]),
                                             std::sin(delta_theta/2)/velocity_values->n_cols());

        const Tensor<1,3> normalized_position = position / position.norm();

        const Tensor<1,3> velocity = (*velocity_values)[idx[0]][idx[1]];
        const Tensor<1,3> old_velocity = (*old_velocity_values)[idx[0]][idx[1]];

        surf_vel += point_weight * time_weight * rotate_grid_velocity(cartesian_point,normalized_position,velocity);
        if (time_weight < 1.0 - 1e-7)
          surf_vel += point_weight * (1-time_weight) * rotate_grid_velocity(cartesian_point,normalized_position,old_velocity);

        return point_weight;
      }

      Tensor<1,3>
      GPlatesLookup::rotate_grid_velocity(const Tensor<1,3> &data_position,
                                          const Tensor<1,3> &point_position,
                                          const Tensor<1,3> &data_velocity) const
      {

        if ((point_position-data_position).norm()/point_position.norm() < 1e-7)
          return data_velocity;
        else
          {
            Tensor<1,3> local_rotation_axis;
            cross_product(local_rotation_axis,data_position,point_position);
            local_rotation_axis /= local_rotation_axis.norm();

            // Calculate the rotation angle from the inner product rule
            const double local_rotation_angle = std::acos(data_position*point_position);

            const Tensor<1,3> point_velocity = rotate(data_velocity,local_rotation_axis,-local_rotation_angle);

            return point_velocity;
          }
      }


      Tensor<1,3>
      GPlatesLookup::get_grid_point_position(const unsigned int theta_index, const unsigned int phi_index, const bool cartesian) const
      {
        Tensor<1,3> spherical_position;
        spherical_position[0] = theta_index * delta_theta;
        spherical_position[1] = phi_index * delta_phi;
        spherical_position[2] = 1.0;

        if (cartesian)
          return cartesian_surface_coordinates(spherical_position);
        else
          return spherical_position;
      }

      double
      GPlatesLookup::arc_distance(const Tensor<1,3> position_1, const Tensor<1,3> position_2) const
      {
        const double cartesian_distance = (position_1-position_2).norm();

        Assert((position_1.norm() - position_2.norm())/position_1.norm() < 1e-14,
               ExcMessage("Error in velocity boundary module interpolation. "
                          "Radius of different surface points is not equal."));

        const double average_radius = (position_1.norm() + position_2.norm())/2;
        const double arc_distance =  average_radius * 2.0 * std::asin(cartesian_distance/(2*average_radius));
        return arc_distance;
      }

      void
      GPlatesLookup::reformat_indices (int idx[2]) const
      {
        int reformatted_idxtheta = idx[0];
        int reformatted_idxphi = idx[1] % velocity_values->n_cols();

        if (reformatted_idxphi < 0)
          reformatted_idxphi += velocity_values->n_cols();

        if (idx[0] >  static_cast<int>(velocity_values->n_rows())-1)
          {
            reformatted_idxtheta = 2*(velocity_values->n_rows()-1) - idx[0];
            reformatted_idxphi += velocity_values->n_cols()/2;
            reformatted_idxphi = reformatted_idxphi % velocity_values->n_cols();
          }
        if (idx[0] < 0)
          {
            reformatted_idxtheta = -1 * idx[0];
            reformatted_idxphi += velocity_values->n_cols()/2;
            reformatted_idxphi = reformatted_idxphi % velocity_values->n_cols();
          }

        idx[0] = reformatted_idxtheta;
        idx[1] = reformatted_idxphi;
      }

      Tensor<1,3>
      GPlatesLookup::rotate (const Tensor<1,3> &position, const Tensor<1,3> &rotation_axis, const double angle) const
      {
        Tensor<1,3> cross;
        cross_product(cross,position,rotation_axis);
        const Tensor<1,3> newpos = (1-std::cos(angle)) * rotation_axis*(rotation_axis*position) +
                                   std::cos(angle) * position + std::sin(angle) * cross;
        return newpos;
      }

      Tensor<1,3>
      GPlatesLookup::spherical_surface_coordinates(const Tensor<1,3> &position) const
      {
        Tensor<1,3> scoord;

        scoord[0] = std::acos(position[2]/position.norm()); // Theta
        scoord[1] = std::atan2(position[1],position[0]); // Phi
        if (scoord[1] < 0.0) scoord[1] = 2*numbers::PI + scoord[1]; // correct phi to [0,2*pi]
        scoord[2] = position.norm(); // R
        return scoord;
      }

      Tensor<1,3>
      GPlatesLookup::cartesian_surface_coordinates(const Tensor<1,3> &sposition) const
      {
        Tensor<1,3> ccoord;

        ccoord[0] = std::sin(sposition[0]) * std::cos(sposition[1]); // X
        ccoord[1] = std::sin(sposition[0]) * std::sin(sposition[1]); // Y
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
      GPlatesLookup::convert_tensor (const Tensor<1,in> &old_tensor) const
      {
        Tensor<1,out> new_tensor;
        for (unsigned int i = 0; i < out; i++)
          if (i < in) new_tensor[i] = old_tensor[i];
          else new_tensor[i] = 0.0;

        return new_tensor;
      }

      void
      GPlatesLookup::calculate_spatial_index(int *index,
                                             const Tensor<1,3> &position) const
      {
        const Tensor<1,3> scoord = spherical_surface_coordinates(position);
        index[0] = lround(scoord[0]/delta_theta);
        index[1] = lround(scoord[1]/delta_phi);
        reformat_indices(index);
      }


      bool
      GPlatesLookup::gplates_1_4_or_higher(const boost::property_tree::ptree &pt) const
      {
        const std::string gpml_version = pt.get<std::string>("gpml:FeatureCollection.<xmlattr>.gpml:version");
        std::vector<std::string> string_versions = dealii::Utilities::split_string_list(gpml_version,'.');
        std::vector<int> int_versions = dealii::Utilities::string_to_int(string_versions);

        const int gplates_1_3_version[3] = {1,6,322};

        for (unsigned int i = 0; i < int_versions.size(); i++)
          {
            if (int_versions[i] > gplates_1_3_version[i]) return true;
            if (int_versions[i] < gplates_1_3_version[i]) return false;
          }

        return false;
      }
    }

    template <int dim>
    GPlates<dim>::GPlates ()
      :
      time_relative_to_vel_file_start_time(0.0),
      current_time_step(-2),
      velocity_file_start_time(0.0),
      time_step(0.0),
      time_weight(0.0),
      time_dependent(true),
      point1("0.0,0.0"),
      point2("0.0,0.0"),
      interpolation_width(0.0),
      lookup()
    {}


    template <int dim>
    void
    GPlates<dim>::initialize ()
    {
      char sep;

      std::stringstream streampoint(point1);
      streampoint >> pointone[0] >> sep >> pointone[1];

      std::stringstream streampoint2(point2);
      streampoint2 >> pointtwo[0] >> sep >> pointtwo[1];

      if (dim == 2)
        Assert (pointone != pointtwo,
                ExcMessage ("To define a plane for the 2D model the two assigned points "
                            "may not be equal."));

      lookup.reset(new internal::GPlatesLookup(pointone,pointtwo,interpolation_width));

      const GeometryModel::Interface<dim> &geometry_model =
          this->get_geometry_model();

      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&geometry_model) != 0,
              ExcMessage ("This boundary condition can only be used if the geometry "
                          "is a spherical shell."));
      Assert((dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&geometry_model))->opening_angle()==360, ExcMessage("gplates velocity model just works for opening angle == 360"));
    }



    template <int dim>
    std::string
    GPlates<dim>::create_filename (const int timestep) const
    {
      std::string templ = data_directory+velocity_file_name;
      const int size = templ.length();
      char *filename = (char *) (malloc ((size + 10) * sizeof(char)));
      snprintf (filename, size + 10, templ.c_str (), timestep);
      std::string str_filename (filename);
      free (filename);
      return str_filename;
    }


    template <int dim>
    void
    GPlates<dim>::update ()
    {
      Interface<dim>::update ();

      time_relative_to_vel_file_start_time = this->get_time() - velocity_file_start_time;

      // If the boundary condition is constant, switch off time_dependence end leave function.
      // This also sets time_weight to 1.0
      if ((create_filename (current_time_step) == create_filename (current_time_step+1)) && time_dependent)
        {
          lookup->screen_output<dim> (pointone, pointtwo,this->get_pcout());

          lookup->load_file (create_filename (current_time_step),
                             this->get_pcout());
          end_time_dependence (current_time_step);
          return;
        }

      if (time_dependent && (time_relative_to_vel_file_start_time >= 0.0))
        {
          if (current_time_step < 0)
            lookup->screen_output<dim> (pointone, pointtwo,this->get_pcout());

          if (static_cast<int> (time_relative_to_vel_file_start_time / time_step) > current_time_step)
            {
              update_velocity_data();
            }

          time_weight = time_relative_to_vel_file_start_time / time_step - current_time_step;

          Assert ((0 <= time_weight) && (time_weight <= 1),
                  ExcMessage (
                    "Error in set_current_time. Time_weight has to be in [0,1]"));
        }
    }

    template <int dim>
    void
    GPlates<dim>::update_velocity_data ()
    {

      const int old_time_step = current_time_step;
      current_time_step =
        static_cast<unsigned int> (time_relative_to_vel_file_start_time / time_step);
      // Load next velocity file for interpolation
      // If the time step was large enough to move forward more
      // then one velocity file, we need to load both current files
      // to stay accurate in interpolation
      if (current_time_step > old_time_step + 1)
        try
          {
            lookup->load_file (create_filename (current_time_step),
                               this->get_pcout());
          }
        catch (...)
          // If loading current_time_step failed, end time dependent part with old_time_step.
          {
            try
              {
                end_time_dependence (old_time_step);
              }
            catch (...)
              {
                // If loading the old file fails (e.g. there was no old file), cancel the model run.
                // We might get here, if the time step is so large that step t is before the whole boundary condition
                // while step t+1 is already behind all files in time.
                AssertThrow (false,
                             ExcMessage (
                               "Loading new and old velocity file did not succeed. "
                               "Maybe the time step was so large we jumped over all files "
                               "or the files were removed during the model run. "
                               "Another possible way here is to restart a model with "
                               "previously time-dependent boundary condition after the "
                               "last file was already read. Aspect has no way to find the "
                               "last readable file from the current model time. Please "
                               "prescribe the last velocity file manually in such a case. "
                               "Cancelling calculation."));
              }
          }

      // Now load the next velocity file. This part is the main purpose of this function.
      try
        {
          lookup->load_file (create_filename (current_time_step + 1),
                             this->get_pcout());
        }

      // If loading current_time_step + 1 failed, end time dependent part with current_time_step.
      // We do not need to check for success here, because current_time_step was guaranteed to be
      // at least tried to be loaded before, and if it fails, it should have done before (except from
      // hard drive errors, in which case the exception is the right thing to be thrown).

      catch (...)
        {
          end_time_dependence (current_time_step);
        }
    }

    template <int dim>
    void
    GPlates<dim>::end_time_dependence (const int timestep)
    {
      // Next velocity file not found --> Constant velocities
      // by simply loading the old file twice
      lookup->load_file (create_filename (timestep), this->get_pcout());
      // no longer consider the problem time dependent from here on out
      // this cancels all attempts to read files at the next time steps
      time_dependent = false;
      // this cancels the time interpolation in lookup
      time_weight = 1.0;
      // Give warning if first processor
      this->get_pcout() << std::endl
                        << "   Loading new velocity file did not succeed." << std::endl
                        << "   Assuming constant boundary conditions for rest of model run."
                        << std::endl << std::endl;
    }

    template <int dim>
    Tensor<1,dim>
    GPlates<dim>::
    boundary_velocity (const Point<dim> &position) const
    {
      if (time_relative_to_vel_file_start_time >= 0.0)
        return scale_factor * lookup->surface_velocity(position,time_weight);
      else
        return Tensor<1,dim> ();
    }

    template <int dim>
    void
    GPlates<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Boundary velocity model");
      {
        prm.enter_subsection ("GPlates model");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/velocity-boundary-conditions/gplates/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the 'data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Velocity file name", "phi.%d",
                             Patterns::Anything (),
                             "The file name of the material data. Provide file in format: (Velocity file name).%d.gpml where %d is any sprintf integer qualifier, specifying the format of the current file number.");
          prm.declare_entry ("Time step", "1e6",
                             Patterns::Double (0),
                             "Time step between following velocity files. "
                             "Depending on the setting of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years. "
                             "The default is one million, i.e., either one million seconds or one million years.");
          prm.declare_entry ("Velocity file start time", "0.0",
                             Patterns::Double (0),
                             "Time at which the velocity file with number 0 shall be loaded. Previous to this "
                             "time, a no-slip boundary condition is assumed. "
                             "Depending on the setting of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years.");
          prm.declare_entry ("Scale factor", "1",
                             Patterns::Double (0),
                             "Scalar factor, which is applied to the boundary velocity. "
                             "You might want to use this to scale the velocities to a "
                             "reference model (e.g. with free-slip boundary) or another "
                             "plate reconstruction.");
          prm.declare_entry ("Point one", "1.570796,0.0",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in. Has to be in the format 'a,b' where a and b are theta (polar angle)  and phi in radians.");
          prm.declare_entry ("Point two", "1.570796,1.570796",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in. Has to be in the format 'a,b' where a and b are theta (polar angle)  and phi in radians.");
          prm.declare_entry ("Interpolation width", "0",
                             Patterns::Double (0),
                             "Determines the width of the velocity interpolation zone around the current point. "
                             "Currently equals the arc distance between evaluation point and velocity data point that "
                             "is still included in the interpolation. The weighting of the points currently only accounts "
                             "for the surface area a single data point is covering ('moving window' interpolation without "
                             "distance weighting).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    GPlates<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Query the unit system for time since we may have to convert below.
      // Note that we can't use this->convert_output_to_years() since this
      // requires the SimulatorAccess base object to have been initialized,
      // but this hasn't happened yet when we get into this function.
      const bool
      use_years_instead_of_seconds
        = prm.get_bool ("Use years in output instead of seconds");

      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("GPlates model");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory        = prm.get ("Data directory");
          {
            const std::string      subst_text = "$ASPECT_SOURCE_DIR";
            std::string::size_type position;
            while (position = data_directory.find (subst_text),  position!=std::string::npos)
              data_directory.replace (data_directory.begin()+position,
                                      data_directory.begin()+position+subst_text.size(),
                                      ASPECT_SOURCE_DIR);
          }

          velocity_file_name    = prm.get ("Velocity file name");
          interpolation_width   = prm.get_double ("Interpolation width");
          scale_factor          = prm.get_double ("Scale factor");
          point1                = prm.get ("Point one");
          point2                = prm.get ("Point two");

          time_step             = prm.get_double ("Time step");
          velocity_file_start_time = prm.get_double ("Velocity file start time");
          if (use_years_instead_of_seconds == true)
            {
              time_step                *= year_in_seconds;
              velocity_file_start_time *= year_in_seconds;
            }
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
                                                 "velocity is derived from files that are generated "
                                                 "by the GPlates program.")
  }
}
