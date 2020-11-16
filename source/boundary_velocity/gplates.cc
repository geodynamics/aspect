/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/boundary_velocity/gplates.h>
#include <aspect/utilities.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>


namespace aspect
{
  namespace BoundaryVelocity
  {
    namespace internal
    {
      template <int dim>
      GPlatesLookup<dim>::GPlatesLookup(const Tensor<1,2> &surface_point_one,
                                        const Tensor<1,2> &surface_point_two)
      {
        // get the Cartesian coordinates of the points the 2D model will lie in
        // this computation is done also for 3D since it is not expensive and the
        // template dim is currently not used here. Could be changed.
        const Tensor<1,3> point_one = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_one));
        const Tensor<1,3> point_two = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_two));

        // Set up the normal vector of an unrotated 2D spherical shell
        // that by default lies in the x-y plane.
        const double normal[3] = {0.0,0.0,1.0};
        const Tensor<1,3> unrotated_normal_vector (normal);

        // Compute the normal vector of the plane that contains
        // the origin and the two user-specified points
        Tensor<1,3> rotated_normal_vector = cross_product_3d(point_one,point_two);

        rotated_normal_vector /= rotated_normal_vector.norm();

        if ((rotated_normal_vector - unrotated_normal_vector).norm() > 1e-3)
          {
            // Calculate the crossing line of the two normals,
            // which will be the rotation axis to transform the one
            // normal into the other
            Tensor<1,3> rotation_axis = cross_product_3d(unrotated_normal_vector,rotated_normal_vector);
            rotation_axis /= rotation_axis.norm();

            // Calculate the rotation angle from the inner product rule
            const double rotation_angle = std::acos(rotated_normal_vector*unrotated_normal_vector);

            rotation_matrix = rotation_matrix_from_axis(rotation_axis,rotation_angle);

            // Now apply the rotation that will project point_one onto the known point
            // (0,1,0).
            const Tensor<1,3> rotated_point_one = transpose(rotation_matrix) * point_one;

            const double point_one_coords[3] = {0.0,1.0,0.0};
            const Tensor<1,3> final_point_one (point_one_coords);

            const double second_rotation_angle = std::acos(rotated_point_one*final_point_one);
            Tensor<1,3> second_rotation_axis = cross_product_3d(final_point_one,rotated_point_one);
            second_rotation_axis /= second_rotation_axis.norm();

            const Tensor<2,3> second_rotation_matrix = rotation_matrix_from_axis(second_rotation_axis,second_rotation_angle);

            // The final rotation used for the model will be the combined
            // rotation of the two operation above. This is achieved by a
            // matrix multiplication of the rotation matrices.
            // This concatenation of rotations is the reason for using a
            // rotation matrix instead of a combined rotation_axis + angle
            rotation_matrix = rotation_matrix * second_rotation_matrix;
          }
        else
          {
            rotation_matrix[0][0] = 1.0;
            rotation_matrix[1][1] = 1.0;
            rotation_matrix[2][2] = 1.0;
          }
      }



      template <int dim>
      std::string
      GPlatesLookup<dim>::screen_output(const Tensor<1,2> &surface_point_one,
                                        const Tensor<1,2> &surface_point_two) const
      {
        const Tensor<1,3> point_one = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_one));
        const Tensor<1,3> point_two = cartesian_surface_coordinates(convert_tensor<2,3>(surface_point_two));

        std::ostringstream output;

        output << std::setprecision (3) << std::setw(3) << std::fixed << std::endl
               << "   Setting up GPlates boundary velocity plugin."  << std::endl
               << std::endl;
        if (dim == 2)
          {
            Tensor<1,3> rotation_axis;
            const double rotation_angle = rotation_axis_from_matrix(rotation_axis,rotation_matrix);

            std::array<double,3> angles = angles_from_matrix(rotation_matrix);
            std::array<double,3> back_angles = angles_from_matrix(transpose(rotation_matrix));

            output << "   Input point 1 spherical coordinates: " << surface_point_one  << std::endl
                   << "   Input point 1 normalized cartesian coordinates: " << point_one  << std::endl
                   << "   Input point 1 rotated model coordinates: " << transpose(rotation_matrix) * point_one  << std::endl
                   << "   Input point 2 spherical coordinates: " << surface_point_two  << std::endl
                   << "   Input point 2 normalized cartesian coordinates: " << point_two  << std::endl
                   << "   Input point 2 rotated model coordinates: " << transpose(rotation_matrix) * point_two << std::endl
                   << std::endl <<  std::setprecision(2)
                   << "   Model will be rotated by " << -rotation_angle * 180.0 / numbers::PI
                   << " degrees around axis " << rotation_axis << std::endl
                   << "   The ParaView rotation angles are: " << angles[0] << " " << angles [1] << " " << angles[2] << std::endl
                   << "   The inverse ParaView rotation angles are: " << back_angles[0] << " " << back_angles [1] << " " << back_angles[2]

                   << std::endl;
          }

        return output.str();
      }



      template <int dim>
      void
      GPlatesLookup<dim>::load_file(const std::string &filename,
                                    const MPI_Comm &comm)
      {
        // Read data from disk and distribute among processes
        std::istringstream filecontent(
          Utilities::read_and_distribute_file_content(filename, comm));

        boost::property_tree::ptree pt;

        // populate tree structure pt
        read_xml(filecontent, pt);

        const unsigned int n_points = pt.get_child("gpml:FeatureCollection.gml:featureMember.gpml:VelocityField.gml:domainSet.gml:MultiPoint").size();

        // These formulas look magic, but they are the proper solution to the equation:
        // n_points = n_theta * n_phi with n_phi = 2 * (n_theta - 1)
        // From the XML information we only know n_points, but need n_theta
        // and n_phi to properly size the arrays and get the grip point positions
        const double dn_theta = 0.5 + std::sqrt(0.25 + n_points/2);
        const unsigned int n_theta = static_cast<unsigned int> (dn_theta);
        const unsigned int n_phi = static_cast<unsigned int> (2 * (dn_theta - 1));

        AssertThrow(dn_theta - n_theta <= 1e-5,
                    ExcMessage("The velocity file has a grid structure that is not readable. Please refer to the manual for a proper grid structure."));

        delta_theta =   numbers::PI / (n_theta-1);
        delta_phi   = 2*numbers::PI / n_phi;

        /**
         * Two tables (one for the theta and one for the phi component)
         * which contain the velocities at every point.
         * velocity_values[0] is the table for the theta component, whereas
         * velocity_values[1] is the table for the phi component.
         */
        Table<2,double> velocity_values[2] = {Table<2,double>(n_theta,n_phi+1), Table<2,double>(n_theta,n_phi+1)};

        std::string velos = pt.get<std::string>("gpml:FeatureCollection.gml:featureMember.gpml:VelocityField.gml:rangeSet.gml:DataBlock.gml:tupleList");
        std::stringstream in(velos, std::ios::in);
        AssertThrow (in,
                     ExcMessage (std::string("Could not find velocities. Is file native gpml format for velocities?")));

        // The lat-lon mesh has changed its starting longitude in gplates1.4
        // correct for this while reading in the velocity data
        unsigned int longitude_correction = 0;
        if (gplates_1_4_or_higher(pt))
          longitude_correction = n_phi/2;

        unsigned int i = 0;
        char sep;
        Tensor<1,2> spherical_velocities;

        while (in >> spherical_velocities[0] >> sep >> spherical_velocities[1])
          {
            const double cmyr_si = 0.01/year_in_seconds;

            const unsigned int idx_theta = i / n_phi;
            const unsigned int idx_phi = (i + longitude_correction) % n_phi;

            velocity_values[0][idx_theta][idx_phi]= spherical_velocities[0] * cmyr_si;
            velocity_values[1][idx_theta][idx_phi]= spherical_velocities[1] * cmyr_si;

            i++;
          }

        // Pad the longitude data with values for phi == 2*pi (== 0),
        // this simplifies interpolation later.
        for (unsigned int i=0; i<n_theta; ++i)
          {
            velocity_values[0][i][n_phi] = velocity_values[0][i][0];
            velocity_values[1][i][n_phi] = velocity_values[1][i][0];
          }

        // number of intervals in the direction of theta and phi
        std::array<unsigned int,2> table_intervals;
        table_intervals[0] = n_theta - 1;
        table_intervals[1] = n_phi;

        // Min and Max coordinates in data file
        std::array<std::pair<double,double>,2> grid_extent;

        // min and max extent of the grid in the direction of theta and phi (whole spheres in GPlates)
        // polar angle theta: from 0° to 180°(PI)
        grid_extent[0].first = 0;
        grid_extent[0].second = numbers::PI;
        // azimuthal angle phi: from 0° to 360°(2*PI)
        grid_extent[1].first = 0;
        grid_extent[1].second = 2 * numbers::PI;

        for (unsigned int i = 0; i < 2; i++)
          {
            velocities[i]
              = std_cxx14::make_unique<Functions::InterpolatedUniformGridData<2>> (grid_extent,
                                                                                   table_intervals,
                                                                                   velocity_values[i]);
          }

        AssertThrow(i == n_points,
                    ExcMessage (std::string("Number of read in points does not match number of points in file. File corrupted?")));
      }



      template <int dim>
      Tensor<1,dim>
      GPlatesLookup<dim>::surface_velocity(const Point<dim> &position) const
      {
        const Point<3> internal_position ((dim == 2)
                                          ?
                                          rotation_matrix * convert_tensor<dim,3>(position)
                                          :
                                          convert_tensor<dim,3>(position));

        // transform internal_position in spherical coordinates
        std::array<double,3> spherical_point =
          Utilities::Coordinates::cartesian_to_spherical_coordinates(internal_position);

        Tensor<1,dim> output_boundary_velocity;
        // Handle all points that are not close to the poles
        if ((spherical_point[2] >= delta_theta) && (spherical_point[2] <= numbers::PI - delta_theta))
          {
            output_boundary_velocity = cartesian_velocity_at_surface_point(spherical_point);
          }

        // The longitude of data points at the poles is set to zero (according to the internal
        // GPlates routine). Because we interpolate velocities before converting them to the
        // cartesian coordinate system, we need to evaluate them twice: First at the poles,
        // and then at the latitude of the first data point that is not at the poles.
        // Afterwards we average the two velocities to the real position.
        else if (spherical_point[2] < delta_theta)
          {
            const double theta = spherical_point[2];
            spherical_point[2] = delta_theta;
            const Tensor<1,dim> first_velocity = cartesian_velocity_at_surface_point(spherical_point);

            spherical_point[1] = 0.0;
            spherical_point[2] = 0.0;
            const Tensor<1,dim> polar_velocity = cartesian_velocity_at_surface_point(spherical_point);

            output_boundary_velocity = (theta / delta_theta) * first_velocity
                                       + (1.0 - theta / delta_theta) * polar_velocity;
          }
        else if (spherical_point[2] > numbers::PI - delta_theta)
          {
            const double theta = spherical_point[2];
            spherical_point[2] = numbers::PI - delta_theta;
            const Tensor<1,dim> first_velocity = cartesian_velocity_at_surface_point(spherical_point);

            spherical_point[1] = 0.0;
            spherical_point[2] = numbers::PI;
            const Tensor<1,dim> polar_velocity = cartesian_velocity_at_surface_point(spherical_point);

            output_boundary_velocity = ((numbers::PI - theta) / delta_theta) * first_velocity
                                       + (1.0 - (numbers::PI - theta) / delta_theta) * polar_velocity;
          }
        else
          Assert(false,ExcInternalError());

        return output_boundary_velocity;
      }



      template <int dim>
      Tensor<1,dim>
      GPlatesLookup<dim>::cartesian_velocity_at_surface_point(const std::array<double,3> &spherical_point) const
      {
        // Re-sort the components of the spherical position from [r,phi,theta] to [theta, phi]
        const Point<2> lookup_coordinates(spherical_point[2],spherical_point[1]);

        // Main work, interpolate velocity at this point
        Tensor<1,2> interpolated_velocity;

        for (unsigned int i = 0; i < 2; i++)
          interpolated_velocity[i] = velocities[i]->value(lookup_coordinates);

        // Transform interpolated_velocity in cartesian coordinates
        const Tensor<1,3> interpolated_velocity_in_cart = sphere_to_cart_velocity(interpolated_velocity,spherical_point);

        // Convert_tensor conveniently also handles the projection to the 2D plane by
        // omitting the z-component of velocity (since the 2D model lies in the x-y plane).
        const Tensor<1,dim> output_boundary_velocity = (dim == 2)
                                                       ?
                                                       convert_tensor<3,dim>(transpose(rotation_matrix) * interpolated_velocity_in_cart)
                                                       :
                                                       convert_tensor<3,dim>(interpolated_velocity_in_cart);

        return output_boundary_velocity;
      }



      template <int dim>
      Tensor<1,3>
      GPlatesLookup<dim>::cartesian_surface_coordinates(const Tensor<1,3> &sposition) const
      {
        Tensor<1,3> ccoord;

        ccoord[0] = std::sin(sposition[0]) * std::cos(sposition[1]); // X
        ccoord[1] = std::sin(sposition[0]) * std::sin(sposition[1]); // Y
        ccoord[2] = std::cos(sposition[0]); // Z
        return ccoord;
      }



      template <int dim>
      Tensor<1,3>
      GPlatesLookup<dim>::sphere_to_cart_velocity(const Tensor<1,2> &s_velocities, const std::array<double,3> &s_position) const
      {
        Tensor<1,3> velocity;

        velocity[0] = std::cos(s_position[2]) * std::cos(s_position[1]) * s_velocities[0]
                      - 1.0 * std::sin(s_position[1]) * s_velocities[1];
        velocity[1] = std::cos(s_position[2]) * std::sin(s_position[1]) * s_velocities[0]
                      + std::cos(s_position[1]) * s_velocities[1];
        velocity[2] = -1.0 * std::sin(s_position[2]) * s_velocities[0];

        return velocity;
      }



      template <int dim>
      Tensor<2,3>
      GPlatesLookup<dim>::rotation_matrix_from_axis (const Tensor<1,3> &rotation_axis,
                                                     const double rotation_angle) const
      {
        Tensor<2,3> rotation_matrix;
        rotation_matrix[0][0] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[0] + std::cos(rotation_angle);
        rotation_matrix[0][1] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[1] - rotation_axis[2] * std::sin(rotation_angle);
        rotation_matrix[0][2] = (1-std::cos(rotation_angle)) * rotation_axis[0]*rotation_axis[2] + rotation_axis[1] * std::sin(rotation_angle);
        rotation_matrix[1][0] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[0] + rotation_axis[2] * std::sin(rotation_angle);
        rotation_matrix[1][1] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[1] + std::cos(rotation_angle);
        rotation_matrix[1][2] = (1-std::cos(rotation_angle)) * rotation_axis[1]*rotation_axis[2] - rotation_axis[0] * std::sin(rotation_angle);
        rotation_matrix[2][0] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[0] - rotation_axis[1] * std::sin(rotation_angle);
        rotation_matrix[2][1] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[1] + rotation_axis[0] * std::sin(rotation_angle);
        rotation_matrix[2][2] = (1-std::cos(rotation_angle)) * rotation_axis[2]*rotation_axis[2] + std::cos(rotation_angle);
        return rotation_matrix;
      }



      template <int dim>
      double
      GPlatesLookup<dim>::rotation_axis_from_matrix (Tensor<1,3> &rotation_axis,
                                                     const Tensor<2,3> &rotation_matrix) const
      {
        double rotation_angle = std::acos(0.5 * (rotation_matrix[0][0] + rotation_matrix[1][1] + rotation_matrix[2][2] - 1));

        if (rotation_angle > std::numeric_limits<double>::min())
          {
            rotation_axis[0] = (rotation_matrix[2][1] - rotation_matrix[1][2]) / (2*std::sin(rotation_angle));
            rotation_axis[1] = (rotation_matrix[0][2] - rotation_matrix[2][0]) / (2*std::sin(rotation_angle));
            rotation_axis[2] = (rotation_matrix[1][0] - rotation_matrix[0][1]) / (2*std::sin(rotation_angle));
          }
        else
          {
            rotation_axis[0] = 0.0;
            rotation_axis[1] = 0.0;
            rotation_axis[2] = 1.0;
          }

        return rotation_angle;
      }



      template <int dim>
      std::array<double,3>
      GPlatesLookup<dim>::angles_from_matrix(const Tensor<2,3> &rotation_matrix) const
      {
        std::array<double,3> orientation;

        /*
         * The following code is part of the VTK project and copied here for
         * compatibility to the paraview rotation formalism. It is protected by
         * the following license:
         *
         *
         * Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
         * All rights reserved.
         *
         * Redistribution and use in source and binary forms, with or without
         * modification, are permitted under certain conditions. See
         * http://www.kitware.com/Copyright.htm for details.

         * This software is distributed WITHOUT ANY WARRANTY; without even
         * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
         * PURPOSE.  See the above copyright notice for more information.
         *
         *
         * The original code in the VTK-6.0.0 source folder is found in:
         * /Common/Transforms/vtkTransform.cxx in the function GetOrientation()
         *
         */

        // first rotate about y axis
        const double x2 = rotation_matrix[2][0];
        const double y2 = rotation_matrix[2][1];
        const double z2 = rotation_matrix[2][2];

        const double x3 = rotation_matrix[1][0];
        const double y3 = rotation_matrix[1][1];
        const double z3 = rotation_matrix[1][2];

        double d1 = sqrt(x2*x2 + z2*z2);

        double cosTheta, sinTheta;
        if (d1 < std::numeric_limits<double>::min())
          {
            cosTheta = 1.0;
            sinTheta = 0.0;
          }
        else
          {
            cosTheta = z2/d1;
            sinTheta = x2/d1;
          }

        double theta = atan2(sinTheta, cosTheta);
        orientation[1] = - theta * 180 / numbers::PI;

        // now rotate about x axis
        double d = sqrt(x2*x2 + y2*y2 + z2*z2);

        double sinPhi, cosPhi;
        if (d < std::numeric_limits<double>::min())
          {
            sinPhi = 0.0;
            cosPhi = 1.0;
          }
        else if (d1 < std::numeric_limits<double>::min())
          {
            sinPhi = y2/d;
            cosPhi = z2/d;
          }
        else
          {
            sinPhi = y2/d;
            cosPhi = (x2*x2 + z2*z2)/(d1*d);
          }

        double phi = atan2(sinPhi, cosPhi);
        orientation[0] = phi * 180 / numbers::PI;

        // finally, rotate about z
        double x3p = x3*cosTheta - z3*sinTheta;
        double y3p = - sinPhi*sinTheta*x3 + cosPhi*y3 - sinPhi*cosTheta*z3;
        double d2 = sqrt(x3p*x3p + y3p*y3p);

        double cosAlpha, sinAlpha;
        if (d2 < std::numeric_limits<double>::min())
          {
            cosAlpha = 1.0;
            sinAlpha = 0.0;
          }
        else
          {
            cosAlpha = y3p/d2;
            sinAlpha = x3p/d2;
          }

        double alpha = atan2(sinAlpha, cosAlpha);
        orientation[2] = alpha * 180 / numbers::PI;
        return orientation;
      }



      template <int dim>
      template <int in, int out>
      Tensor<1,out>
      GPlatesLookup<dim>::convert_tensor (const Tensor<1,in> &old_tensor) const
      {
        Tensor<1,out> new_tensor;
        for (unsigned int i = 0; i < out; i++)
          if (i < in) new_tensor[i] = old_tensor[i];
          else new_tensor[i] = 0.0;

        return new_tensor;
      }



      template <int dim>
      bool
      GPlatesLookup<dim>::gplates_1_4_or_higher(const boost::property_tree::ptree &pt) const
      {
        const std::string gpml_version = pt.get<std::string>("gpml:FeatureCollection.<xmlattr>.gpml:version");
        const std::vector<std::string> string_versions = dealii::Utilities::split_string_list(gpml_version,'.');
        const std::vector<int> int_versions = dealii::Utilities::string_to_int(string_versions);

        const int gplates_1_3_version[3] = {1,6,322};

        for (unsigned int i = 0; i < int_versions.size(); i++)
          {
            if (int_versions[i] > gplates_1_3_version[i])
              return true;
            if (int_versions[i] < gplates_1_3_version[i])
              return false;
          }

        return false;
      }
    }



    template <int dim>
    GPlates<dim>::GPlates ()
      :
      current_file_number(0),
      first_data_file_model_time(0.0),
      first_data_file_number(0),
      decreasing_file_order(false),
      data_file_time_step(0.0),
      time_weight(0.0),
      time_dependent(true),
      point1("0.0,0.0"),
      point2("0.0,0.0"),
      lithosphere_thickness(0.0),
      lookup(),
      old_lookup()
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

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
                   ExcMessage ("This gplates plugin can only be used when the "
                               "preferred coordinate system of the geometry model is spherical "
                               "(e.g. spherical shell, chunk, sphere)."));

      lookup = std_cxx14::make_unique<internal::GPlatesLookup<dim>>(pointone, pointtwo);
      old_lookup = std_cxx14::make_unique<internal::GPlatesLookup<dim>>(pointone, pointtwo);

      // display the GPlates module information at model start.
      this->get_pcout() << lookup->screen_output(pointone, pointtwo);

      // Set the first file number and load the first files
      current_file_number = first_data_file_number;

      const int next_file_number =
        (decreasing_file_order)
        ?
        (current_file_number - 1)
        :
        (current_file_number + 1);

      this->get_pcout() << std::endl << "   Loading GPlates data boundary file "
                        << create_filename (current_file_number) << "." << std::endl << std::endl;

      const std::string filename (create_filename (current_file_number));
      if (Utilities::fexists(filename))
        lookup->load_file(filename,this->get_mpi_communicator());
      else
        AssertThrow(false,
                    ExcMessage (std::string("GPlates data file <")
                                +
                                filename
                                +
                                "> not found!"));

      // If the boundary condition is constant, switch off time_dependence
      // immediately. If not, also load the second file for interpolation.
      // This catches the case that many files are present, but the
      // parameter file requests a single file.
      if (create_filename (current_file_number) == create_filename (current_file_number+1))
        {
          end_time_dependence ();
        }
      else
        {
          const std::string filename (create_filename (next_file_number));
          this->get_pcout() << std::endl << "   Loading GPlates data boundary file "
                            << filename << "." << std::endl << std::endl;
          if (Utilities::fexists(filename))
            {
              lookup.swap(old_lookup);
              lookup->load_file(filename,this->get_mpi_communicator());
            }
          else
            end_time_dependence ();
        }
    }



    template <int dim>
    std::string
    GPlates<dim>::create_filename (const int timestep) const
    {
      std::string templ = data_directory+velocity_file_name;
      const int size = templ.length();
      std::vector<char> buffer(size+10);
      snprintf (buffer.data(), size + 10, templ.c_str(), timestep);
      std::string str_filename (buffer.data());
      return str_filename;
    }



    template <int dim>
    void
    GPlates<dim>::update ()
    {
      const double time_since_start = this->get_time() - first_data_file_model_time;

      if (time_dependent && (time_since_start >= 0.0))
        {
          const int time_steps_since_start = static_cast<int> (time_since_start / data_file_time_step);

          // whether we need to update our data files. This looks so complicated
          // because we need to catch increasing and decreasing file orders and all
          // possible first_data_file_model_times and first_data_file_numbers.
          const bool need_update = time_steps_since_start
                                   > std::abs(current_file_number - first_data_file_number);

          if (need_update)
            {
              // The last file, which was tried to be loaded was
              // number current_file_number +/- 1, because current_file_number
              // is the file older than the current model time
              const int old_file_number =
                (decreasing_file_order)
                ?
                (current_file_number - 1)
                :
                (current_file_number + 1);

              // Calculate new file_number
              current_file_number =
                (decreasing_file_order)
                ?
                (first_data_file_number - time_steps_since_start)
                :
                (first_data_file_number + time_steps_since_start);

              const bool load_both_files = std::abs(current_file_number - old_file_number) >= 1;

              update_data(load_both_files);
            }

          time_weight = (time_since_start / data_file_time_step)
                        - std::abs(current_file_number - first_data_file_number);

          Assert ((0 <= time_weight) && (time_weight <= 1),
                  ExcMessage (
                    "Error in set_current_time. Time_weight has to be in [0,1]"));
        }
    }



    template <int dim>
    void
    GPlates<dim>::update_data (const bool load_both_files)
    {
      // If the time step was large enough to move forward more
      // then one data file we need to load both current files
      // to stay accurate in interpolation
      if (load_both_files)
        {
          const std::string filename (create_filename (current_file_number));
          this->get_pcout() << std::endl << "   Loading GPlates data boundary file "
                            << filename << "." << std::endl << std::endl;
          if (Utilities::fexists(filename))
            {
              lookup.swap(old_lookup);
              lookup->load_file(filename,this->get_mpi_communicator());
            }

          // If loading current_time_step failed, end time dependent part with old_file_number.
          else
            end_time_dependence ();
        }

      // Now load the next data file. This part is the main purpose of this function.
      const int next_file_number =
        (decreasing_file_order) ?
        current_file_number - 1
        :
        current_file_number + 1;

      const std::string filename (create_filename (next_file_number));
      this->get_pcout() << std::endl << "   Loading GPlates data boundary file "
                        << filename << "." << std::endl << std::endl;
      if (Utilities::fexists(filename))
        {
          lookup.swap(old_lookup);
          lookup->load_file(filename,this->get_mpi_communicator());
        }

      // If next file does not exist, end time dependent part with current_time_step.
      else
        end_time_dependence ();
    }



    template <int dim>
    void
    GPlates<dim>::end_time_dependence ()
    {
      // no longer consider the problem time dependent from here on out
      // this cancels all attempts to read files at the next time steps
      time_dependent = false;
      // Give warning if first processor
      this->get_pcout() << std::endl
                        << "   Loading new gplates velocity file did not succeed." << std::endl
                        << "   Assuming constant boundary conditions for rest of model run."
                        << std::endl << std::endl;
    }



    template <int dim>
    Tensor<1,dim>
    GPlates<dim>::
    boundary_velocity (const types::boundary_id /*boundary_indicator*/,
                       const Point<dim> &position) const
    {
      // We compare the depth of the current point to the lithosphere thickness.
      // The depth is calculated using squares, sums, square-roots and differences
      // of large numbers, and we possibly end up with a very small number close to
      // the surface, thus rounding errors can get quite large. We therefore
      // compare the depth to lithosphere_thickness plus a magic number, which we
      // choose as 1e-7 times the maximal model depth, because we safely assume no
      // model will have more than 1e7 quadrature points in depth direction.
      // Without the magic number it may unintentionally happen that the GPlates
      // velocities are not prescribed at every point on the surface.
      const double magic_number = 1e-7 * this->get_geometry_model().maximal_depth();

      if ((this->get_time() - first_data_file_model_time >= 0.0) && (this->get_geometry_model().depth(position) <= lithosphere_thickness + magic_number))
        {
          const Tensor<1,dim> data = lookup->surface_velocity(position);

          if (!time_dependent)
            return data;

          const Tensor<1,dim> old_data = old_lookup->surface_velocity(position);

          return time_weight * data + (1 - time_weight) * old_data;
        }
      else
        return Tensor<1,dim>();
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
                             "$ASPECT_SOURCE_DIR/data/boundary-velocity/gplates/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the model data. This path "
                             "may either be absolute (if starting with a '/') or relative to "
                             "the current directory. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Velocity file name", "phi.%d",
                             Patterns::Anything (),
                             "The file name of the material data. Provide file in format: "
                             "(Velocity file name).\\%d.gpml where \\%d is any sprintf integer "
                             "qualifier, specifying the format of the current file number.");
          prm.declare_entry ("First data file model time", "0.",
                             Patterns::Double (0.),
                             "Time from which on the velocity file with number 'First velocity "
                             "file number' is used as boundary condition. Previous to this "
                             "time, a no-slip boundary condition is assumed. Depending on the setting "
                             "of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years.");
          prm.declare_entry ("First data file number", "0",
                             Patterns::Integer (),
                             "Number of the first velocity file to be loaded when the model time "
                             "is larger than 'First velocity file model time'.");
          prm.declare_entry ("Decreasing file order", "false",
                             Patterns::Bool (),
                             "In some cases the boundary files are not numbered in increasing "
                             "but in decreasing order (e.g. 'Ma BP'). If this flag is set to "
                             "'True' the plugin will first load the file with the number "
                             "'First velocity file number' and decrease the file number during "
                             "the model run.");
          prm.declare_entry ("Data file time step", "1e6",
                             Patterns::Double (0.),
                             "Time step between following velocity files. "
                             "Depending on the setting of the global 'Use years in output instead of seconds' flag "
                             "in the input file, this number is either interpreted as seconds or as years. "
                             "The default is one million, i.e., either one million seconds or one million years.");
          prm.declare_entry ("Scale factor", "1.",
                             Patterns::Double (),
                             "Scalar factor, which is applied to the boundary velocity. "
                             "You might want to use this to scale the velocities to a "
                             "reference model (e.g. with free-slip boundary) or another "
                             "plate reconstruction.");
          prm.declare_entry ("Point one", "1.570796,0.0",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in. Has to be in the format `a,b' where a and b are theta (polar angle) and "
                             "phi in radians. This value is not utilized in 3D geometries, and can therefore be set to the default or any user-defined quantity.");
          prm.declare_entry ("Point two", "1.570796,1.570796",
                             Patterns::Anything (),
                             "Point that determines the plane in which a 2D model lies in. Has to be in the format `a,b' where a and b are theta (polar angle) and "
                             "phi in radians. This value is not utilized in 3D geometries, and can therefore be set to the default or any user-defined quantity.");
          prm.declare_entry ("Lithosphere thickness", "100000.",
                             Patterns::Double (0.),
                             "Determines the depth of the lithosphere, so that the GPlates velocities can be applied at the sides of the model "
                             "as well as at the surface.");
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
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));

          velocity_file_name         = prm.get        ("Velocity file name");
          data_file_time_step        = prm.get_double ("Data file time step");
          first_data_file_model_time = prm.get_double ("First data file model time");
          first_data_file_number     = prm.get_integer("First data file number");
          decreasing_file_order      = prm.get_bool   ("Decreasing file order");
          scale_factor               = prm.get_double ("Scale factor");
          point1                     = prm.get        ("Point one");
          point2                     = prm.get        ("Point two");
          lithosphere_thickness      = prm.get_double ("Lithosphere thickness");

          if (this->convert_output_to_years())
            {
              data_file_time_step        *= year_in_seconds;
              first_data_file_model_time *= year_in_seconds;
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
  namespace BoundaryVelocity
  {
    namespace internal
    {
#define INSTANTIATE(dim) \
  template class GPlatesLookup<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(GPlates,
                                            "gplates",
                                            "Implementation of a model in which the boundary "
                                            "velocity is derived from files that are generated "
                                            "by the GPlates program.")
  }
}
