/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_velocity_gplates_h
#define _aspect_boundary_velocity_gplates_h

#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>

#include <array>
#include <deal.II/base/function_lib.h>


namespace aspect
{
  namespace BoundaryVelocity
  {
    using namespace dealii;

    namespace internal
    {
      /**
       * GPlatesLookup handles all kinds of tasks around looking up a certain
       * velocity boundary condition from a gplates .gpml file.
       */
      template <int dim>
      class GPlatesLookup
      {
        public:

          /**
           * Initialize all members and calculates any necessary rotation
           * parameters for a 2D model.
           */
          GPlatesLookup(const Tensor<1,2> &surface_point_one,
                        const Tensor<1,2> &surface_point_two);

          /**
           * Outputs the GPlates module information at model start.
           */
          std::string
          screen_output(const Tensor<1,2> &surface_point_one,
                        const Tensor<1,2> &surface_point_two) const;

          /**
           * Loads a gplates .gpml velocity file. Throws an exception if the
           * file does not exist.
           */
          void load_file(const std::string &filename,
                         const MPI_Comm &comm);

          /**
           * Returns the computed surface velocity in cartesian coordinates.
           * Takes as input the position. Actual velocity interpolation is
           * performed in spherical coordinates.
           *
           * @param position The current position to compute velocity
           */
          Tensor<1,dim> surface_velocity(const Point<dim> &position) const;

        private:
          /**
           * Interpolation functions to access the velocities.
           */
          std::array<std::unique_ptr<Functions::InterpolatedUniformGridData<2>>, 2> velocities;

          /**
           * Distances between adjacent point in the Lat/Long grid
           */
          double delta_phi, delta_theta;

          /**
           * The matrix, which describes the rotation by which a 2D model
           * needs to be transformed to a plane that contains the origin and
           * the two user prescribed points. Is not necessary and therefore
           * not used for 3D models.
           */
          Tensor<2,3> rotation_matrix;

          /**
           * A function that returns the corresponding paraview angles for a
           * rotation described by a rotation matrix. These differ from the
           * usually used euler angles by assuming a rotation around the
           * coordinate axes in the order y-x-z (instead of the often used
           * z-x-z)
           */
          std::array<double,3>
          angles_from_matrix (const Tensor<2,3> &rotation_matrix) const;

          /**
           * A function that returns the corresponding rotation axis/angle for
           * a rotation described by a rotation matrix.
           */
          double
          rotation_axis_from_matrix (Tensor<1,3> &rotation_axis,
                                     const Tensor<2,3> &rotation_matrix) const;

          /**
           * A function that returns the corresponding euler angles for a
           * rotation described by rotation axis and angle.
           */
          Tensor<2,3>
          rotation_matrix_from_axis (const Tensor<1,3> &rotation_axis,
                                     const double rotation_angle) const;

          /**
           * Convert a tensor of rank 1 and dimension in to rank 1 and
           * dimension out. If $out < in$ the last elements will be discarded,
           * if $out > in$ zeroes will be appended to fill the tensor.
           */
          template <int in, int out>
          Tensor<1,out> convert_tensor (const Tensor<1,in> &old_tensor) const;

          /**
           * Return the cartesian coordinates of a spherical surface position
           * defined by theta (polar angle, not geographical latitude) and
           * phi.
           */
          Tensor<1,3>
          cartesian_surface_coordinates(const Tensor<1,3> &sposition) const;

          /**
           * This function looks up the north- and east-velocities at a given
           * position and converts them to cartesian velocities.
           */
          Tensor<1,dim>
          cartesian_velocity_at_surface_point(const std::array<double,3> &spherical_point) const;

          /**
           * Returns cartesian velocities calculated from surface velocities
           * and position in spherical coordinates
           *
           * @param s_velocities Surface velocities in spherical coordinates
           * (theta, phi)
           * @param s_position Position in spherical coordinates
           * (radius,phi,theta)
           */
          Tensor<1,3> sphere_to_cart_velocity(const Tensor<1,2> &s_velocities,
                                              const std::array<double,3> &s_position) const;

          /**
           * Check whether the gpml file was created by GPlates1.4 or later.
           * We need to know this, because the mesh has changed its longitude
           * origin from 0 to -180 degrees and we need to correct for this.
           */
          bool
          gplates_1_4_or_higher(const boost::property_tree::ptree &pt) const;
      };
    }

    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from GPlates input files. The interpolation in time is
     * performed between two objects of the GPlatesLookup class.
     *
     * @ingroup BoundaryVelocities
     */
    template <int dim>
    class GPlates : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        GPlates ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns value from gplates.
         */
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id boundary_indicator,
                           const Point<dim> &position) const override;

        // avoid -Woverloaded-virtual warning until the deprecated function
        // is removed from the interface:
        using Interface<dim>::boundary_velocity;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        /**
         * A function that is called at the beginning of each time step. For
         * the current plugin, this function loads the next velocity files if
         * necessary and outputs a warning if the end of the set of velocity
         * files is reached.
         */
        void
        update () override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * A variable that stores the currently used data file of a series. It
         * gets updated if necessary by update().
         */
        int current_file_number;

        /**
         * Time from which on the data file with number 'First data file
         * number' is used as boundary condition. Previous to this time, 0 is
         * returned for every field. Depending on the setting of the global
         * 'Use years in output instead of seconds' flag in the input file,
         * this number is either interpreted as seconds or as years."
         */
        double first_data_file_model_time;

        /**
         * Number of the first data file to be loaded when the model time is
         * larger than 'First data file model time'.
         */
        int first_data_file_number;

        /**
         * In some cases the boundary files are not numbered in increasing but
         * in decreasing order (e.g. 'Ma BP'). If this flag is set to 'True'
         * the plugin will first load the file with the number 'First data
         * file number' and decrease the file number during the model run.
         */
        bool decreasing_file_order;

        /**
         * Time in model units (depends on other model inputs) between two
         * velocity files.
         */
        double data_file_time_step;

        /**
         * Weight between velocity file n and n+1 while the current time is
         * between the two values t(n) and t(n+1).
         */
        double time_weight;

        /**
         * State whether we have time_dependent boundary conditions. Switched
         * off after finding no more velocity files to suppress attempts to read
         * in new files.
         */
        bool time_dependent;

        /**
         * Directory in which the gplates velocity files are present.
         */
        std::string data_directory;

        /**
         * First part of filename of velocity files. The files have to have
         * the pattern velocity_file_name.n.gpml where n is the number of the
         * current timestep (starts from 0).
         */
        std::string velocity_file_name;

        /**
         * Scale the velocity boundary condition by a scalar factor.
         */
        double scale_factor;

        /**
         * Two user defined points that prescribe the plane from which the 2D
         * model takes the velocity boundary condition. One can think of this,
         * as if the model is lying in this plane although no actual model
         * coordinate is changed. The strings need to have the format "a,b"
         * where a and b are doubles and define theta and phi on a sphere.
         */
        std::string point1;
        std::string point2;

        /**
         * Parsed user input of point1 and point2
         */
        Tensor<1,2> pointone;
        Tensor<1,2> pointtwo;

        /**
         * Determines the depth of the lithosphere. The user might want to apply
         * the GPlates velocities not only at the surface of the model, but also
         * in the whole lithosphere. At every side boundary point with a depth
         * smaller than this value (and thus being located in the lithosphere),
         * the surface velocity will be described.
         */
        double lithosphere_thickness;

        /**
         * Pointer to an object that reads and processes data we get from
         * gplates files.
         */
        std::unique_ptr<internal::GPlatesLookup<dim> > lookup;

        /**
         * Pointer to an object that reads and processes data we get from
         * gplates files. This saves the previous data time step.
         */
        std::unique_ptr<internal::GPlatesLookup<dim> > old_lookup;

        /**
         * Handles the update of the velocity data in lookup. The input
         * parameter makes sure that both velocity files (n and n+1) can be
         * reloaded if the model time step is larger than the velocity file
         * time step.
         */
        void
        update_data (const bool load_both_files);

        /**
         * Handles settings and user notification in case the time-dependent
         * part of the boundary condition is over.
         */
        void
        end_time_dependence ();

        /**
         * Create a filename out of the name template.
         */
        std::string
        create_filename (const int timestep) const;
    };
  }
}


#endif
