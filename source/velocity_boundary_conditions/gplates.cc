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
/*  $Id: zero_velocity.cc 1071 2012-08-01 16:50:02Z bangerth $  */


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
      /**
       * GPlatesLookup handles all kinds of tasks around looking up a certain velocity boundary
       * condition from a gplates .gpml file. This class keeps around the contents of two sets
       * of files, corresponding to two instances in time where GPlates provides us with data;
       * the boundary values at one particular time are interpolated between the two currently
       * loaded data sets.
       */
      class GPlatesLookup
      {
        public:
          /**
           * Initialize all members and the two pointers referring to the actual velocities
           */
          GPlatesLookup()
          {
            velocity_vals.reinit(0,0);
            old_velocity_vals.reinit(0,0);

            velocity_values = &velocity_vals;
            old_velocity_values = &old_velocity_vals;

            delta_phi = 0.0;
            delta_theta = 0.0;
          }

          /**
           * Check whether a file named filename exists.
           */
          bool fexists(const std::string &filename) const
          {
            std::ifstream ifile(filename.c_str());
            return ifile;
          }

          /**
           * Loads a gplates .gpml velocity file. Throws an exception if
           * the file does not exist.
           */
          void load_file(const std::string &filename)
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
            const int n_phi = std::sqrt(2*n_points);
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


          /**
           * Returns the computed surface velocity in cartesian coordinates. Takes
           * as input the position and current time weight.
           *
           * @param position The current position to compute velocity @param time_weight A weighting between
           * the two current timesteps n and n+1
           */
          template <int dim>
          Tensor<1,dim> surface_velocity(const Point<dim> &position, const double time_weight) const
          {
            const Point<dim> scoord = spherical_surface_coordinates(position);

            const double idtheta = get_idtheta(scoord(0));
            const unsigned int idxtheta = static_cast<unsigned int>(idtheta);
            const unsigned int nexttheta = idxtheta + 1;


            const double idphi = get_idphi(scoord(1));
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


            return sphere_to_cart_velocities(surf_vel,scoord);
          }


        private:

          /**
           * Tables which contain the velocities
           */
          dealii::Table<2,Tensor<1,2> > velocity_vals;
          dealii::Table<2,Tensor<1,2> > old_velocity_vals;

          /**
           * Pointers to the actual tables.
           * Used to avoid unnecessary copying
           * of values.
           */
          dealii::Table<2,Tensor<1,2> > * velocity_values;
          dealii::Table<2,Tensor<1,2> > * old_velocity_values;

          /**
           * Distances between adjacent point in the Lat/Lon grid
           */
          double delta_phi,delta_theta;

          /**
           * Returns spherical coordinates of a cartesian position.
           */
          template <int dim>
          Point<dim> spherical_surface_coordinates(const Point<dim> &position) const
          {
            double coord[dim];
            Point<dim> scoord;

            if (dim == 3)
              {
                coord[0] = std::acos(position(2)/std::sqrt(position.square())); // Theta
                coord[1] = std::atan2(position(1),position(0)); // Phi
                if (coord[1] < 0.0) coord[1] = 2*numbers::PI + coord[1]; // correct phi to [0,2*pi]
                coord[2] = std::sqrt(position.square()); // R

                return Point<dim> (coord[0],coord[1],coord[2]);
              }
            else
              return Point<dim>(); // TODO: define for 2 Dimensions
          }

          /**
           * Returns cartesian velocities calculated from surface velocities and position in spherical coordinates
           *
           * @param s_velocities Surface velocities in spherical coordinates (theta, phi) @param s_position Position
           * in spherical coordinates (theta,phi,radius)
           */
          template <int dim>
          Tensor<1,dim> sphere_to_cart_velocities(const Tensor<1,2> &s_velocities, const Point<dim> &s_position) const
          {
            Tensor<1,dim> velocities;
            if (dim == 3)
              {
                velocities[0] = -1.0 * s_velocities[1] * std::sin(s_position[1]) + s_velocities[0]*std::cos(s_position[0])*std::cos(s_position[1]);
                velocities[1] = s_velocities[1]*std::cos(s_position[1])+s_velocities[0]*std::cos(s_position[0])*std::sin(s_position[1]);
                velocities[2] = -1.0*s_velocities[0]*std::sin(s_position[0]);
                return velocities;
              }
            else
              return Tensor<1,dim>();  // TODO: define for 2 Dimensions
          }

          /**
           * calculates the phi-index given a certain phi
           */
          double get_idphi(const double phi_) const
          {
            const double phi = std::max(std::min(phi_,2*numbers::PI-1e-7),0.0);

            Assert(phi>=0, ExcMessage("not in range"));
            Assert(phi<=2*numbers::PI, ExcMessage("not in range"));
            return phi/delta_phi;
          }

          /**
           * calculates the theta-index given a certain polar angle
           */
          double get_idtheta(const double theta_) const
          {
            double theta = std::max(std::min(theta_,numbers::PI-1e-7),0.0);

            Assert(theta>=0, ExcMessage("not in range"));
            Assert(theta<=numbers::PI, ExcMessage("not in range"));
            return theta/delta_theta;
          }


      };
    }


    template <int dim>
    GPlates<dim>::GPlates ()
      :
      lookup (new internal::GPlatesLookup)
    {}



    template <int dim>
    void
    GPlates<dim>::initialize (const GeometryModel::Interface<dim> &geometry_model)
    {
      int loaded;

      Assert(dim==3, ExcMessage("gplates velocity model just works for dim = 3"));
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&geometry_model) != 0,
              ExcMessage ("This boundary condition can only be used if the geometry "
                          "is a spherical shell."));

      // Initialize variables
      current_time_step = 0;
      time_weight = 0.0;

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
      current_time = time;
      if (time_dependent)
        {
          Assert ((0 <= time_weight) && (time_weight <= 1),
                  ExcMessage("Error in set_current_time. Time_weight has to be in [0,1]"));

          if (static_cast<unsigned int>(current_time / time_step) > current_time_step)
            {
              current_time_step = static_cast<unsigned int>(current_time / time_step);

              // Load next velocity file for interpolation
              try
                {
                  lookup->load_file(create_filename(current_time_step+1));
                }
              catch (...)
                // of loading failed, simply take the previous file a second time
                {
                  // Next velocity file not found --> Constant velocities
                  lookup->load_file(create_filename(current_time_step));
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

        }
    }


    template <int dim>
    Tensor<1,dim>
    GPlates<dim>::
    boundary_velocity (const Point<dim> &position) const
    {
      return lookup->surface_velocity(position,time_weight);
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
          velocity_file_name   = prm.get ("Velocity file name");
          time_step             = prm.get_double ("Time step");
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
                                                 "velocity is derived from GPlates.");
  }
}
