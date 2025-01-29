/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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


#include <aspect/boundary_traction/ascii_data.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace BoundaryTraction
  {
    /**
     * A class that implements prescribed traction boundary conditions determined
     * from pressures given in an AsciiData input file.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class InfillAsciiData : public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        InfillAsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataBoundary<dim>::initialize;

        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        Tensor<1,dim>
        boundary_traction (const types::boundary_id surface_boundary_id,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        void
        update () override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        types::boundary_id surface_boundary_id;
        double rock_density;
        double sediment_density;
        double crustal_density;
        double rock_infill_height;
    };

    template <int dim>
    InfillAsciiData<dim>::InfillAsciiData ()
      :
      surface_boundary_id(numbers::invalid_unsigned_int)
    {}


    template <int dim>
    void
    InfillAsciiData<dim>::initialize ()
    {
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(surface_boundary_id);

      // The input ascii table contains one data column (Traction [Pa]) due to
      // the load in addition to the coordinate columns.
      Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,
                                                    1);
    }


    template <int dim>
    Tensor<1,dim>
    InfillAsciiData<dim>::
    boundary_traction (const types::boundary_id surface_boundary_id,
                       const Point<dim> &position,
                       const Tensor<1,dim> &normal_vector) const
    {
      const double load = Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id,
                                                                                position,
                                                                                0);
      Tensor<1, dim> traction = normal_vector;
      const double gravity_norm = this->get_gravity_model().gravity_vector(position).norm();
      const double elevation = this->get_geometry_model().height_above_reference_surface(position);

      // Check if the load at the current point has a height which exceeds the
      // 'rock_infill height'. The assumption is that the values taken from the
      // ASCII file are specifying the part of the load that is above the undeformed
      // reference surface. If we consider a seamount as the load, the ASCII
      // file takes the bathymetry of the seamount and provides it as a traction.
      // However, the volcanic rock comprising the seamount also fills in the
      // flexural moat, providing an additional traction. Since it is much harder to
      // specify this traction, 'rock_infill_height' determines where the rock_density
      // that the seamount is made of will infill the flexural moat, and where
      // sediment_density will infill the flexural moat. For this test, this
      // variable is not that important given the idealized load, but this is useful
      // when using more complicated bathymetry maps where small scale seafloor
      // features makes it unrealistic to set rock_infill_height=0.
      if (load >= rock_infill_height*gravity_norm*rock_density)
        {
          traction = (-load + elevation*gravity_norm*rock_density) * normal_vector;
        }
      else
        {
          traction = (-load + elevation*gravity_norm*sediment_density) * normal_vector;
        }
      return traction;
    }


    template <int dim>
    void
    InfillAsciiData<dim>::update()
    {
      Interface<dim>::update ();
      Utilities::AsciiDataBoundary<dim>::update();
    }


    template <int dim>
    void
    InfillAsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        // Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
        //                                                       "$ASPECT_SOURCE_DIR/tests/infill_density/",
        //                                                       "box_2d_%s.%d.txt");
        prm.enter_subsection("Infill ascii data");
        {
          Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                                "$ASPECT_SOURCE_DIR/tests/infill_density/",
                                                                "box_2d_%s.%d.txt");
          prm.declare_entry ("Rock density", "2800",
                             Patterns::Double(0.),
                             "Density of the volcanic edifice that infills the flexural moat.");
          prm.declare_entry ("Sediment density", "2300",
                             Patterns::Double(0.),
                             "Density of the sediment that infills the flexural moat.");
          prm.declare_entry ("Height for specifying rock infill", "500",
                             Patterns::Double(0.),
                             "If the load defined in the ASCII file has a height equal to or greater than "
                             "the Height for specifying rock infill, then the infill material will be rock, "
                             "otherwise it will be sediment");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    InfillAsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
      {
        // Utilities::AsciiDataBoundary<dim>::parse_parameters(prm);
        prm.enter_subsection("Infill ascii data");
        {
          Utilities::AsciiDataBoundary<dim>::parse_parameters(prm);
          rock_density = prm.get_double("Rock density");
          sediment_density = prm.get_double("Sediment density");
          rock_infill_height = prm.get_double("Height for specifying rock infill");
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
  namespace BoundaryTraction
  {
    ASPECT_REGISTER_BOUNDARY_TRACTION_MODEL(InfillAsciiData,
                                            "infill ascii data",
                                            "Implementation of a model in which the boundary "
                                            "traction is derived from files containing pressure data "
                                            "in ascii format. This model differs from the normal ascii data "
                                            "model as the pressure given in the data file is summed with another "
                                            "traction that depends on the current state of flexure in the model. "
                                            "The traction sum is applied normal to the surface of a given boundary, "
                                            "pointing inwards. Note the required format of the "
                                            "input data: The first lines may contain any number of comments "
                                            "if they begin with `#', but one of these lines needs to "
                                            "contain the number of grid points in each dimension as "
                                            "for example `# POINTS: 3 3'. "
                                            "The order of the data columns "
                                            "has to be `x', `Pressure [Pa]' in a 2d model and "
                                            " `x', `y', `Pressure [Pa]' in a 3d model, which means that "
                                            "there has to be a single column "
                                            "containing the pressure. "
                                            "Note that the data in the input "
                                            "files need to be sorted in a specific order: "
                                            "the first coordinate needs to ascend first, "
                                            "followed by the second in order to "
                                            "assign the correct data to the prescribed coordinates. "
                                            "If you use a spherical model, "
                                            "then the data will still be handled as Cartesian, "
                                            "however the assumed grid changes. `x' will be replaced by "
                                            "the radial distance of the point to the bottom of the model, "
                                            "`y' by the azimuth angle and `z' by the polar angle measured "
                                            "positive from the north pole. The grid will be assumed to be "
                                            "a latitude-longitude grid. Note that the order "
                                            "of spherical coordinates is `r', `phi', `theta' "
                                            "and not `r', `theta', `phi', since this allows "
                                            "for dimension independent expressions.")
  }
}
