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

#include <aspect/material_model/replace_lithosphere_viscosity.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <limits>

#include <array>
#include <utility>



namespace aspect
{
  namespace MaterialModel
  {
  template <int dim>
  ReplaceLithosphereViscosity<dim>::ReplaceLithosphereViscosity ()
    :
    lab_depths(1, 1.0)
  {}

    template <int dim>
    void
    ReplaceLithosphereViscosity<dim>::initialize()
    {
      base_model->initialize();

      //Get LAB depths
      const std::string filename = data_directory+LAB_file_name;
      std::cout << "   Loading Ascii data lookup file " << filename << "." << std::endl;

      lab_depths.load_file(filename,this->get_mpi_communicator());
    }

    //Read in ascii data - two dimensions so the third column is treated as data
    template <int dim>
    double
	ReplaceLithosphereViscosity<dim>::ascii_lab (const Point<2> &position) const
    {
      const double lab = lab_depths.get_data(position,0)*1000; //In km
      return lab;
    }

    template <int dim>
    void
	ReplaceLithosphereViscosity<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                                  typename Interface<dim>::MaterialModelOutputs &out) const
    {
      base_model->evaluate(in,out);

      for (unsigned int i=0; i < in.position.size(); ++i)
      {
	  double depth = this->SimulatorAccess<dim>::get_geometry_model().depth(in.position[i]);

      //Get spherical coordinates for model
      std::array<double,dim> scoord      = Utilities::Coordinates::cartesian_to_spherical_coordinates(in.position[i]);
      const double phi = scoord[1];
      const double theta = scoord[2];
      const Point<2> phi_theta (phi, theta);

      //Get lab depth for specific phi and theta
      const double lab_depth = ascii_lab(phi_theta);

      if (depth <= lab_depth)
        out.viscosities[i] = lithosphere_viscosity;
      else
    	out.viscosities[i] *= 1;
      }
    }



    template <int dim>
    void
	ReplaceLithosphereViscosity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Replace lithosphere viscosity");
        {
          prm.declare_entry("Base model","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by a depth "
                            "dependent viscosity. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.declare_entry ("Lithosphere viscosity", "1600",
                             Patterns::Double (0),
                             "The viscosity within lithosphere, applied above"
                             "the maximum lithosphere depth.");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/replace-lithosphere-viscosity/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text `$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("LAB depth filename",
                             "LAB_CAM2016.txt",
                             Patterns::FileName (),
                             "File from which the lithosphere-asthenosphere boundary depth data is read.");
        }
          prm.leave_subsection();
        }
        prm.leave_subsection();
    }

    template <int dim>
    void
	ReplaceLithosphereViscosity<dim>::parse_parameters (ParameterHandler &prm)
    {
        AssertThrow (dim == 3,
                     ExcMessage ("The 'Replace lithosphere viscosity' material model "
                                 "is only available for 3d computations."));

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Replace lithosphere viscosity");
        {
          AssertThrow( prm.get("Base model") != "replace lithosphere viscosity",
                       ExcMessage("You may not use ``replace lithosphere viscosity'' as the base model for "
                                  "a replace lithosphere viscosity model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          lithosphere_viscosity   = prm.get_double ("Lithosphere viscosity");

          data_directory = Utilities::expand_ASPECT_SOURCE_DIR (prm.get("Data directory"));
          LAB_file_name   = prm.get("LAB depth filename");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for replace lithosphere viscosity, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }

    template <int dim>
    bool
	ReplaceLithosphereViscosity<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }

    template <int dim>
       double
	   ReplaceLithosphereViscosity<dim>::
       reference_viscosity() const
       {
         /* Return reference viscosity from base model*/
         return base_model->reference_viscosity();
       }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ReplaceLithosphereViscosity,
                                   "replace lithosphere viscosity",
                                   "The ``depth dependent'' Material model applies a depth-dependent scaling "
                                   "to any of the other available material models. In other words, it "
                                   "is a ``compositing material model''."
                                   "\n\n"
                                   "Parameters related to the depth dependent model are read from a subsection "
                                   "``Material model/Depth dependent model''. "
                                   "The user must specify a ``Base model'' from which material properties are "
                                   "derived. Currently the depth dependent model only allows depth dependence of "
                                   "viscosity - other material properties are taken from the ``Base model''. "
                                   "Viscosity $\\eta$ at depth $z$ is calculated according to:"
                                   "\\begin{equation}"
                                   "\\eta(z,p,T,X,...) = \\eta(z) \\eta_b(p,T,X,..)/\\eta_{rb}"
                                   "\\end{equation}"
                                   "where $\\eta(z)$ is the depth-dependence specified by the depth dependent "
                                   "model, $\\eta_b(p,T,X,...)$ is the viscosity calculated from the base model, "
                                   "and $\\eta_{rb}$ is the reference viscosity of the ``Base model''. "
                                   "In addition to the specification of the ``Base model'', the user must specify "
                                   "the method to be used to calculate the depth-dependent viscosity $\\eta(z)$ as "
                                   "``Material model/Depth dependent model/Depth dependence method'', which can be "
                                   "chosen among ``None|Function|File|List''. Each method and the associated parameters "
                                   "are as follows:"
                                   "\n"
                                   "\n"
                                   "``Function'': read a user-specified parsed function from the input file in a "
                                   "subsection ``Material model/Depth dependent model/Viscosity depth function''. "
                                   "By default, this function is uniformly equal to 1.0e21. Specifying a function "
                                   "that returns a value less than or equal to 0.0 anywhere in the model domain will "
                                   "produce an error. "
                                   "\n"
                                   "\n"
                                   "``File'': read a user-specified file containing viscosity values at specified "
                                   "depths. The file containing depth-dependent viscosities is read from a "
                                   "directory specified by the user as "
                                   "``Material model/Depth dependent model/Data directory'', from a file with name "
                                   "specified as ``Material model/Depth dependent model/Viscosity depth file''. "
                                   "The format of this file is ascii text and contains two columns with one header line:"
                                   "\n"
                                   "\n"
                                   "example Viscosity depth file:\\\\"
                                   "Depth (m)    Viscosity (Pa-s)\\\\"
                                   "0.0000000e+00     1.0000000e+21\\\\"
                                   "6.7000000e+05     1.0000000e+22\\\\"
                                   "\n"
                                   "\n"
                                   "Viscosity is interpolated from this file using linear interpolation. "
                                   "``None'': no depth-dependence. Viscosity is taken directly from ``Base model''"
                                   "\n"
                                   "\n"
                                   "``List:'': read a comma-separated list of depth values corresponding to the maximum "
                                   "depths of layers having constant depth-dependence $\\eta(z)$. The layers must be "
                                   "specified in order of increasing depth, and the last layer in the list must have a depth "
                                   "greater than or equal to the maximal depth of the model. The list of layer depths is "
                                   "specified as ``Material model/Depth dependent model/Depth list'' and the corresponding "
                                   "list of layer viscosities is specified as "
                                   "``Material model/Depth dependent model/Viscosity list''"
                                  )
  }
}
