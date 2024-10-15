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

#include <aspect/material_model/depth_dependent.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <boost/lexical_cast.hpp>

#include <limits>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    DepthDependent<dim>::initialize()
    {
      base_model->initialize();
    }



    template <int dim>
    void
    DepthDependent<dim>::update()
    {
      base_model->update();

      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        viscosity_function.set_time (this->get_time() / year_in_seconds);
      else
        viscosity_function.set_time (this->get_time());
    }



    template <int dim>
    bool
    DepthDependent<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }



    template <int dim>
    double
    DepthDependent<dim>::calculate_depth_dependent_prefactor(const double &depth) const
    {
      if (viscosity_source == file)
        {
          return depth_dependent_rheology->compute_viscosity(depth) / reference_viscosity;
        }
      else if (viscosity_source == function)
        {
          const Point<1> dpoint(depth);
          const double viscosity = viscosity_function.value(dpoint);
          Assert (viscosity > 0.0, ExcMessage("Viscosity depth function should be larger than zero"));
          return viscosity / reference_viscosity;
        }
      else if (viscosity_source == list)
        {
          const unsigned int nlayers = depth_values.size()-1;
          unsigned int i=0;
          /* find the layer containing the specified depth and return the corresponding viscosity */
          while (depth > depth_values[i] && i<nlayers)
            {
              /* increment i until depth is above base of layer i */
              ++i;
            }

          return viscosity_values[i] / reference_viscosity;
        }
      else if (viscosity_source == none)
        {
          return 1.0;
        }
      else
        {
          Assert(false, ExcMessage("Invalid method for viscosity depth dependence"));
          return 0.0;
        }
    }



    template <int dim>
    void
    DepthDependent<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                                  typename Interface<dim>::MaterialModelOutputs &out) const
    {
      base_model->evaluate(in,out);
      if (in.requests_property(MaterialProperties::viscosity))
        {
          // Scale the base model viscosity value by the depth dependent prefactor
          for (unsigned int i=0; i < out.n_evaluation_points(); ++i)
            {
              const double depth = this->get_geometry_model().depth(in.position[i]);
              out.viscosities[i] *= calculate_depth_dependent_prefactor(depth);
            }
        }
    }



    template <int dim>
    void
    DepthDependent<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        // Depth-dependent parameters from the rheology plugin
        Rheology::AsciiDepthProfile<dim>::declare_parameters(prm,
                                                             "Depth dependent model");

        prm.enter_subsection("Depth dependent model");
        {
          prm.declare_entry("Base model","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by a depth "
                            "dependent viscosity. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.declare_entry ("Depth dependence method", "None",
                             Patterns::Selection("Function|File|List|None"),
                             "Method that is used to specify how the viscosity should vary with depth.");
          prm.declare_entry("Depth list", "", Patterns::List(Patterns::Double ()),
                            "A comma-separated list of depth values for use with the ``List'' "
                            "``Depth dependence method''. The list must be provided in order of "
                            "increasing depth, and the last value must be greater than or equal to "
                            "the maximal depth of the model. The depth list is interpreted as a layered "
                            "viscosity structure and the depth values specify the maximum depths of each "
                            "layer.");
          prm.declare_entry("Viscosity list", "", Patterns::List(Patterns::Double ()),
                            "A comma-separated list of viscosity values, corresponding to the depth values "
                            "provided in ``Depth list''. The number of viscosity values specified here must "
                            "be the same as the number of depths provided in ``Depth list''.");

          prm.declare_entry ("Reference viscosity",
                             boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                             Patterns::Double (0.),
                             "The value of the constant reference viscosity $\\eta_r$ that is used to scale "
                             "the non-dimensional depth-dependent viscosity prefactor. "
                             "Units: \\si{\\pascal\\second}.");

          prm.declare_alias ("Data file name","Viscosity depth file");

          prm.enter_subsection("Viscosity depth function");
          {
            Functions::ParsedFunction<1>::declare_parameters(prm,1);
            prm.declare_entry("Function expression","1.0e21");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DepthDependent<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Depth dependent model");
        {
          AssertThrow( prm.get("Base model") != "depth dependent",
                       ExcMessage("You may not use ``depth dependent'' as the base model for "
                                  "a depth-dependent model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model = create_material_model<dim>(prm.get("Base model"));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          if (prm.get("Depth dependence method") == "Function")
            viscosity_source = function;
          else if (prm.get("Depth dependence method") == "File")
            viscosity_source = file;
          else if (prm.get("Depth dependence method") == "List")
            viscosity_source = list;
          else if (prm.get("Depth dependence method") == "None")
            viscosity_source = none;
          else
            {
              AssertThrow(false, ExcMessage("Unknown method for depth dependence."));
            }

          depth_values     = Utilities::string_to_double(Utilities::split_string_list(prm.get("Depth list")));
          viscosity_values = Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosity list")));
          /*
           * check sanity of viscosity list values and depth list values input
           */
          if (viscosity_source == list)
            {
              /* check that length of depth values and viscosity values are compatible */
              AssertThrow( depth_values.size() == viscosity_values.size() ,
                           ExcMessage("Depth list must be same size as Viscosity list"));
              /* check that list is in ascending order */
              for (unsigned int i=1; i<depth_values.size(); ++i)
                AssertThrow(depth_values[i] > depth_values[i-1],
                            ExcMessage("Viscosity depth values must be strictly ascending"));
              /* check that last layer includes base of model */
              AssertThrow( *(depth_values.end()-1) >= this->get_geometry_model().maximal_depth(),
                           ExcMessage("Last value in Depth list must be greater than or equal to maximal depth of domain"));
            }

          prm.enter_subsection("Viscosity depth function");
          {
            try
              {
                viscosity_function.parse_parameters(prm);
              }
            catch (...)
              {
                std::cerr << "FunctionParser failed to parse\n"
                          << "\t Viscosity depth function\n"
                          << "with expression \n"
                          << "\t' " << prm.get("Function expression") << "'";
                throw;
              }
          }
          prm.leave_subsection();

          reference_viscosity = prm.get_double("Reference viscosity");

          AssertThrow(reference_viscosity != std::numeric_limits<double>::max(),
                      ExcMessage("You have to set a reference viscosity for the depth dependent model."));
        }
        prm.leave_subsection();

        if (viscosity_source == file)
          {
            depth_dependent_rheology = std::make_unique<Rheology::AsciiDepthProfile<dim>>();
            depth_dependent_rheology->initialize_simulator (this->get_simulator());
            depth_dependent_rheology->parse_parameters(prm, "Depth dependent model");
            depth_dependent_rheology->initialize();
          }
      }
      prm.leave_subsection();

      /* After parsing the parameters for depth dependent, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }



    template <int dim>
    void
    DepthDependent<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      base_model->create_additional_named_outputs(out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(DepthDependent,
                                   "depth dependent",
                                   "The ``depth dependent'' Material model applies a depth-dependent scaling "
                                   "to the viscosity of any other available material models. In other words, it "
                                   "is a ``compositing material model''."
                                   "\n\n"
                                   "Parameters related to the depth dependent model are read from a subsection "
                                   "``Material model/Depth dependent model''. "
                                   "The user must specify a ``Base model'' from which material properties are "
                                   "derived. Currently the depth dependent model only allows depth dependence of "
                                   "viscosity - other material properties are taken from the ``Base model''. "
                                   "Viscosity $\\eta$ at depth $z$ is calculated according to:"
                                   "$ \\eta(z,p,T,X,...) = \\eta(z) \\eta_b(p,T,X,..)/\\eta_{r}$ "
                                   "where $\\eta(z)$ is the depth-dependence specified by the depth dependent "
                                   "model, $\\eta_b(p,T,X,...)$ is the viscosity calculated from the base model, "
                                   "and $\\eta_{r}$ is the reference viscosity. "
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
