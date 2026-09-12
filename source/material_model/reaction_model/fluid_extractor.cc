/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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


#include <aspect/material_model/reaction_model/fluid_extractor.h>
#include <aspect/material_model/reaction_model/katz2003_mantle_melting.h>
#include <aspect/material_model/reaction_model/tian2019_solubility.h>
#include <aspect/material_model/reactive_fluid_transport.h>
#include <aspect/utilities.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      template <int dim>
      void
      FluidExtractor<dim>::set_parameters (const std::string &scheme_name,
                                           const double depth,
                                           const std::string &method)
      {
        reaction_scheme_name = scheme_name;
        extraction_depth = depth;
        extraction_method = method;

        if (reaction_scheme_name == "katz2003")
          if (this->include_melt_transport())
            AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                        ExcMessage("The fluid extractor reaction model, with the Katz2003 reaction model and "
                                   "melt transport enabled can only be used if there is a compositional field "
                                   "named 'peridotite'."));
      }



      template <int dim>
      void
      FluidExtractor<dim>::calculate_reaction_rate_outputs (const typename Interface<dim>::MaterialModelInputs  &in,
                                                            typename Interface<dim>::MaterialModelOutputs       &out) const
      {
        std::shared_ptr<ReactionRateOutputs<dim>> reaction_rate_out
          = out.template get_additional_output_object<ReactionRateOutputs<dim>>();

        // Depending on what the parent reaction model is, we need to determine which compositional field is used to
        // react with the fluid. If the reaction scheme is katz2003, then this compositional field is "peridotite",
        // and if the reaction scheme is tian2019, then this compositional field is "bound_fluid". A caveat for the
        // katz2003 reaction scheme is that the porosity will only react with the solid if melt transport is enabled,
        // otherwise no "peridotite" compositional field is required. Handle this edge case by making the reaction
        // index negative and catching this farther down.
        int reaction_index = std::numeric_limits<unsigned int>::max();
        if (reaction_scheme_name == "katz2003")
          reaction_index = this->include_melt_transport() ? this->introspection().compositional_index_for_name("peridotite") : -1;
        else if (reaction_scheme_name == "tian approximation")
          reaction_index = this->introspection().compositional_index_for_name("bound_fluid");
        else
          AssertThrow(false, ExcMessage("The reaction scheme " + reaction_scheme_name + " is not recognized as a valid reaction scheme for extracting fluids."));

        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          {
            if (reaction_rate_out != nullptr && in.requests_property(MaterialProperties::reaction_rates))
              {
                const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");
                const double depth = this->get_geometry_model().depth(in.position[i]);
                const bool above_extraction_depth = depth < extraction_depth;

                if (above_extraction_depth)
                  {
                    const double porosity = in.composition[i][porosity_index];
                    double porosity_change = 0.0;

                    if (extraction_method == "linear")
                      porosity_change = -porosity * (extraction_depth - depth) / extraction_depth;
                    else if (extraction_method == "constant")
                      porosity_change = -porosity;
                    else
                      AssertThrow(false, ExcMessage("The extraction method " + extraction_method + " is not recognized as a valid extraction method for extracting fluids."));

                    // Prevent negative porosity from developing
                    porosity_change = std::max(-porosity, porosity_change);
                    reaction_rate_out->reaction_rates[i][porosity_index] = porosity_change / (100.0 * year_in_seconds);

                    // Set the reaction rates for the compositional field that reacts with the fluid to be 0. This is to prevent the "porosity" from being extracted while
                    // simultaneously reacting with the solid phase. Do not override the reaction rate if it is negative, because we still want to allow the solid phase
                    // to create porosity, which will then be extracted since it is above the extraction depth.
                    if (reaction_index >= 0)
                      if (reaction_rate_out->reaction_rates[i][reaction_index] > 0)
                        reaction_rate_out->reaction_rates[i][reaction_index] -= reaction_rate_out->reaction_rates[i][reaction_index];
                  }
              }
          }
      }



      template <int dim>
      void
      FluidExtractor<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Melt extraction depth", "1000.0",
                           Patterns::Double (),
                           "Depth below the surface at which melt is extracted."
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Extraction method", "linear",
                           Patterns::Selection("linear|constant"),
                           "Method for extracting melt. Options are 'linear' or 'constant'. "
                           "All methods are proportional to the volume fraction of the porosity "
                           "above the extraction depth. 'linear' increases the extraction rate "
                           "from a minimum at the extraction depth to a maximum at the model surface, "
                           "'constant' is equal to the volume fraction of the porosity above the "
                           "extraction depth.");
      }



      template <int dim>
      void
      FluidExtractor<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("The reaction model 'melt extractor' "
                               "can only be used if there is a compositional field named "
                               "'porosity'."));

        extraction_depth  = prm.get_double ("Melt extraction depth");
        extraction_method = prm.get ("Extraction method");
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
#define INSTANTIATE(dim) \
  template class FluidExtractor<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)
#undef INSTANTIATE
    }
  }
}
