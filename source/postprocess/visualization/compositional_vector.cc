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

#include <aspect/postprocess/visualization/compositional_vector.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      CompositionalVector<dim>::
      CompositionalVector ()
        :
        DataPostprocessor<dim> ()
      {}


      template <int dim>
      std::vector<std::string>
      CompositionalVector<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;
        for (unsigned int i=0; i<vector_names.size(); ++i)
          for (unsigned int j=0; j<dim; ++j)
            solution_names.push_back(vector_names[i]);
        return solution_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      CompositionalVector<dim>::
      get_data_component_interpretation () const
      {
        return std::vector<DataComponentInterpretation::DataComponentInterpretation>
               (vector_names.size()*dim,
                DataComponentInterpretation::component_is_part_of_vector);
      }


      template <int dim>
      UpdateFlags
      CompositionalVector<dim>::
      get_needed_update_flags () const
      {
        return update_values;
      }


      template <int dim>
      void
      CompositionalVector<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> > &,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError ());
        Assert (uh[0].size() == this->introspection().n_components, ExcInternalError ());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            for (unsigned int i=0; i<vector_names.size(); ++i)
              for (unsigned int j=0; j<dim; ++j)
                computed_quantities[q][i*dim+j] =
                  uh[q][this->introspection().component_indices.compositional_fields[sets[i][j]]];
          }
      }


      template <int dim>
      void
      CompositionalVector<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Compositional fields as vectors");
            {
              prm.declare_entry("Names of vectors", "",
                                Patterns::List (Patterns::Anything()),
                                "Names of vectors as they will appear in the output.");

              prm.declare_entry("Names of fields", "",
                                Patterns::Anything (),
                                "A list of sets of compositional fields which should be output "
                                "as vectors. Sets are separated from each other by semicolons "
                                "and vector components within each set are separated by commas "
                                "(e.g. $vec1_x$, $vec1_y$ ; $vec2_x$, $vec2_y$) where each name must be "
                                "a defined named compositional field. If only one name is given "
                                "in a set, it is interpreted as the first in a sequence of dim "
                                "consecutive compositional fields.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      CompositionalVector<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Compositional fields as vectors");
            {
              vector_names = Utilities::split_string_list(prm.get("Names of vectors"), ',');

              const std::vector<std::string> sets_list = Utilities::split_string_list(prm.get("Names of fields"), ';');
              for (unsigned int i=0; i<sets_list.size(); ++i)
                {
                  const std::vector<std::string> set = Utilities::split_string_list(sets_list[i], ',');
                  AssertThrow((set.size() == dim || set.size() == 1),
                              ExcMessage("Sets of compositional fields to be output as vectors must have dim components, "
                                         "or one component (i.e. the first field in a sequence of dim consecutive fields)."));

                  std::vector<unsigned int> set_idx;
                  for (unsigned int j=0; j<set.size(); ++j)
                    {
                      AssertThrow(this->introspection().compositional_name_exists (set[j]),
                                  ExcMessage("All fields in compositional vector sets must match names of compositional "
                                             "fields as assigned in the \"Compositional fields/Names of fields\" parameter."));
                      set_idx.push_back(this->introspection().compositional_index_for_name(set[j]));
                    }

                  while (set_idx.size() < dim)
                    {
                      AssertThrow(set_idx[0]+set_idx.size() < this->n_compositional_fields(),
                                  ExcMessage(set[0] + " must be the first in a sequence of dim compositional fields to be "
                                             "interpreted as a vector. You have too few compositional fields."));
                      set_idx.push_back(set_idx[0]+set_idx.size());
                    }
                  sets.push_back(set_idx);
                }

              AssertThrow(vector_names.size() == sets.size(),
                          ExcMessage("You must define names for each of the sets of compositional fields which will "
                                     "be output as vectors, and the length of that list must match the number of vectors."));
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(CompositionalVector,
                                                  "compositional vector",
                                                  "A visualization output object that outputs vectors whose "
                                                  "components are derived from compositional fields. Input "
                                                  "parameters for this postprocessor are defined in section "
                                                  "Postprocess/Visualization/Compositional fields as vectors")
    }
  }
}
