/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_elastic_h
#define _aspect_elastic_h

#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/assembly.h>
#include <aspect/material_model/interface.h>

namespace aspect
{
  using namespace dealii;

  namespace MaterialModel
  {
    template <int dim>
    class ElasticOutputs : public AdditionalMaterialOutputs<dim>
    {
      public:
        ElasticOutputs (const unsigned int n_points,
                        const unsigned int /*n_comp*/)
        {
          elastic_viscosities.resize(n_points);
          elastic_evolutions.resize(n_points);
        }

        /**
         * Viscoelasticity value $\eta_ve$ at the given positions.
         * Units: Pa s
         */
        std::vector<double> elastic_viscosities;

        /**
         * Elastic stress evolution values $\chi_t$ at the given positions.
         * Units: none
         */
        std::vector<double> elastic_evolutions;


        /**
         * Do the requested averaging operation for the elastic outputs.
         * The projection matrix argument is only used if the operation
         * chosen is project_to_Q1.
         */
        /** void average (const MaterialAveraging::AveragingOperation operation,
         *              const FullMatrix<double>  &projection_matrix,
         *             const FullMatrix<double>  &expansion_matrix);
         */
    };

  }

}

#endif
