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

// A postprocessor that simply outputs which of the output variables
// of the material model depends on the input variables

#include <aspect/postprocess/interface.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class MaterialModelDependencies : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler & /*statistics*/) override
        {
          using namespace MaterialModel;
          using namespace NonlinearDependence;

          const aspect::MaterialModel::Interface<dim> &model
            = this->get_material_model();
          std::cout << "viscosity depends on"
                    << ((model.get_model_dependence().viscosity & MaterialModel::NonlinearDependence::pressure)
                        ?
                        " pressure" : "")
                    << ((model.get_model_dependence().viscosity & MaterialModel::NonlinearDependence::temperature)
                        ?
                        " temperature" : "")
                    << ((model.get_model_dependence().viscosity & MaterialModel::NonlinearDependence::strain_rate)
                        ?
                        " strainrate" : "")
                    << ((model.get_model_dependence().viscosity & MaterialModel::NonlinearDependence::compositional_fields)
                        ?
                        " composition" : "")
                    << std::endl;
          std::cout << "density depends on"
                    << ((model.get_model_dependence().density & MaterialModel::NonlinearDependence::pressure)
                        ?
                        " pressure" : "")
                    << ((model.get_model_dependence().density & MaterialModel::NonlinearDependence::temperature)
                        ?
                        " temperature" : "")
                    << ((model.get_model_dependence().density & MaterialModel::NonlinearDependence::strain_rate)
                        ?
                        " strainrate" : "")
                    << ((model.get_model_dependence().density & MaterialModel::NonlinearDependence::compositional_fields)
                        ?
                        " composition" : "")
                    << std::endl;
          std::cout << "compressibility depends on"
                    << ((model.get_model_dependence().compressibility & MaterialModel::NonlinearDependence::pressure)
                        ?
                        " pressure" : "")
                    << ((model.get_model_dependence().compressibility & MaterialModel::NonlinearDependence::temperature)
                        ?
                        " temperature" : "")
                    << ((model.get_model_dependence().compressibility & MaterialModel::NonlinearDependence::strain_rate)
                        ?
                        " strainrate" : "")
                    << ((model.get_model_dependence().compressibility & MaterialModel::NonlinearDependence::compositional_fields)
                        ?
                        " composition" : "")
                    << std::endl;
          std::cout << "specific_heat depends on"
                    << ((model.get_model_dependence().specific_heat & MaterialModel::NonlinearDependence::pressure)
                        ?
                        " pressure" : "")
                    << ((model.get_model_dependence().specific_heat & MaterialModel::NonlinearDependence::temperature)
                        ?
                        " temperature" : "")
                    << ((model.get_model_dependence().specific_heat & MaterialModel::NonlinearDependence::strain_rate)
                        ?
                        " strainrate" : "")
                    << ((model.get_model_dependence().specific_heat & MaterialModel::NonlinearDependence::compositional_fields)
                        ?
                        " composition" : "")
                    << std::endl;
          std::cout << "thermal_conductivity depends on"
                    << ((model.get_model_dependence().thermal_conductivity & MaterialModel::NonlinearDependence::pressure)
                        ?
                        " pressure" : "")
                    << ((model.get_model_dependence().thermal_conductivity & MaterialModel::NonlinearDependence::temperature)
                        ?
                        " temperature" : "")
                    << ((model.get_model_dependence().thermal_conductivity & MaterialModel::NonlinearDependence::strain_rate)
                        ?
                        " strainrate" : "")
                    << ((model.get_model_dependence().thermal_conductivity & MaterialModel::NonlinearDependence::compositional_fields)
                        ?
                        " composition" : "")
                    << std::endl;
          std::cout << "model is "
                    << (model.is_compressible() ? "compressible" : "incompressible")
                    << std::endl;
          return std::pair<std::string,std::string>();
        }
    };
  }
}


namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MaterialModelDependencies,
                                  "material model dependencies",
                                  ".")
  }
}
