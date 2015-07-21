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
       execute (TableHandler &statistics) 
       {
	 using namespace MaterialModel;
	 using namespace NonlinearDependence;

	 const aspect::MaterialModel::Interface<dim> &model
	   = this->get_material_model();
	 std::cout << "viscosity depends on "
		   << (model.viscosity_depends_on (pressure)
		       ?
		       "pressure " : "")
		   << (model.viscosity_depends_on (temperature)
		       ?
		       "temperature " : "")
		   << (model.viscosity_depends_on (strain_rate)
		       ?
		       "strainrate " : "")
		   << (model.viscosity_depends_on (compositional_fields)
		       ?
		       "composition " : "")
		   << std::endl;
	 std::cout << "density depends on "
		   << (model.density_depends_on (pressure)
		       ?
		       "pressure " : "")
		   << (model.density_depends_on (temperature)
		       ?
		       "temperature " : "")
		   << (model.density_depends_on (strain_rate)
		       ?
		       "strainrate " : "")
		   << (model.density_depends_on (compositional_fields)
		       ?
		       "composition " : "")
		   << std::endl;
	 std::cout << "compressibility depends on "
		   << (model.compressibility_depends_on (pressure)
		       ?
		       "pressure " : "")
		   << (model.compressibility_depends_on (temperature)
		       ?
		       "temperature " : "")
		   << (model.compressibility_depends_on (strain_rate)
		       ?
		       "strainrate " : "")
		   << (model.compressibility_depends_on (compositional_fields)
		       ?
		       "composition " : "")
		   << std::endl;
	 std::cout << "specific_heat depends on "
		   << (model.specific_heat_depends_on (pressure)
		       ?
		       "pressure " : "")
		   << (model.specific_heat_depends_on (temperature)
		       ?
		       "temperature " : "")
		   << (model.specific_heat_depends_on (strain_rate)
		       ?
		       "strainrate " : "")
		   << (model.specific_heat_depends_on (compositional_fields)
		       ?
		       "composition " : "")
		   << std::endl;
	 std::cout << "thermal_conductivity depends on "
		   << (model.thermal_conductivity_depends_on (pressure)
		       ?
		       "pressure " : "")
		   << (model.thermal_conductivity_depends_on (temperature)
		       ?
		       "temperature " : "")
		   << (model.thermal_conductivity_depends_on (strain_rate)
		       ?
		       "strainrate " : "")
		   << (model.thermal_conductivity_depends_on (compositional_fields)
		       ?
		       "composition " : "")
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

