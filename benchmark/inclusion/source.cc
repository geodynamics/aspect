#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class NoAdiabaticHeating :  public MaterialModel::Simple<dim>
    {
      public:
        virtual double thermal_expansion_coefficient (const double,
                                                      const double,
                                                      const std::vector<double> &,
                                                      const Point<dim> &) const
	{
	  std::cout << "hi" << std::endl;
	  
	  return 0;
	}
    };

  }
}








// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(NoAdiabaticHeating,
                                   "local_Inclusion",
                                   "As described in the .prm file.")
  }
}
