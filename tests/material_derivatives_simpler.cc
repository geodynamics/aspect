// make sure we use the additional material model outputs
#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

#include <iostream>
#include <typeinfo>


int f()
{
  const int dim=2;

  using namespace aspect::MaterialModel;
  MaterialModelInputs<dim> in_base(3,3);
  in_base.composition[0][0] = 1;
  in_base.composition[0][1] = 0;
  in_base.composition[0][2] = 0;
  in_base.composition[1][0] = 0.5;
  in_base.composition[1][1] = 0.25;
  in_base.composition[1][2] = 0.25;
  in_base.composition[2][0] = 0;
  in_base.composition[2][1] = 0.5;
  in_base.composition[2][2] = 0.5;

  in_base.pressure[0] = 10;
  in_base.pressure[1] = 1000;
  in_base.pressure[2] = 1e6;

  in_base.strain_rate[0] = SymmetricTensor<2,dim>();
  in_base.strain_rate[1][0][0] = 1e-12;
  in_base.strain_rate[1][0][1] = 2e-12;
  in_base.strain_rate[1][1][0] = 5e-13;
  in_base.strain_rate[1][1][1] = 5e-12;
  in_base.strain_rate[1] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[1][0][1] = 1e-11;
  in_base.strain_rate[2] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[2][1][1] = 1e-15;

  in_base.temperature[0] = 293;
  in_base.temperature[1] = 1000;
  in_base.temperature[2] = 2000;

  MaterialModelInputs<dim> in_dviscositydpressure(3,3);
  in_dviscositydpressure = in_base;
  in_dviscositydpressure.pressure[0] *= 1.0000000001;
  in_dviscositydpressure.pressure[1] *= 1.0000000001;
  in_dviscositydpressure.pressure[2] *= 1.0000000001;
  MaterialModelInputs<dim> in_dviscositydstrainrate_0(3,3);
  MaterialModelInputs<dim> in_dviscositydstrainrate_1(3,3);
  MaterialModelInputs<dim> in_dviscositydstrainrate_2(3,3);
  in_dviscositydstrainrate_0 = in_base;
  in_dviscositydstrainrate_1 = in_base;
  in_dviscositydstrainrate_2 = in_base;
  in_dviscositydstrainrate_0.strain_rate[0][0][0] *= 1.0000000001;
  in_dviscositydstrainrate_0.strain_rate[1][0][0] *= 1.0000000001;
  in_dviscositydstrainrate_0.strain_rate[2][0][0] *= 1.0000000001;
  in_dviscositydstrainrate_1.strain_rate[0][1][0] *= 1.0000000001;
  in_dviscositydstrainrate_1.strain_rate[1][1][0] *= 1.0000000001;
  in_dviscositydstrainrate_1.strain_rate[2][1][0] *= 1.0000000001;
  in_dviscositydstrainrate_2.strain_rate[0][1][1] *= 1.0000000001;
  in_dviscositydstrainrate_2.strain_rate[1][1][1] *= 1.0000000001;
  in_dviscositydstrainrate_2.strain_rate[2][1][1] *= 1.0000000001;

  MaterialModelInputs<dim> in_dviscositydtemperature(3,3);
  in_dviscositydtemperature = in_base;
  in_dviscositydtemperature.temperature[0] *= 1.0000000001;
  in_dviscositydtemperature.temperature[1] *= 1.0000000001;
  in_dviscositydtemperature.temperature[2] *= 1.0000000001;

  MaterialModelOutputs<dim> out_base(3,3);

  MaterialModelOutputs<dim> out_dviscositydpressure(3,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_0(3,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_1(3,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_2(3,3);
  MaterialModelOutputs<dim> out_dviscositydtemperature(3,3);

  if (out_base.get_additional_output<MaterialModelDerivatives<dim> >() != NULL)
    throw "error";

  out_base.additional_outputs.push_back(std::make_shared<MaterialModelDerivatives<dim> > (3, 3));

  Simple<dim> mat;
  ParameterHandler prm;
  mat.declare_parameters(prm);
  mat.parse_parameters(prm);
  mat.evaluate(in_base, out_base);

  mat.evaluate(in_dviscositydpressure, out_dviscositydpressure);
  mat.evaluate(in_dviscositydstrainrate_0, out_dviscositydstrainrate_0);
  mat.evaluate(in_dviscositydstrainrate_1, out_dviscositydstrainrate_1);
  mat.evaluate(in_dviscositydstrainrate_2, out_dviscositydstrainrate_2);
  mat.evaluate(in_dviscositydtemperature, out_dviscositydtemperature);

  //set up additional output for the derivatives
  MaterialModelDerivatives<dim> *derivatives;
  derivatives = out_base.get_additional_output<MaterialModelDerivatives<dim> >();

  double temp;
  for(unsigned int i = 0; i < 3; i++)
  {
	  // prevent division by zero. If it is zero, the test has passed, because or
	  // the finite difference and the analytical result match perfectly, or (more
	  // likely) the material model in independent of this variable.
      temp = (out_dviscositydpressure.viscosities[i] - out_base.viscosities[i]);
	  std::cout << "out_base.viscosities[" << i << "] = " << out_base.viscosities[i] << ", out_dviscositydpressure.viscosity[" << i << "] = " << out_dviscositydpressure.viscosities[i] << ". Difference = " << temp << std::endl;
      if(temp != 0)
      {
    	  temp /= (in_base.pressure[i] * 1.0000000001);
    	  if(temp == 0 || temp / (in_base.pressure[i] * 1.0000000001) > derivatives->dviscosities_dpressure[i] * 1.0000000001 || temp < derivatives->dviscosities_dpressure[i] * 0.9999999999)
    		 throw "error: The derivative of the viscosity to the pressure is too different from the analitical value.";
      }

  }


  for(unsigned int i = 0; i < 3; i++)
  {
	  // prevent division by zero. If it is zero, the test has passed, because or
	  // the finite difference and the analytical result match perfectly, or (more
	  // likely) the material model in independent of this variable.
	  temp = out_dviscositydstrainrate_0.viscosities[i] - out_base.viscosities[i];
	  std::cout << "out_base.viscosities[" << i << "] = " << out_base.viscosities[i] << ", out_dviscositydstrainrate_0.viscosity[" << i << "] = " << out_dviscositydstrainrate_0.viscosities[i] << ". Difference = " << temp << std::endl;
	  if(temp != 0)
	  {
		  temp /= (std::sqrt(std::fabs(second_invariant(deviator(in_base.strain_rate[i])))) * 1.0000000001);
		  if(temp > std::sqrt(std::fabs(second_invariant(deviator(derivatives->dviscosities_dstrain_rate[i])))) * 1.0000000001 ||
			 temp < std::sqrt(std::fabs(second_invariant(deviator(derivatives->dviscosities_dstrain_rate[i])))) * 0.9999999999)
				throw "error: The derivative of the viscosity to the strain rate is too different from the analitical value.";
	  }
  }

  for(unsigned int i = 0; i < 3; i++)
  {
	  // prevent division by zero. If it is zero, the test has passed, because or
	  // the finite difference and the analytical result match perfectly, or (more
	  // likely) the material model in independent of this variable.
	  temp = out_dviscositydstrainrate_1.viscosities[i] - out_base.viscosities[i];
	  std::cout << "out_base.viscosities[" << i << "] = " << out_base.viscosities[i] << ", out_dviscositydstrainrate_1.viscosity[" << i << "] = " << out_dviscositydstrainrate_1.viscosities[i] << ". Difference = " << temp << std::endl;
	  if(temp != 0)
	  {
		  temp /= (std::sqrt(std::fabs(second_invariant(deviator(in_base.strain_rate[i])))) * 1.0000000001);
		  if(temp > std::sqrt(std::fabs(second_invariant(deviator(derivatives->dviscosities_dstrain_rate[i])))) * 1.0000000001 ||
			 temp < std::sqrt(std::fabs(second_invariant(deviator(derivatives->dviscosities_dstrain_rate[i])))) * 0.9999999999)
				throw "error: The derivative of the viscosity to the strain rate is too different from the analitical value.";
	  }
  }

  for(unsigned int i = 0; i < 3; i++)
  {
	  // prevent division by zero. If it is zero, the test has passed, because or
	  // the finite difference and the analytical result match perfectly, or (more
	  // likely) the material model in independent of this variable.
	  temp = out_dviscositydstrainrate_2.viscosities[i] - out_base.viscosities[i];
	  std::cout << "out_base.viscosities[" << i << "] = " << out_base.viscosities[i] << ", out_dviscositydstrainrate_2.viscosity[" << i << "] = " << out_dviscositydstrainrate_2.viscosities[i] << ". Difference = " << temp << std::endl;
	  	  if(temp != 0)
	  {
		  temp /= (std::sqrt(std::fabs(second_invariant(deviator(in_base.strain_rate[i])))) * 1.0000000001);
		  std::cout << "out_base.viscosities[" << i << "] = " << out_base.viscosities[i] << ", out_dviscositydstrainrate_2.viscosity[" << i << "] = " << out_dviscositydstrainrate_2.viscosities[i] << ". Difference = " << temp << std::endl;
		  if(temp > std::sqrt(std::fabs(second_invariant(deviator(derivatives->dviscosities_dstrain_rate[i])))) * 1.0000000001 ||
			 temp < std::sqrt(std::fabs(second_invariant(deviator(derivatives->dviscosities_dstrain_rate[i])))) * 0.9999999999)
				throw "error: The derivative of the viscosity to the strain rate is too different from the analitical value.";
	  }
  }
  for(unsigned int i = 0; i < 3; i++)
  {
	  // prevent division by zero. If it is zero, the test has passed, because or
	  // the finite difference and the analytical result match perfectly, or (more
	  // likely) the material model in independent of this variable.
	  temp = out_dviscositydtemperature.viscosities[i] - out_base.viscosities[i];
	  std::cout << "out_base.viscosities[" << i << "] = " << out_base.viscosities[i] << ", out_dviscositydtemperature.viscosity[" << i << "] = " << out_dviscositydtemperature.viscosities[i] << ". Difference = " << temp << std::endl;
	  if(temp != 0)
	  {
		  temp /= in_base.temperature[i] * 1.0000000001;
		  std::cout << "out_base.temperature[" << i << "] = " << out_base.viscosities[i] << ", out_ddensitydpressure.temperature[" << i << "] = " << out_dviscositydtemperature.viscosities[i] << ". Difference = " << temp << std::endl;
		  if(temp > derivatives->dviscosities_dtemperature[i] * 1.0000000001 ||
			 temp < derivatives->dviscosities_dtemperature[i] * 0.9999999999)
				throw "error: The derivative of the viscosity to the strain rate is too different from the analitical value.";
	  }
  }
  std::cout << "OK" << std::endl;
  exit(0);
  return 42;
}

// run this function by initializing a global variable by it
int i = f();
