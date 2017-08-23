#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/drucker_prager.h>
#include <aspect/simulator_access.h>
#include <aspect/newton.h>

#include <iostream>

int f(double parameter)
{

  std::cout << std::endl << "Test for p = " << parameter << std::endl;

  const int dim=2;
  using namespace aspect::MaterialModel;
  MaterialModelInputs<dim> in_base(5,3);
  in_base.composition[0][0] = 0;
  in_base.composition[0][1] = 0;
  in_base.composition[0][2] = 0;
  in_base.composition[1][0] = 0.75;
  in_base.composition[1][1] = 0.15;
  in_base.composition[1][2] = 0.10;
  in_base.composition[2][0] = 0;
  in_base.composition[2][1] = 0.2;
  in_base.composition[2][2] = 0.4;
  in_base.composition[3][0] = 0;
  in_base.composition[3][1] = 0.2;
  in_base.composition[3][2] = 0.4;
  in_base.composition[4][0] = 1;
  in_base.composition[4][1] = 0;
  in_base.composition[4][2] = 0;

  in_base.pressure[0] = 1e9;
  in_base.pressure[1] = 5e9;
  in_base.pressure[2] = 2e10;
  in_base.pressure[3] = 2e11;
  in_base.pressure[4] = 2e12;

  /**
   * We can't take to small strain-rates, because then the difference in the
   * visocisty will be too small for the double accuracy which stores
   * the visocity solutions and the finite diference solution.
   */
  in_base.strain_rate[0] = SymmetricTensor<2,dim>();
  in_base.strain_rate[0][0][0] = 1e-12;
  in_base.strain_rate[0][0][1] = 1e-12;
  in_base.strain_rate[0][1][1] = 1e-11;
  in_base.strain_rate[1] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[1][0][0] = -1.71266e-13;
  in_base.strain_rate[1][0][1] = -5.82647e-12;
  in_base.strain_rate[1][1][1] = 4.21668e-14;
  in_base.strain_rate[2] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[2][1][1] = 1e-13;
  in_base.strain_rate[2][0][1] = 1e-11;
  in_base.strain_rate[2][0][0] = -1e-12;
  in_base.strain_rate[3] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[3][1][1] = 4.9e-21;
  in_base.strain_rate[3][0][1] = 4.9e-21;
  in_base.strain_rate[3][0][0] = 4.9e-21;
  in_base.strain_rate[4] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
  in_base.strain_rate[4][1][1] = 1e-11;
  in_base.strain_rate[4][0][1] = 1e-11;
  in_base.strain_rate[4][0][0] = 1e-11;

  in_base.temperature[0] = 293;
  in_base.temperature[1] = 1600;
  in_base.temperature[2] = 2000;
  in_base.temperature[3] = 2100;
  in_base.temperature[4] = 2200;

  SymmetricTensor<2,dim> zerozero = SymmetricTensor<2,dim>();
  SymmetricTensor<2,dim> onezero = SymmetricTensor<2,dim>();
  SymmetricTensor<2,dim> oneone = SymmetricTensor<2,dim>();

  zerozero[0][0] = 1;
  onezero[1][0]  = 0.5; // because symmetry doubles this entry
  oneone[1][1]   = 1;

  double finite_difference_accuracy = 1e-7;
  double finite_difference_factor = 1+finite_difference_accuracy;

  bool Error = false;

  MaterialModelInputs<dim> in_dviscositydpressure(in_base);
  in_dviscositydpressure.pressure[0] *= finite_difference_factor;
  in_dviscositydpressure.pressure[1] *= finite_difference_factor;
  in_dviscositydpressure.pressure[2] *= finite_difference_factor;
  in_dviscositydpressure.pressure[3] *= finite_difference_factor;
  in_dviscositydpressure.pressure[4] *= finite_difference_factor;

  MaterialModelInputs<dim> in_dviscositydstrainrate_zerozero(in_base);
  MaterialModelInputs<dim> in_dviscositydstrainrate_onezero(in_base);
  MaterialModelInputs<dim> in_dviscositydstrainrate_oneone(in_base);

  in_dviscositydstrainrate_zerozero.strain_rate[0] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[0][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[1] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[1][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[2] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[2][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[3] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[3][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_zerozero.strain_rate[4] += std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[4][0][0]) * finite_difference_accuracy * zerozero;
  in_dviscositydstrainrate_onezero.strain_rate[0]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[0][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[1]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[1][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[2]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[2][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[3]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[3][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_onezero.strain_rate[4]  += std::fabs(in_dviscositydstrainrate_onezero.strain_rate[4][1][0]) * finite_difference_accuracy * onezero;
  in_dviscositydstrainrate_oneone.strain_rate[0]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[0][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[1]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[1][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[2]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[2][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[3]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[3][1][1]) * finite_difference_accuracy * oneone;
  in_dviscositydstrainrate_oneone.strain_rate[4]   += std::fabs(in_dviscositydstrainrate_oneone.strain_rate[4][1][1]) * finite_difference_accuracy * oneone;

  MaterialModelInputs<dim> in_dviscositydtemperature(in_base);
  in_dviscositydtemperature.temperature[0] *= 1.0000000001;
  in_dviscositydtemperature.temperature[1] *= 1.0000000001;
  in_dviscositydtemperature.temperature[2] *= 1.0000000001;
  in_dviscositydtemperature.temperature[3] *= 1.0000000001;
  in_dviscositydtemperature.temperature[4] *= 1.0000000001;


  MaterialModelOutputs<dim> out_base(5,3);

  MaterialModelOutputs<dim> out_dviscositydpressure(5,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_zerozero(5,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_onezero(5,3);
  MaterialModelOutputs<dim> out_dviscositydstrainrate_oneone(5,3);
  MaterialModelOutputs<dim> out_dviscositydtemperature(5,3);

  if (out_base.get_additional_output<MaterialModelDerivatives<dim> >() != NULL)
    throw "error";

  out_base.additional_outputs.push_back(std::make_shared<MaterialModelDerivatives<dim> > (5));

  DruckerPrager<dim> mat;
  ParameterHandler prm;
  mat.declare_parameters(prm);

  prm.enter_subsection("Material model");
  {
    prm.enter_subsection ("Drucker Prager");
    {
      prm.enter_subsection ("Viscosity");
      {
        prm.set ("Reference strain rate", "1e-20");
        prm.set ("Angle of internal friction", "30");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }
  prm.leave_subsection();

  mat.parse_parameters(prm);

  mat.evaluate(in_base, out_base);
  mat.evaluate(in_dviscositydpressure, out_dviscositydpressure);
  mat.evaluate(in_dviscositydstrainrate_zerozero, out_dviscositydstrainrate_zerozero);
  mat.evaluate(in_dviscositydstrainrate_onezero, out_dviscositydstrainrate_onezero);
  mat.evaluate(in_dviscositydstrainrate_oneone, out_dviscositydstrainrate_oneone);
  mat.evaluate(in_dviscositydtemperature, out_dviscositydtemperature);

  // set up additional output for the derivatives
  MaterialModelDerivatives<dim> *derivatives;
  derivatives = out_base.get_additional_output<MaterialModelDerivatives<dim> >();

  double temp;
  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = (out_dviscositydpressure.viscosities[i] - out_base.viscosities[i]);
      if (in_base.pressure[i] != 0)
        {
          temp /= (in_base.pressure[i] * finite_difference_accuracy);
        }
      std::cout << "pressure " << i << ": Finite difference = " << temp << ". Analytical derivative = " << derivatives->viscosity_derivative_wrt_pressure[i]  << std::endl;
      if (std::fabs(temp - derivatives->viscosity_derivative_wrt_pressure[i]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_pressure[i])))
        // if (temp > derivatives->viscosity_derivative_wrt_pressure[i] * finite_difference_factor || temp < derivatives->viscosity_derivative_wrt_pressure[i] * (2-finite_difference_factor))
        {
          std::cout << "   Error: The derivative of the viscosity to the pressure is too different from the analitical value." << std::endl;
          Error = true;
        }

    }

  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = out_dviscositydstrainrate_zerozero.viscosities[i] - out_base.viscosities[i];
      if (temp != 0)
        {
          temp /= std::fabs(in_dviscositydstrainrate_zerozero.strain_rate[i][0][0]) * finite_difference_accuracy;
        }
      std::cout << "zerozero " << i << ": Finite difference = " << temp << ". Analytical derivative = " << derivatives->viscosity_derivative_wrt_strain_rate[i][0][0]  << std::endl;
      if (std::fabs(temp - derivatives->viscosity_derivative_wrt_strain_rate[i][0][0]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_strain_rate[i][0][0])))
        {
          std::cout << "   Error: The derivative of the viscosity to the strain rate is too different from the analitical value." << std::endl;
          Error = true;
        }



    }

  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = out_dviscositydstrainrate_onezero.viscosities[i] - out_base.viscosities[i];
      if (temp != 0)
        {
          temp /= std::fabs(in_dviscositydstrainrate_onezero.strain_rate[i][1][0]) * finite_difference_accuracy;
        }
      std::cout << "onezero " << i << ": Finite difference = " << temp << ". Analytical derivative = " << derivatives->viscosity_derivative_wrt_strain_rate[i][1][0]   << std::endl;
      if (std::fabs(temp - derivatives->viscosity_derivative_wrt_strain_rate[i][1][0]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_strain_rate[i][1][0])) )
        {
          std::cout << "   Error: The derivative of the viscosity to the strain rate is too different from the analitical value." << std::endl;
          Error = true;
        }
    }

  for (unsigned int i = 0; i < 5; i++)
    {
      // prevent division by zero. If it is zero, the test has passed, because or
      // the finite difference and the analytical result match perfectly, or (more
      // likely) the material model in independent of this variable.
      temp = out_dviscositydstrainrate_oneone.viscosities[i] - out_base.viscosities[i];
      if (temp != 0)
        {
          temp /= std::fabs(in_dviscositydstrainrate_oneone.strain_rate[i][1][1]) * finite_difference_accuracy;
        }
      std::cout << "oneone " << i << ": Finite difference = " << temp << ". Analytical derivative = " << derivatives->viscosity_derivative_wrt_strain_rate[i][1][1]  << std::endl;
      if (std::fabs(temp - derivatives->viscosity_derivative_wrt_strain_rate[i][1][1]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_strain_rate[i][1][1])) )
        {
          std::cout << "   Error: The derivative of the viscosity to the strain rate is too different from the analitical value." << std::endl;
          Error = true;
        }

    }

  if (Error)
    {
      std::cout << "Some parts of the test where not succesful." << std::endl;
    }
  else
    {
      std::cout << "OK" << std::endl;
    }

  return 42;
}

int exit_function()
{
  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int ii = f(-1000); // Testing min function
int iz = f(-2); // Testing generalized p norm mean with negative p
int ij = f(-1.5); // Testing generalized p norm mean with negative, non int p
int ik = f(-1); // Testing harmonic mean
int ji = f(0); // Testing geometric mean
int jj = f(1); // Testing arithmetic mean
int jk = f(2); // Testing generalized p norm mean with positive p
int kj = f(1000); // Testing max function
int kl = exit_function();


