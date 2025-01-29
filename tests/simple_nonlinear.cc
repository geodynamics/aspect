/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <deal.II/grid/tria.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/newton.h>

#include <iostream>

#include <aspect/utilities.h>

#include "../benchmarks/newton_solver_benchmark_set/nonlinear_channel_flow/simple_nonlinear.cc"

namespace aspect
{
  template <int dim>
  int f(double parameter)
  {

    std::cout << std::endl << "Test for p = " << parameter << " with dimension " << dim << std::endl;

    using namespace aspect::MaterialModel;

    // first set all material model values
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

    in_base.temperature[0] = 293;
    in_base.temperature[1] = 1600;
    in_base.temperature[2] = 2000;
    in_base.temperature[3] = 2100;
    in_base.temperature[4] = 2200;

    in_base.pressure[0] = 1e9;
    in_base.pressure[1] = 5e9;
    in_base.pressure[2] = 2e10;
    in_base.pressure[3] = 2e11;
    in_base.pressure[4] = 2e12;

    /**
     * We can't take to small strain-rates, because then the difference in the
     * viscosity will be too small for the double accuracy which stores
     * the viscosity solutions and the finite difference solution.
     */
    in_base.strain_rate[0] = SymmetricTensor<2,dim>();
    in_base.strain_rate[0][0][0] = 1e-12;
    in_base.strain_rate[0][0][1] = 1e-12;
    in_base.strain_rate[0][1][1] = 1e-11;
    if (dim == 3)
      {
        in_base.strain_rate[0][2][0] = 1e-12;
        in_base.strain_rate[0][2][1] = 1e-12;
        in_base.strain_rate[0][2][2] = 1e-11;
      }

    in_base.strain_rate[1] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
    in_base.strain_rate[1][0][0] = -1.71266e-13;
    in_base.strain_rate[1][0][1] = -5.82647e-12;
    in_base.strain_rate[1][1][1] = 4.21668e-14;
    if (dim == 3)
      {
        in_base.strain_rate[1][2][0] = -5.42647e-12;
        in_base.strain_rate[1][2][1] = -5.22647e-12;
        in_base.strain_rate[1][2][2] = 4.21668e-14;
      }
    in_base.strain_rate[2] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
    in_base.strain_rate[2][1][1] = 1e-13;
    in_base.strain_rate[2][0][1] = 1e-11;
    in_base.strain_rate[2][0][0] = -1e-12;
    if (dim == 3)
      {
        in_base.strain_rate[2][2][0] = 1e-11;
        in_base.strain_rate[2][2][1] = 1e-11;
        in_base.strain_rate[2][2][2] = -1e-12;
      }
    in_base.strain_rate[3] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
    in_base.strain_rate[3][1][1] = 4.9e-21;
    in_base.strain_rate[3][0][1] = 4.9e-21;
    in_base.strain_rate[3][0][0] = 4.9e-21;
    if (dim == 3)
      {
        in_base.strain_rate[3][2][0] = 4.9e-21;
        in_base.strain_rate[3][2][1] = 4.9e-21;
        in_base.strain_rate[3][2][2] = 4.9e-21;
      }
    in_base.strain_rate[4] = SymmetricTensor<2,dim>(in_base.strain_rate[0]);
    in_base.strain_rate[4][1][1] = 1e-11;
    in_base.strain_rate[4][0][1] = 1e-11;
    in_base.strain_rate[4][0][0] = 1e-11;
    if (dim == 3)
      {
        in_base.strain_rate[4][2][0] = 1e-11;
        in_base.strain_rate[4][2][1] = 1e-11;
        in_base.strain_rate[4][2][2] = 1e-11;
      }

    // initialize some variables we will need later.
    double finite_difference_accuracy = 1e-7;
    double finite_difference_factor = 1+finite_difference_accuracy;


    MaterialModelInputs<dim> in_dviscositydstrainrate(in_base);

    MaterialModelOutputs<dim> out_base(5,3);
    MaterialModelOutputs<dim> out_dviscositydstrainrate(5,3);

    if (out_base.template get_additional_output<MaterialModelDerivatives<dim>>() != nullptr)
      throw "error";

    out_base.additional_outputs.push_back(std::make_unique<MaterialModelDerivatives<dim>> (5));

    // initialize the material we want to test.
    SimpleNonlinear<dim> mat;
    ParameterHandler prm;
    mat.declare_parameters(prm);

    prm.enter_subsection("Compositional fields");
    {
      prm.set("Number of fields","3");
    }
    prm.leave_subsection();
    prm.enter_subsection("Material model");
    {
      prm.enter_subsection ("Simple nonlinear");
      {
        prm.set ("Viscosity prefactor", "1e-37,1e-36,1e-35,5e-36");
        prm.set ("Viscosity averaging p", std::to_string(parameter));
        prm.set ("Minimum strain rate", 1.4e-20);
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();

    mat.parse_parameters(prm);

    mat.evaluate(in_base, out_base);

    // set up additional output for the derivatives
    MaterialModelDerivatives<dim> *derivatives;
    derivatives = out_base.template get_additional_output<MaterialModelDerivatives<dim>>();
    double temp;

    // have a bool so we know whether the test has succeed or not.
    bool Error = false;

    // this material is not pressure dependent, so we do not test it.

    // test the strain-rate derivative.
    for (unsigned int component = 0; component < SymmetricTensor<2,dim>::n_independent_components; ++component)
      {
        const TableIndices<2> strain_rate_indices = SymmetricTensor<2,dim>::unrolled_to_component_indices (component);

        for (unsigned int i = 0; i < 5; i++)
          {
            in_dviscositydstrainrate.strain_rate[i] = in_base.strain_rate[i]
                                                      + std::fabs(in_base.strain_rate[i][strain_rate_indices])
                                                      * finite_difference_accuracy
                                                      * aspect::Utilities::nth_basis_for_symmetric_tensors<dim>(component);
          }


        mat.evaluate(in_dviscositydstrainrate, out_dviscositydstrainrate);

        for (unsigned int i = 0; i < 5; i++)
          {
            // prevent division by zero. If it is zero, the test has passed, because or
            // the finite difference and the analytical result match perfectly, or (more
            // likely) the material model in independent of this variable.
            temp = out_dviscositydstrainrate.viscosities[i] - out_base.viscosities[i];
            if (temp != 0)
              {
                temp /= std::fabs(in_dviscositydstrainrate.strain_rate[i][strain_rate_indices]) * finite_difference_accuracy;
              }
            std::cout << "strain-rate: component = " << component << ", point = " << i << ", Finite difference = " << temp << ", Analytical derivative = " << derivatives->viscosity_derivative_wrt_strain_rate[i][strain_rate_indices]  << std::endl;
            if (std::fabs(temp - derivatives->viscosity_derivative_wrt_strain_rate[i][strain_rate_indices]) > 1e-3 * (std::fabs(temp) + std::fabs(derivatives->viscosity_derivative_wrt_strain_rate[i][strain_rate_indices])))
              {
                std::cout << "   Error: The derivative of the viscosity to the strain rate is too different from the analytical value." << std::endl;
                Error = true;
              }

          }

      }

    if (Error)
      {
        std::cout << "Some parts of the test were not successful." << std::endl;
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
// test 2D
  int ii2 = f<2>(-1000); // Testing min function
  int iz2 = f<2>(-2); // Testing generalized p norm mean with negative p
  int ij2 = f<2>(-1.5); // Testing generalized p norm mean with negative, non int p
  int ik2 = f<2>(-1); // Testing harmonic mean
  int ji2 = f<2>(0); // Testing geometric mean
  int jj2 = f<2>(1); // Testing arithmetic mean
  int jk2 = f<2>(2); // Testing generalized p norm mean with positive p
  int kj2 = f<2>(1000); // Testing max function
// test 3D
  int ii3 = f<3>(-1000); // Testing min function
  int iz3 = f<3>(-2); // Testing generalized p norm mean with negative p
  int ij3 = f<3>(-1.5); // Testing generalized p norm mean with negative, non int p
  int ik3 = f<3>(-1); // Testing harmonic mean
  int ji3 = f<3>(0); // Testing geometric mean
  int jj3 = f<3>(1); // Testing arithmetic mean
  int jk3 = f<3>(2); // Testing generalized p norm mean with positive p
  int kj3 = f<3>(1000); // Testing max function
// exit
  int kl2 = exit_function();

}
