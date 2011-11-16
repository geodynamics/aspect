/* $Id: step-32.cc 234 2011-10-19 18:07:35Z bangerth $ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University,
     Timo Heister, University of Goettingen, 2008-2011 */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <aspect/simulator.h>
#include <aspect/equation_data.h>
#include <aspect/postprocess_visualization.h>
#include <aspect/model_base.h>
#include <aspect/model_simple.h>
#include <aspect/model_table.h>
#include <aspect/solver.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/matrices.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/base/index_set.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/split_member.hpp>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <locale>
#include <string>

using namespace dealii;

// In the following namespace, we define the
// various pieces of equation data. All of
// these are exhaustively discussed in the
// description of the testcase in the
// introduction:
namespace EquationData
{
  double kappa                 = 1e-6;
  const double reference_specific_heat = 1250;    /* J / K / kg */  //??
  double radiogenic_heating    = 7.4e-12; /* W / kg     */  //??
  double thermal_conductivity = 4.7;
  double reference_gravity    = 30;

  // scale not by R1-R0, but by a
  // typical length scale, say 10km,
  // of variation ("plume
  // diameter"). this choice also
  // roughly equilibrates the sizes
  // of the velocity and pressure
  // components of the solution
  // vectors
  const double length_scale = 10000;


  double R0      = 6371000.-2890000.;     /* m          */
  double R1      = 6371000.-  35000.;     /* m          */
  double apperture_angle = numbers::PI;     /* m          */

  double T0      = 4000+273;              /* K          */
  double T1      =  700+273;              /* K          */

  const double year_in_seconds  = 60*60*24*365.2425;

  struct Perturbation
  {
    double Angle;
    double depth;
    double Amplitude;
    double Sigma;
    double Sign;
    bool GaussianPerturbation;
  };

  Perturbation perturbation;






  template <int dim>
  Tensor<1,dim> gravity_vector (const Point<dim> &p)
  {
// interpolate the following values with a physically realistic model:
//    const double g0      = 10.7;                  /* m / s^2    */
//    const double g1      = 9.81;                  /* m / s^2    */

    const double r = p.norm();
    return -reference_gravity*p/r;
    /*  for now we use a constant gravity */
    /*    return -(1.245e-6 * r + 7.714e13/r/r) * p / r;*/
  }




  template <int dim>
  AdiabaticConditions<dim>::AdiabaticConditions(const aspect::MaterialModel::Interface<dim> * material_model)
    :
    n_points(1000),
    temperatures(n_points, -1),
    pressures(n_points, -1)
  {
    const double delta_z = (R1-R0)/(n_points-1);
    //TODO: look up real value!
    const double dTdp = 2.5e-8;

    // start with these values: 1200K, 1MPa
    temperatures[0] = 1200;
    pressures[0] = 1e6;

    // now integrate downward using the explicit Euler method for simplicity
    //
    // note: p'(z) = rho(p,T) * g
    //       T'(z) = dT/dp|s dp/dz = dT/dp|S rho(p,T) * g
    double z = delta_z;
    for (unsigned int i=1; i<n_points; ++i, z+=delta_z)
      {
        Assert (i < pressures.size(), ExcInternalError());
        Assert (i < temperatures.size(), ExcInternalError());

        //TODO: use the real gravity model here as a function of z
        const Point<dim> representative_point
          = Point<dim>::unit_vector(0) * (R1-z);

        const double density = material_model->density(temperatures[i-1], pressures[i-1], representative_point);

        pressures[i] = (pressures[i-1]
                        + pressures[i-1] * 2/z
                        - density *
                        (gravity_vector(representative_point)*Point<dim>::unit_vector(0)) * delta_z);
        temperatures[i] = (temperatures[i-1] -
                           dTdp * density *
                           (gravity_vector(representative_point)*Point<dim>::unit_vector(0)) * delta_z);
      }

    Assert (*min_element (pressures.begin(), pressures.end()) >= 0, ExcInternalError());
    Assert (*min_element (temperatures.begin(), temperatures.end()) >= 0, ExcInternalError());
  }

  template <int dim>
  double AdiabaticConditions<dim>::pressure (const Point<dim> &p) const
  {
    const double delta_z = (R1-R0)/(n_points-1);

    // clamp the depth to be positive, can happen due to rounding errors on the mesh
    const double z = std::min(std::max(R1 - p.norm(), 0.0), R1-R0 - delta_z);

    const unsigned int i = z/delta_z;
    Assert (z >= 0, ExcInternalError());

    const double d=1.0+i-z/delta_z;
    return d*pressures[i]+(1-d)*pressures[i+1];
  }

  template <int dim>
  double AdiabaticConditions<dim>::temperature (const Point<dim> &p) const
  {
    // clamp the depth to be positive, can happen due to rounding errors on the mesh
    const double z = std::max(R1 - p.norm(), 0.0);
    const double delta_z = (R1-R0)/(n_points-1);

    const unsigned int i = z/delta_z;
    Assert (i < pressures.size(), ExcInternalError());

    //TODO: interpolate linearly
    return temperatures[i];
  }


//instantiation:
  template class AdiabaticConditions<deal_II_dimension>;


  template <int dim>
  class AdiabaticPressure : public Function<dim>
  {
    public:
      AdiabaticPressure(const AdiabaticConditions<dim> & ad_c)
        : adiabatic_conditions(ad_c)
      {}
      double value (const Point<dim> &p,
                    const unsigned int = 0) const
      {
        return adiabatic_conditions.pressure (p);
      }
    private:
      const AdiabaticConditions<dim> & adiabatic_conditions;
  };


  template <int dim>
  class TemperatureInitialValues : public Function<dim>
  {
    public:
      TemperatureInitialValues () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
                                 Vector<double>   &value) const;
  };



  template <int dim>
  double
  TemperatureInitialValues<dim>::value (const Point<dim>  &p,
                                        const unsigned int) const
  {
    const double r = p.norm();
    const double h = R1-R0;

    // s = fraction of the way from
    // the inner to the outer
    // boundary; 0<=s<=1
    const double s = (r-R0)/h;
    double Perturbation = 0e0;
    double InterpolVal = 0e0;
    double depth[4];
    double geotherm[4];
    double x, y;
    double scale=R1/(R1 - R0);
    float eps = 1e-4;
    if (!EquationData::perturbation.GaussianPerturbation)
      {

        /* now compute an angular variation of the linear temperature field by
           stretching the variable s appropriately. note that the following
           formula leaves the end points s=0 and s=1 fixed, but stretches the
           region in between depending on the angle phi=atan2(x,y).

           For a plot, see
           http://www.wolframalpha.com/input/?i=plot+%28%282*sqrt%28x^2%2By^2%29-1%29%2B0.2*%282*sqrt%28x^2%2By^2%29-1%29*%281-%282*sqrt%28x^2%2By^2%29-1%29%29*sin%286*atan2%28x%2Cy%29%29%29%2C+x%3D-1+to+1%2C+y%3D-1+to+1
        */
        const double scale = (dim==3)?std::max(0.0,cos(3.14159*abs(p(2)/R1))):1.0;
        const double phi   = std::atan2(p(0),p(1));
        const double s_mod = s
                             +
                             0.2 * s * (1-s) * std::sin(6*phi) * scale;

        return T0*(1.0-s_mod) + T1*s_mod;
      }
    else
      {
        geotherm[3]=T1;
        geotherm[2]=T1 + 1200e0;
        geotherm[1]=T0 - 1300e0;
        geotherm[0]=T0;
        depth[0]=R0-1e-2*R0;
        depth[1]=R0+500e3;
        depth[2]=R1-500e3;
        depth[3]=R1+1e-2*R1;

        int indx = -1;
        for (unsigned int i=0; i<3; ++i)
          {
            if ((depth[i] - r) < eps && (depth[i+1] - r) > eps)
              {
                indx = i;
                break;
              }
          }
        Assert (indx >= 0,                  ExcInternalError());
        Assert (indx < 3,                  ExcInternalError());
        int indx1 = indx + 1;
        float dx = depth[indx1] - depth[indx];
        float dy = geotherm[indx1] - geotherm[indx];

        if ( dx > 0.5*eps)
          {
            // linear interpolation
            InterpolVal = std::max(geotherm[3],geotherm[indx] + (r-depth[indx]) * (dy/dx));
          }
        else
          {
            // eval.point in discontinuity
            InterpolVal = 0.5*( geotherm[indx] + geotherm[indx1] );
          }
        x = (scale - EquationData::perturbation.depth)*std::cos(EquationData::perturbation.Angle);
        y = (scale - EquationData::perturbation.depth)*std::sin(EquationData::perturbation.Angle);
        Perturbation = EquationData::perturbation.Sign*EquationData::perturbation.Amplitude*T0*std::exp( -( std::pow((p(0)*scale/R1-x),2) +std::pow((p(1)*scale/R1-y),2) ) / EquationData::perturbation.Sigma) ;
        if (r > R1 - 1e-2*R1)
          {
            return geotherm[3];
          }
        else if (r < R0 + 1e-2*R0)
          {
            return geotherm[0];
          }
        else
          {
            return InterpolVal + Perturbation;
          }
      }
  }



  template <int dim>
  void
  TemperatureInitialValues<dim>::vector_value (const Point<dim> &p,
                                               Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = TemperatureInitialValues<dim>::value (p, c);
  }
}




namespace aspect
{
  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        struct StokesPreconditioner
        {
          StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
                                const Quadrature<dim>    &stokes_quadrature,
                                const Mapping<dim>       &mapping,
                                const UpdateFlags         stokes_update_flags,
                                const FiniteElement<dim> &temperature_fe,
                                const UpdateFlags         temperature_update_flags);
          StokesPreconditioner (const StokesPreconditioner &data);

          FEValues<dim>               stokes_fe_values;
          FEValues<dim>               temperature_fe_values;

          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  phi_p;

          std::vector<double>                  old_temperature_values;
          std::vector<double>                  old_pressure_values;
        };

        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
                              const Quadrature<dim>    &stokes_quadrature,
                              const Mapping<dim>       &mapping,
                              const UpdateFlags         stokes_update_flags,
                              const FiniteElement<dim> &temperature_fe,
                              const UpdateFlags         temperature_update_flags)
          :
          stokes_fe_values (mapping, stokes_fe, stokes_quadrature,
                            stokes_update_flags),
          temperature_fe_values (mapping,
                                 temperature_fe, stokes_quadrature,
                                 temperature_update_flags),
          grads_phi_u (stokes_fe.dofs_per_cell),
          phi_p (stokes_fe.dofs_per_cell),
          old_temperature_values (stokes_quadrature.size()),
          old_pressure_values (stokes_quadrature.size())
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const StokesPreconditioner &scratch)
          :
          stokes_fe_values (scratch.stokes_fe_values.get_mapping(),
                            scratch.stokes_fe_values.get_fe(),
                            scratch.stokes_fe_values.get_quadrature(),
                            scratch.stokes_fe_values.get_update_flags()),
          temperature_fe_values (scratch.temperature_fe_values.get_mapping(),
                                 scratch.temperature_fe_values.get_fe(),
                                 scratch.temperature_fe_values.get_quadrature(),
                                 scratch.temperature_fe_values.get_update_flags()),
          grads_phi_u (scratch.grads_phi_u),
          phi_p (scratch.phi_p),
          old_temperature_values (scratch.old_temperature_values),
          old_pressure_values (scratch.old_pressure_values)
        {}



        // Observe that we derive the
        // StokesSystem scratch array from the
        // StokesPreconditioner array. We do this
        // because all the objects that are
        // necessary for the assembly of the
        // preconditioner are also needed for the
        // actual matrix system and right hand
        // side, plus some extra data. This makes
        // the program more compact. Note also
        // that the assembly of the Stokes system
        // and the temperature right hand side
        // further down requires data from
        // temperature and velocity,
        // respectively, so we actually need two
        // FEValues objects for those two cases.
        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &stokes_fe,
                        const Mapping<dim>       &mapping,
                        const Quadrature<dim>    &stokes_quadrature,
                        const UpdateFlags         stokes_update_flags,
                        const FiniteElement<dim> &temperature_fe,
                        const UpdateFlags         temperature_update_flags);

          StokesSystem (const StokesSystem<dim> &data);

          std::vector<Tensor<1,dim> >          phi_u;
          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  div_phi_u;
          std::vector<Tensor<1,dim> > old_velocity_values;
        };


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const FiniteElement<dim> &stokes_fe,
                      const Mapping<dim>       &mapping,
                      const Quadrature<dim>    &stokes_quadrature,
                      const UpdateFlags         stokes_update_flags,
                      const FiniteElement<dim> &temperature_fe,
                      const UpdateFlags         temperature_update_flags)
          :
          StokesPreconditioner<dim> (stokes_fe, stokes_quadrature,
                                     mapping,
                                     stokes_update_flags,
                                     temperature_fe,
                                     temperature_update_flags),
          phi_u (stokes_fe.dofs_per_cell),
          grads_phi_u (stokes_fe.dofs_per_cell),
          div_phi_u (stokes_fe.dofs_per_cell),
          old_velocity_values (stokes_quadrature.size())
        {}


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &scratch)
          :
          StokesPreconditioner<dim> (scratch),
          phi_u (scratch.phi_u),
          grads_phi_u (scratch.grads_phi_u),
          div_phi_u (scratch.div_phi_u),
          old_velocity_values (scratch.old_velocity_values)
        {}



        template <int dim>
        struct TemperatureSystem
        {
          TemperatureSystem (const FiniteElement<dim> &temperature_fe,
                             const FiniteElement<dim> &stokes_fe,
                             const Mapping<dim>       &mapping,
                             const Quadrature<dim>    &quadrature);
          TemperatureSystem (const TemperatureSystem &data);

          FEValues<dim>               temperature_fe_values;
          FEValues<dim>               stokes_fe_values;

          std::vector<double>         phi_T;
          std::vector<Tensor<1,dim> > grad_phi_T;

          std::vector<Tensor<1,dim> > old_velocity_values;
          std::vector<Tensor<1,dim> > old_old_velocity_values;

          std::vector<double>         old_pressure;
          std::vector<double>         old_old_pressure;

          std::vector<SymmetricTensor<2,dim> > old_strain_rates;
          std::vector<SymmetricTensor<2,dim> > old_old_strain_rates;

          std::vector<double>         old_temperature_values;
          std::vector<double>         old_old_temperature_values;
          std::vector<Tensor<1,dim> > old_temperature_grads;
          std::vector<Tensor<1,dim> > old_old_temperature_grads;
          std::vector<double>         old_temperature_laplacians;
          std::vector<double>         old_old_temperature_laplacians;
        };

        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const FiniteElement<dim> &temperature_fe,
                           const FiniteElement<dim> &stokes_fe,
                           const Mapping<dim>       &mapping,
                           const Quadrature<dim>    &quadrature)
          :
          temperature_fe_values (mapping,
                                 temperature_fe, quadrature,
                                 update_values    |
                                 update_gradients |
                                 update_hessians  |
                                 update_quadrature_points |
                                 update_JxW_values),
          stokes_fe_values(mapping, stokes_fe, quadrature,
                           update_values    | update_gradients |
                           update_JxW_values),
          phi_T (temperature_fe.dofs_per_cell),
          grad_phi_T (temperature_fe.dofs_per_cell),
          old_velocity_values (quadrature.size()),
          old_old_velocity_values (quadrature.size()),
          old_pressure (quadrature.size()),
          old_old_pressure (quadrature.size()),
          old_strain_rates (quadrature.size()),
          old_old_strain_rates (quadrature.size()),
          old_temperature_values (quadrature.size()),
          old_old_temperature_values(quadrature.size()),
          old_temperature_grads(quadrature.size()),
          old_old_temperature_grads(quadrature.size()),
          old_temperature_laplacians(quadrature.size()),
          old_old_temperature_laplacians(quadrature.size())
        {}


        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const TemperatureSystem &scratch)
          :
          temperature_fe_values (scratch.temperature_fe_values.get_mapping(),
                                 scratch.temperature_fe_values.get_fe(),
                                 scratch.temperature_fe_values.get_quadrature(),
                                 scratch.temperature_fe_values.get_update_flags()),
          stokes_fe_values (scratch.stokes_fe_values.get_mapping(),
                            scratch.stokes_fe_values.get_fe(),
                            scratch.stokes_fe_values.get_quadrature(),
                            scratch.stokes_fe_values.get_update_flags()),
          phi_T (scratch.phi_T),
          grad_phi_T (scratch.grad_phi_T),
          old_velocity_values (scratch.old_velocity_values),
          old_old_velocity_values (scratch.old_old_velocity_values),
          old_pressure (scratch.old_pressure),
          old_old_pressure (scratch.old_old_pressure),
          old_strain_rates (scratch.old_strain_rates),
          old_old_strain_rates (scratch.old_old_strain_rates),
          old_temperature_values (scratch.old_temperature_values),
          old_old_temperature_values (scratch.old_old_temperature_values),
          old_temperature_grads (scratch.old_temperature_grads),
          old_old_temperature_grads (scratch.old_old_temperature_grads),
          old_temperature_laplacians (scratch.old_temperature_laplacians),
          old_old_temperature_laplacians (scratch.old_old_temperature_laplacians)
        {}
      }


      // The CopyData arrays are similar to the
      // Scratch arrays. They provide a
      // constructor, a copy operation, and
      // some arrays for local matrix, local
      // vectors and the relation between local
      // and global degrees of freedom (a.k.a.
      // <code>local_dof_indices</code>).
      namespace CopyData
      {
        template <int dim>
        struct StokesPreconditioner
        {
          StokesPreconditioner (const FiniteElement<dim> &stokes_fe);
          StokesPreconditioner (const StokesPreconditioner &data);

          FullMatrix<double>          local_matrix;
          std::vector<unsigned int>   local_dof_indices;
        };

        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const FiniteElement<dim> &stokes_fe)
          :
          local_matrix (stokes_fe.dofs_per_cell,
                        stokes_fe.dofs_per_cell),
          local_dof_indices (stokes_fe.dofs_per_cell)
        {}



        template <int dim>
        StokesPreconditioner<dim>::
        StokesPreconditioner (const StokesPreconditioner &data)
          :
          local_matrix (data.local_matrix),
          local_dof_indices (data.local_dof_indices)
        {}



        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &stokes_fe);
          StokesSystem (const StokesSystem<dim> &data);

          Vector<double> local_rhs;
          Vector<double> local_rhs_helper;

        };


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const FiniteElement<dim> &stokes_fe)
          :
          StokesPreconditioner<dim> (stokes_fe),
          local_rhs (stokes_fe.dofs_per_cell),
          local_rhs_helper (stokes_fe.dofs_per_cell)
        {}


        template <int dim>
        StokesSystem<dim>::
        StokesSystem (const StokesSystem<dim> &data)
          :
          StokesPreconditioner<dim> (data),
          local_rhs (data.local_rhs),
          local_rhs_helper (data.local_rhs_helper)
        {}


        template <int dim>
        struct TemperatureSystem
        {
          TemperatureSystem (const FiniteElement<dim> &temperature_fe);
          TemperatureSystem (const TemperatureSystem &data);

          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;
          std::vector<unsigned int>   local_dof_indices;
        };

        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const FiniteElement<dim> &temperature_fe)
          :
          local_matrix (temperature_fe.dofs_per_cell,
                        temperature_fe.dofs_per_cell),
          local_rhs (temperature_fe.dofs_per_cell),
          local_dof_indices (temperature_fe.dofs_per_cell)
        {}


        template <int dim>
        TemperatureSystem<dim>::
        TemperatureSystem (const TemperatureSystem &data)
          :
          local_matrix (data.local_matrix),
          local_rhs (data.local_rhs),
          local_dof_indices (data.local_dof_indices)
        {}
      }
    }



    template <class T>
    inline T bdf2_extrapolate(
      bool use_bdf_scheme,
      double old_time_step, double time_step,
      const T &old_data, const T &new_data)
    {
      return (use_bdf_scheme) ?
             (new_data * (1 + time_step/old_time_step)
              - old_data * time_step/old_time_step)
             : new_data;
    }
  }
}




namespace aspect
{

  template <int dim>
  Simulator<dim>::Parameters::Parameters (ParameterHandler &prm)
  {
    parse_parameters (prm);
  }


  template <int dim>
  void
  Simulator<dim>::Parameters::
  declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Resume computation", "false",
                       Patterns::Bool (),
                       "should the last computation be resumed?");

    prm.declare_entry ("End time", "1e8",
                       Patterns::Double (0),
                       "The end time of the simulation in years.");

    prm.enter_subsection ("Mesh refinement");
    {
      prm.declare_entry ("Initial global refinement", "2",
                         Patterns::Integer (0),
                         "The number of global refinement steps performed on "
                         "the initial coarse mesh, before the problem is first "
                         "solved there.");
      prm.declare_entry ("Initial adaptive refinement", "2",
                         Patterns::Integer (0),
                         "The number of adaptive refinement steps performed after "
                         "initial global refinement.");
      prm.declare_entry ("Time steps between mesh refinement", "10",
                         Patterns::Integer (1),
                         "The number of time steps after which the mesh is to be "
                         "adapted based on computed error indicators.");
      prm.declare_entry ("Refinement fraction", "0.3",
                         Patterns::Double(0,1),
                         "The fraction of cells with the largest error that "
                         "should be flagged for refinement.");
      prm.declare_entry ("Coarsening fraction", "0.05",
                         Patterns::Double(0,1),
                         "The fraction of cells with the smallest error that "
                         "should be flagged for coarsening.");
      prm.declare_entry ("Additional refinement times", "",
                         Patterns::List (Patterns::Double(0)),
                         "A list of times (in years) so that if the end time of a time step "
                         "is beyond this time, an additional round of mesh refinement "
                         "is triggered. This is mostly useful to make sure we "
                         "can get through the initial transient phase of a simulation "
                         "on a relatively coarse mesh, and then refine again when we "
                         "are in a time range that we are interested in and where "
                         "we would like to use a finer mesh.");
    }
    prm.leave_subsection();

    prm.enter_subsection ("Stabilization parameters");
    {
      prm.declare_entry ("alpha", "2",
                         Patterns::Double (1, 2),
                         "The exponent in the entropy viscosity stabilization.");
      prm.declare_entry ("c_R", "0.11",
                         Patterns::Double (0),
                         "The c_R factor in the entropy viscosity "
                         "stabilization.");
      prm.declare_entry ("beta", "0.078",
                         Patterns::Double (0),
                         "The beta factor in the artificial viscosity "
                         "stabilization. An appropriate value for 2d is 0.052 "
                         "and 0.078 for 3d.");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Thermal perturbation");
    {
      prm.declare_entry ("Angle", "0e0",
                         Patterns::Double (0),
                         "The angle where the center of the perturbation is placed.");
      prm.declare_entry ("non-dim depth", "0.7",
                         Patterns::Double (0),
                         "The radial distance where the center of the perturbation is placed.");
      prm.declare_entry ("Amplitude", "0.01",
                         Patterns::Double (0),
                         "The amplitude of the perturbation.");
      prm.declare_entry ("Sigma", "0.2",
                         Patterns::Double (0),
                         "The standard deviation of the Gaussian perturbation.");
      prm.declare_entry ("Sign", "1",
                         Patterns::Double (),
                         "The sign of the perturbation.");
      prm.declare_entry ("gaussian perturbation", "true",
                         Patterns::Bool (),
                         "The sign of the perturbation.");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("ModelSettings");
    {
      prm.declare_entry ("Include shear heating", "true",
                         Patterns::Bool (),
                         "Whether to include shear heating into the model or not. From a "
                         "physical viewpoint, shear heating should always be used but may "
                         "be undesirable when comparing results with known benchmarks that "
                         "do not include this term in the temperature equation.");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Material model");
    {
      prm.declare_entry ("model", "simple",
                         Patterns::Selection ("simple|table"),
                         "select the active model: simple|table");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("ModelParameters");
    {
      prm.declare_entry ("kappa", "4.548e-7",
                         Patterns::Double (),
                         "thermal diffusivity (k/(rho*cp)");
      prm.declare_entry ("radiogenic_heating", "0e0",
                         Patterns::Double (),
                         "H0");
      prm.declare_entry ("thermal_conductivity", "4.7",
                         Patterns::Double (),
                         "k");
      prm.declare_entry ("R1", "10415e3",
                         Patterns::Double (),
                         "Outer radius");
      prm.declare_entry ("R0", "4717e3",
                         Patterns::Double (),
                         "Inner radius");
      prm.declare_entry ("apperture_angle", "180",
                         Patterns::Double (),
                         "apperture angle (opening angle)");
      prm.declare_entry ("T1", "0",
                         Patterns::Double (),
                         "temperature at outer boundary (lythosphere water/air)");
      prm.declare_entry ("T0", "6000",
                         Patterns::Double (),
                         "temperature at inner boundary (core mantle boundary");
      prm.declare_entry ("reference_gravity", "30",
                         Patterns::Double (),
                         "g0");
    }
    prm.leave_subsection ();


    prm.enter_subsection ("Discretization");
    {
      prm.declare_entry ("Stokes velocity polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the velocity variables "
                         "in the Stokes system.");
      prm.declare_entry ("Temperature polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the temperature variable.");
      prm.declare_entry ("Use locally conservative discretization", "true",
                         Patterns::Bool (),
                         "Whether to use a Stokes discretization that is locally "
                         "conservative at the expense of a larger number of degrees "
                         "of freedom, or to go with a cheaper discretization "
                         "that does not locally conserve mass (although it is "
                         "globally conservative.");
    }
    prm.leave_subsection ();
  }



  template <int dim>
  void
  Simulator<dim>::Parameters::
  parse_parameters (ParameterHandler &prm)
  {
    resume_computation      = prm.get_bool ("Resume computation");
    end_time                = prm.get_double ("End time");

    prm.enter_subsection ("Mesh refinement");
    {
      initial_global_refinement   = prm.get_integer ("Initial global refinement");
      initial_adaptive_refinement = prm.get_integer ("Initial adaptive refinement");

      adaptive_refinement_interval= prm.get_integer ("Time steps between mesh refinement");
      refinement_fraction         = prm.get_double ("Refinement fraction");
      coarsening_fraction         = prm.get_double ("Coarsening fraction");

      // extract the list of times at which additional refinement is requested
      // then sort it and convert it to seconds
      additional_refinement_times
        = Utilities::string_to_double
          (Utilities::split_string_list(prm.get ("Additional refinement times")));
      std::sort (additional_refinement_times.begin(),
                 additional_refinement_times.end());
      for (unsigned int i=0; i<additional_refinement_times.size(); ++i)
        additional_refinement_times[i] *= EquationData::year_in_seconds;
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Stabilization parameters");
    {
      stabilization_alpha = prm.get_double ("alpha");
      stabilization_c_R   = prm.get_double ("c_R");
      stabilization_beta  = prm.get_double ("beta");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Thermal perturbation");
    {
      EquationData::perturbation.Angle = prm.get_double ("Angle");
      EquationData::perturbation.depth = prm.get_double ("non-dim depth");
      //EquationData::perturbation.depth = EquationData::R1 - (EquationData::R1 - EquationData::R0)*prm.get_double ("non-dim depth");
      EquationData::perturbation.Amplitude  = prm.get_double ("Amplitude");
      EquationData::perturbation.Sigma  = prm.get_double ("Sigma");
      EquationData::perturbation.Sign  = prm.get_double ("Sign");
      EquationData::perturbation.GaussianPerturbation  = prm.get_bool ("gaussian perturbation");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("ModelSettings");
    {
      include_shear_heating = prm.get_bool ("Include shear heating");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Material model");
    {
      model = prm.get ("model");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("ModelParameters");
    {
      EquationData::kappa = prm.get_double ("kappa");
      EquationData::radiogenic_heating = prm.get_double ("radiogenic_heating");
      EquationData::thermal_conductivity = prm.get_double ("thermal_conductivity");
      EquationData::R1 = prm.get_double ("R1");
      EquationData::R0 = prm.get_double ("R0");
      EquationData::apperture_angle = (std::acos(-1e0)/180e0)*prm.get_double ("apperture_angle");
      EquationData::T0 = prm.get_double ("T0");
      EquationData::T1 = prm.get_double ("T1");
      EquationData::reference_gravity = prm.get_double ("reference_gravity");
    }

    prm.leave_subsection ();

    prm.enter_subsection ("Discretization");
    {
      stokes_velocity_degree = prm.get_integer ("Stokes velocity polynomial degree");
      temperature_degree     = prm.get_integer ("Temperature polynomial degree");
      use_locally_conservative_discretization
        = prm.get_bool ("Use locally conservative discretization");
    }
    prm.leave_subsection ();
  }



  template <int dim>
  Simulator<dim>::Simulator (ParameterHandler &prm)
    :
    parameters (prm),
    pcout (std::cout,
           (Utilities::MPI::
            this_mpi_process(MPI_COMM_WORLD)
            == 0)),

    material_model (MaterialModel::create<dim>(parameters.model)),
    adiabatic_conditions(material_model.get()),

    triangulation (MPI_COMM_WORLD,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening),
                   parallel::distributed::Triangulation<dim>::mesh_reconstruction_after_repartitioning
                  ),

    mapping (4),

    stokes_fe (FE_Q<dim>(parameters.stokes_velocity_degree),
               dim,
               (parameters.use_locally_conservative_discretization
                ?
                static_cast<const FiniteElement<dim> &>
                (FE_DGP<dim>(parameters.stokes_velocity_degree-1))
                :
                static_cast<const FiniteElement<dim> &>
                (FE_Q<dim>(parameters.stokes_velocity_degree-1))),
               1),

    stokes_dof_handler (triangulation),

    temperature_fe (parameters.temperature_degree),
    temperature_dof_handler (triangulation),

    time (0),
    time_step (0),
    old_time_step (0),
    timestep_number (0),
    rebuild_stokes_matrix (true),
    rebuild_stokes_preconditioner (true),

    computing_timer (pcout, TimerOutput::summary,
                     TimerOutput::wall_times)
  {
    postprocess_manager.parse_parameters (prm);
    postprocess_manager.initialize (*this);

    material_model->parse_parameters(prm);

    pressure_scaling = material_model->reference_viscosity() / EquationData::length_scale;
    
    // make sure that we don't have to fill every column of the statistics
    // object in each time step.
    statistics.set_auto_fill_mode(true);
  }



  template <int dim>
  void Simulator<dim>::declare_parameters (ParameterHandler &prm)
  {
    Parameters::declare_parameters (prm);
    Postprocess::Manager<dim>::declare_parameters (prm);
    MaterialModel::Simple<dim>::declare_parameters (prm);
    MaterialModel::Table<dim>::declare_parameters (prm);
  }



  template <int dim>
  double Simulator<dim>::get_maximal_velocity () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, stokes_fe, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_velocity = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (stokes_solution,
                                                     velocity_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
        }

    return Utilities::MPI::max (max_local_velocity, MPI_COMM_WORLD);
  }



  template <int dim>
  double
  Simulator<dim>::get_entropy_variation (const double average_temperature) const
  {
    // only do this if we really need entropy
    // variation
    if (parameters.stabilization_alpha != 2)
      return 1.;

    // record maximal entropy on Gauss quadrature
    // points
    const QGauss<dim> quadrature_formula (parameters.temperature_degree+1);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (temperature_fe, quadrature_formula,
                             update_values | update_JxW_values);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<double> old_old_temperature_values(n_q_points);

    double min_entropy = std::numeric_limits<double>::max(),
           max_entropy = -std::numeric_limits<double>::max(),
           area = 0,
           entropy_integrated = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values.get_function_values (old_temperature_solution,
                                         old_temperature_values);
          fe_values.get_function_values (old_old_temperature_solution,
                                         old_old_temperature_values);
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const double T = (old_temperature_values[q] +
                                old_old_temperature_values[q]) / 2;
              const double entropy = ((T-average_temperature) *
                                      (T-average_temperature));

              min_entropy = std::min (min_entropy, entropy);
              max_entropy = std::max (max_entropy, entropy);
              area += fe_values.JxW(q);
              entropy_integrated += fe_values.JxW(q) * entropy;
            }
        }

    // do MPI data exchange: we need to sum over
    // the two integrals (area,
    // entropy_integrated), and get the extrema
    // for maximum and minimum. combine
    // MPI_Allreduce for two values since that is
    // an expensive operation
    const double local_for_sum[2] = { entropy_integrated, area },
                                    local_for_max[2] = { -min_entropy, max_entropy };
    double global_for_sum[2], global_for_max[2];

    Utilities::MPI::sum (local_for_sum, MPI_COMM_WORLD, global_for_sum);
    Utilities::MPI::max (local_for_max, MPI_COMM_WORLD, global_for_max);

    const double average_entropy = global_for_sum[0] / global_for_sum[1];
    const double entropy_diff = std::max(global_for_max[1] - average_entropy,
                                         average_entropy - (-global_for_max[0]));
    return entropy_diff;
  }



  // Similar function to before, but we now
  // compute the cfl number, i.e., maximal
  // velocity on a cell divided by the cell
  // diameter
  template <int dim>
  double Simulator<dim>::get_cfl_number () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, stokes_fe, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_cfl = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (stokes_solution,
                                                     velocity_values);

          double max_local_velocity = 1e-10;
          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
          max_local_cfl = std::max(max_local_cfl,
                                   max_local_velocity / cell->diameter());
        }

    return Utilities::MPI::max (max_local_cfl, MPI_COMM_WORLD);
  }



  template <int dim>
  std::pair<double,double>
  Simulator<dim>::get_extrapolated_temperature_range () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.temperature_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, temperature_fe, quadrature_formula,
                             update_values);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<double> old_old_temperature_values(n_q_points);

    // This presets the minimum with a bigger
    // and the maximum with a smaller number
    // than one that is going to appear. Will
    // be overwritten in the cell loop or in
    // the communication step at the
    // latest.
    double min_local_temperature = std::numeric_limits<double>::max(),
           max_local_temperature = -std::numeric_limits<double>::max();

    if (timestep_number != 0)
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = temperature_dof_handler.begin_active(),
        endc = temperature_dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values.get_function_values (old_temperature_solution,
                                             old_temperature_values);
              fe_values.get_function_values (old_old_temperature_solution,
                                             old_old_temperature_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double temperature =
                    (1. + time_step/old_time_step) * old_temperature_values[q]-
                    time_step/old_time_step * old_old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    temperature);
                }
            }
      }
    else
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = temperature_dof_handler.begin_active(),
        endc = temperature_dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values.get_function_values (old_temperature_solution,
                                             old_temperature_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double temperature = old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    temperature);
                }
            }
      }

    return std::make_pair(-Utilities::MPI::max (-min_local_temperature,
                                                MPI_COMM_WORLD),
                          Utilities::MPI::max (max_local_temperature,
                                               MPI_COMM_WORLD));
  }



  // The function that calculates the
  // viscosity is purely local
  template <int dim>
  double
  Simulator<dim>::
  compute_viscosity (const std::vector<double>          &old_temperature,
                     const std::vector<double>          &old_old_temperature,
                     const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                     const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                     const std::vector<double>          &old_temperature_laplacians,
                     const std::vector<double>          &old_old_temperature_laplacians,
                     const std::vector<Tensor<1,dim> >  &old_velocity_values,
                     const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                     const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
                     const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
                     const std::vector<double>          &old_pressure,
                     const std::vector<double>          &old_old_pressure,
                     const double                        global_u_infty,
                     const double                        global_T_variation,
                     const double                        average_temperature,
                     const double                        global_entropy_variation,
                     const std::vector<Point<dim> >     &evaluation_points,
                     const double                        cell_diameter) const
  {
    if (global_u_infty == 0)
      return 5e-3 * cell_diameter;

    const unsigned int n_q_points = old_temperature.size();

    double max_residual = 0;
    double max_velocity = 0;

    for (unsigned int q=0; q < n_q_points; ++q)
      {
        const Tensor<1,dim> u = (old_velocity_values[q] +
                                 old_old_velocity_values[q]) / 2;

        const SymmetricTensor<2,dim> strain_rate = (old_strain_rates[q] +
                                                    old_old_strain_rates[q]) / 2;

        const double T = (old_temperature[q] + old_old_temperature[q]) / 2;
        const double p = (old_pressure[q] + old_old_pressure[q]) / 2;
        const double dT_dt = (old_temperature[q] - old_old_temperature[q])
                             / old_time_step;
        const double u_grad_T = u * (old_temperature_grads[q] +
                                     old_old_temperature_grads[q]) / 2;

        const double kappa_Delta_T = EquationData::kappa
                                     * (old_temperature_laplacians[q] +
                                        old_old_temperature_laplacians[q]) / 2;

        const double density = material_model->density(T, p, evaluation_points[q]);

        const double gamma
          = (EquationData::radiogenic_heating * density
             +
             (parameters.include_shear_heating
              ?
              2 * material_model->viscosity(T, p, evaluation_points[q]) *
              strain_rate * strain_rate  /
              (density * material_model->specific_heat(T, p, evaluation_points[q]))
              :
              0));

        double residual
          = std::abs(dT_dt + u_grad_T - kappa_Delta_T - gamma);
        if (parameters.stabilization_alpha == 2)
          residual *= std::abs(T - average_temperature);

        max_residual = std::max (residual,        max_residual);
        max_velocity = std::max (std::sqrt (u*u), max_velocity);
      }

    const double max_viscosity = (parameters.stabilization_beta *
                                  max_velocity * cell_diameter);
    if (timestep_number == 0)
      return max_viscosity;
    else
      {
        Assert (old_time_step > 0, ExcInternalError());

        double entropy_viscosity;
        if (parameters.stabilization_alpha == 2)
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * cell_diameter *
                               max_residual /
                               global_entropy_variation);
        else
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * global_Omega_diameter *
                               max_velocity * max_residual /
                               (global_u_infty * global_T_variation));

        return std::min (max_viscosity, entropy_viscosity);
      }
  }



  template <int dim>
  void Simulator<dim>::set_initial_temperature_field ()
  {
    // create a fully distributed vector since
    // the VectorTools::interpolate function
    // needs to write into it and we can not
    // write into vectors with ghost elements
    TrilinosWrappers::MPI::Vector
    solution (temperature_matrix.row_partitioner());

    // interpolate the initial values
    VectorTools::interpolate (mapping,
                              temperature_dof_handler,
                              EquationData::TemperatureInitialValues<dim>(),
                              solution);

    // then apply constraints and copy the
    // result into vectors with ghost elements
    temperature_constraints.distribute (solution);
    temperature_solution = solution;
    old_temperature_solution = solution;
    old_old_temperature_solution = solution;
  }



  template <int dim>
  void Simulator<dim>::compute_initial_pressure_field ()
  {
    // we'd like to interpolate the initial pressure onto the pressure
    // variable but that's a bit involved because the pressure may either
    // be an FE_Q (for which we can interpolate) or an FE_DGP (for which
    // we can't since the element has no nodal basis.
    //
    // fortunately, in the latter case, the element is discontinuous and
    // we can compute a local projection onto the pressure space
    if (parameters.use_locally_conservative_discretization == false)
      {
        class InitialConditions : public Function<dim>
        {
          public:
            InitialConditions (const EquationData::AdiabaticConditions<dim> &ad_c)
              : Function<dim> (dim+1), adiabatic_conditions (ad_c)
            {}

            double value (const Point<dim> &p,
                          const unsigned int component) const
            {
              switch (component)
                {
                  case dim:
                    return adiabatic_conditions.pressure (p);
                  default:
                    return 0;
                }
            }

            const EquationData::AdiabaticConditions<dim> & adiabatic_conditions;
        };

        TrilinosWrappers::MPI::BlockVector stokes_tmp;
        stokes_tmp.reinit (stokes_rhs);
        VectorTools::interpolate (mapping, stokes_dof_handler,
                                  InitialConditions(adiabatic_conditions),
                                  stokes_tmp);
        old_stokes_solution = stokes_tmp;
      }
    else
      {
        // implement the local projection for the discontinuous pressure
        // element. this is only going to work if, indeed, the element
        // is discontinuous
        const FiniteElement<dim> &pressure_fe = stokes_fe.base_element(1);
        Assert (pressure_fe.dofs_per_face == 0,
                ExcNotImplemented());

        QGauss<dim> quadrature(parameters.stokes_velocity_degree+1);
        UpdateFlags update_flags = UpdateFlags(update_values   |
                                               update_quadrature_points |
                                               update_JxW_values);
        FEValues<dim> fe_values (mapping, stokes_fe, quadrature, update_flags);
        const FEValuesExtractors::Scalar pressure (dim);

        const unsigned int
        dofs_per_cell = fe_values.dofs_per_cell,
        n_q_points    = fe_values.n_quadrature_points;

        std::vector<unsigned int> local_dof_indices (dofs_per_cell);
        Vector<double> cell_vector (dofs_per_cell);
        Vector<double> local_projection (dofs_per_cell);
        FullMatrix<double> local_mass_matrix (dofs_per_cell, dofs_per_cell);

        std::vector<double> rhs_values(n_q_points);

        EquationData::AdiabaticPressure<dim> ad_pressure (adiabatic_conditions);


        typename DoFHandler<dim>::active_cell_iterator
        cell = stokes_dof_handler.begin_active(),
        endc = stokes_dof_handler.end();

        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              cell->get_dof_indices (local_dof_indices);
              fe_values.reinit(cell);

              ad_pressure.value_list
              (fe_values.get_quadrature_points(), rhs_values);

              cell_vector = 0;
              local_mass_matrix = 0;
              for (unsigned int point=0; point<n_q_points; ++point)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    if (stokes_fe.system_to_component_index(i).first == dim)
                      cell_vector(i)
                      +=
                        rhs_values[point] *
                        fe_values[pressure].value(i,point) *
                        fe_values.JxW(point);

                    // populate the local matrix; create the pressure mass matrix
                    // in the pressure pressure block and the identity matrix
                    // for all other variables so that the whole thing remains
                    // invertible
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      if ((stokes_fe.system_to_component_index(i).first == dim)
                          &&
                          (stokes_fe.system_to_component_index(j).first == dim))
                        local_mass_matrix(j,i) += (fe_values[pressure].value(i,point) *
                                                   fe_values[pressure].value(j,point) *
                                                   fe_values.JxW(point));
                      else if (i == j)
                        local_mass_matrix(i,j) = 1;
                  }

              // now invert the local mass matrix and multiply it with the rhs
              local_mass_matrix.gauss_jordan();
              local_mass_matrix.vmult (local_projection, cell_vector);

              // then set the global solution vector to the values just computed
              cell->set_dof_values (local_projection, old_stokes_solution);
            }
      }

    // normalize the pressure in such a way that the surface pressure
    // equals a known and desired value
    normalize_pressure(old_stokes_solution);

    // set the current solution to the same value as the previous solution
    stokes_solution = old_stokes_solution;
  }


  template <int dim>
  void Simulator<dim>::
  setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning)
  {
    stokes_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioning,
                                               MPI_COMM_WORLD);

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);

    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (! ((c==dim) && (d==dim)))
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                     coupling, sp,
                                     stokes_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    stokes_matrix.reinit (sp);
  }



  template <int dim>
  void Simulator<dim>::
  setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning)
  {
    Amg_preconditioner.reset ();
    Mp_preconditioner.reset ();

    stokes_preconditioner_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp (stokes_partitioning,
                                               MPI_COMM_WORLD);

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (c == d)
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                     coupling, sp,
                                     stokes_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    stokes_preconditioner_matrix.reinit (sp);
  }


  template <int dim>
  void Simulator<dim>::
  setup_temperature_matrix (const IndexSet &temperature_partitioner)
  {
    T_preconditioner.reset ();
    temperature_matrix.clear ();

    TrilinosWrappers::SparsityPattern sp (temperature_partitioner,
                                          MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern (temperature_dof_handler, sp,
                                     temperature_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    temperature_matrix.reinit (sp);
  }



  template <int dim>
  void Simulator<dim>::setup_dofs ()
  {
    computing_timer.enter_section("Setup dof systems");

    stokes_dof_handler.distribute_dofs (stokes_fe);

    // Renumber the DoFs hierarchical so that we get the
    // same numbering if we resume the computation. This
    // is because the numbering depends on the order the
    // cells are created.
    DoFRenumbering::hierarchical (stokes_dof_handler);
    std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
    stokes_sub_blocks[dim] = 1;
    DoFRenumbering::component_wise (stokes_dof_handler, stokes_sub_blocks);

    temperature_dof_handler.distribute_dofs (temperature_fe);

    DoFRenumbering::hierarchical (temperature_dof_handler);
    std::vector<unsigned int> stokes_dofs_per_block (2);
    DoFTools::count_dofs_per_block (stokes_dof_handler, stokes_dofs_per_block,
                                    stokes_sub_blocks);

    const unsigned int n_u = stokes_dofs_per_block[0],
                       n_p = stokes_dofs_per_block[1],
                       n_T = temperature_dof_handler.n_dofs();

    // print dof numbers with 1000s
    // separator since they are frequently
    // large
    std::locale s = pcout.get_stream().getloc();
    // Creating std::locale with an empty string causes problems
    // on some platforms, so catch the exception and ignore
    try
      {
        pcout.get_stream().imbue(std::locale(""));
      }
    catch (std::runtime_error e)
      {
        // If the locale doesn't work, just give up
      }
    pcout << "Number of active cells: "
          << triangulation.n_global_active_cells()
          << " (on "
          << triangulation.n_levels()
          << " levels)"
          << std::endl
          << "Number of degrees of freedom: "
          << n_u + n_p + n_T
          << " (" << n_u << '+' << n_p << '+'<< n_T <<')'
          << std::endl
          << std::endl;
    pcout.get_stream().imbue(s);



    std::vector<IndexSet> stokes_partitioning, stokes_relevant_partitioning;
    IndexSet temperature_partitioning (n_T), temperature_relevant_partitioning (n_T);
    IndexSet stokes_relevant_set;
    {
      IndexSet stokes_index_set = stokes_dof_handler.locally_owned_dofs();
      stokes_partitioning.push_back(stokes_index_set.get_view(0,n_u));
      stokes_partitioning.push_back(stokes_index_set.get_view(n_u,n_u+n_p));

      DoFTools::extract_locally_relevant_dofs (stokes_dof_handler,
                                               stokes_relevant_set);
      stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(0,n_u));
      stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(n_u,n_u+n_p));

      temperature_partitioning = temperature_dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (temperature_dof_handler,
                                               temperature_relevant_partitioning);
    }

    {

      stokes_constraints.clear ();
//    IndexSet stokes_la;
//    DoFTools::extract_locally_active_dofs (stokes_dof_handler,
//             stokes_la);
      stokes_constraints.reinit (stokes_relevant_set);

      DoFTools::make_hanging_node_constraints (stokes_dof_handler,
                                               stokes_constraints);

      /*    std::vector<bool> velocity_mask (dim+1, true);
      velocity_mask[dim] = false;
      VectorTools::interpolate_boundary_values (stokes_dof_handler,
                  0,
                  ZeroFunction<dim>(dim+1),
                  stokes_constraints,
                  velocity_mask);
      */
      std::set<unsigned char> no_normal_flux_boundaries;
      no_normal_flux_boundaries.insert (0);
      no_normal_flux_boundaries.insert (1);
      no_normal_flux_boundaries.insert (2);
      no_normal_flux_boundaries.insert (3);
      VectorTools::compute_no_normal_flux_constraints (stokes_dof_handler, 0,
                                                       no_normal_flux_boundaries,
                                                       stokes_constraints,
                                                       mapping);
      stokes_constraints.close ();
    }
    {
      temperature_constraints.clear ();
      temperature_constraints.reinit (temperature_relevant_partitioning);//temp_locally_active);

      DoFTools::make_hanging_node_constraints (temperature_dof_handler,
                                               temperature_constraints);
      VectorTools::interpolate_boundary_values (temperature_dof_handler,
                                                0,
                                                EquationData::TemperatureInitialValues<dim>(),
                                                temperature_constraints);
      VectorTools::interpolate_boundary_values (temperature_dof_handler,
                                                1,
                                                EquationData::TemperatureInitialValues<dim>(),
                                                temperature_constraints);
      temperature_constraints.close ();
    }

    setup_stokes_matrix (stokes_partitioning);
    setup_stokes_preconditioner (stokes_partitioning);
    setup_temperature_matrix (temperature_partitioning);

    stokes_rhs.reinit (stokes_partitioning, MPI_COMM_WORLD);
    stokes_rhs_helper.reinit (stokes_partitioning, MPI_COMM_WORLD);
    stokes_solution.reinit (stokes_relevant_partitioning, MPI_COMM_WORLD);
    old_stokes_solution.reinit (stokes_solution);

    temperature_rhs.reinit (temperature_partitioning, MPI_COMM_WORLD);
    temperature_solution.reinit (temperature_relevant_partitioning, MPI_COMM_WORLD);
    old_temperature_solution.reinit (temperature_solution);
    old_old_temperature_solution.reinit (temperature_solution);

    rebuild_stokes_matrix              = true;
    rebuild_stokes_preconditioner      = true;

    computing_timer.exit_section();
  }



  template <int dim>
  void
  Simulator<dim>::
  local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                        internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
    const unsigned int   n_q_points      = scratch.stokes_fe_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    scratch.stokes_fe_values.reinit (cell);

    typename DoFHandler<dim>::active_cell_iterator
    temperature_cell (&triangulation,
                      cell->level(),
                      cell->index(),
                      &temperature_dof_handler);
    scratch.temperature_fe_values.reinit (temperature_cell);
    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);
    scratch.stokes_fe_values[pressure].get_function_values(old_stokes_solution,
                                                           scratch.old_pressure_values);

    data.local_matrix = 0;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
//TODO: make the temperature be something useful here
        const double old_temperature = scratch.old_temperature_values[q];
        const double old_pressure = scratch.old_pressure_values[q];

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grads_phi_u[k] = scratch.stokes_fe_values[velocities].symmetric_gradient(k,q);
            scratch.phi_p[k]       = scratch.stokes_fe_values[pressure].value (k, q);
          }

        double eta = material_model->viscosity(old_temperature,
                                               old_pressure,
                                               scratch.stokes_fe_values.quadrature_point(q) );

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            if (stokes_fe.system_to_component_index(i).first
                ==
                stokes_fe.system_to_component_index(j).first)
              data.local_matrix(i,j) += (eta *
                                         (scratch.grads_phi_u[i] *
                                          scratch.grads_phi_u[j])
                                         +
                                         (1./eta) *
                                         pressure_scaling *
                                         pressure_scaling *
                                         (scratch.phi_p[i] * scratch.phi_p[j]))
                                        * scratch.stokes_fe_values.JxW(q);
      }

    cell->get_dof_indices (data.local_dof_indices);
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_preconditioner (const internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    stokes_constraints.distribute_local_to_global (data.local_matrix,
                                                   data.local_dof_indices,
                                                   stokes_preconditioner_matrix);
  }



  template <int dim>
  void
  Simulator<dim>::assemble_stokes_preconditioner ()
  {
    stokes_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_stokes_preconditioner,
                          this,
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_stokes_preconditioner,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         StokesPreconditioner<dim> (stokes_fe, quadrature_formula,
                                    mapping,
                                    update_JxW_values |
                                    update_values |
                                    update_gradients |
                                    update_quadrature_points,
                                    temperature_fe,
                                    update_values),
         internal::Assembly::CopyData::
         StokesPreconditioner<dim> (stokes_fe));

    stokes_preconditioner_matrix.compress();
  }



  template <int dim>
  void
  Simulator<dim>::build_stokes_preconditioner ()
  {
    if (rebuild_stokes_preconditioner == false)
      return;

    computing_timer.enter_section ("   Build Stokes preconditioner");
    pcout << "   Rebuilding Stokes preconditioner..." << std::flush;

    assemble_stokes_preconditioner ();

    std::vector<std::vector<bool> > constant_modes;
    std::vector<bool>  velocity_components (dim+1,true);
    velocity_components[dim] = false;
    DoFTools::extract_constant_modes (stokes_dof_handler, velocity_components,
                                      constant_modes);

    Mp_preconditioner.reset (new TrilinosWrappers::PreconditionILU());
    Amg_preconditioner.reset (new TrilinosWrappers::PreconditionAMG());

    TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = true;
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.02;

    Mp_preconditioner->initialize (stokes_preconditioner_matrix.block(1,1));
    Amg_preconditioner->initialize (stokes_preconditioner_matrix.block(0,0),
                                    Amg_data);

    rebuild_stokes_preconditioner = false;

    pcout << std::endl;
    computing_timer.exit_section();
  }


  template <int dim>
  void
  Simulator<dim>::
  local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                internal::Assembly::Scratch::StokesSystem<dim> &scratch,
                                internal::Assembly::CopyData::StokesSystem<dim> &data)
  {
    const unsigned int dofs_per_cell = scratch.stokes_fe_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.stokes_fe_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    scratch.stokes_fe_values.reinit (cell);

    typename DoFHandler<dim>::active_cell_iterator
    temperature_cell (&triangulation,
                      cell->level(),
                      cell->index(),
                      &temperature_dof_handler);
    scratch.temperature_fe_values.reinit (temperature_cell);
    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);
    scratch.stokes_fe_values[pressure].get_function_values(old_stokes_solution,
                                                           scratch.old_pressure_values);
    scratch.stokes_fe_values[velocities].get_function_values(old_stokes_solution,
                                                             scratch.old_velocity_values);

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;
    data.local_rhs_helper = 0;

    // cache whether the model is compressible or not
    const bool is_compressible = material_model->is_compressible ();


    for (unsigned int q=0; q<n_q_points; ++q)
      {
        const double old_temperature = scratch.old_temperature_values[q];
        const double old_pressure = scratch.old_pressure_values[q];

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_u[k] = scratch.stokes_fe_values[velocities].value (k,q);
            scratch.phi_p[k]       = scratch.stokes_fe_values[pressure].value (k, q);
            if (rebuild_stokes_matrix)
              {
                scratch.grads_phi_u[k] = scratch.stokes_fe_values[velocities].symmetric_gradient(k,q);
                scratch.div_phi_u[k]   = scratch.stokes_fe_values[velocities].divergence (k, q);
              }
          }

        const double eta = material_model->viscosity(old_temperature,
                                                     old_pressure,
                                                     scratch.stokes_fe_values.quadrature_point(q));

        const Tensor<1,dim>
        gravity = EquationData::gravity_vector (scratch.stokes_fe_values.quadrature_point(q));


        const double compressibility = material_model->compressibility(old_temperature,
                                                                       old_pressure,
                                                                       scratch.stokes_fe_values
                                                                       .quadrature_point(q));
        const double density = material_model->density(old_temperature,
                                                       old_pressure,
                                                       scratch.stokes_fe_values
                                                       .quadrature_point(q));

        if (rebuild_stokes_matrix)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              data.local_matrix(i,j) += ( eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j])
                                          - (is_compressible
                                             ?
                                             eta * 2.0/3.0 * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                             :
                                             0)
                                          - (pressure_scaling *
                                             scratch.div_phi_u[i] * scratch.phi_p[j])
                                          - (pressure_scaling *
                                             scratch.phi_p[i] * scratch.div_phi_u[j]))
                                        * scratch.stokes_fe_values.JxW(q);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          data.local_rhs(i) += (  // TODO: extrapolation von old_velocity
                                 (density * gravity * scratch.phi_u[i])
                                 + (is_compressible
                                    ?
                                    (pressure_scaling *
                                     compressibility * density *
                                     (scratch.old_velocity_values[q] * gravity) *
                                     scratch.phi_p[i])
                                    :
                                    0)
                               )
                               * scratch.stokes_fe_values.JxW(q);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          data.local_rhs_helper(i) += scratch.phi_p[i] * scratch.stokes_fe_values.JxW(q);


      }

    cell->get_dof_indices (data.local_dof_indices);
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data)
  {
    if (rebuild_stokes_matrix == true)
      stokes_constraints.distribute_local_to_global (data.local_matrix,
                                                     data.local_rhs,
                                                     data.local_dof_indices,
                                                     stokes_matrix,
                                                     stokes_rhs);
    else
      stokes_constraints.distribute_local_to_global (data.local_rhs,
                                                     data.local_dof_indices,
                                                     stokes_rhs);

    stokes_constraints.distribute_local_to_global (data.local_rhs_helper,
                                                   data.local_dof_indices,
                                                   stokes_rhs_helper);
  }



  template <int dim>
  void Simulator<dim>::assemble_stokes_system ()
  {
    computing_timer.enter_section ("   Assemble Stokes system");

    if (rebuild_stokes_matrix == true)
      stokes_matrix=0;

    stokes_rhs=0;
    stokes_rhs_helper=0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_stokes_system,
                          this,
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_stokes_system,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         StokesSystem<dim> (stokes_fe, mapping, quadrature_formula,
                            (update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             (rebuild_stokes_matrix == true
                              ?
                              update_gradients
                              :
                              UpdateFlags(0))),
                            temperature_fe,
                            update_values),
         internal::Assembly::CopyData::
         StokesSystem<dim> (stokes_fe));

    stokes_matrix.compress();
    stokes_rhs.compress(Add);
    stokes_rhs_helper.compress(Add);

    rebuild_stokes_matrix = false;

    pcout << std::endl;
    computing_timer.exit_section();
  }


  template <int dim>
  void Simulator<dim>::
  local_assemble_temperature_system (const std::pair<double,double> global_T_range,
                                     const double                   global_max_velocity,
                                     const double                   global_entropy_variation,
                                     const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     internal::Assembly::Scratch::TemperatureSystem<dim> &scratch,
                                     internal::Assembly::CopyData::TemperatureSystem<dim> &data)
  {
    const bool use_bdf2_scheme = (timestep_number != 0);



    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    const unsigned int dofs_per_cell = scratch.temperature_fe_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.temperature_fe_values.n_quadrature_points;

    scratch.temperature_fe_values.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);

    data.local_matrix = 0;
    data.local_rhs = 0;

    typename DoFHandler<dim>::active_cell_iterator
    stokes_cell (&triangulation,
                 cell->level(),
                 cell->index(),
                 &stokes_dof_handler);
    scratch.stokes_fe_values.reinit (stokes_cell);

    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);
    scratch.temperature_fe_values.get_function_values (old_old_temperature_solution,
                                                       scratch.old_old_temperature_values);

    scratch.temperature_fe_values.get_function_gradients (old_temperature_solution,
                                                          scratch.old_temperature_grads);
    scratch.temperature_fe_values.get_function_gradients (old_old_temperature_solution,
                                                          scratch.old_old_temperature_grads);

    scratch.temperature_fe_values.get_function_laplacians (old_temperature_solution,
                                                           scratch.old_temperature_laplacians);
    scratch.temperature_fe_values.get_function_laplacians (old_old_temperature_solution,
                                                           scratch.old_old_temperature_laplacians);

    scratch.stokes_fe_values[velocities].get_function_values (stokes_solution,
                                                              scratch.old_velocity_values);
    scratch.stokes_fe_values[velocities].get_function_values (old_stokes_solution,
                                                              scratch.old_old_velocity_values);
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (stokes_solution,
                                                                           scratch.old_strain_rates);
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (old_stokes_solution,
                                                                           scratch.old_old_strain_rates);

    scratch.stokes_fe_values[pressure].get_function_values (stokes_solution,
                                                            scratch.old_pressure);
    scratch.stokes_fe_values[pressure].get_function_values (old_stokes_solution,
                                                            scratch.old_old_pressure);

    const double nu
      = compute_viscosity (scratch.old_temperature_values,
                           scratch.old_old_temperature_values,
                           scratch.old_temperature_grads,
                           scratch.old_old_temperature_grads,
                           scratch.old_temperature_laplacians,
                           scratch.old_old_temperature_laplacians,
                           scratch.old_velocity_values,
                           scratch.old_old_velocity_values,
                           scratch.old_strain_rates,
                           scratch.old_old_strain_rates,
                           scratch.old_pressure,
                           scratch.old_old_pressure,
                           global_max_velocity,
                           global_T_range.second - global_T_range.first,
                           0.5 * (global_T_range.second + global_T_range.first),
                           global_entropy_variation,
                           scratch.temperature_fe_values.get_quadrature_points(),
                           cell->diameter());

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grad_phi_T[k] = scratch.temperature_fe_values.shape_grad (k,q);
            scratch.phi_T[k]      = scratch.temperature_fe_values.shape_value (k, q);
          }

        const double T_term_for_rhs
          = (use_bdf2_scheme ?
             (scratch.old_temperature_values[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_temperature_values[q] *
              (time_step * time_step) /
              (old_time_step * (time_step + old_time_step)))
             :
             scratch.old_temperature_values[q]);

        const double ext_T
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_temperature_values[q], scratch.old_temperature_values[q]);

        const Tensor<1,dim> extrapolated_u
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_velocity_values[q], scratch.old_velocity_values[q]);

        const SymmetricTensor<2,dim> extrapolated_strain_rate
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_strain_rates[q], scratch.old_strain_rates[q]);

        const double ext_pressure
          = aspect::internal::bdf2_extrapolate(
              use_bdf2_scheme, old_time_step, time_step,
              scratch.old_old_pressure[q], scratch.old_pressure[q]);

        const double density
          = material_model->density(ext_T,
                                    ext_pressure,
                                    scratch.temperature_fe_values.quadrature_point(q));
        const double gamma
          = (EquationData::radiogenic_heating * density
             +
             (parameters.include_shear_heating
              ?
              2 * material_model->viscosity(ext_T,
                                            ext_pressure,
                                            scratch.temperature_fe_values.quadrature_point(q)) *
              extrapolated_strain_rate * extrapolated_strain_rate /
              (density * material_model->specific_heat(ext_T, ext_pressure,
                                                       scratch.temperature_fe_values.quadrature_point(q)))
              :
              0)
            );

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            data.local_rhs(i) += (T_term_for_rhs * scratch.phi_T[i]
                                  +
                                  time_step *
                                  gamma * scratch.phi_T[i])
                                 *
                                 scratch.temperature_fe_values.JxW(q);

            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const double factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                          (time_step + old_time_step)) : 1.0;
                data.local_matrix(i,j)
                += (
                     (time_step * (EquationData::kappa+nu) * scratch.grad_phi_T[i] * scratch.grad_phi_T[j])
                     + (time_step * (extrapolated_u * scratch.grad_phi_T[j] * scratch.phi_T[i]))
                     + (factor * scratch.phi_T[i] * scratch.phi_T[j])
                   )
                   * scratch.temperature_fe_values.JxW(q);

              }
          }
      }


  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_temperature_system (const internal::Assembly::CopyData::TemperatureSystem<dim> &data)
  {
    temperature_constraints.distribute_local_to_global (data.local_matrix,
                                                        data.local_rhs,
                                                        data.local_dof_indices,
                                                        temperature_matrix,
                                                        temperature_rhs
                                                       );
  }


  template <int dim>
  void Simulator<dim>::assemble_temperature_system ()
  {
    computing_timer.enter_section ("   Assemble temperature system");
    temperature_matrix = 0;
    temperature_rhs = 0;

    const QGauss<dim> quadrature_formula(parameters.temperature_degree+2);
    const std::pair<double,double>
    global_T_range = get_extrapolated_temperature_range();

    // use midpoint between maximum and minimum
    // temperature for definition of average
    // temperature in entropy viscosity. Could
    // also use the integral average, but the
    // results are not very sensitive to this
    // choice.
    const double average_temperature = 0.5 * (global_T_range.first +
                                              global_T_range.second);
    const double global_entropy_variation =
      get_entropy_variation (average_temperature);

    const double maximal_velocity = get_maximal_velocity();

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.end()),
         std_cxx1x::bind (&Simulator<dim>::
                          local_assemble_temperature_system,
                          this,
                          global_T_range,
                          maximal_velocity,
                          global_entropy_variation,
                          std_cxx1x::_1,
                          std_cxx1x::_2,
                          std_cxx1x::_3),
         std_cxx1x::bind (&Simulator<dim>::
                          copy_local_to_global_temperature_system,
                          this,
                          std_cxx1x::_1),
         internal::Assembly::Scratch::
         TemperatureSystem<dim> (temperature_fe, stokes_fe, mapping, quadrature_formula),
         internal::Assembly::CopyData::
         TemperatureSystem<dim> (temperature_fe));

    temperature_matrix.compress();
    temperature_rhs.compress(Add);

    T_preconditioner.reset (new TrilinosWrappers::PreconditionILU());
    T_preconditioner->initialize (temperature_matrix);

    computing_timer.exit_section();
  }



  template <int dim>
  void Simulator<dim>::solve ()
  {
    computing_timer.enter_section ("   Solve Stokes system");

    {
      pcout << "   Solving Stokes system... " << std::flush;

      TrilinosWrappers::MPI::BlockVector
      distributed_stokes_solution (stokes_rhs);
      distributed_stokes_solution = stokes_solution;

      // before solving we scale the initial solution to the right dimensions
      distributed_stokes_solution.block(1) /= pressure_scaling;

      const unsigned int
      start = (distributed_stokes_solution.block(0).size() +
               distributed_stokes_solution.block(1).local_range().first),
              end   = (distributed_stokes_solution.block(0).size() +
                       distributed_stokes_solution.block(1).local_range().second);
      for (unsigned int i=start; i<end; ++i)
        if (stokes_constraints.is_constrained (i))
          distributed_stokes_solution(i) = 0;

      // if the model is compressible then we need to adjust the right hand
      // side of the equation to make it compatible with the matrix on the
      // left
      if (material_model->is_compressible ())
        make_pressure_rhs_compatible(stokes_rhs, stokes_rhs_helper);

      PrimitiveVectorMemory< TrilinosWrappers::MPI::BlockVector > mem;

      // step 1: try if the simple and fast solver
      // succeeds in 30 steps or less.
      const double solver_tolerance = 1e-7 * stokes_rhs.l2_norm();
      SolverControl solver_control_cheap (30, solver_tolerance);
      SolverControl solver_control_expensive (stokes_matrix.m(), solver_tolerance);

      try
        {
          const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                TrilinosWrappers::PreconditionILU>
                preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                                *Mp_preconditioner, *Amg_preconditioner,
                                false);

          SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
          solver(solver_control_cheap, mem,
                 SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
                 AdditionalData(30, true));
          solver.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
                       preconditioner);
        }

      // step 2: take the stronger solver in case
      // the simple solver failed
      catch (SolverControl::NoConvergence)
        {
          const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                TrilinosWrappers::PreconditionILU>
                preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                                *Mp_preconditioner, *Amg_preconditioner,
                                true);

          SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
          solver(solver_control_expensive, mem,
                 SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
                 AdditionalData(50, true));
          solver.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
                       preconditioner);
        }


      stokes_constraints.distribute (distributed_stokes_solution);

      // now rescale the pressure back to real physical units
      distributed_stokes_solution.block(1) *= pressure_scaling;

      stokes_solution = distributed_stokes_solution;

      normalize_pressure(stokes_solution);

      // print the number of iterations to screen and record it in the
      // statistics file
      if (solver_control_expensive.last_step() == 0)
        pcout << solver_control_cheap.last_step()  << " iterations.";
      else
        pcout << solver_control_cheap.last_step() << '+'
              << solver_control_expensive.last_step() << " iterations.";
      pcout << std::endl;

      statistics.add_value("Iterations for Stokes solver",
                           solver_control_cheap.last_step() + solver_control_expensive.last_step());
    }
    computing_timer.exit_section();


    {
      old_time_step = time_step;
      const double cfl_number = get_cfl_number();

      // we found out that we need
      // approximately a quarter the time step
      // size in 3d
      double scaling = (dim==3)?0.25:1.0;
      time_step = (scaling/(2.1*dim*std::sqrt(1.*dim)) /
                   (parameters.temperature_degree *
                    cfl_number));

      statistics.add_value("Time step size (year)", time_step / EquationData::year_in_seconds);

      temperature_solution = old_temperature_solution;
      assemble_temperature_system ();
    }

    computing_timer.enter_section ("   Solve temperature system");
    {
      pcout << "   Solving temperature system... " << std::flush;

      SolverControl solver_control (temperature_matrix.m(),
                                    1e-12*temperature_rhs.l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector>   solver (solver_control,
                                                           SolverGMRES<TrilinosWrappers::MPI::Vector>::AdditionalData(30,true));

      TrilinosWrappers::MPI::Vector
      distributed_temperature_solution (temperature_rhs);
      distributed_temperature_solution = temperature_solution;

      solver.solve (temperature_matrix, distributed_temperature_solution,
                    temperature_rhs, *T_preconditioner);

      temperature_constraints.distribute (distributed_temperature_solution);
      temperature_solution = distributed_temperature_solution;

      // print number of iterations and also record it in the
      // statistics file
      pcout << solver_control.last_step()
            << " iterations." << std::endl;

      statistics.add_value("Iterations for temperature solver",
                           solver_control.last_step());
    }
    computing_timer.exit_section();
  }

  /*
   * normalize the pressure by calculating the surface integral of the pressure on the outer
   * shell and subtracting this from all pressure nodes.
   */
  template <int dim>
  void Simulator<dim>::normalize_pressure(TrilinosWrappers::MPI::BlockVector &vector)
  {
    double my_pressure = 0.0;
    double my_area = 0.0;
    {
      QGauss < dim - 1 > quadrature (parameters.stokes_velocity_degree + 1);

      const unsigned int n_q_points = quadrature.size();
      FEFaceValues<dim> fe_face_values (mapping, stokes_fe,  quadrature,
                                        update_JxW_values | update_values);
      const FEValuesExtractors::Scalar pressure (dim);

      std::vector<double> pressure_values(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator
      cell = stokes_dof_handler.begin_active(),
      endc = stokes_dof_handler.end();
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
              {
                const typename DoFHandler<dim>::face_iterator face = cell->face (face_no);
                if (face->at_boundary() && face->boundary_indicator() == 1) // outer shell boundary
                  {
                    fe_face_values.reinit (cell, face_no);
                    fe_face_values[pressure].get_function_values(vector,
                                                                 pressure_values);

                    for (unsigned int q = 0; q < n_q_points; ++q)
                      {
                        my_pressure += pressure_values[q]
                                       * fe_face_values.JxW (q);
                        my_area += fe_face_values.JxW (q);
                      }
                  }
              }
          }
    }

    double surf_pressure = 0;
    // sum up the surface integrals from each processor
    {
      const double my_temp[2] = {my_pressure, my_area};
      double temp[2];
      Utilities::MPI::sum (my_temp, MPI_COMM_WORLD, temp);
      surf_pressure = temp[0]/temp[1];
    }

    const double adjust = -surf_pressure + 1e7;
    if (parameters.use_locally_conservative_discretization == false)
      vector.block(1).add(adjust);
    else
      {
        // this case is a bit more complicated: if the condition above is false
        // then we use the FE_DGP element for which the shape functions do not
        // add up to one; consequently, adding a constant to all degrees of
        // freedom does not alter the overall function by that constant, but
        // by something different
        //
        // we can work around this by using the documented property of the
        // FE_DGP element that the first shape function is constant.
        // consequently, adding the adjustment to the global function is
        // achieved by adding the adjustment to the first pressure degree
        // of freedom on each cell.
        Assert (dynamic_cast<const FE_DGP<dim>*>(&stokes_fe.base_element(1)) != 0,
                ExcInternalError());
        std::vector<unsigned int> local_dof_indices (stokes_fe.dofs_per_cell);
        typename DoFHandler<dim>::active_cell_iterator
        cell = stokes_dof_handler.begin_active(),
        endc = stokes_dof_handler.end();
        for (; cell != endc; ++cell)
          if (cell->is_locally_owned())
            {
              // identify the first pressure dof
              cell->get_dof_indices (local_dof_indices);
              const unsigned int first_pressure_dof
                = stokes_fe.component_to_system_index (dim, 0);

              // make sure that this DoF is really owned by the current processor
              // and that it is in fact a pressure dof
              Assert (stokes_dof_handler.locally_owned_dofs().is_element(first_pressure_dof),
                      ExcInternalError());
              Assert (local_dof_indices[first_pressure_dof] >= vector.block(0).size(),
                      ExcInternalError());

              // then adjust its value
              vector(local_dof_indices[first_pressure_dof]) += adjust;
            }
      }
  }



// This routine adjusts the second block of the right hand side of the
// system containing the compressibility, so that the system becomes
// compatible: 0=\int div u = \int g
// the helper vector h contains h_i=(q_i,1) with the pressure functions q_i
// and we adjust the right hand side g by h_i \int g / |\Omega|
  template <int dim>
  void Simulator<dim>::make_pressure_rhs_compatible(TrilinosWrappers::MPI::BlockVector &vector,
                                                    const TrilinosWrappers::MPI::BlockVector &helper)
  {
    if (parameters.use_locally_conservative_discretization)
      throw ExcNotImplemented();

    double mean = vector.block(1).mean_value();
    double correct = -mean*vector.block(1).size()/global_volume;

//  pcout << "    pressure correction: " << correct << std::endl;
    vector.block(1).add(correct, helper.block(1));
  }



  template <int dim>
  void Simulator<dim>::create_snapshot()
  {
    unsigned int myid = Utilities::System::get_this_mpi_process (MPI_COMM_WORLD);

    if (myid == 0)
      {
        // keep the last snapshot in case this one fails to save
        system ("mv bin/mesh bin/mesh.old");
        system ("mv bin/mesh.info bin/mesh.info.old");
        system ("mv bin/resume.txt bin/resume.txt.old");
      }

    //save Triangulation and Solution vectors:
    {
      std::vector<const TrilinosWrappers::MPI::Vector *> x_temperature (3);
      x_temperature[0] = &temperature_solution;
      x_temperature[1] = &old_temperature_solution;
      x_temperature[2] = &old_old_temperature_solution;
      std::vector<const TrilinosWrappers::MPI::BlockVector *> x_stokes (2);
      x_stokes[0] = &stokes_solution;
      x_stokes[1] = &old_stokes_solution;

      parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
      temperature_trans (temperature_dof_handler);
      parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector>
      stokes_trans (stokes_dof_handler);

      temperature_trans.prepare_serialization (x_temperature);
      stokes_trans.prepare_serialization (x_stokes);

      const char *filename = "bin/mesh";
      triangulation.save (filename);
    }

    //save general information
    if (myid == 0)
      {
        std::ofstream ofs ("bin/resume.txt");
        boost::archive::text_oarchive oa (ofs);
        oa << (*this);
      }
    pcout << "*** Snapshot created!" << std::endl;
  }

  template <int dim>
  void Simulator<dim>::resume_from_snapshot()
  {
    triangulation.load ("bin/mesh");
    global_volume = GridTools::volume (triangulation, mapping);
    setup_dofs();

    TrilinosWrappers::MPI::Vector
    distributed_temp1 (temperature_rhs);
    TrilinosWrappers::MPI::Vector
    distributed_temp2 (temperature_rhs);
    TrilinosWrappers::MPI::Vector
    distributed_temp3 (temperature_rhs);

    std::vector<TrilinosWrappers::MPI::Vector *> x_temperature (3);
    x_temperature[0] = & (distributed_temp1);
    x_temperature[1] = & (distributed_temp2);
    x_temperature[2] = & (distributed_temp3);

    TrilinosWrappers::MPI::BlockVector
    distributed_stokes (stokes_rhs);
    TrilinosWrappers::MPI::BlockVector
    old_distributed_stokes (stokes_rhs);
    std::vector<TrilinosWrappers::MPI::BlockVector *> x_stokes (2);
    x_stokes[0] = & (distributed_stokes);
    x_stokes[1] = & (old_distributed_stokes);

    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector>
    temperature_trans (temperature_dof_handler);
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector>
    stokes_trans (stokes_dof_handler);

    temperature_trans.deserialize (x_temperature);
    stokes_trans.deserialize (x_stokes);


    temperature_solution = distributed_temp1;
    old_temperature_solution = distributed_temp2;
    old_old_temperature_solution = distributed_temp3;

    stokes_solution = distributed_stokes;
    old_stokes_solution = old_distributed_stokes;

    std::ifstream ifs ("bin/resume.txt");
    boost::archive::text_iarchive ia (ifs);
    ia >> (*this);

    // re-initialize the postprocessors with the current object
    postprocess_manager.initialize (*this);

    pcout << "*** resuming from Snapshot!" << std::endl;
  }

}

//why do we need this?!
BOOST_CLASS_TRACKING (aspect::Simulator<2>, boost::serialization::track_never)
BOOST_CLASS_TRACKING (aspect::Simulator<3>, boost::serialization::track_never)


namespace aspect
{

  template <int dim>
  template<class Archive>
  void Simulator<dim>::serialize (Archive &ar, const unsigned int)
  {
    ar &time;
    ar &time_step;
    ar &old_time_step;
    ar &timestep_number;

    ar &postprocess_manager &statistics;

// how about global_volume, global_Omega_diameter
  }


  template <int dim>
  void Simulator<dim>::postprocess ()
  {
    computing_timer.enter_section ("Postprocessing");
    pcout << "   Postprocessing:" << std::endl;

    // run all the postprocessing routines and then write
    // the current state of the statistics table to a file
    std::list<std::pair<std::string,std::string> >
    output_list = postprocess_manager.execute (statistics);

    std::ofstream stat_file ("bin/statistics");
    statistics.set_scientific("Time (years)", true);
    statistics.set_scientific("Time step size (year)", true);
    statistics.write_text (stat_file,
                           TableHandler::table_with_separate_column_description);

    // determine the width of the first column of text so that
    // everything gets nicely aligned; then output everything
    {
      unsigned int width = 0;
      for (std::list<std::pair<std::string,std::string> >::const_iterator
           p = output_list.begin();
           p != output_list.end(); ++p)
        width = std::max<unsigned int> (width, p->first.size());

      for (std::list<std::pair<std::string,std::string> >::const_iterator
           p = output_list.begin();
           p != output_list.end(); ++p)
        pcout << "     "
              << std::left
              << std::setw(width)
              << p->first
              << " "
              << p->second
              << std::endl;
    }

    pcout << std::endl;
    computing_timer.exit_section ();
  }



// Contrary to step-32, we have found that just refining by the temperature
// works well in 2d, but only leads to refinement in the boundary
// layer at the core-mantle boundary in 3d. consequently, we estimate
// the error based both on the temperature and on the velocity; the
// vectors with the resulting error indicators are then both normalized
// to a maximal value of one, and we take the maximum of the two indicators
// to decide whether we want to refine or not. this ensures that we
// also refine into plumes where maybe the temperature gradients aren't
// as strong as in the boundary layer but where nevertheless the gradients
// in the velocity are large
  template <int dim>
  void Simulator<dim>::refine_mesh (const unsigned int max_grid_level)
  {
    computing_timer.enter_section ("Refine mesh structure, part 1");

    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    // compute the errors for
    // temperature and stokes solution,
    // then scale them and find the
    // maximum between the two
    {
      Vector<float> estimated_error_per_cell_T (triangulation.n_active_cells());

      KellyErrorEstimator<dim>::estimate (temperature_dof_handler,
                                          QGauss<dim-1>(parameters.temperature_degree+1),
                                          typename FunctionMap<dim>::type(),
                                          temperature_solution,
                                          estimated_error_per_cell_T,
                                          std::vector<bool>(),
                                          0,
                                          0,
                                          triangulation.locally_owned_subdomain());
      estimated_error_per_cell_T /= Utilities::MPI::max (estimated_error_per_cell_T.linfty_norm(),
                                                         MPI_COMM_WORLD);

      Vector<float> estimated_error_per_cell_u (triangulation.n_active_cells());
      std::vector<bool> velocity_mask (dim+1, true);
      velocity_mask[dim] = false;
      KellyErrorEstimator<dim>::estimate (stokes_dof_handler,
                                          QGauss<dim-1>(parameters.stokes_velocity_degree+1),
                                          typename FunctionMap<dim>::type(),
                                          stokes_solution,
                                          estimated_error_per_cell_u,
                                          velocity_mask,
                                          0,
                                          0,
                                          triangulation.locally_owned_subdomain());
      estimated_error_per_cell_u /= Utilities::MPI::max (estimated_error_per_cell_u.linfty_norm(),
                                                         MPI_COMM_WORLD);

      for (unsigned int i=0; i<estimated_error_per_cell.size(); ++i)
        estimated_error_per_cell(i) = std::max (estimated_error_per_cell_T(i),
                                                estimated_error_per_cell_u(i));
    }

    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_fraction (triangulation,
                                       estimated_error_per_cell,
                                       parameters.refinement_fraction,
                                       parameters.coarsening_fraction);

    // limit maximum refinement level
    if (triangulation.n_levels() > max_grid_level)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active(max_grid_level);
           cell != triangulation.end(); ++cell)
        cell->clear_refine_flag ();

    std::vector<const TrilinosWrappers::MPI::Vector *> x_temperature (2);
    x_temperature[0] = &temperature_solution;
    x_temperature[1] = &old_temperature_solution;
    std::vector<const TrilinosWrappers::MPI::BlockVector *> x_stokes (2);
    x_stokes[0] = &stokes_solution;
    x_stokes[1] = &old_stokes_solution;

    parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::Vector>
    temperature_trans(temperature_dof_handler);
    parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::BlockVector>
    stokes_trans(stokes_dof_handler);

    triangulation.prepare_coarsening_and_refinement();
    temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
    stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);

    triangulation.execute_coarsening_and_refinement ();
    global_volume = GridTools::volume (triangulation, mapping);
    computing_timer.exit_section();

    setup_dofs ();

    computing_timer.enter_section ("Refine mesh structure, part 2");

    {
      TrilinosWrappers::MPI::Vector
      distributed_temp1 (temperature_rhs);
      TrilinosWrappers::MPI::Vector
      distributed_temp2 (temperature_rhs);

      std::vector<TrilinosWrappers::MPI::Vector *> tmp (2);
      tmp[0] = &(distributed_temp1);
      tmp[1] = &(distributed_temp2);
      temperature_trans.interpolate(tmp);

      temperature_solution     = distributed_temp1;
      old_temperature_solution = distributed_temp2;
    }

    {
      TrilinosWrappers::MPI::BlockVector
      distributed_stokes (stokes_rhs);
      TrilinosWrappers::MPI::BlockVector
      old_distributed_stokes (stokes_rhs);
      std::vector<TrilinosWrappers::MPI::BlockVector *> stokes_tmp (2);
      stokes_tmp[0] = &(distributed_stokes);
      stokes_tmp[1] = &(old_distributed_stokes);

      stokes_trans.interpolate (stokes_tmp);
      stokes_solution     = distributed_stokes;
      old_stokes_solution = old_distributed_stokes;
    }

    computing_timer.exit_section();
  }



  template <int dim>
  void Simulator<dim>::run ()
  {
    if (EquationData::apperture_angle == std::acos(-1e0)*2e0)
      {

        GridGenerator::hyper_shell (triangulation,
                                    Point<dim>(),
                                    EquationData::R0,
                                    EquationData::R1,
                                    (dim==3) ? 96 : 12,
                                    true);
      }
    else if (EquationData::apperture_angle == std::acos(-1e0)/2e0)
      {
        GridGenerator::quarter_hyper_shell (triangulation,
                                            Point<dim>(),
                                            EquationData::R0,
                                            EquationData::R1,0,
                                            true);
      }
    else if (EquationData::apperture_angle == std::acos(-1e0))
      {
        GridGenerator::half_hyper_shell (triangulation,
                                         Point<dim>(),
                                         EquationData::R0,
                                         EquationData::R1,0,
                                         true);
      }
    else
      {
        Assert (false, ExcInternalError());
      }

    static HyperShellBoundary<dim> boundary;
    triangulation.set_boundary (0, boundary);
    triangulation.set_boundary (1, boundary);

    global_Omega_diameter = GridTools::diameter (triangulation);

    if (parameters.resume_computation == true)
      {
        resume_from_snapshot();
      }
    else
      {
        triangulation.refine_global (parameters.initial_global_refinement);
        global_volume = GridTools::volume (triangulation, mapping);

        setup_dofs();
      }

    unsigned int max_refinement_level = parameters.initial_global_refinement +
                                        parameters.initial_adaptive_refinement;

    unsigned int pre_refinement_step = 0;

  start_time_iteration:

    if (parameters.resume_computation == false)
      {
        set_initial_temperature_field ();
        compute_initial_pressure_field ();

        time                      = 0;
        timestep_number           = 0;
        time_step = old_time_step = 0;
      }

    do
      {
        pcout << "*** Timestep " << timestep_number
              << ":  t=" << time/EquationData::year_in_seconds
              << " years"
              << std::endl;

        // set global statistics about this time step
        statistics.add_value("Time step number", timestep_number);
        statistics.add_value("Time (years)", time / EquationData::year_in_seconds);

        assemble_stokes_system ();
        build_stokes_preconditioner ();

        solve ();

        pcout << std::endl;

        // see if we have to start over with a new refinement cycle
        // at the beginning of the simulation
        if ((timestep_number == 0) &&
            (pre_refinement_step < parameters.initial_adaptive_refinement))
          {
            refine_mesh (max_refinement_level);
            ++pre_refinement_step;
            goto start_time_iteration;
          }

        postprocess ();

        // see if this is a time step where additional refinement is requested
        // if so, then loop over as many times as this is necessary
        if ((parameters.additional_refinement_times.size() > 0)
            &&
            (parameters.additional_refinement_times.front () < time+time_step))
          {
            while ((parameters.additional_refinement_times.size() > 0)
                   &&
                   (parameters.additional_refinement_times.front () < time+time_step))
              {
                ++max_refinement_level;
                refine_mesh (max_refinement_level);

                parameters.additional_refinement_times
                .erase (parameters.additional_refinement_times.begin());
              }
          }
        else
          // see if this is a time step where regular refinement is necessary, but only
          // if the previous rule wasn't triggered
          if ((timestep_number > 0)
              &&
              (timestep_number % parameters.adaptive_refinement_interval == 0))
            refine_mesh (max_refinement_level);

        // prepare for the next time
        // step
        TrilinosWrappers::MPI::BlockVector old_old_stokes_solution;
        old_old_stokes_solution      = old_stokes_solution;
        old_stokes_solution          = stokes_solution;
        old_old_temperature_solution = old_temperature_solution;
        old_temperature_solution     = temperature_solution;
        if (old_time_step > 0)
          {
            stokes_solution.sadd (1.+time_step/old_time_step, -time_step/old_time_step,
                                  old_old_stokes_solution);
            temperature_solution.sadd (1.+time_step/old_time_step,
                                       -time_step/old_time_step,
                                       old_old_temperature_solution);
          }

        // every 100 time steps output
        // a summary of the current
        // timing information
        if ((timestep_number > 0) && (timestep_number % 100 == 0))
          computing_timer.print_summary ();

        time += time_step;
        ++timestep_number;

        if (timestep_number % 50 == 0)
          {
            create_snapshot();
            // matrices will be regenerated after a resume, so do that here too
            // to be consistent.
            rebuild_stokes_matrix =
              rebuild_stokes_preconditioner = true;
          }

        // if we are at the end of
        // time, stop now
        if (time > parameters.end_time * EquationData::year_in_seconds)
          break;
      }
    while (true);
  }


}



int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  try
    {
      deallog.depth_console (0);

      // see which parameter file to use
      std::string parameter_filename;
      if (argc>=2)
        parameter_filename = argv[1];
      else
        parameter_filename = "aspect.prm";

      // declare parameters so that we can create a default file
      // if there is no parameter file
      ParameterHandler prm;
      aspect::Simulator<deal_II_dimension>::declare_parameters(prm);

      std::ifstream parameter_file (parameter_filename.c_str());
      if (!parameter_file)
        {
          parameter_file.close ();

          std::ostringstream message;
          message << "Input parameter file <"
                  << parameter_filename << "> not found. Creating a"
                  << std::endl
                  << "template file of the same name."
                  << std::endl;

          std::ofstream parameter_out (parameter_filename.c_str());
          prm.print_parameters (parameter_out,
                                ParameterHandler::Text);

          AssertThrow (false, ExcMessage (message.str().c_str()));
        }

      const bool success = prm.read_input (parameter_file);
      AssertThrow (success, ExcMessage ("Invalid input parameter file."));

      aspect::Simulator<deal_II_dimension> flow_problem (prm);
      flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
