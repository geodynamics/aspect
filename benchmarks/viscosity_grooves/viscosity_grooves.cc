/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
#include <aspect/material_model/simple.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace ViscosityGroovesBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {


      Tensor<1,2>
      ViscosityGrooves_velocity (const Point<2> &pos)
      {
        const double x = pos[0];
        const double y = pos[1];
        const double v_x = x*x*x*y+x*x+x*y+x;
        const double v_y = -1.5*x*x*y*y-2*x*y-0.5*y*y-y;
        return Point<2> (v_x,v_y);
      }

      double
      ViscosityGrooves_pressure (const Point<2> &pos,
                                 const GeometryModel::Interface<2> &geometry_model)
      {
        const double x = pos[0];
        const double y = pos[1];
        const GeometryModel::Box<2> *geometry
          = dynamic_cast<const GeometryModel::Box<2>*> (&geometry_model);
        const double L=geometry->get_extents()[0];
        return x*x*y*y+x*y+5. - pow(L,4.)/9.-pow(L,2.)/4.-5.;
      }

      double
      ViscosityGrooves_pressure (const Point<2> &, const GeometryModel::Interface<3> &)
      {
        Assert (false, ExcNotImplemented());
        return 0.;
      }

      template <int dim>
      class FunctionViscosityGrooves : public Function<dim>
      {
        public:
          FunctionViscosityGrooves (const GeometryModel::Interface<dim> &geometry_model)
            :
            Function<dim>(dim+2),
            geometry_model (geometry_model)
          {}

          virtual void vector_value (const Point< dim >   &pos,
                                     Vector< double >   &values) const
          {
            Assert (dim == 2, ExcNotImplemented());
            Assert (values.size() >= 4, ExcInternalError());

            const Point<2> p (pos[0], pos[1]);

            const Tensor<1,2> v = AnalyticSolutions::ViscosityGrooves_velocity (p);
            values[0] = v[0];
            values[1] = v[1];
            values[2] = AnalyticSolutions::ViscosityGrooves_pressure (p, geometry_model);
          }

        private:
          const GeometryModel::Interface<dim> &geometry_model;
      };
    }


    template <int dim>
    class ViscosityGroovesBoundary : public BoundaryVelocity::Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        ViscosityGroovesBoundary();

        /**
         * Return the boundary velocity as a function of position.
         */
        virtual
        Tensor<1,dim>
        boundary_velocity (const types::boundary_id ,
                           const Point<dim> &position) const;

      private:
        const double epsilon;

    };


    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class ViscosityGroovesMaterial : public MaterialModel::Interface<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const Point<dim> &pos = in.position[i];
              const double x = pos[0];
              const double y = pos[1];
              out.viscosities[i] = -std::sin(x*x*y*y+x*y+5.)+1.+epsilon;
              out.densities[i] = 1.;
              out.compressibilities[i] = 0.;
              out.specific_heat[i] = 0.;
              out.thermal_expansion_coefficients[i] = 0.;
              out.thermal_conductivities[i] = 0.;
            }
        }
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the continuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;
        /**
         * @}
         */

        double get_epsilon() const;

        double epsilon;

    };


    template <int dim>
    double
    ViscosityGroovesMaterial<dim>::
    reference_viscosity () const
    {
      return 1.;
    }


    template <int dim>
    bool
    ViscosityGroovesMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    ViscosityGroovesMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      //create a global section in the parameter file for parameters
      //that describe this benchmark. note that we declare them here
      //in the material model, but other kinds of plugins (e.g., the gravity
      //model below) may also read these parameters even though they do not
      //declare them
      prm.enter_subsection("ViscosityGrooves benchmark");
      {
        prm.declare_entry("Viscosity parameter", "0.1",
                          Patterns::Double (0.),
                          "Viscosity parameter epsilon in the ViscosityGrooves benchmark.");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    ViscosityGroovesMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("ViscosityGrooves benchmark");
      {
        epsilon = prm.get_double ("Viscosity parameter");
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
      this->model_dependence.density = MaterialModel::NonlinearDependence::none;
      this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
      this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
    }

    template <int dim>
    double
    ViscosityGroovesMaterial<dim>::get_epsilon() const
    {
      return epsilon;
    }

    template <int dim>
    ViscosityGroovesBoundary<dim>::ViscosityGroovesBoundary ()
      :
      epsilon (0.)
    {}

    template <>
    Tensor<1,2>
    ViscosityGroovesBoundary<2>::
    boundary_velocity (const types::boundary_id ,
                       const Point<2> &p) const
    {
      return AnalyticSolutions::ViscosityGrooves_velocity (p);
    }


    template <>
    Tensor<1,3>
    ViscosityGroovesBoundary<3>::
    boundary_velocity (const types::boundary_id ,
                       const Point<3> &) const
    {
      Assert (false, ExcNotImplemented());
      return Tensor<1,3>();
    }

    /**
     * Gravity model for the ViscosityGrooves benchmark
    */

    template <int dim>
    class ViscosityGroovesGravity : public aspect::GravityModel::Interface<dim>
    {
      public:
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &pos) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
        */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        double epsilon;

    };

    template <int dim>
    Tensor<1,dim>
    ViscosityGroovesGravity<dim>::
    gravity_vector(const Point<dim> &pos) const
    {
      const double x=pos[0];
      const double y=pos[1];
      Tensor<1,dim> g;

      const double eta=-std::sin(x*x*y*y+x*y+5.)+1.+epsilon;

      const double costerm=std::cos(x*x*y*y+x*y+5.);

      const double deta_dx=-y*(2.*x*y+1.)*costerm;
      const double deta_dy=-x*(2.*x*y+1.)*costerm;

      const double dpdx=2. * x  *y*y +y ;
      const double dpdy=2. * x*x*y   +x ;

      const double exx= 3.*x*x * y +2.*x +y +1.;
      const double eyy=-3.*x*x * y -2.*x -y -1.;

      const double exy=0.5*(x*x*x + x -3.*x*y*y -2.*y);
      const double eyx=0.5*(x*x*x + x -3.*x*y*y -2.*y);

      const double dexxdx= 6.*x*y+2.;
      const double deyxdy=-3.*x*y-1.;

      const double dexydx= 0.5*(3.*x*x +1. -3.*y*y);
      const double deyydy= -3.*x*x -1.;

      const double gx =-dpdx + 2.*eta*dexxdx + 2.*deta_dx*exx + 2.*eta*deyxdy + 2.*deta_dy*eyx;
      const double gy =-dpdy + 2.*eta*dexydx + 2.*deta_dx*exy + 2.*eta*deyydy + 2.*deta_dy*eyy;

      g[0]=-gx;
      g[1]=-gy;

      return g;
    }

    template <int dim>
    void
    ViscosityGroovesGravity<dim>::declare_parameters (ParameterHandler &)
    {
      //nothing to declare here. This plugin will however, read parameters
      //declared by the material model in the "ViscosityGrooves benchmark" section
    }

    template <int dim>
    void
    ViscosityGroovesGravity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("ViscosityGrooves benchmark");
      {
        epsilon = prm.get_double ("Viscosity parameter");
      }
      prm.leave_subsection();
    }


    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the Donea and Huerta FEM book (see manual).
      */
    template <int dim>
    class ViscosityGroovesPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };

    template <int dim>
    std::pair<std::string,std::string>
    ViscosityGroovesPostprocessor<dim>::execute (TableHandler &)
    {
      std::shared_ptr<Function<dim> > ref_func;
      {

        ref_func.reset (new AnalyticSolutions::FunctionViscosityGrooves<dim>(this->get_geometry_model()));
      }

      const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities+2);

      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                          dim+2);
      ComponentSelectFunction<dim> comp_p(dim, dim+2);

      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_u,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_p);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_ul2,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_pl2,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p);

      const double u_l1 =  VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u, VectorTools::L1_norm);
      const double p_l1 =  VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L1_norm);
      const double u_l2 =  VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_ul2, VectorTools::L2_norm);
      const double p_l2 =  VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_pl2, VectorTools::L2_norm);

      std::ostringstream os;

      os << std::scientific <<  u_l1
         << ", " << p_l1
         << ", " << u_l2
         << ", " << p_l2;

      return std::make_pair("Errors u_L1, p_L1, u_L2, p_L2:", os.str());
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace ViscosityGroovesBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscosityGroovesMaterial,
                                   "ViscosityGroovesMaterial",
                                   "A material model that corresponds to the `ViscosityGrooves' benchmark. "
                                   "See the manual for more information.")

    ASPECT_REGISTER_BOUNDARY_VELOCITY_MODEL(ViscosityGroovesBoundary,
                                            "ViscosityGroovesBoundary",
                                            "Implementation of the velocity boundary conditions for the "
                                            "`ViscosityGrooves' benchmark. See the manual for more information about this "
                                            "benchmark.")

    ASPECT_REGISTER_POSTPROCESSOR(ViscosityGroovesPostprocessor,
                                  "ViscosityGroovesPostprocessor",
                                  "A postprocessor that compares the solution of the `ViscosityGrooves' benchmark "
                                  "with the one computed by ASPECT "
                                  "and reports the error. See the manual for more information.")

    ASPECT_REGISTER_GRAVITY_MODEL(ViscosityGroovesGravity,
                                  "ViscosityGroovesGravity",
                                  "A gravity model in corresponding to the `ViscosityGrooves' benchmark. "
                                  "See the manual for more information.")
  }
}

