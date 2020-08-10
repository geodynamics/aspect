#include <aspect/material_model/interface.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/mesh_refinement/interface.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/derivative_approximation.h>


namespace aspect
{
  template <int dim>
  class MeltVelocityRefinement : public MeshRefinement::Interface<dim>,
    public SimulatorAccess<dim>
  {
    public:
      /**
       * Execute this mesh refinement criterion.
       *
       * @param[out] error_indicators A vector that for every active cell of
       * the current mesh (which may be a partition of a distributed mesh)
       * provides an error indicator. This vector will already have the
       * correct size when the function is called.
       */
      virtual
      void
      execute (Vector<float> &indicators) const
      {
        indicators = 0;

        KellyErrorEstimator<dim>::estimate (this->get_mapping(),
                                            this->get_dof_handler(),
                                            QGauss<dim-1>(this->introspection().polynomial_degree.velocities +1),
                                            std::map<types::boundary_id,const Function<dim>*>(),
                                            this->get_solution(),
                                            indicators,
                                            this->introspection().variable("fluid velocity").component_mask,
                                            0,
                                            0,
                                            this->get_triangulation().locally_owned_subdomain());
      }

  };

  template <int dim>
  class PCRefinement : public MeshRefinement::Interface<dim>,
    public SimulatorAccess<dim>
  {
    public:
      /**
       * Execute this mesh refinement criterion.
       *
       * @param[out] error_indicators A vector that for every active cell of
       * the current mesh (which may be a partition of a distributed mesh)
       * provides an error indicator. This vector will already have the
       * correct size when the function is called.
       */
      virtual
      void
      execute (Vector<float> &indicators) const
      {
        indicators = 0;

        // create a vector in which we set the compaction pressure block
        // (that contains the scaled compaction pressure p_c_bar) to
        // be a finite element interpolation of the (real) compaction pressure.
        // we do so by setting up a quadrature formula and looping over these
        // points, rescaling the compaction pressure, and writing
        // the result into the output vector. As the scaling factor between
        // the pressures is cell-wise constant (and the error indicator is
        // a cellwise quantity too), we use a quadrature formula with only
        // one point.
        LinearAlgebra::BlockVector vec_distributed (this->introspection().index_sets.system_partitioning,
                                                    this->get_mpi_communicator());

        const unsigned int pc_component_index = this->introspection().variable("compaction pressure").first_component_index;
        const unsigned int pc_base_index = this->introspection().variable("compaction pressure").base_index;

        // Use a quadrature formula with only one point, as the scaling factor we need is cell-wise constant
        const QMidpoint<dim> quadrature;

        std::vector<types::global_dof_index> local_dof_indices (this->get_fe().dofs_per_cell);
        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature,
                                 update_quadrature_points | update_values | update_gradients);

        MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
        MeltHandler<dim>::create_material_model_outputs(out);

        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);
              // Set use_strain_rates to false since we don't need viscosity
              in.reinit(fe_values, cell, this->introspection(), this->get_solution(), false);

              this->get_material_model().evaluate(in, out);

              MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
              AssertThrow(melt_out != nullptr,
                          ExcMessage("Need MeltOutputs from the material model for computing the melt properties."));

              const double p_c_scale = dynamic_cast<const MaterialModel::MeltInterface<dim>*>(&this->get_material_model())->p_c_scale(in, out, this->get_melt_handler(), true);

              cell->get_dof_indices (local_dof_indices);

              // for each compaction pressure dof, write into the output
              // vector the scaled compaction pressure. note that quadrature points and
              // dofs are enumerated in the same order
              for (unsigned int i=0; i<this->get_fe().base_element(pc_base_index).dofs_per_cell; ++i)
                {
                  const unsigned int system_local_dof
                    = this->get_fe().component_to_system_index(pc_component_index,
                                                               /*dof index within component=*/i);

                  vec_distributed(local_dof_indices[system_local_dof])
                    = p_c_scale * this->get_solution()[local_dof_indices[system_local_dof]];
                }
            }

        vec_distributed.compress(VectorOperation::insert);

        // now create a vector with the requisite ghost elements
        // and use it for estimating the gradients
        LinearAlgebra::BlockVector vec (this->introspection().index_sets.system_partitioning,
                                        this->introspection().index_sets.system_relevant_partitioning,
                                        this->get_mpi_communicator());
        vec = vec_distributed;

        DerivativeApproximation::approximate_gradient  (this->get_mapping(),
                                                        this->get_dof_handler(),
                                                        vec,
                                                        indicators,
                                                        pc_component_index);

        // Scale gradient in each cell with the correct power of h. Otherwise,
        // error indicators do not reduce when refined if there is a density
        // jump. We need at least order 1 for the error not to grow when
        // refining, so anything >1 should work. (note that the gradient
        // itself scales like 1/h, so multiplying it with any factor h^s, s>1
        // will yield convergence of the error indicators to zero as h->0)
        const double power = 1.0 + dim/2.0;
        {
          unsigned int i=0;
          for (cell = this->get_dof_handler().begin_active(); cell!=endc; ++cell, ++i)
            if (cell->is_locally_owned())
              indicators(i) *= std::pow(cell->diameter(), power);
        }

      }

  };





  template <int dim>
  class TestMeltMaterial:
    public MaterialModel::MeltInterface<dim>, public SimulatorAccess<dim>
  {
    public:
      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1.0;
      }

      virtual double reference_darcy_coefficient () const
      {
        const double porosity = 0.01;
        const double permeability = porosity;
        return permeability / 1.0;
      }

      virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                            typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const
      {
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const double x = in.position[i](0);
            const double z = in.position[i](1);
            out.viscosities[i] = 1.0;
            out.thermal_expansion_coefficients[i] = 0.0;
            out.specific_heat[i] = 1.0;
            out.thermal_conductivities[i] = 1.0;
            out.compressibilities[i] = 0.0;
            out.densities[i] = 1.0;
            // This is the RHS we need to use for the manufactured solution.
            // We calculated it by subtracting the first term of the RHS (\Nabla (K_D \rho_f g)) from the LHS
            // we computed using our analytical solution.
            double reactionterm;
            reactionterm = (0.8544665250e4 * exp(-0.60e2 * x * x - 0.160e3 * x * z - 0.180e3 * z * z + z + 0.1e1) - 0.4950000000e1 * exp(-0.100e3 * x * x - 0.320e3 * x * z - 0.340e3 * z * z + z + 0.1e1) + 0.3245696667e4 * exp(-0.60e2 * x * x - 0.240e3 * x * z - 0.240e3 * z * z + z) - 0.4687000001e3 * exp(-0.80e2 * x * x - 0.320e3 * x * z - 0.320e3 * z * z + z) + 0.160e3 * exp(-0.120e3 * x * x - 0.400e3 * x * z - 0.420e3 * z * z + z + 0.1e1) * x * x + 0.3992000000e3 * exp(-0.100e3 * x * x - 0.320e3 * x * z - 0.340e3 * z * z + z + 0.1e1) * z + 0.1600000000e1 * exp(-0.100e3 * x * x - 0.320e3 * x * z - 0.340e3 * z * z + z + 0.1e1) * x - 0.1399880000e4 * exp(-0.80e2 * x * x - 0.240e3 * x * z - 0.260e3 * z * z + z + 0.1e1) * x + 0.8826800400e5 * exp(-0.40e2 * x * x - 0.160e3 * x * z - 0.160e3 * z * z + z) * x + 0.1765360080e6 * exp(-0.40e2 * x * x - 0.160e3 * x * z - 0.160e3 * z * z + z) * z - 0.7840998000e5 * exp(-0.60e2 * x * x - 0.160e3 * x * z - 0.180e3 * z * z + z + 0.1e1) * x * x + 0.4943188800e6 * exp(-0.60e2 * x * x - 0.160e3 * x * z - 0.180e3 * z * z + z + 0.1e1) * z * z - 0.1925014410e6 * exp(-0.400e2 * pow(x + 0.2e1 * z, 0.2e1)) * z + 0.3698838000e5 * exp(-0.600e2 * pow(x + 0.2e1 * z, 0.2e1)) * x + 0.7397676000e5 * exp(-0.600e2 * pow(x + 0.2e1 * z, 0.2e1)) * z - 0.16e2 * exp(-0.1200e3 * pow(x + 0.2e1 * z, 0.2e1)) * x - 0.32e2 * exp(-0.1200e3 * pow(x + 0.2e1 * z, 0.2e1)) * z - 0.6688400000e4 * exp(-0.800e2 * pow(x + 0.2e1 * z, 0.2e1)) * x - 0.1337680000e5 * exp(-0.800e2 * pow(x + 0.2e1 * z, 0.2e1)) * z + 0.556e3 * exp(-0.1000e3 * pow(x + 0.2e1 * z, 0.2e1)) * x + 0.1112e4 * exp(-0.1000e3 * pow(x + 0.2e1 * z, 0.2e1)) * z - 0.9625072050e5 * exp(-0.400e2 * pow(x + 0.2e1 * z, 0.2e1)) * x + 0.8e1 * exp(-0.120e3 * x * x - 0.400e3 * x * z - 0.420e3 * z * z + z + 0.1e1) * x - 0.5763380040e5 * exp(-0.40e2 * x * x - 0.80e2 * x * z - 0.100e3 * z * z + z + 0.1e1) * z + 0.1971823986e7 * exp(-0.40e2 * x * x - 0.80e2 * x * z - 0.100e3 * z * z + z + 0.1e1) * z * z - 0.1950349995e5 * exp(-0.40e2 * x * x - 0.80e2 * x * z - 0.100e3 * z * z + z + 0.1e1) * x + 0.7723580040e6 * exp(-0.40e2 * x * x - 0.80e2 * x * z - 0.100e3 * z * z + z + 0.1e1) * x * x - 0.552024e6 * exp(-0.60e2 * x * x - 0.240e3 * x * z - 0.240e3 * z * z + z) * z * z - 0.138006e6 * exp(-0.60e2 * x * x - 0.240e3 * x * z - 0.240e3 * z * z + z) * x * x + 0.39792e5 * exp(-0.100e3 * x * x - 0.320e3 * x * z - 0.340e3 * z * z + z + 0.1e1) * z * z + 0.9860004000e4 * exp(-0.60e2 * x * x - 0.160e3 * x * z - 0.180e3 * z * z + z + 0.1e1) * x + 0.3318598800e5 * exp(-0.60e2 * x * x - 0.160e3 * x * z - 0.180e3 * z * z + z + 0.1e1) * z + 0.8e1 * exp(-0.100e3 * x * x - 0.400e3 * x * z - 0.400e3 * z * z + z) * x + 0.16e2 * exp(-0.100e3 * x * x - 0.400e3 * x * z - 0.400e3 * z * z + z) * z + 0.1599287976e7 * exp(-0.40e2 * x * x - 0.80e2 * x * z - 0.100e3 * z * z + z + 0.1e1) * x * z - 0.3447976000e5 * exp(-0.60e2 * x * x - 0.240e3 * x * z - 0.240e3 * z * z + z) * z + 0.1068266666e4 * exp(-0.80e2 * x * x - 0.320e3 * x * z - 0.320e3 * z * z + z) * x - 0.320e3 * exp(-0.120e3 * x * x - 0.400e3 * x * z - 0.420e3 * z * z + z + 0.1e1) * z * z + 0.4008e4 * exp(-0.100e3 * x * x - 0.320e3 * x * z - 0.340e3 * z * z + z + 0.1e1) * x * x + 0.1930895010e5 * z * z * exp(-0.20e2 * x * x - 0.20e2 * z * z + z + 0.1e1) + 0.1930895010e5 * x * x * exp(-0.20e2 * x * x - 0.20e2 * z * z + z + 0.1e1) - 0.9654475050e3 * z * exp(-0.20e2 * x * x - 0.20e2 * z * z + z + 0.1e1) + 0.9751995000e4 * exp(-0.20e2 * x * x - 0.80e2 * x * z - 0.80e2 * z * z + z) * x * x + 0.3900798000e5 * exp(-0.20e2 * x * x - 0.80e2 * x * z - 0.80e2 * z * z + z) * z * z + 0.3881097000e6 * exp(-0.40e2 * x * x - 0.160e3 * x * z - 0.160e3 * z * z + z) * x * x + 0.1552438800e7 * exp(-0.40e2 * x * x - 0.160e3 * x * z - 0.160e3 * z * z + z) * z * z + 0.1552438800e7 * exp(-0.40e2 * x * x - 0.160e3 * x * z - 0.160e3 * z * z + z) * x * z + 0.3900798000e5 * exp(-0.20e2 * x * x - 0.80e2 * x * z - 0.80e2 * z * z + z) * x * z - 0.9097753500e4 * exp(-0.40e2 * x * x - 0.160e3 * x * z - 0.160e3 * z * z + z) - 0.640e3 * exp(-0.120e3 * x * x - 0.400e3 * x * z - 0.420e3 * z * z + z + 0.1e1) * x * z + 0.4200000000e1 * exp(-0.120e3 * x * x - 0.400e3 * x * z - 0.420e3 * z * z + z + 0.1e1) - 0.9533794112e3 * exp(-0.20e2 * x * x - 0.20e2 * z * z + z + 0.1e1) - 0.1036095000e4 * exp(-0.80e2 * x * x - 0.240e3 * x * z - 0.260e3 * z * z + z + 0.1e1) + 0.2253333334e2 * exp(-0.100e3 * x * x - 0.400e3 * x * z - 0.400e3 * z * z + z) + 0.1123231354e5 * exp(-0.20e2 * x * x - 0.80e2 * x * z - 0.80e2 * z * z + z) - 0.1863776536e5 * exp(-0.40e2 * x * x - 0.80e2 * x * z - 0.100e3 * z * z + z + 0.1e1) - 0.3960040000e5 * exp(-0.80e2 * x * x - 0.240e3 * x * z - 0.260e3 * z * z + z + 0.1e1) * x * x - 0.3696256000e6 * exp(-0.80e2 * x * x - 0.240e3 * x * z - 0.260e3 * z * z + z + 0.1e1) * z * z + 0.7636384800e6 * exp(-0.60e2 * x * x - 0.160e3 * x * z - 0.180e3 * z * z + z + 0.1e1) * x * z - 0.4400336000e6 * x * exp(-0.80e2 * x * x - 0.240e3 * x * z - 0.260e3 * z * z + z + 0.1e1) * z + 0.47840e5 * exp(-0.80e2 * x * x - 0.320e3 * x * z - 0.320e3 * z * z + z) * x * z + 0.47840e5 * exp(-0.80e2 * x * x - 0.320e3 * x * z - 0.320e3 * z * z + z) * z * z + 0.11960e5 * exp(-0.80e2 * x * x - 0.320e3 * x * z - 0.320e3 * z * z + z) * x * x - 0.6320160000e4 * exp(-0.80e2 * x * x - 0.240e3 * x * z - 0.260e3 * z * z + z + 0.1e1) * z + 0.9509415350e5 * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)) * x + 0.1901883070e6 * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)) * z - 0.1191523129e5 * exp(z) - 0.552024e6 * exp(-0.60e2 * x * x - 0.240e3 * x * z - 0.240e3 * z * z + z) * x * z + 0.47712e5 * exp(-0.100e3 * x * x - 0.320e3 * x * z - 0.340e3 * z * z + z + 0.1e1) * x * z - 0.1723988000e5 * exp(-0.60e2 * x * x - 0.240e3 * x * z - 0.240e3 * z * z + z) * x - 0.1488766999e6 * exp(-0.20e2 * x * x - 0.80e2 * x * z - 0.80e2 * z * z + z) * x - 0.2977534000e6 * exp(-0.20e2 * x * x - 0.80e2 * x * z - 0.80e2 * z * z + z) * z + 0.2136533332e4 * exp(-0.80e2 * x * x - 0.320e3 * x * z - 0.320e3 * z * z + z) * z) * pow(-0.4950000000e1 + exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)), -0.3e1) * pow(-0.9950000000e1 + exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1)), -0.2e1);

            out.reaction_terms[i][porosity_idx] = reactionterm;
          }

        // fill melt outputs if they exist
        aspect::MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<aspect::MaterialModel::MeltOutputs<dim> >();

        if (melt_out != nullptr)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                double porosity = in.composition[i][porosity_idx];

                const double x = in.position[i](0);
                const double z = in.position[i](1);
                // porosity = 0.1000000000e-1 + 0.1000000000e0 * exp(-0.40e1 * pow(x + 0.20e1 * z, 0.2e1));
                // porosity = 0.1000000000e-1 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.2e1 * z, 0.2e1));
                melt_out->compaction_viscosities[i] = 0.1e0 + 0.1e0 * exp(-0.20e2 * x * x - 0.20e2 * z * z + 0.1e1);
                melt_out->fluid_viscosities[i] = 1.0;
                melt_out->permeabilities[i] = porosity;// K_D
                melt_out->fluid_density_gradients[i] = Tensor<1,dim>();
                melt_out->fluid_densities[i] = 0.5;
              }
          }
      }

  };


  template <int dim>
  class Gravity : public aspect::GravityModel::Interface<dim>
  {
    public:
      virtual Tensor<1,dim> gravity_vector (const Point<dim> &pos) const
      {
        const double x=pos[0];
        const double z=pos[1];
        Tensor<1,dim> gravity;
        gravity[0] = (0.4000000000e1 * x * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.400e2 * x - 0.8000e2 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)) + 0.4000000000e1 * x * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z)) / (0.9950000000e0 - 0.1000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)));
        gravity[1] = (-0.1333333333e1 * exp(z) + 0.4000000000e1 * z * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) - (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.8000e2 * x - 0.160000e3 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)) - 0.1e1 + 0.4000000000e1 * z * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) - (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z)) / (0.9950000000e0 - 0.1000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)));


        return gravity;
      }


  };


  template <int dim>
  class RefFunction : public Function<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      RefFunction () : Function<dim>(2*dim+3+2) {}
      virtual void vector_value (const Point< dim >   &p,
                                 Vector< double >   &values) const
      {
        double x = p(0);
        double z = p(1);

        values[0] = x;
        values[1] = -z + exp(z);
        values[2] = -(0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + 0.1e1 - z;
        values[3] = -(0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z);
        values[4] = x - 0.4000000000e1 * x * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) - 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.400e2 * x - 0.8000e2 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)) + 0.5e0 * (0.4000000000e1 * x * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.400e2 * x - 0.8000e2 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)) + 0.4000000000e1 * x * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z)) / (0.9950000000e0 - 0.1000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)));
        values[5] = -z + exp(z) - 0.4000000000e1 * z * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) - 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.8000e2 * x - 0.160000e3 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)) + 0.1e1 + 0.5e0 * (-0.1333333333e1 * exp(z) + 0.4000000000e1 * z * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) - (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.8000e2 * x - 0.160000e3 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)) - 0.1e1 + 0.4000000000e1 * z * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) - (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z)) / (0.9950000000e0 - 0.1000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)));
        values[6] = 0.1e1 - z;
        values[7] = 0;
        values[8] = 0.1000000000e-1 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1));

        // We have to scale the compaction pressure solution to p_c_bar using sqrt(K_D / ref_K_D).
        // K_D is equal to the porosity (as defined in the material model).
        const double K_D = values[8];
        const double ref_K_D = 0.01;
        const double p_c_scale = std::sqrt(K_D / ref_K_D);

        if (p_c_scale > 0)
          values[3] /= p_c_scale;
      }
  };



  /**
    * A postprocessor that evaluates the accuracy of the solution
    * by using the L2 norm.
    */
  template <int dim>
  class ConvergenceMeltPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
  ConvergenceMeltPostprocessor<dim>::execute (TableHandler &statistics)
  {
    RefFunction<dim> ref_func;
    const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +2);

    const unsigned int n_total_comp = this->introspection().n_components;

    Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p_c_bar (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_u_f (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_porosity (this->get_triangulation().n_active_cells());

    ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                        n_total_comp);
    ComponentSelectFunction<dim> comp_p_f(dim, n_total_comp);
    ComponentSelectFunction<dim> comp_p_c(dim+1, n_total_comp);
    ComponentSelectFunction<dim> comp_u_f(std::pair<unsigned int, unsigned int>(dim+2,dim+2+
                                                                                dim),
                                          n_total_comp);
    ComponentSelectFunction<dim> comp_p(dim+2+dim, n_total_comp);
    ComponentSelectFunction<dim> comp_porosity(dim+2+dim+2, n_total_comp);

    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_f);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_p_c_bar,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_p_c);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_porosity,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_porosity);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_u_f,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_u_f);


    // Loop over all cells to compute the error for p_c from p_c_bar
    const QGauss<dim> quadrature(this->get_parameters().stokes_velocity_degree+1);
    FEValues<dim> fe_values (this->get_mapping(),
                             this->get_fe(),
                             quadrature,
                             update_quadrature_points | update_values | update_gradients | update_JxW_values);

    MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
    MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());

    MeltHandler<dim>::create_material_model_outputs(out);

    typename DoFHandler<dim>::active_cell_iterator
    cell = this->get_dof_handler().begin_active(),
    endc = this->get_dof_handler().end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          in.reinit(fe_values, cell, this->introspection(), this->get_solution());

          this->get_material_model().evaluate(in, out);

          const double p_c_scale = dynamic_cast<const MaterialModel::MeltInterface<dim>*>(&this->get_material_model())->p_c_scale(in, out, this->get_melt_handler(), true);

          const unsigned int i = cell->active_cell_index();
          cellwise_errors_p_c[i] = cellwise_errors_p_c_bar[i] * p_c_scale;
        }

    const double u_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u, VectorTools::L2_norm);
    const double p_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p, VectorTools::L2_norm);
    const double p_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_f, VectorTools::L2_norm);
    const double p_c_bar_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_c_bar, VectorTools::L2_norm);
    const double p_c_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_p_c, VectorTools::L2_norm);
    const double phi_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_porosity, VectorTools::L2_norm);
    const double u_f_l2 = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_u_f, VectorTools::L2_norm);

    std::ostringstream os;
    os << std::scientific
       << "ndofs= " << this->get_solution().size()
       << " u_L2= " << u_l2
       << " p_L2= "  << p_l2
       << " p_f_L2= " << p_f_l2
       << " p_c_bar_L2= " << p_c_bar_l2
       << " p_c_L2= " << p_c_l2
       << " phi_L2= " << phi_l2
       << " u_f_L2= " << u_f_l2
       ;

    return std::make_pair("Errors", os.str());
  }


  template <int dim>
  class PressureBdry:

    public BoundaryFluidPressure::Interface<dim>
  {
    public:
      virtual
      void fluid_pressure_gradient (
        const types::boundary_id boundary_indicator,
        const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
        const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
        const std::vector<Tensor<1,dim> > &normal_vectors,
        std::vector<double> &output
      ) const
      {
        for (unsigned int q=0; q<output.size(); ++q)
          {
            const double x = material_model_inputs.position[q][0];
            const double z = material_model_inputs.position[q][1];
            Tensor<1,dim> gradient;
            gradient[0] = 0.4000000000e1 * x * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.400e2 * x - 0.8000e2 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1));
            gradient[1] = 0.4000000000e1 * z * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) - (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) / (-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1))) + 0.2000000000e0 * (0.1e0 + 0.1e0 * exp(0.1e1 - 0.2000000000e2 * x * x - 0.2000000000e2 * z * z)) * exp(z) * pow(-0.9900000000e0 + 0.2000000000e0 * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)), -0.2e1) * (-0.8000e2 * x - 0.160000e3 * z) * exp(-0.200e2 * pow(x + 0.20e1 * z, 0.2e1)) - 0.1e1;
            output[q] = gradient * normal_vectors[q];
          }
      }



  };

}

// explicit instantiations
namespace aspect
{
  ASPECT_REGISTER_GRAVITY_MODEL(Gravity,
                                "MyGravity",
                                "")

  ASPECT_REGISTER_MATERIAL_MODEL(TestMeltMaterial,
                                 "test melt material",
                                 "")
  ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(MeltVelocityRefinement,
                                            "melt velocity",
                                            "")
  ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(PCRefinement,
                                            "pc",
                                            "")

  ASPECT_REGISTER_POSTPROCESSOR(ConvergenceMeltPostprocessor,
                                "melt error calculation",
                                "A postprocessor that compares the numerical solution to the analytical "
                                "solution derived for incompressible melt transport in a 2D box as described "
                                "in the manuscript and reports the error.")

  ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(PressureBdry,
                                                "PressureBdry",
                                                "A fluid pressure boundary condition that prescribes the "
                                                "gradient of the fluid pressure at the boundaries as "
                                                "calculated in the analytical solution. ")

}
