#include <aspect/material_model/interface.h>
#include <aspect/simulator.h>
#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

#include <iostream>

using namespace aspect;


template <int dim>
void post_constraints_creation (const SimulatorAccess<dim> &sim,
                                AffineConstraints<double> &)
{
  static bool first = true;
  if (first)
    {
      std::stringstream ss;
      sim.get_simulator().write_plugin_graph(ss);
      std::string content = ss.str();

      const char *checks[] =
      {
        "digraph Plugins",
        "label=\"TestMaterial\"",
        "N6aspect12TestMaterialILi2EEE -> N6aspect13MaterialModel9InterfaceILi2EEE",
        "label=\"simpler\"",
        "Simulator -> SimulatorAccess",
        "N6aspect19BoundaryComposition9AsciiDataILi2EEE [label=\"ascii data\""
      };

      for (unsigned int i = 0; i < sizeof(checks)/sizeof(checks[0]); ++i)
        {
          bool found = content.find(checks[i])!=std::string::npos;

          if (found)
            std::cout << "found " << checks[i] << std::endl;
          Assert(found, ExcInternalError());
        }

    }
  first = false;
}


template <int dim>
void signal_connector (SimulatorSignals<dim> &signals)
{
  signals.post_constraints_creation.connect (&post_constraints_creation<dim>);
}


ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)


// make a dummy material model and instantiate and register it
namespace aspect
{
  template <int dim>
  class TestMaterial : public MaterialModel::Interface<dim>
  {
    public:
      virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                            MaterialModel::MaterialModelOutputs<dim> &out) const
      {}

      virtual bool is_compressible () const
      {
        return false;
      }

      virtual double reference_viscosity () const
      {
        return 1;
      }
  };


  ASPECT_REGISTER_MATERIAL_MODEL(TestMaterial,
                                 "TestMaterial",
                                 "A material model.")
}
