#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/particle/property/interface.h>

#include <iostream>

unsigned int counter_without = 0, counter_with = 0;
bool quiet = true;


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      class PostInitializeParticleProperty : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:

          virtual
          void
          initialize ()
          {
            std::cout << "initialize" << std::endl;
          }

          virtual
          void
          post_initialize ()
          {
            const Particle::Property::Manager<dim> &manager = this->get_particle_world().get_property_manager();
            post_initialized_info = manager.get_data_info().get_field_index_by_name("initial position");
            std::cout << "initial position: post_initialized_info = " << post_initialized_info << std::endl;

            post_initialized_info = manager.get_data_info().get_field_index_by_name("initial C_1");
            std::cout << "initial composition C_1: post_initialized_info = " << post_initialized_info << std::endl;

            post_initialized_info = manager.get_data_info().get_field_index_by_name("velocity");
            std::cout << "velocity: post_initialized_info = " << post_initialized_info << std::endl;
            exit(0);
          }

          virtual
          std::vector<std::pair<std::string, unsigned int> >
          get_property_information() const
          {
            return std::vector<std::pair<std::string, unsigned int> >();
          }

        private:
          unsigned int post_initialized_info;

      };
    }
  }
}


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(PostInitializeParticleProperty,
                                        "PostInitializeParticleProperty",
                                        "")
    }
  }
}
