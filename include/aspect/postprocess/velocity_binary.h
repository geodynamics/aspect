#ifndef __aspect__postprocess_velocity_binary_h 
#define __aspect__postprocess_velocity_binary_h 


#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect {
    namespace Postprocess {

        template<int dim>
        class VelocityBinary : public Interface<dim>, public ::aspect::SimulatorAccess<dim> {
        public:
            VelocityBinary();

            ~VelocityBinary();

            void initialize();

            std::pair<std::string, std::string> execute(TableHandler &statistics);

           /* static void declare_parameters(ParameterHandler &prm);

            virtual void parse_parameters(ParameterHandler &prm);

            void save();

            void load();
*/
/*            template<class Archive>
            void serialize(Archive &ar, const unsigned int version);
*/
            // private:
        };
    }
}

#endif
