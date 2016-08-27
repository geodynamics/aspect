#ifndef __aspect__postprocess_binary_data_h
#define __aspect__postprocess_binary_data_h


#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <boost/serialization/serialization.hpp>

namespace aspect
{
  namespace Postprocess
  {

//    struct Field
//    {
//        double time;
//        double time_step;
//        double old_time_step;
//        unsigned int timestep_number;
//
//      public:
//        template <class Archive>
//        void serialize(Archive &ar, const unsigned int version)
//        {
//          ar &time;
//          ar &time_step;
//          ar &old_time_step;
//          ar &timestep_number;
//        }
//    };

    template<int dim>
    class BinaryData : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        BinaryData();

        ~BinaryData();

        void initialize();

        void update_time();

        std::pair<std::string, std::string> execute(TableHandler &statistics);

        static void declare_parameters(ParameterHandler &prm);

        virtual void parse_parameters(ParameterHandler &prm);
        /*
                    void save();

                    void load();
        */
      private:
        unsigned int my_id;
        std::string filename_prefix;

        aspect::Utilities::BinaryInputFields attributes;
    };
  }
}

#endif
