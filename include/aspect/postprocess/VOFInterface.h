/*
 Copyright (C) 2016 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#ifndef __aspect__vofinterface_VOFInterface_h
#define __aspect__vofinterface_VOFInterface_h

// Aspect includes
#include <aspect/global.h>
#include <aspect/postprocess/interface.h>

// Deal II includes
#include <deal.II/base/parsed_function.h>

// Local includes

#include <aspect/vofinterface/VOFEngine.h>
#include <aspect/vofinterface/VOFOutput.h>

// Definition

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class VOFInterface : public Interface<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        //Constructor
        VOFInterface ();
        //Destructor
        virtual ~VOFInterface ();

        // Execution
        virtual std::pair<std::string, std::string> execute (TableHandler &statistics);

        // Utility
        double get_next_t (double time, double interval);

        // Parameter declarations.
        static void declare_parameters (ParameterHandler &prm);
        virtual void parse_parameters (ParameterHandler &prm);

      private:
        // Config vars
        double voleps;
        double out_interval;
        double next_out_t;
        unsigned int n_i_samp, n_e_samp;
        std::string err_out_fn;

        bool mms;
        double err_interval;
        double next_err_t;

        Functions::ParsedFunction<dim> initFunc;
        ::aspect::InterfaceTracker::VOFEngine<dim> engine;

        // Status vars
        bool initialized;
        ::aspect::InterfaceTracker::Output::VOFOutput<dim> *output;
    };
  }
}

#endif
