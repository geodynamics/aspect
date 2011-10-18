//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__post_processor_h
#define __aspect__post_processor_h


namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    class Base
    {
      public:
        double get_time () const;

        double get_timestep_number () const;

        const TrilinosWrappers::Vector &
        get_solution () const;

        const TrilinosWrappers::Vector &
        get_old_solution () const;


        const TrilinosWrappers::Vector &
        get_old_old_solution () const;

        const parallel::distributed::DoFhandler<dim> &
        get_dof_handler () const;
    };
  }
}


#endif
