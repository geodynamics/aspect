#ifndef _aspect_boundary_traction_infill_ascii_data_h
#define _aspect_boundary_traction_infill_ascii_data_h

#include <aspect/boundary_traction/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace BoundaryTraction
  {
    using namespace dealii;

    /**
     * A class that implements prescribed traction boundary conditions determined
     * from pressures given in an AsciiData input file.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class InfillAsciiData : public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        InfillAsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataBoundary<dim>::initialize;

        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        Tensor<1,dim>
        boundary_traction (const types::boundary_id surface_boundary_id,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        void
        update () override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        types::boundary_id surface_boundary_id;
        double rock_density;
        double sediment_density;
        double crustal_density;
        double rock_infill_height;
    };
  }
}
#endif
