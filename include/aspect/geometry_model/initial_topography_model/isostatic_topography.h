#ifndef _aspect_geometry_model_initial_topography_model_isostatic_topography_h
#define _aspect_geometry_model_initial_topography_model_isostatic_topography_h

#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace InitialTopographyModel
  {
    /**
    * An initial topography model that computes isostatic topography
    * from the initial temperature and compositional fields by
    * balancing the mass of vertical columns.
    */
    template <int dim>
    class IsostaticTopography : public SimulatorAccess<dim>, public Interface<dim>
    {
      public:
        /**
        * Constructor.
        */
        IsostaticTopography();

        /**
        * Perform basic initialization and verify that the selected
        * geometry model is supported.
        */
        void initialize() override;

        /**
        * Perform an additional initialization step after the material model,
        * initial temperature model, initial composition model, gravity model,
        * and adiabatic conditions have been initialized.
        */
        void required_initialize() override;

        /**
         * Return the initial topography at the given surface point.
         */
        double
        value (const Point<dim-1> &surface_point) const override;


        /**
         * Return the maximum value of the initial topography.
         */
        double
        max_topography () const override;

        /**
         * Declare the parameters this class takes.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares.
         */
        void
        parse_parameters (ParameterHandler &prm) override;


      private:
        /**
         * Number of equally spaced lateral sampling points used to
         * construct the isostatic topography profile.
         */
        double n_lateral_points;

        /**
         * Number of equally spaced vertical sampling points used to
         * integrate the mass of each column.
         */
        unsigned n_vertical_points;

        /**
         * The depth above which the model is assumed to be in
         * isostatic equilibrium, and where the reference density is evaluated
         */
        double compensation_depth;

        /**
         * Maximum absolute value of the computed isostatic topography.
         */
        double max_isostatic_topography;

        /**
         * The characteristic horizontal length scale (in m) over which
         * topographic loads are assumed to attain isostatic compensation.
         */
        double isostatic_length_scale;

        /**
         * Maximum computed initial topography.
         */
        double maximal_topography;

        /**
         * Isostatic topography evaluated at the lateral sampling points.
         */
        std::vector<double> topography;
    };
  }
}

#endif
