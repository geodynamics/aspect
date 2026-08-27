/*
  A geometry model plugin that extends ASPECT's built-in
  "box with lithosphere boundary indicators" (GeometryModel::TwoMergedBoxes)
  with two additional boundary indicators, and renames the resulting three
  vertically stacked zones on the left/right boundaries (from the surface
  down) to:
    - "left/right lithosphere": a new, thin strip at the very top, whose
      thickness is set by the "Oceanic plate thickness" parameter,
    - "left/right upper mantle": the remainder of "Lithospheric thickness"
      below that strip (this is what the base class itself calls
      "left/right lithosphere"), and
    - "left"/"right": the mantle below "Lithospheric thickness", unchanged
      from the base class.

  This makes it possible to prescribe, on the left (or right) side of the
  box, three independent boundary conditions stacked vertically:
    - a traction (stress) condition on "left/right lithosphere" (the thin
      plate strip at the top), e.g. to drive inflow of an oceanic plate,
    - a different condition (e.g. tangential/free-slip or zero/no-slip
      velocity) on "left/right upper mantle" below it, and
    - a traction condition on the underlying "left"/"right" mantle
      boundary below that (unchanged from the base class), e.g. to allow
      outflow.

  This class derives from TwoMergedBoxes (instead of reimplementing
  everything from scratch) so that it is still recognized wherever ASPECT
  specifically checks for a GeometryModel::TwoMergedBoxes geometry via
  Plugins::plugin_type_matches (e.g. the free surface diffusion mesh
  deformation plugin, the adiabatic/harmonic/random initial temperature
  perturbation plugins, the initial lithostatic pressure boundary
  traction plugin, lateral averaging, ...). All base class functionality
  (mesh generation for the mantle/lithosphere split, depth(),
  get_extents(), maximal_depth(), natural coordinates, etc.) is left
  completely unchanged; this class only adds the two extra boundary
  indicators and renames the affected symbolic boundary names.

  Notes / limitations:
   - The split between "upper mantle" and "lithosphere" is NOT guaranteed
     to fall exactly on a coarse mesh cell face (unlike the
     mantle/upper-mantle split, which by construction of the two merged
     triangulations always does). The actual boundary between the two
     indicators will therefore follow the nearest existing cell edge to
     the requested "Oceanic plate thickness".
   - In 3D, the additional split is only applied to the "left"/"right"
     boundaries, not to "front"/"back": those keep exactly the behavior
     (and naming) of the base class.
*/
#ifndef _aspect_geometry_model_box_with_oceanic_plate_h
#define _aspect_geometry_model_box_with_oceanic_plate_h

#include <aspect/geometry_model/two_merged_boxes.h>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    class BoxWithOceanicPlate : public TwoMergedBoxes<dim>
    {
      public:
        /**
         * Generate a coarse mesh for the geometry described by this class,
         * and assign the additional "left/right lithosphere" (thin plate
         * strip) boundary indicators on top of what the base class sets up.
         */
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        /**
         * Return the set of boundary indicators used by this model: the
         * ones from TwoMergedBoxes, plus two more for the thin plate strip.
         */
        std::set<types::boundary_id>
        get_used_boundary_indicators () const override;

        /**
         * Return the mapping from symbolic boundary names to boundary
         * indicators. The base class' own "left/right lithosphere" names
         * are renamed to "left/right upper mantle", and "left/right
         * lithosphere" is used instead for the new, thinner strip at the
         * very top.
         */
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const override;

        /**
         * Declare the parameters this class takes through input files
         * (those of TwoMergedBoxes, plus one more).
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Set the "left/right lithosphere" boundary indicators on the top
         * part of the base class' own "left/right lithosphere" boundaries
         * (renamed "left/right upper mantle" by get_symbolic_boundary_names_map()),
         * based on the vertical position of each boundary face relative to
         * height_oceanic. Called once after the coarse mesh is built, and
         * again after every mesh refinement (child faces would otherwise
         * simply inherit their parent's indicator, which is fine, but we
         * reapply this in the same way the base class reapplies its own
         * split for robustness).
         */
        void
        set_oceanic_plate_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const;

        /**
         * The height (measured in the same vertical coordinate and from
         * the same origin as the rest of the box geometry) above which
         * the left/right boundary gets the "lithosphere" (thin plate
         * strip) indicator instead of the "upper mantle" indicator. Set in
         * parse_parameters() from the "Oceanic plate thickness" input
         * parameter.
         */
        double height_oceanic;
    };
  }
}

#endif
