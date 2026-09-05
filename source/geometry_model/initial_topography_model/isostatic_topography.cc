#include <aspect/geometry_model/initial_topography_model/isostatic_topography.h>
#include <aspect/geometry_model/box.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace InitialTopographyModel
  {

    template <int dim>
    IsostaticTopography<dim>::
    IsostaticTopography():
      maximal_topography(0.0)
    {
      this->set_required_initialized(false);
    }



    template <int dim>
    void
    IsostaticTopography<dim>::
    initialize()
    {
      AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("The isostatic topography model currently only supports the box geometry model."));
    }

    template <int dim>
    void
    IsostaticTopography<dim>::
    required_initialize()
    {
      AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("The isostatic topography model currently only supports the box geometry model."));

      const auto &geometry_model = dynamic_cast<const GeometryModel::Box<dim> &>(this->get_geometry_model());

      AssertThrow(compensation_depth >= 0.0 &&
                  compensation_depth <= geometry_model.maximal_depth(),
                  ExcMessage("The compensation depth must be between zero and the maximal model depth."));

      const double x_min = 0.0;
      const double x_max = geometry_model.get_extents()[0];

      const double y_top = geometry_model.get_extents()[1];

      // Parameters for the sampling mesh
      const double grid_spacing = (x_max - x_min) / (n_lateral_points - 1.0);
      const double delta_z = this->get_geometry_model().maximal_depth() / (n_vertical_points-1);

      // Parameters for the gaussian filter
      const bool use_gaussian_filter = (isostatic_length_scale > 0.0);
      const double sigma = use_gaussian_filter?
                           isostatic_length_scale * std::sqrt(2.0 * std::log(2.0)) / (2.0 * numbers::PI):
                           0.0;
      AssertThrow(!use_gaussian_filter || sigma > grid_spacing,
                  ExcMessage("The isostatic length scale is too small compared "
                             "to the lateral sampling spacing. Increase the "
                             "isostatic length scale or disable the Gaussian "
                             "filter by setting it to zero."));

      const unsigned int kernel_radius = use_gaussian_filter?
                                         std::ceil(3.0 * sigma / grid_spacing):
                                         0.0;
      std::vector<double> kernel = use_gaussian_filter ?
                                   std::vector<double>(2 * kernel_radius + 1) :
                                   std::vector<double>();
      if (use_gaussian_filter)
        {
          for (unsigned int k = 0; k < kernel.size(); ++k)
            {
              const double distance = std::abs(static_cast<int>(k) - static_cast<int>(kernel_radius))
                                      * grid_spacing;

              kernel[k] = std::exp(-0.5 * Utilities::fixed_power<2>(distance / sigma));
            }
        }

      // Determine the reference density at the compensation depth beneath
      // the center of the model domain.
      double ref_density = std::numeric_limits<double>::lowest();
      const unsigned int ref_column_index = static_cast<unsigned int>(std::round((x_max - x_min) / 2.0 / grid_spacing));
      const unsigned int ref_depth_index = static_cast<unsigned int>(std::round(compensation_depth/delta_z));

      // Compute the mass of each vertical column and accumulate the
      // total column mass.
      std::vector<double> column_masses(n_lateral_points);
      double total_column_mass = 0.0;

      for (unsigned int i = 0; i < n_lateral_points; ++i)
        {
          const double x = x_min + i * grid_spacing;
          double column_mass = 0.0;

          // Compute pressure by vertically integrating the density,
          // following the same approach used by the adiabatic conditions.
          double pressure;
          double density;

          MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
          MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());

          in.requested_properties = MaterialModel::MaterialProperties::equation_of_state_properties;
          in.velocity[0] = Tensor <1,dim> ();

          const Tensor <1,dim> g = this->get_gravity_model().gravity_vector(this->get_geometry_model().representative_point(0));
          const Point<dim> point_surf = this->get_geometry_model().representative_point(0);
          const Point<dim> point_bot = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
          const int gravity_direction =  (g * (point_bot - point_surf) >= 0) ?
                                         1 :
                                         -1;

          for (unsigned int j=0; j<n_vertical_points; ++j)
            {
              Point<dim> position;
              position[0] = x;
              position[1] = y_top - j * delta_z;

              if (j==0)
                {
                  // Use the surface pressure of the adiabatic condition
                  pressure = this->get_adiabatic_conditions().pressure(position);
                }
              else
                {
                  // Update the pressure using the density evaluated at the
                  // previous vertical sampling point.
                  const double gravity = gravity_direction * this->get_gravity_model().gravity_vector(position).norm();

                  pressure = pressure + density * gravity * delta_z;
                }

              // Use the initial temperature and composition at the local position
              std::vector<double> composition(this->n_compositional_fields());

              for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
                composition[c] =
                  this->get_initial_composition_manager().initial_composition(position, c);

              in.position[0] = position;
              in.temperature[0] =  this->get_initial_temperature_manager().initial_temperature(position);
              in.pressure[0] = pressure;
              in.composition[0] = composition;

              this->get_material_model().evaluate(in, out);
              density = out.densities[0];

              // Store the reference density at the compensation depth beneath
              // the center of the model domain.
              if (i==ref_column_index and j==ref_depth_index)
                ref_density = density;

              const double integration_weight = (j == 0 || j == n_vertical_points-1) ? 0.5 : 1.0;
              column_mass += integration_weight * density * delta_z;
            }

          column_masses[i] = column_mass;
          total_column_mass += column_mass;
        }

      AssertThrow(ref_density > 0.0,
                  ExcMessage("The reference density is not set properly by the isostatic topography model."));

      const double average_column_mass = total_column_mass / n_lateral_points;

      std::vector<double> raw_topography(n_lateral_points, std::numeric_limits<double>::max());
      for (unsigned int i = 0; i < n_lateral_points; ++i)
        {
          raw_topography[i] = (average_column_mass - column_masses[i]) / ref_density;
          raw_topography[i] = std::clamp(raw_topography[i], -max_isostatic_topography, max_isostatic_topography);
        }

      topography.resize(n_lateral_points);
      if (use_gaussian_filter)
        {
          for (unsigned int i=0; i<n_lateral_points; ++i)
            {
              double sum = 0;
              const int begin = std::max<int>(0, i-kernel_radius);
              const int end = std::min<int>(n_lateral_points-1, i+kernel_radius);

              double weight_sum = 0.0;

              for (int j=begin; j<=end; ++j)
                {
                  const double w = kernel[j-i+kernel_radius];
                  weight_sum += w;
                  sum += w * raw_topography[j];
                }
              topography[i] = sum/weight_sum;
            }
        }
      else
        {
          for (unsigned int i = 0; i < n_lateral_points; ++i)
            topography[i] = raw_topography[i];
        }

      maximal_topography = *std::max_element(topography.begin(), topography.end());

      this->set_required_initialized(true);
    }




    template <int dim>
    double
    IsostaticTopography<dim>::
    value(const Point<dim-1> &surface_point) const
    {
      if (!this->is_required_initialized())
        return 0.0;

      AssertThrow(!topography.empty(),
                  ExcMessage("The isostatic topography profile has not been properly initialized."));

      const auto &geometry_model = dynamic_cast<const GeometryModel::Box<dim> &>(this->get_geometry_model());

      const double x_min = 0.0;
      const double x_max = geometry_model.get_extents()[0];
      const double x = surface_point[0];

      const double grid_spacing = (x_max - x_min) / (n_lateral_points - 1.0);
      const double index = (x - x_min) / grid_spacing;
      const unsigned int i = static_cast<unsigned int>(std::floor(index));
      const double alpha = index - i;

      return (1.0 - alpha) * topography[i] + alpha * topography[i+1];
    }



    template <int dim>
    double
    IsostaticTopography<dim>::
    max_topography() const
    {
      return maximal_topography;
    }



    template <int dim>
    void
    IsostaticTopography<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Initial topography model");
        {
          prm.enter_subsection("Isostatic topography");
          {
            prm.declare_entry("Number of lateral points",
                              "2000",
                              Patterns::Integer(2),
                              "Number of equally spaced lateral points used to "
                              "calculate and store the initial isostatic "
                              "topography profile.");

            prm.declare_entry("Number of vertical points",
                              "2000",
                              Patterns::Integer(2),
                              "Number of equally spaced vertical points used to "
                              "integrate the mass of each model column.");


            prm.declare_entry("Compensation depth",
                              "500000",
                              Patterns::Double(0),
                              "The depth above which the model is assumed to be in "
                              "isostatic equilibrium. The implementation assumes that "
                              "material flows horizontally at this depth to balance "
                              "column masses. Consequently, the reference density used "
                              "for the isostatic compensation is taken as the density of "
                              "the material at the compensation depth. Units: m.");

            prm.declare_entry("Maximum isostatic topography",
                              "10e3",
                              Patterns::Double(0.0),
                              "Maximum absolute value of the calculated isostatic "
                              "topography. The calculated topography is limited "
                              "to the interval between the negative and positive "
                              "values of this parameter. Units: m.");

            prm.declare_entry("Isostatic length scale",
                              "500e3",
                              Patterns::Double(0.0),
                              "The characteristic horizontal length scale (in m) over which "
                              "topographic loads are assumed to attain isostatic compensation. "
                              "The computed isostatic topography is smoothed with a Gaussian filter "
                              "whose half-power wavelength is given by this value. "
                              "Set to 0 to disable Gaussian filtering. Units: m.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    IsostaticTopography<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Initial topography model");
        {
          prm.enter_subsection("Isostatic topography");
          {
            // todo_df
            n_lateral_points = prm.get_integer("Number of lateral points");
            n_vertical_points = prm.get_integer("Number of vertical points");
            compensation_depth = prm.get_double("Compensation depth");
            max_isostatic_topography = prm.get_double("Maximum isostatic topography");
            isostatic_length_scale = prm.get_double("Isostatic length scale");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      AssertThrow(n_lateral_points >= 2 && n_vertical_points >= 2,
                  ExcMessage("The number of lateral and vertical sampling points must be at least two."));

    }
  }
}



namespace aspect
{
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(
      IsostaticTopography,
      "isostatic topography",
      "An initial topography model that computes an isostatic "
      "topography profile by balancing the mass of vertical columns "
      "using the initial temperature and compositional fields.")
  }
}
