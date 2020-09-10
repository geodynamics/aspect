#include <aspect/material_model/simple.h>
#include <aspect/global.h>

#include <deal.II/base/point.h>
#include <deal.II/base/exceptions.h>


namespace aspect
{
  /**
   * This is the "NSinker" benchmark defined in \cite May2015496 in the implementation
   * of \cite rudi2017weighted. It creates a number of spherical high-viscosity, high-density
   * sinking spheres in a box geometry that provide a challenge for the Stokes preconditioner.
   * The difficulty of the problem is determined by the number of sinkers and the viscosity
   * contrast between sinkers and background.
   */
  namespace NSinkerBenchmark
  {
    using namespace dealii;

    /**
     * @note This benchmark only talks about the flow field, not about a
     * temperature field. All quantities related to the temperature are
     * therefore set to zero in the implementation of this class.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class NSinkerMaterial : public MaterialModel::Interface<dim>
    {
      public:
        /**
         * Constructor
         */
        NSinkerMaterial ();

        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          const double sqrt_dynamic_viscosity_ratio = std::sqrt(dynamic_viscosity_ratio);

          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {
              const double chi = inside_outside_factor(in.position[i]);
              out.viscosities[i] = (sqrt_dynamic_viscosity_ratio - 1./sqrt_dynamic_viscosity_ratio)*(1-chi) + 1./sqrt_dynamic_viscosity_ratio;
              out.densities[i] = sinker_density * (1.0 - chi);
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;
            }
        }

        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the continuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        virtual bool is_compressible () const
        {
          return false;
        }
        /**
         * @}
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("NSinker");
            {
              prm.declare_entry ("Dynamic viscosity ratio", "1e3",
                                 Patterns::Double (0),
                                 "Viscosity in the sinkers.");
              prm.declare_entry ("Sinker density", "10",
                                 Patterns::Double (0),
                                 "Density in the sinkers.");
              prm.declare_entry ("Number of sinkers", "4",
                                 Patterns::Integer (0,75),
                                 "Number of sinking spheres.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("NSinker");
            {
              dynamic_viscosity_ratio = prm.get_double ("Dynamic viscosity ratio");
              sinker_density = prm.get_double ("Sinker density");
              n_sinkers = prm.get_integer ("Number of sinkers");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();

          // Declare dependencies on solution variables
          this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
          this->model_dependence.density = MaterialModel::NonlinearDependence::none;
          this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
          this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
          this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
        }



        /**
         * The reference viscosity was chosen to coincide with the reference length scale of
         * a box of size 1 (0.01) so that the resulting pressure scaling is equal to one, and
         * therefore does not influence the scaling of the equations.
         */
        virtual double reference_viscosity () const
        {
          return 0.01;
        }

      private:
        /**
         * Ratio of viscosities between sinkers and background material.
         */
        double dynamic_viscosity_ratio;

        /**
         * Sinker density. For a gravity of 1 this corresponds to
         * the size of the right-hand-side forcing term.
         */
        double sinker_density;

        /**
         * Number of sinking spheres.
         */
        unsigned int n_sinkers;

        /**
         * Centers for the sinkers provided by Cedric Thielot (pers. comm. from Dave May)
         */
        std::vector<Point<3> > centers;

        /**
         * Parameters for evaluating viscosity
         */
        double delta;
        double omega;

        /**
         * Return a factor that indicates whether a point is inside or
         * outside the sinker spheres. The function can be thought of
         * as a 1 (outside) vs 0 (inside) factor, but is in reality a
         * smoothed version of this that takes into account the
         * distance from the centers of the sinkers.
         */
        double inside_outside_factor (const Point<dim> &p) const;
    };


    template<int dim>
    NSinkerMaterial<dim>::NSinkerMaterial ()
    {
      delta = 200.0;
      omega = 0.1;

      centers.resize(75);
      centers[0] = Point<3>(2.4257829890e-01, 1.3469574514e-02, 3.8313885004e-01);
      centers[1] = Point<3>(4.1465269048e-01, 6.7768972864e-02, 9.9312692973e-01);
      centers[2] = Point<3>(4.8430804651e-01, 7.6533776604e-01, 3.1833815403e-02);
      centers[3] = Point<3>(3.0935481671e-02, 9.3264044027e-01, 8.8787953411e-01);
      centers[4] = Point<3>(5.9132973039e-01, 4.7877868473e-01, 8.3335433660e-01);
      centers[5] = Point<3>(1.8633519681e-01, 7.3565270739e-01, 1.1505317181e-01);
      centers[6] = Point<3>(6.9865863058e-01, 3.5560411138e-01, 6.3830000658e-01);
      centers[7] = Point<3>(9.0821050755e-01, 2.9400041480e-01, 2.6497158886e-01);
      centers[8] = Point<3>(3.7749399775e-01, 5.4162011554e-01, 9.2818150340e-03);
      centers[9] = Point<3>(8.5247022139e-01, 4.6701098395e-01, 5.3607231962e-02);
      centers[10] = Point<3>(9.7674759057e-01, 1.9675474344e-01, 8.5697294067e-01);
      centers[11] = Point<3>(1.4421375987e-01, 8.0066218823e-01, 7.2939761948e-01);
      centers[12] = Point<3>(9.8579064709e-01, 1.8340570954e-01, 4.9976021075e-01);
      centers[13] = Point<3>(4.6986202126e-01, 9.7099129947e-01, 4.5077026191e-01);
      centers[14] = Point<3>(9.5791877292e-02, 9.7408164664e-01, 3.9023506101e-01);
      centers[15] = Point<3>(6.8067035576e-01, 2.6669318800e-02, 2.3124107450e-01);
      centers[16] = Point<3>(4.6873909443e-01, 9.7960100555e-02, 4.1541002524e-01);
      centers[17] = Point<3>(7.9629418710e-01, 3.1640260216e-01, 7.7853444953e-01);
      centers[18] = Point<3>(8.2849331472e-01, 4.8714042059e-01, 3.6904878000e-01);
      centers[19] = Point<3>(6.0284549678e-01, 2.4264360789e-02, 8.1111178631e-01);
      centers[20] = Point<3>(3.5579259291e-01, 8.0610905439e-01, 2.7487712366e-01);
      centers[21] = Point<3>(8.5981739865e-01, 9.5101905612e-01, 7.7727618477e-01);
      centers[22] = Point<3>(6.8083745971e-01, 8.3518540665e-01, 9.6112961413e-01);
      centers[23] = Point<3>(7.0542474869e-01, 7.3751226102e-02, 5.3685709440e-01);
      centers[24] = Point<3>(9.5718558131e-01, 4.1806501915e-01, 4.1877679639e-01);
      centers[25] = Point<3>(3.8161700050e-01, 8.3692747440e-01, 2.4006224854e-01);
      centers[26] = Point<3>(7.2621119848e-01, 4.3161282150e-01, 1.1669089744e-01);
      centers[27] = Point<3>(2.2391322592e-01, 3.0958795748e-01, 2.4480139429e-01);
      centers[28] = Point<3>(3.7703382754e-01, 8.0753940242e-01, 3.1473643301e-01);
      centers[29] = Point<3>(7.7522956709e-01, 2.8333410774e-01, 9.9634871585e-01);
      centers[30] = Point<3>(6.3286731189e-01, 6.0091089904e-01, 5.0948022423e-01);
      centers[31] = Point<3>(8.3412860373e-01, 1.9944285005e-01, 3.5980841627e-02);
      centers[32] = Point<3>(7.3000523063e-01, 1.9791117972e-01, 2.9319749786e-01);
      centers[33] = Point<3>(7.7034656693e-01, 2.1475035521e-01, 3.0922000730e-01);
      centers[34] = Point<3>(6.0662675677e-02, 5.5759010630e-01, 4.1691651960e-01);
      centers[35] = Point<3>(1.1594487686e-01, 6.8554530558e-01, 9.5995079957e-01);
      centers[36] = Point<3>(2.7973348288e-02, 1.4806467395e-01, 5.2297503060e-01);
      centers[37] = Point<3>(6.4133927209e-01, 9.8914607800e-01, 5.7813295237e-01);
      centers[38] = Point<3>(6.8053043246e-01, 6.7497840462e-01, 3.6204645148e-01);
      centers[39] = Point<3>(9.1470996426e-01, 5.3036934674e-01, 9.1761070439e-01);
      centers[40] = Point<3>(2.8310876353e-01, 2.0898862472e-01, 4.7181570645e-01);
      centers[41] = Point<3>(8.0657831198e-01, 1.6168943288e-01, 5.1429839456e-01);
      centers[42] = Point<3>(8.1311740159e-01, 6.4168478858e-02, 4.7962416312e-01);
      centers[43] = Point<3>(4.3309508843e-02, 9.0291512474e-01, 2.9450144167e-01);
      centers[44] = Point<3>(6.8573011443e-01, 6.6033273035e-02, 8.2121989495e-01);
      centers[45] = Point<3>(2.4277445452e-01, 3.1025718772e-01, 4.9255406554e-01);
      centers[46] = Point<3>(3.5617944848e-01, 3.0799053857e-01, 3.9698166931e-01);
      centers[47] = Point<3>(7.0916077621e-02, 8.8651657239e-01, 6.8403214295e-01);
      centers[48] = Point<3>(5.2822650202e-01, 9.0281945043e-01, 6.8650344000e-01);
      centers[49] = Point<3>(6.3316007640e-02, 1.5214040370e-01, 2.3765034985e-02);
      centers[50] = Point<3>(4.1894298765e-01, 1.7479340461e-01, 7.5275125343e-01);
      centers[51] = Point<3>(4.9031640053e-01, 7.4774375406e-01, 3.2927456281e-01);
      centers[52] = Point<3>(1.1757708859e-01, 1.1812786251e-01, 3.7498524244e-01);
      centers[53] = Point<3>(3.7696964032e-01, 7.2874483733e-01, 1.4480990830e-02);
      centers[54] = Point<3>(3.8201288152e-01, 4.9049964756e-01, 8.2757658503e-01);
      centers[55] = Point<3>(7.9664661586e-02, 9.2396727806e-01, 1.1804237828e-01);
      centers[56] = Point<3>(9.3825167927e-01, 1.9597347043e-01, 7.2611756191e-01);
      centers[57] = Point<3>(8.5786301170e-01, 1.0363770514e-01, 8.3891028205e-01);
      centers[58] = Point<3>(5.6511039453e-01, 8.1040084307e-01, 4.0696941614e-01);
      centers[59] = Point<3>(9.3497714490e-01, 1.6087440083e-01, 8.1605472361e-01);
      centers[60] = Point<3>(4.3173963829e-01, 2.4810082244e-01, 8.3052277138e-01);
      centers[61] = Point<3>(5.9621858625e-01, 6.4577903070e-01, 6.0816894547e-01);
      centers[62] = Point<3>(4.9546643556e-01, 3.0438243752e-01, 7.5562733447e-01);
      centers[63] = Point<3>(8.2861043319e-01, 4.5555055302e-01, 4.3814466774e-01);
      centers[64] = Point<3>(8.9743076959e-01, 1.1894442752e-01, 9.8993320995e-02);
      centers[65] = Point<3>(6.9884936497e-01, 5.6127713367e-01, 3.8478565932e-01);
      centers[66] = Point<3>(9.2576270966e-02, 9.2938612771e-01, 1.9264837596e-01);
      centers[67] = Point<3>(8.4125479722e-01, 9.6937695284e-01, 3.1844636161e-01);
      centers[68] = Point<3>(1.2799954700e-01, 2.8838638276e-01, 9.0993508972e-01);
      centers[69] = Point<3>(2.7905288352e-01, 4.1813262758e-02, 7.5550716964e-01);
      centers[70] = Point<3>(8.0900019305e-01, 8.6624463269e-01, 9.7354159503e-01);
      centers[71] = Point<3>(3.1358765965e-01, 4.6779574243e-01, 2.4304298462e-01);
      centers[72] = Point<3>(8.2344259034e-01, 5.9961585635e-01, 7.4369772512e-01);
      centers[73] = Point<3>(3.2766604253e-01, 8.3176720460e-02, 9.5114077951e-01);
      centers[74] = Point<3>(8.2308128282e-01, 5.2712029523e-01, 3.1080186614e-01);
    }



    template <>
    double
    NSinkerMaterial<2>::inside_outside_factor (const Point<2> &p) const
    {
      double chi = 1.0;

      for (unsigned int s=0; s<n_sinkers; ++s)
        {
          double dist = p.distance(Point<2>(centers[s](0), centers[s](1)));
          double temp = 1-std::exp(-delta*
                                   std::pow(std::max(0.0,dist-omega/2.0),2));
          chi *= temp;
        }
      return chi;
    }



    template <>
    double
    NSinkerMaterial<3>::inside_outside_factor (const Point<3> &p) const
    {
      double chi = 1.0;

      for (unsigned int s=0; s<n_sinkers; ++s)
        {
          double dist = p.distance(centers[s]);
          double temp = 1-std::exp(-delta*
                                   std::pow(std::max(0.0,dist-omega/2.0),2));
          chi *= temp;
        }
      return chi;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace NSinkerBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(NSinkerMaterial,
                                   "nsinker",
                                   "A material model that corresponds to the 'NSinker' benchmark "
                                   "defined in May et al., G-Cubed, 2015. We here implement the "
                                   "version defined in Rudi et al., SIAM Journal on Scientific "
                                   "Computing, 2017. Number of sinkers and viscosity ratio can be chosen "
                                   "as input parameters.")
  }
}
