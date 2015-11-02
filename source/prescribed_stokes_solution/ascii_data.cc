#include <aspect/global.h>
#include <aspect/prescribed_stokes_solution/ascii_data.h>



namespace aspect
{
  namespace PrescribedStokesSolution
  {
    template <int dim>
    AsciiData<dim>::AsciiData ()
    {}


    template <int dim>
    void
    AsciiData<dim>::initialize ()
    {
      Utilities::AsciiDataInitial<dim>::initialize(dim);
    }


    template <int dim>
    void
    AsciiData<dim>::
    stokes_solution (const Point<dim> &position, Vector<double> &value) const
    {
      value(0) = Utilities::AsciiDataInitial<dim>::get_data_component(position,0);
      value(1) = Utilities::AsciiDataInitial<dim>::get_data_component(position,1);
      if (dim == 3)
        value(2) = Utilities::AsciiDataInitial<dim>::get_data_component(position,2);
      value(dim) = 0;  // makes pressure 0, must set pressure
    }


    template <int dim>
    void
    AsciiData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed Stokes solution");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/prescribed-stokes-solution/",
                                                          "box_2d.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AsciiData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed Stokes solution");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace PrescribedStokesSolution
  {
    ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION(AsciiData,
                                               "ascii data",
                                               "Implementation of a model in which the velocity"
                                               "is derived from files containing data "
                                               "in ascii format. Note the required format of the "
                                               "input data: The first lines may contain any number of comments"
                                               "if they begin with '#', but one of these lines needs to"
                                               "contain the number of grid points in each dimension as"
                                               "for example '# POINTS: 3 3'."
                                               "The order of the data columns "
                                               "has to be 'x', 'y', 'v$_x$' , 'v$_y$' in a 2d model and "
                                               " 'x', 'y', 'z', 'v$_x$' , 'v$_y$' , 'v$_z$' in a 3d model."
                                               "Note that the data in the input "
                                               "files need to be sorted in a specific order: "
                                               "the first coordinate needs to ascend first, "
                                               "followed by the second and the third at last in order to "
                                               "assign the correct data to the prescribed coordinates."
                                               "If you use a spherical model, "
                                               "then the data will still be handled as cartesian,"
                                               "however the assumed grid changes. 'x' will be replaced by "
                                               "the radial distance of the point to the bottom of the model, "
                                               "'y' by the azimuth angle and 'z' by the polar angle measured "
                                               "positive from the north pole. The grid will be assumed to be "
                                               "a latitude-longitude grid. Note that the order "
                                               "of spherical coordinates is 'r', 'phi', 'theta' "
                                               "and not 'r', 'theta', 'phi', since this allows "
                                               "for dimension independent expressions. ")
  }
}
