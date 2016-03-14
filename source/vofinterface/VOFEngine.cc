/*
 * FEM Interface tracking plugin for ASPECT
 *
 * Copyright (C) 2015 Jonathan M Robey
 *
 */

#include <aspect/vofinterface/VOFEngine.h>

namespace aspect
{

  namespace InterfaceTracker
  {
    using namespace dealii;

    template <int dim>
    std::vector<std::string> VOFEngine<dim>::component_names ()
    {
      std::vector<std::string> names (1, "VOF1");
      names.push_back ("vofNormal");
      names.push_back ("vofNormal");
      names.push_back ("vofD");
      return names;
    }

    template <int dim>
    std::vector<DataComponentInterpretation::DataComponentInterpretation> VOFEngine<
    dim>::component_interpretation ()
    {
      std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation (
        1, DataComponentInterpretation::component_is_scalar);
      data_component_interpretation.push_back (
        DataComponentInterpretation::component_is_part_of_vector);
      data_component_interpretation.push_back (
        DataComponentInterpretation::component_is_part_of_vector);
      data_component_interpretation.push_back (
        DataComponentInterpretation::component_is_scalar);
      return data_component_interpretation;
    }

    template <int dim>
    VOFEngine<dim>::VOFEngine ()
      : voleps (std::numeric_limits<double>::quiet_NaN ()),
        interfaceFE (FE_DGP<dim> (0), VOFEngine<dim>::n_components),
        dof_handler (NULL),
        normals_calced (false),
        old_vel_set (false),
        dir_order (0),
        qorder (3),
        triangulation (NULL),
        parent_sim_handler (NULL),
        mapping (NULL),
        communicator (NULL),
        world_size (0),
        self_rank (0)
    {
    }

    template <int dim>
    VOFEngine<dim>::~VOFEngine ()
    {
      clear_dof_handler ();
    }

    template <int dim>
    void VOFEngine<dim>::clear_dof_handler ()
    {
      triangulation = NULL;
      parent_sim_handler = NULL;
      if (dof_handler != NULL)
        {
          delete dof_handler;
          dof_handler = NULL;
        }
    }

    template <int dim>
    void VOFEngine<dim>::init (Function<dim> &initFunc,
                               unsigned int n_samp,
                               double time)
    {
      // Check for required info
      AssertThrow(triangulation!=NULL,
                  ExcMessage ("vofinterface requires triangulation."));
      AssertThrow(parent_sim_handler!=NULL,
                  ExcMessage ("vofinterface requires parent simulation dofs."));

      const QCMid<dim> quadrature (n_samp);
      FEValues<dim> fe_init (*mapping, interfaceFE, quadrature,
                             update_JxW_values | update_quadrature_points);

      // Create location to store required data
      dof_handler = new DoFHandler<dim> (*triangulation);
      dof_handler->distribute_dofs (interfaceFE);
      triangulation->signals.clear.connect (
        std_cxx11::bind (&VOFEngine<dim>::clear_dof_handler,
                         std_cxx11::ref (*this)));

      state.reinit (dof_handler->n_dofs ());
      deltaState.reinit (dof_handler->n_dofs ());

      // Required references
      const unsigned int dofs_per_cell = interfaceFE.dofs_per_cell;
      std::vector<types::global_dof_index> local_dof_indicies (dofs_per_cell);

      AssertThrow(dofs_per_cell == VOFEngine<dim>::n_components,
                  ExcMessage ("VOF must be single valued per cell."));

      double h = 1.0/n_samp;

      initFunc.set_time (time);

      vol_init = 0.0;

      // Initialize state based on provided function
      for (auto cell : dof_handler->active_cell_iterators ())
        {
          // Calculate approximation for volume
          double cell_vol, cell_diam, d_func;
          cell->get_dof_indices (local_dof_indicies);
          cell_vol = cell->measure ();
          cell_diam = cell->diameter();
          d_func = initFunc.value(cell->barycenter());
          fe_init.reinit (cell);

          if (d_func <=-0.5*cell_diam)
            {
              continue;
            }
          if (d_func >= 0.5*cell_diam)
            {
              state (local_dof_indicies[VOFEngine<dim>::vof_ind]) = 1.0;
              vol_init += cell->measure ();
              continue;
            }

          double val = 0.0;
          for (unsigned int i = 0; i < fe_init.n_quadrature_points; ++i)
            {
              double d = 0.0;
              Point<dim> grad;
              Point<dim> xU = quadrature.point (i);
              for (unsigned int di = 0; di < dim; ++di)
                {
                  Point<dim> xH, xL;
                  xH = xU;
                  xL = xU;
                  xH[di] += 0.5*h;
                  xL[di] -= 0.5*h;
                  double dH = initFunc.value(cell->intermediate_point(xH));
                  double dL = initFunc.value(cell->intermediate_point(xL));
                  grad[di] = (dL-dH);
                  d += (0.5/dim)*(dH+dL);
                }
              double ptvof = VOFEngine<dim>::vol_from_d (grad, d);
              val += ptvof * (fe_init.JxW (i) / cell->measure ());
            }
          state (local_dof_indicies[VOFEngine<dim>::vof_ind]) = val;
          vol_init += val*cell->measure ();
        }

      normals_calced = false;
    }

    template <int dim>
    void VOFEngine<dim>::set_voleps (double new_voleps)
    {
      voleps = new_voleps;
    }

    template <int dim>
    void VOFEngine<dim>::set_triangulation (const parallel::distributed::Triangulation<
                                            dim> *new_tria)
    {
      triangulation = new_tria;
    }

    template <int dim>
    void VOFEngine<dim>::set_mapping (const Mapping<dim> *new_mapping)
    {
      mapping = new_mapping;

    }

    template <int dim>
    void VOFEngine<dim>::set_parent_dofs (const DoFHandler<dim> *new_par_dofs)
    {
      parent_sim_handler = new_par_dofs;
    }

    template <int dim>
    void VOFEngine<dim>::set_comm (MPI_Comm new_comm_world)
    {
      communicator = new_comm_world;

      world_size = Utilities::MPI::n_mpi_processes (communicator);
      self_rank = Utilities::MPI::this_mpi_process (communicator);
    }

    template <int dim>
    DoFHandler<dim> *VOFEngine<dim>::get_dof_handler ()
    {
      return dof_handler;
    }

    template <int dim>
    Vector<double> &VOFEngine<dim>::get_state ()
    {
      return state;
    }

    template <>
    double VOFEngine<2>::vol_from_d (Point<2> normal,
                                     double d)
    {
      //Handle basic recursion case
      const int dim = 2;
      double norm1, mpos, dtest;

      //Get 1-Norm
      norm1 = 0.0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          norm1 += numbers::NumberTraits<double>::abs (normal (i));
        }

      //Find min component
      mpos = 0.6;
      for (unsigned int i = 0; i < dim; ++i)
        {
          double mcand = numbers::NumberTraits<double>::abs (normal (i))
                         / norm1;
          mpos = (mcand < mpos) ? mcand : mpos;
        }

      //Obtain volume
      dtest = d / norm1;
      if (dtest <= -0.5)
        {
          return 0.0;
        }
      if (dtest >= 0.5)
        {
          return 1.0;
        }
      if (dtest < mpos - 0.5)
        {
          return (dtest + 0.5) * (dtest + 0.5) / (2.0*mpos * (1.0 - mpos));
        }
      if (dtest > 0.5 - mpos)
        {
          return 1.0 - (dtest - 0.5) * (dtest - 0.5) / (2.0*mpos * (1.0 - mpos));
        }
      return 0.5 + dtest / (1.0 - mpos);
    }

    template <>
    double VOFEngine<2>::d_from_vol (Point<2> normal,
                                     double vol)
    {
      //Handle basic recursion case
      const int dim = 2;
      double norm1, mpos;

      //Get 1-Norm
      norm1 = 0.0;
      for (unsigned int i = 0; i < dim; ++i)
        {
          norm1 += numbers::NumberTraits<double>::abs (normal (i));
        }

      //Find min component
      mpos = 0.6;
      for (unsigned int i = 0; i < dim; ++i)
        {
          double mcand = numbers::NumberTraits<double>::abs (normal (i))
                         / norm1;
          mpos = (mcand < mpos) ? mcand : mpos;
        }

      //Obtain const
      if (vol <= 0.0)
        {
          return -0.5 * norm1;
        }
      if (vol >= 1.0)
        {
          return 0.5 * norm1;
        }
      if (vol < 0.5 * mpos / (1 - mpos))
        {
          return norm1 * (-0.5 + sqrt (2.0*vol * mpos * (1 - mpos)));
        }
      if (vol > 1.0 - 0.5 * mpos / (1 - mpos))
        {
          return norm1 * (0.5 - sqrt (2.0*(1.0 - vol) * mpos * (1 - mpos)));
        }
      return norm1 * (1 - mpos) * (vol - 0.5);
    }

    template <>
    double VOFEngine<3>::vol_from_d (Point<3> normal,
                                     double d)
    {
      // 3D vol calculation not yet implemented
      return 0.0;
    }

    template <>
    double VOFEngine<3>::d_from_vol (Point<3> normal,
                                     double vol)
    {
      // 3D interface location calculation not yet implemented
      return 0.0;
    }

    template <int dim>
    void VOFEngine<dim>::do_step (const LinearAlgebra::BlockVector &par_new_soln,
                                  const LinearAlgebra::BlockVector &par_old_soln,
                                  double timestep)
    {
      if (dir_order == 0)
        {
          for (unsigned int dir = 0; dir < dim; ++dir)
            {
              calc_flux (par_new_soln, par_old_soln, timestep, dir);
              update_vof (par_new_soln, par_old_soln, timestep, dir);
            }
        }
      else
        {
          for (unsigned int dir = 0; dir < dim; ++dir)
            {
              calc_flux (par_new_soln, par_old_soln, timestep, dim - dir - 1);
              update_vof (par_new_soln, par_old_soln, timestep,
                          dim - dir - 1);
            }
        }
      dir_order = (dir_order + 1) % 2;
      old_vel_set = true;
    }

    template <>
    void VOFEngine<2>::calc_normals ()
    {
      const int dim = 2;
      if (normals_calced)
        return;

      // Boundary reference
      typename DoFHandler<dim>::active_cell_iterator endc =
        dof_handler->end ();

      // Interface Reconstruction vars
      const unsigned int n_local = 9;

      Vector<double> local_vofs (n_local);
      std::vector<Point<dim>> resc_cell_centers (n_local);

      const unsigned int n_sums = 3;
      std::vector<double> strip_sums (dim * n_sums);

      const unsigned int n_normals = 6;
      std::vector<Point<dim>> normals (n_normals);
      std::vector<double> errs (2 * n_normals);

      // Normal holding vars
      Point<dim> normal;
      double d;

      //Iterate over cells
      DoFHandler<dim>::active_cell_iterator par_cell =
        parent_sim_handler->begin_active ();
      for (auto cell : dof_handler->active_cell_iterators ())
        {
          std::vector<types::global_dof_index> cell_dof_indicies (
            interfaceFE.dofs_per_cell);
          types::global_dof_index cell_vof_index, cell_normal_x_index,
                cell_normal_y_index, cell_d_index;
          double cell_vof, cell_vol;

          // Obtain data for this cell and neighbors
          cell->get_dof_indices (cell_dof_indicies);
          cell_vof_index = cell_dof_indicies[VOFEngine<dim>::vof_ind];
          cell_normal_x_index =
            cell_dof_indicies[VOFEngine<dim>::first_normal_ind];
          cell_normal_y_index =
            cell_dof_indicies[VOFEngine<dim>::first_normal_ind + 1];
          cell_d_index = cell_dof_indicies[VOFEngine<dim>::d_ind];
          cell_vof = state (cell_vof_index);
          cell_vol = cell->measure ();

          normal[0] = 0.0;
          normal[1] = 0.0;
          d = -0.5;

          if (cell_vof > 1.0 - voleps)
            {
              d = 0.5;
            }
          else if (cell_vof > voleps)
            {
              //Identify best normal
              // Get references to neighboring cells
              for (unsigned int i = 0; i < 3; ++i)
                {
                  typename DoFHandler<dim>::active_cell_iterator cen;
                  if (i == 0)
                    cen = cell->neighbor (0);
                  if (i == 1)
                    cen = cell;
                  if (i == 2)
                    cen = cell->neighbor (1);
                  for (unsigned int j = 0; j < 3; ++j)
                    {
                      typename DoFHandler<dim>::active_cell_iterator curr;
                      if (cen == endc)
                        {
                          curr = endc;
                        }
                      else
                        {
                          if (j == 0)
                            curr = cen->neighbor (2);
                          if (j == 1)
                            curr = cen;
                          if (j == 2)
                            curr = cen->neighbor (3);
                        }
                      if (curr != endc)
                        {
                          curr->get_dof_indices (cell_dof_indicies);
                          resc_cell_centers[3 * j + i] = Point<dim> (-1.0 + i,
                                                                     -1.0 + j);
                        }
                      else
                        {
                          cell->get_dof_indices (cell_dof_indicies);
                          resc_cell_centers[3 * j + i] = Point<dim> (0.0,
                                                                     0.0);
                        }
                      local_vofs (3 * j + i) = state (
                                                 cell_dof_indicies[VOFEngine<dim>::vof_ind]);
                    }
                }
              // Gather cell strip sums
              for (unsigned int i = 0; i < dim * n_sums; ++i)
                strip_sums[i] = 0.0;

              for (unsigned int i = 0; i < 3; ++i)
                {
                  for (unsigned int j = 0; j < 3; ++j)
                    {
                      strip_sums[3 * 0 + i] += local_vofs (3 * j + i);
                      strip_sums[3 * 1 + j] += local_vofs (3 * j + i);
                    }
                }

              // std::cout << "  Strip sums: " << std::endl;
              // for (unsigned int i = 0; i < dim * n_sums; ++i)
              //   std::cout << "    " << strip_sums[i] << std::endl;

              for (unsigned int di = 0; di < dim; ++di)
                {
                  unsigned int di2 = (di + 1) % dim;
                  for (unsigned int i = 0; i < 3; ++i)
                    {
                      normals[3 * di + i][di] = 0.0;
                      normals[3 * di + i][di2] = 0.0;
                      if (i % 2 == 0)
                        {
                          //use low sum
                          normals[3 * di + i][di] += strip_sums[3 * di + 0];
                          normals[3 * di + i][di2] += 1.0;
                        }
                      else
                        {
                          //use high sum
                          normals[3 * di + i][di] += strip_sums[3 * di + 1];
                          normals[3 * di + i][di2] += 0.0;
                        }
                      if (i == 0)
                        {
                          //use low sum
                          normals[3 * di + i][di] -= strip_sums[3 * di + 1];
                          normals[3 * di + i][di2] += 0.0;
                        }
                      else
                        {
                          //use high sum
                          normals[3 * di + i][di] -= strip_sums[3 * di + 2];
                          normals[3 * di + i][di2] += 1.0;
                        }
                      if (strip_sums[3 * di2 + 2] > strip_sums[3 * di2 + 0])
                        normals[3 * di + i][di2] *= -1.0;
                    }
                }

              unsigned int mn_ind = 0;
              for (unsigned int nind = 0; nind < n_normals; ++nind)
                {
                  errs[nind] = 0.0;
                  d = VOFEngine<dim>::d_from_vol (normals[nind], cell_vof);
                  for (unsigned int i = 0; i < n_local; ++i)
                    {
                      double dot = 0.0;
                      for (unsigned int di = 0; di < dim; di++)
                        dot += normals[nind][di] * resc_cell_centers[i][di];
                      double val = local_vofs (i)
                                   - VOFEngine<dim>::vol_from_d (normals[nind],
                                                                 d - dot);
                      errs[nind] += val * val;
                    }
                  if (errs[mn_ind] >= errs[nind])
                    mn_ind = nind;
                  // std::cout << "   " << normals[nind] << " e ";
                  // std::cout  << errs[nind] << " " << mn_ind << std::endl;
                }

              normal = normals[mn_ind];
              d = VOFEngine<dim>::d_from_vol (normal, cell_vof);
            }

          double n2 = sqrt (normal * normal);
          if (n2 > voleps)
            {
              normal = (normal / n2);
              d = VOFEngine<dim>::d_from_vol (normal, cell_vof);
            }
          else
            {
              normal[0] = 0.0;
              normal[1] = 0.0;
            }
          state (cell_normal_x_index) = normal[0];
          state (cell_normal_y_index) = normal[1];
          state (cell_d_index) = d;
        }

      normals_calced = true;
    }

    template <>
    void VOFEngine<3>::calc_normals ()
    {
        // 3D interface reconstruction not yet implemented
    }


    template <>
    void VOFEngine<2>::calc_flux (const LinearAlgebra::BlockVector &par_new_soln,
                                  const LinearAlgebra::BlockVector &par_old_soln,
                                  double timestep,
                                  unsigned int dir)
    {
      const int dim = 2;

      // Boundary reference
      typename DoFHandler<dim>::active_cell_iterator endc =
        dof_handler->end ();

      // Interface holding vars

      calc_normals ();
      Point<dim> normal;
      double d;

      // Flux Computation Vars

      const unsigned int n_neighbor = 2;
      Vector<double> dvofs (n_neighbor);

      Functions::FEFieldFunction<dim, DoFHandler<dim>,
                LinearAlgebra::BlockVector> par_new_state (*parent_sim_handler,
                                                           par_new_soln, *mapping);
      Functions::FEFieldFunction<dim, DoFHandler<dim>,
                LinearAlgebra::BlockVector> par_old_state (*parent_sim_handler,
                                                           par_old_soln, *mapping);
      const QGauss<dim - 1> quadrature (qorder);
      FEFaceValues<dim> fe_face (*mapping, interfaceFE, quadrature,
                                 update_JxW_values | update_quadrature_points
                                 | update_normal_vectors);
      std::vector<Vector<double>> nvels (fe_face.n_quadrature_points);
      std::vector<Vector<double>> ovels (fe_face.n_quadrature_points);

      Vector<double> flux (n_neighbor);
      Vector<double> flux_vof (n_neighbor);

      deltaState = 0.0;

      // std::cout << "Dir: " << dir << std::endl;

      //Iterate over cells
      DoFHandler<dim>::active_cell_iterator par_cell =
        parent_sim_handler->begin_active ();
      for (auto cell : dof_handler->active_cell_iterators ())
        {
          std::vector<types::global_dof_index> cell_dof_indicies (
            interfaceFE.dofs_per_cell);
          types::global_dof_index cell_vof_index, cell_normal_x_index,
                cell_normal_y_index, cell_d_index;
          double cell_vof, cell_vol;

          // Obtain data for this cell and neighbors
          cell->get_dof_indices (cell_dof_indicies);
          cell_vof_index = cell_dof_indicies[VOFEngine<dim>::vof_ind];
          cell_normal_x_index =
            cell_dof_indicies[VOFEngine<dim>::first_normal_ind];
          cell_normal_y_index =
            cell_dof_indicies[VOFEngine<dim>::first_normal_ind + 1];
          cell_d_index = cell_dof_indicies[VOFEngine<dim>::d_ind];
          cell_vof = state (cell_vof_index);
          cell_vol = cell->measure ();

          if (cell_vof > voleps)
            {
              // std::cout << "Cen: " << cell->center () << std::endl;
              // std::cout << "  Ind: " << cell_vof_index << std::endl;
              // std::cout << "  VOF: " << cell_vof << std::endl;

              normal[0] = state (cell_normal_x_index);
              normal[1] = state (cell_normal_y_index);
              d = state (cell_d_index);

              if (cell_vof > 1.0 - voleps)
                {
                  normal[0] = 0.0;
                  normal[1] = 1.0;
                }

              // std::cout << "  N,d: " << normal << ", " << d << std::endl;
              // std::cout << "  Vol Fluxes: " << std::endl;

              par_new_state.set_active_cell (par_cell);
              par_old_state.set_active_cell (par_cell);
              for (unsigned int i = 0; i < 2; ++i)
                {
                  fe_face.reinit (cell, 2 * dir + i);
                  flux[i] = 0.0;
                  const std::vector<double> &jxw = fe_face.get_JxW_values ();
                  const std::vector<Tensor<1, dim>> &normals =
                                                   fe_face.get_all_normal_vectors ();
                  par_new_state.vector_value_list (
                    fe_face.get_quadrature_points (), nvels);
                  par_old_state.vector_value_list (
                    fe_face.get_quadrature_points (), ovels);
                  for (unsigned int j = 0; j < fe_face.n_quadrature_points;
                       ++j)
                    {
                      Point<dim> vel;
                      Point<dim> vGrad;
                      for (unsigned int vdi = 0; vdi < dim; ++vdi)
                        {
                          if (old_vel_set)
                            {
                              vel[vdi] = 0.5 * (nvels[j][vdi] + ovels[j][vdi]);
                            }
                          else
                            {
                              vel[vdi] = nvels[j][vdi];
                            }
                        }
                      flux[i] += (jxw[j]/cell_vol) * (normals[j] * (timestep*vel));
                    }
                  // std::cout << "    " << i << " " << flux[i] << std::endl;
                }

              // std::cout << "  Flux Fracs: " << std::endl;

              for (unsigned int i = 0; i < 2; ++i)
                {
                  double adj = flux[i];
                  if (adj > voleps)
                    {
                      double dot = normal[dir] * (1.0 - adj);
                      dot *= 0.5 * (i == 0 ? -1.0 : 1.0);

                      Point<dim> adj_norm = normal;
                      adj_norm[dir] = normal[dir] * adj;

                      // std::cout << "     " << adj_norm << "*x=" << d-dot << std::endl;

                      flux_vof[i] = VOFEngine<dim>::vol_from_d (adj_norm,
                                                                d - dot);
                    }
                  else
                    {
                      flux_vof[i] = 0.0;
                    }
                  // std::cout << "    " << i << " " << flux_vof[i] << std::endl;
                }

              // std::cout << "  Fluxes:" << std::endl;

              for (unsigned int i = 0; i < n_neighbor; ++i)
                {
                  dvofs[i] = flux_vof[i] * flux[i];

                  // std::cout << "    " << i << " " << dvofs[i] << std::endl;
                }

              // Set valid fluxes
              for (unsigned int i = 0; i < n_neighbor; ++i)
                {
                  deltaState (cell_vof_index) -= dvofs[i];
                  typename DoFHandler<dim>::active_cell_iterator curr;
                  curr = cell->neighbor (2 * dir + i);
                  if (curr != endc)
                    {
                      curr->get_dof_indices (cell_dof_indicies);
                      double curr_vol = curr->measure ();
                      types::global_dof_index curr_vof_ind;
                      curr_vof_ind = cell_dof_indicies[VOFEngine<dim>::vof_ind];
                      deltaState (curr_vof_ind) += dvofs[i] * (cell_vol / curr_vol);
                    }
                }
            }
          ++par_cell;
        }
    }

    template <>
    void VOFEngine<3>::calc_flux (const LinearAlgebra::BlockVector &par_new_soln,
                                  const LinearAlgebra::BlockVector &par_old_soln,
                                  double timestep,
                                  unsigned int dir)
    {
        // 3D flux calculation not yet implemented, may be later handled with
        // more general interface code
    }

    template <int dim>
    void VOFEngine<dim>::update_vof (const LinearAlgebra::BlockVector &par_new_soln,
                                     const LinearAlgebra::BlockVector &par_old_soln,
                                     double timestep,
                                     unsigned int dir)
    {
      bool last_update = false;

      Functions::FEFieldFunction<dim, DoFHandler<dim>,
                LinearAlgebra::BlockVector> par_new_state (*parent_sim_handler,
                                                           par_new_soln, *mapping);
      Functions::FEFieldFunction<dim, DoFHandler<dim>,
                LinearAlgebra::BlockVector> par_old_state (*parent_sim_handler,
                                                           par_old_soln, *mapping);
      const QGauss<dim - 1> quadrature (qorder);
      FEFaceValues<dim> fe_face (*mapping, interfaceFE, quadrature,
                                 update_JxW_values | update_quadrature_points
                                 | update_normal_vectors);
      std::vector<Vector<double>> nvels (fe_face.n_quadrature_points);
      std::vector<Vector<double>> ovels (fe_face.n_quadrature_points);

      last_update = (dir == (dim - 1) * dir_order);

      typename DoFHandler<dim>::active_cell_iterator par_cell =
        parent_sim_handler->begin_active ();
      // double term, sum;
      // sum = 0.0;
      // vofs += deltaVOF;

      for (auto cell : dof_handler->active_cell_iterators ())
        {
          double cell_vol = cell->measure ();

          std::vector<types::global_dof_index> cell_dof_indicies (
            interfaceFE.dofs_per_cell);

          cell->get_dof_indices (cell_dof_indicies);
          types::global_dof_index cell_vof_ind = cell_dof_indicies[VOFEngine<dim>::vof_ind];

          par_new_state.set_active_cell (par_cell);
          par_old_state.set_active_cell (par_cell);
          double dflux = 0.0;
          for (unsigned int i = 0; i < 2; ++i)
            {
              fe_face.reinit (cell, 2 * dir + i);
              const std::vector<double> &jxw = fe_face.get_JxW_values ();
              const std::vector<Tensor<1, dim>> &normals =
                                               fe_face.get_all_normal_vectors ();
              par_new_state.vector_value_list (
                fe_face.get_quadrature_points (), nvels);
              par_old_state.vector_value_list (
                fe_face.get_quadrature_points (), ovels);
              for (unsigned int j = 0; j < fe_face.n_quadrature_points; ++j)
                {
                  Point<dim> vel;
                  Point<dim> vGrad;
                  for (unsigned int vdi = 0; vdi < dim; ++vdi)
                    {
                      vel[vdi] = 0.5 * (nvels[j][vdi] + ovels[j][vdi]);
                    }
                  dflux += jxw[j] * normals[j] * vel;
                }
            }
          dflux *= timestep / cell_vol;

          state(cell_vof_ind) += deltaState(cell_vof_ind);
          //state (cell_vof_ind) = (state (cell_vof_ind) * (1 + 0.5 * dflux)
          //    + deltaState (cell_vof_ind))
          //                       / (1 - 0.5 * dflux);

          ++par_cell;
        }
      normals_calced = false;
    }

    template <int dim>
    std::vector<std::string> VOFEngine<dim>::error_names ()
    {
      std::vector<std::string> names (1, "L1 Interface Error");
      names.push_back ("L1 Field Error");
      names.push_back ("Volume delta");
      return names;
    }


    template <int dim>
    std::vector<std::string> VOFEngine<dim>::error_abrev ()
    {
      std::vector<std::string> names (1, "I_EL1");
      names.push_back ("F_EL1");
      names.push_back ("D_vol");
      return names;
    }

    template <int dim>
    std::vector<double> VOFEngine<dim>::calc_error (Function<dim> &func,
                                                    unsigned int n_samp,
                                                    double time)
    {
      double I_err_est = 0.0;
      double F_err_est = 0.0;
      double curr_vol = 0.0;
      double h = 1.0/n_samp;

      Point<dim> xU;
      Point<dim> normal;
      double d;

      const unsigned int dofs_per_cell = interfaceFE.dofs_per_cell;
      std::vector<types::global_dof_index> local_dof_indicies (dofs_per_cell);

      const QCMid<dim> quadrature (n_samp);
      FEValues<dim> fe_err (*mapping, interfaceFE, quadrature,
                            update_JxW_values | update_quadrature_points);

      calc_normals ();

      func.set_time (time);

      // Initialize state based on provided function
      for (auto cell : dof_handler->active_cell_iterators ())
        {

          std::vector<types::global_dof_index> cell_dof_indicies (
            interfaceFE.dofs_per_cell);
          types::global_dof_index cell_vof_index, cell_normal_x_index,
                cell_normal_y_index, cell_d_index;
          double cell_vof, cell_vol;
          double cell_diam;
          double d_func;

          // Obtain data for this cell and neighbors
          cell->get_dof_indices (local_dof_indicies);
          cell_vof_index = local_dof_indicies[VOFEngine<dim>::vof_ind];
          cell_normal_x_index =
            local_dof_indicies[VOFEngine<dim>::first_normal_ind];
          cell_normal_y_index =
            local_dof_indicies[VOFEngine<dim>::first_normal_ind + 1];
          cell_d_index = local_dof_indicies[VOFEngine<dim>::d_ind];

          // Calculate approximation for volume
          fe_err.reinit (cell);

          cell_vol = cell->measure ();
          cell_diam = cell->diameter();
          d_func = func.value(cell->barycenter());
          cell_vof = state (cell_vof_index);
          normal[0] = state (cell_normal_x_index);
          normal[1] = state (cell_normal_y_index);
          d = state (cell_d_index);

          if (numbers::NumberTraits<double>::abs (d) >= 0.5 &&
              numbers::NumberTraits<double>::abs (d_func) >= 0.5*cell_diam &&
              d*d_func>=0.0)
            {
              if (d >= 0.5)
                {
                  curr_vol += cell->measure ();
                }
              continue;
            }

          Point<dim> grad_t;
          double d_t;
          double val = 0.0;
          double vof_reinit = 0.0;
          for (unsigned int i = 0; i < fe_err.n_quadrature_points; ++i)
            {
              d_t = 0.0;
              xU = quadrature.point (i);
              for (unsigned int di = 0; di < dim; ++di)
                {
                  Point<dim> xH, xL;
                  xH = xU;
                  xL = xU;
                  xH[di] += 0.5*h;
                  xL[di] -= 0.5*h;
                  double dH = func.value(cell->intermediate_point(xH));
                  double dL = func.value(cell->intermediate_point(xL));
                  grad_t[di] = (dL-dH);
                  d_t += (0.5/dim)*(dH+dL);
                }
              double ptvof_t = VOFEngine<dim>::vol_from_d (grad_t, d_t);
              vof_reinit += ptvof_t *(fe_err.JxW (i)/cell->measure ());
              for (unsigned int di = 0; di < dim; ++di)
                {
                  xU[di] -= 0.5;
                }
              double dot = normal * xU;
              double ptvof_e = VOFEngine<dim>::vol_from_d (h*normal,
                                                           (d - dot));
              double diff = numbers::NumberTraits<double>::abs (ptvof_t - ptvof_e);
              val += diff * fe_err.JxW (i);
            }
          I_err_est += val;
          F_err_est += numbers::NumberTraits<double>::abs (cell_vof-vof_reinit)
                       *cell->measure ();
          curr_vol += vof_reinit*cell->measure ();
        }

      std::vector<double> errors(1, I_err_est);
      errors.push_back (F_err_est);
      errors.push_back (curr_vol-vol_init);
      return errors;
    }

    template class VOFEngine<2> ;
    template class VOFEngine<3> ;
  }
}
