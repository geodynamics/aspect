/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/compat.h>

// deal.II versions up to 9.5 had a poorly designed interface of the
// SphericalManifold class that made it impossible for us to use.
// This file thus contains a copy of it.
#if !DEAL_II_VERSION_GTE(9,6,0)

namespace aspect
{
  namespace
  {
    template <int dim, int spacedim>
    bool
    spherical_face_is_horizontal(
      const typename Triangulation<dim, spacedim>::face_iterator &face,
      const Point<spacedim>                                      &manifold_center)
    {
      // We test whether a face is horizontal by checking that the vertices
      // all have roughly the same distance from the center: If the
      // maximum deviation for the distances from the vertices to the
      // center is less than 1.e-5 of the distance between vertices (as
      // measured by the minimum distance from any of the other vertices
      // to the first vertex), then we call this a horizontal face.
      constexpr unsigned int n_vertices =
        GeometryInfo<spacedim>::vertices_per_face;
      std::array<double, n_vertices>     sqr_distances_to_center;
      std::array<double, n_vertices - 1> sqr_distances_to_first_vertex;
      sqr_distances_to_center[0] =
        (face->vertex(0) - manifold_center).norm_square();
      for (unsigned int i = 1; i < n_vertices; ++i)
        {
          sqr_distances_to_center[i] =
            (face->vertex(i) - manifold_center).norm_square();
          sqr_distances_to_first_vertex[i - 1] =
            (face->vertex(i) - face->vertex(0)).norm_square();
        }
      const auto minmax_sqr_distance =
        std::minmax_element(sqr_distances_to_center.begin(),
                            sqr_distances_to_center.end());
      const auto min_sqr_distance_to_first_vertex =
        std::min_element(sqr_distances_to_first_vertex.begin(),
                         sqr_distances_to_first_vertex.end());

      return (*minmax_sqr_distance.second - *minmax_sqr_distance.first <
              1.e-10 * *min_sqr_distance_to_first_vertex);
    }
  } // namespace


  // ============================================================
  // SphericalManifold
  // ============================================================

  template <int dim, int spacedim>
  SphericalManifold<dim, spacedim>::SphericalManifold(
    const Point<spacedim> center)
    : center(center)
    , polar_manifold(center)
  {}



  template <int dim, int spacedim>
  std::unique_ptr<Manifold<dim, spacedim>>
  SphericalManifold<dim, spacedim>::clone() const
  {
    return std::make_unique<SphericalManifold<dim, spacedim>>(*this);
  }



  template <int dim, int spacedim>
  Point<spacedim>
  SphericalManifold<dim, spacedim>::get_intermediate_point(
    const Point<spacedim> &p1,
    const Point<spacedim> &p2,
    const double           w) const
  {
    const double tol = 1e-10;

    if ((p1 - p2).norm_square() < tol * tol || std::abs(w) < tol)
      return p1;
    else if (std::abs(w - 1.0) < tol)
      return p2;

    // If the points are one dimensional then there is no need for anything but
    // a linear combination.
    if (spacedim == 1)
      return Point<spacedim>(w * p2 + (1 - w) * p1);

    const Tensor<1, spacedim> v1 = p1 - center;
    const Tensor<1, spacedim> v2 = p2 - center;
    const double              r1 = v1.norm();
    const double              r2 = v2.norm();

    Assert(r1 > tol && r2 > tol,
           ExcMessage("p1 and p2 cannot coincide with the center."));

    const Tensor<1, spacedim> e1 = v1 / r1;
    const Tensor<1, spacedim> e2 = v2 / r2;

    // Find the cosine of the angle gamma described by v1 and v2.
    const double cosgamma = e1 * e2;

    // Points are collinear with the center (allow for 8*eps as a tolerance)
    if (cosgamma < -1 + 8. * std::numeric_limits<double>::epsilon())
      return center;

    // Points are along a line, in which case e1 and e2 are essentially the same.
    if (cosgamma > 1 - 8. * std::numeric_limits<double>::epsilon())
      return Point<spacedim>(center + w * v2 + (1 - w) * v1);

    // Find the angle sigma that corresponds to arclength equal to w. acos
    // should never be undefined because we have ruled out the two special cases
    // above.
    const double sigma = w * std::acos(cosgamma);

    // Normal to v1 in the plane described by v1,v2,and the origin.
    // Since p1 and p2 do not coincide n is not zero and well defined.
    Tensor<1, spacedim> n      = v2 - (v2 * e1) * e1;
    const double        n_norm = n.norm();
    Assert(n_norm > 0,
           ExcInternalError("n should be different from the null vector. "
                            "Probably, this means v1==v2 or v2==0."));

    n /= n_norm;

    // Find the point Q along O,v1 such that
    // P1,V,P2 has measure sigma.
    const Tensor<1, spacedim> P = std::cos(sigma) * e1 + std::sin(sigma) * n;

    // Project this point on the manifold.
    return Point<spacedim>(center + (w * r2 + (1.0 - w) * r1) * P);
  }



  template <int dim, int spacedim>
  Tensor<1, spacedim>
  SphericalManifold<dim, spacedim>::get_tangent_vector(
    const Point<spacedim> &p1,
    const Point<spacedim> &p2) const
  {
    const double tol = 1e-10;
    (void)tol;

    Assert(p1 != p2, ExcMessage("p1 and p2 should not concide."));

    const Tensor<1, spacedim> v1 = p1 - center;
    const Tensor<1, spacedim> v2 = p2 - center;
    const double              r1 = v1.norm();
    const double              r2 = v2.norm();

    Assert(r1 > tol, ExcMessage("p1 cannot coincide with the center."));

    Assert(r2 > tol, ExcMessage("p2 cannot coincide with the center."));

    const Tensor<1, spacedim> e1 = v1 / r1;
    const Tensor<1, spacedim> e2 = v2 / r2;

    // Find the cosine of the angle gamma described by v1 and v2.
    const double cosgamma = e1 * e2;

    Assert(cosgamma > -1 + 8. * std::numeric_limits<double>::epsilon(),
           ExcMessage("p1 and p2 cannot lie on the same diameter and be opposite "
                      "respect to the center."));

    if (cosgamma > 1 - 8. * std::numeric_limits<double>::epsilon())
      return v2 - v1;

    // Normal to v1 in the plane described by v1,v2,and the origin.
    // Since p1 and p2 do not coincide n is not zero and well defined.
    Tensor<1, spacedim> n      = v2 - (v2 * e1) * e1;
    const double        n_norm = n.norm();
    Assert(n_norm > 0,
           ExcInternalError("n should be different from the null vector. "
                            "Probably, this means v1==v2 or v2==0."));

    n /= n_norm;

    // this is the derivative of the geodesic in get_intermediate_point
    // derived with respect to w and inserting w=0.
    const double gamma = std::acos(cosgamma);
    return (r2 - r1) * e1 + r1 * gamma * n;
  }



  template <int dim, int spacedim>
  Tensor<1, spacedim>
  SphericalManifold<dim, spacedim>::normal_vector(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    const Point<spacedim>                                      &p) const
  {
    // Let us first test whether we are on a "horizontal" face
    // (tangential to the sphere).  In this case, the normal vector is
    // easy to compute since it is proportional to the vector from the
    // center to the point 'p'.
    if (spherical_face_is_horizontal<dim, spacedim>(face, center))
      {
        // So, if this is a "horizontal" face, then just compute the normal
        // vector as the one from the center to the point 'p', adequately
        // scaled.
        const Tensor<1, spacedim> unnormalized_spherical_normal = p - center;
        const Tensor<1, spacedim> normalized_spherical_normal =
          unnormalized_spherical_normal / unnormalized_spherical_normal.norm();
        return normalized_spherical_normal;
      }
    else
      // If it is not a horizontal face, just use the machinery of the
      // base class.
      return Manifold<dim, spacedim>::normal_vector(face, p);

    return Tensor<1, spacedim>();
  }



  template <>
  void
  SphericalManifold<1, 1>::get_normals_at_vertices(
    const Triangulation<1>::face_iterator &,
    Manifold<1, 1>::FaceVertexNormals &) const
  {
    Assert(false, ExcImpossibleInDim(1));
  }



  template <>
  void
  SphericalManifold<1, 2>::get_normals_at_vertices(
    const Triangulation<1, 2>::face_iterator &,
    Manifold<1, 2>::FaceVertexNormals &) const
  {
    Assert(false, ExcImpossibleInDim(1));
  }



  template <int dim, int spacedim>
  void
  SphericalManifold<dim, spacedim>::get_normals_at_vertices(
    const typename Triangulation<dim, spacedim>::face_iterator &face,
    typename Manifold<dim, spacedim>::FaceVertexNormals &face_vertex_normals)
  const
  {
    // Let us first test whether we are on a "horizontal" face
    // (tangential to the sphere).  In this case, the normal vector is
    // easy to compute since it is proportional to the vector from the
    // center to the point 'p'.
    if (spherical_face_is_horizontal<dim, spacedim>(face, center))
      {
        // So, if this is a "horizontal" face, then just compute the normal
        // vector as the one from the center to the point 'p', adequately
        // scaled.
        for (unsigned int vertex = 0;
             vertex < GeometryInfo<spacedim>::vertices_per_face;
             ++vertex)
          face_vertex_normals[vertex] = face->vertex(vertex) - center;
      }
    else
      Manifold<dim, spacedim>::get_normals_at_vertices(face, face_vertex_normals);
  }



  template <int dim, int spacedim>
  void
  SphericalManifold<dim, spacedim>::get_new_points(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    const Table<2, double>                 &weights,
    ArrayView<Point<spacedim>>              new_points) const
  {
    AssertDimension(new_points.size(), weights.size(0));
    AssertDimension(surrounding_points.size(), weights.size(1));

    do_get_new_points(surrounding_points, make_array_view(weights), new_points);

    return;
  }



  template <int dim, int spacedim>
  Point<spacedim>
  SphericalManifold<dim, spacedim>::get_new_point(
    const ArrayView<const Point<spacedim>> &vertices,
    const ArrayView<const double>          &weights) const
  {
    // To avoid duplicating all of the logic in get_new_points, simply call it
    // for one position.
    Point<spacedim> new_point;
    do_get_new_points(vertices,
                      weights,
                      make_array_view(&new_point, &new_point + 1));

    return new_point;
  }



  namespace internal
  {
    // Rotate a given unit vector u around the axis dir
    // where the angle is given by the length of dir.
    // This is the exponential map for a sphere.
    Tensor<1, 3>
    apply_exponential_map(const Tensor<1, 3> &u, const Tensor<1, 3> &dir)
    {
      const double theta = dir.norm();
      if (theta < 1.e-10)
        {
          return u;
        }
      else
        {
          const Tensor<1, 3> dirUnit = dir / theta;
          const Tensor<1, 3> tmp =
            std::cos(theta) * u + std::sin(theta) * dirUnit;
          return tmp / tmp.norm();
        }
    }


    // Returns the direction to go from v to u
    // projected to the plane perpendicular to the unit vector v.
    // This one is more stable when u and v are nearly equal.
    Tensor<1, 3>
    projected_direction(const Tensor<1, 3> &u, const Tensor<1, 3> &v)
    {
      Tensor<1, 3> ans = u - v;
      ans -= (ans * v) * v;
      return ans; // ans = (u-v) - ((u-v)*v)*v
    }


    // helper function to compute a vector orthogonal to a given one.
    // does nothing unless spacedim == 3.
    template <int spacedim>
    Point<spacedim>
    compute_normal(const Tensor<1, spacedim> & /*vector*/,
                   bool /*normalize*/ = false)
    {
      return {};
    }

    Point<3>
    compute_normal(const Tensor<1, 3> &vector, bool normalize = false)
    {
      Assert(vector.norm_square() != 0.,
             ExcMessage("The direction parameter must not be zero!"));
      Point<3> normal;
      if (std::abs(vector[0]) >= std::abs(vector[1]) &&
          std::abs(vector[0]) >= std::abs(vector[2]))
        {
          normal[1] = -1.;
          normal[2] = -1.;
          normal[0] = (vector[1] + vector[2]) / vector[0];
        }
      else if (std::abs(vector[1]) >= std::abs(vector[0]) &&
               std::abs(vector[1]) >= std::abs(vector[2]))
        {
          normal[0] = -1.;
          normal[2] = -1.;
          normal[1] = (vector[0] + vector[2]) / vector[1];
        }
      else
        {
          normal[0] = -1.;
          normal[1] = -1.;
          normal[2] = (vector[0] + vector[1]) / vector[2];
        }
      if (normalize)
        normal /= normal.norm();
      return normal;
    }


    namespace SphericalManifold
    {
      namespace
      {
        template <int spacedim>
        Point<spacedim>
        do_get_new_point(
          const ArrayView<const Tensor<1, spacedim>> & /*directions*/,
          const ArrayView<const double> & /*distances*/,
          const ArrayView<const double> & /*weights*/,
          const Point<spacedim> & /*candidate_point*/)
        {
          Assert(false, ExcNotImplemented());
          return {};
        }

        template <>
        Point<3>
        do_get_new_point(const ArrayView<const Tensor<1, 3>> &directions,
                         const ArrayView<const double>       &distances,
                         const ArrayView<const double>       &weights,
                         const Point<3>                      &candidate_point)
        {
          (void)distances;

          AssertDimension(directions.size(), distances.size());
          AssertDimension(directions.size(), weights.size());

          Point<3>           candidate       = candidate_point;
          const unsigned int n_merged_points = directions.size();
          const double       tolerance       = 1e-10;
          const int          max_iterations  = 10;

          {
            // If the candidate happens to coincide with a normalized
            // direction, we return it. Otherwise, the Hessian would be singular.
            for (unsigned int i = 0; i < n_merged_points; ++i)
              {
                const double squared_distance =
                  (candidate - directions[i]).norm_square();
                if (squared_distance < tolerance * tolerance)
                  return candidate;
              }

            // check if we only have two points now, in which case we can use the
            // get_intermediate_point function
            if (n_merged_points == 2)
              {
                static const dealii::SphericalManifold<3, 3> unit_manifold;
                Assert(std::abs(weights[0] + weights[1] - 1.0) < 1e-13,
                       ExcMessage("Weights do not sum up to 1"));
                const Point<3> intermediate =
                  unit_manifold.get_intermediate_point(Point<3>(directions[0]),
                                                       Point<3>(directions[1]),
                                                       weights[1]);
                return intermediate;
              }

            Tensor<1, 3> vPerp;
            Tensor<2, 2> Hessian;
            Tensor<1, 2> gradient;
            Tensor<1, 2> gradlocal;

            // On success we exit the loop early.
            // Otherwise, we just take the result after max_iterations steps.
            for (unsigned int i = 0; i < max_iterations; ++i)
              {
                // Step 2a: Find new descent direction

                // Get local basis for the estimate candidate
                const Tensor<1, 3> Clocalx = internal::compute_normal(candidate);
                const Tensor<1, 3> Clocaly = cross_product_3d(candidate, Clocalx);

                // For each vertices vector, compute the tangent vector from
                // candidate towards the vertices vector -- its length is the
                // spherical length from candidate to the vertices vector. Then
                // compute its contribution to the Hessian.
                gradient = 0.;
                Hessian  = 0.;
                for (unsigned int i = 0; i < n_merged_points; ++i)
                  if (std::abs(weights[i]) > 1.e-15)
                    {
                      vPerp =
                        internal::projected_direction(directions[i], candidate);
                      const double sinthetaSq = vPerp.norm_square();
                      const double sintheta   = std::sqrt(sinthetaSq);
                      if (sintheta < tolerance)
                        {
                          Hessian[0][0] += weights[i];
                          Hessian[1][1] += weights[i];
                        }
                      else
                        {
                          const double costheta = (directions[i]) * candidate;
                          const double theta    = std::atan2(sintheta, costheta);
                          const double sincthetaInv = theta / sintheta;

                          const double cosphi = vPerp * Clocalx;
                          const double sinphi = vPerp * Clocaly;

                          gradlocal[0] = cosphi;
                          gradlocal[1] = sinphi;
                          gradient += (weights[i] * sincthetaInv) * gradlocal;

                          const double wt       = weights[i] / sinthetaSq;
                          const double sinphiSq = sinphi * sinphi;
                          const double cosphiSq = cosphi * cosphi;
                          const double tt       = sincthetaInv * costheta;
                          const double offdiag =
                            cosphi * sinphi * wt * (1.0 - tt);
                          Hessian[0][0] += wt * (cosphiSq + tt * sinphiSq);
                          Hessian[0][1] += offdiag;
                          Hessian[1][0] += offdiag;
                          Hessian[1][1] += wt * (sinphiSq + tt * cosphiSq);
                        }
                    }

                Assert(determinant(Hessian) > tolerance, ExcInternalError());

                const Tensor<2, 2> inverse_Hessian = invert(Hessian);

                const Tensor<1, 2> xDisplocal = inverse_Hessian * gradient;
                const Tensor<1, 3> xDisp =
                  xDisplocal[0] * Clocalx + xDisplocal[1] * Clocaly;

                // Step 2b: rotate candidate in direction xDisp for a new
                // candidate.
                const Point<3> candidateOld = candidate;
                candidate =
                  Point<3>(internal::apply_exponential_map(candidate, xDisp));

                // Step 2c: return the new candidate if we didn't move
                if ((candidate - candidateOld).norm_square() <
                    tolerance * tolerance)
                  break;
              }
          }
          return candidate;
        }
      } // namespace
    }   // namespace SphericalManifold
  } // namespace internal



  template <int dim, int spacedim>
  void
  SphericalManifold<dim, spacedim>::do_get_new_points(
    const ArrayView<const Point<spacedim>> &surrounding_points,
    const ArrayView<const double>          &weights,
    ArrayView<Point<spacedim>>              new_points) const
  {
    AssertDimension(weights.size(),
                    new_points.size() * surrounding_points.size());
    const unsigned int weight_rows    = new_points.size();
    const unsigned int weight_columns = surrounding_points.size();

    if (surrounding_points.size() == 2)
      {
        for (unsigned int row = 0; row < weight_rows; ++row)
          new_points[row] =
            SphericalManifold<dim, spacedim>::get_intermediate_point(surrounding_points[0],
                                                                     surrounding_points[1],
                                                                     weights[row * weight_columns + 1]);
        return;
      }

    small_vector<std::pair<double, Tensor<1, spacedim>>>
    new_candidates(new_points.size());
    small_vector<Tensor<1, spacedim>> directions(
      surrounding_points.size(), Point<spacedim>());
    small_vector<double> distances(
      surrounding_points.size(), 0.0);
    double max_distance = 0.;
    for (unsigned int i = 0; i < surrounding_points.size(); ++i)
      {
        directions[i] = surrounding_points[i] - center;
        distances[i]  = directions[i].norm();

        if (distances[i] != 0.)
          directions[i] /= distances[i];
        else
          Assert(false,
                 ExcMessage("One of the vertices coincides with the center. "
                            "This is not allowed!"));

        // Check if an estimate is good enough,
        // this is often the case for sufficiently refined meshes.
        for (unsigned int k = 0; k < i; ++k)
          {
            const double squared_distance =
              (directions[i] - directions[k]).norm_square();
            max_distance = std::max(max_distance, squared_distance);
          }
      }

    // Step 1: Check for some special cases, create simple linear guesses
    // otherwise.
    const double                              tolerance = 1e-10;
    small_vector<bool> accurate_point_was_found(
      new_points.size(), false);
    const ArrayView<const Tensor<1, spacedim>> array_directions =
      make_array_view(directions.begin(), directions.end());
    const ArrayView<const double> array_distances =
      make_array_view(distances.begin(), distances.end());
    for (unsigned int row = 0; row < weight_rows; ++row)
      {
        new_candidates[row] =
          guess_new_point(array_directions,
                          array_distances,
                          ArrayView<const double>(&weights[row * weight_columns],
                                                  weight_columns));

        // If the candidate is the center, mark it as found to avoid entering
        // the Newton iteration in step 2, which would crash.
        if (new_candidates[row].first == 0.0)
          {
            new_points[row]               = center;
            accurate_point_was_found[row] = true;
            continue;
          }

        // If not in 3d, just use the implementation from PolarManifold
        // after we verified that the candidate is not the center.
        if (spacedim < 3)
          new_points[row] = polar_manifold.get_new_point(
                              surrounding_points,
                              ArrayView<const double>(&weights[row * weight_columns],
                                                      weight_columns));
      }

    // In this case, we treated the case that the candidate is the center and
    // obtained the new locations from the PolarManifold object otherwise.
    if (spacedim < 3)
      return;
    else
      {
        // If all the points are close to each other, we expect the estimate to
        // be good enough. This tolerance was chosen such that the first iteration
        // for a at least three time refined HyperShell mesh with radii .5 and 1.
        // doesn't already succeed.
        if (max_distance < 2e-2)
          {
            for (unsigned int row = 0; row < weight_rows; ++row)
              new_points[row] =
                center + new_candidates[row].first * new_candidates[row].second;

            return;
          }

        // Step 2:
        // Do more expensive Newton-style iterations to improve the estimate.

        // Search for duplicate directions and merge them to minimize the cost of
        // the get_new_point function call below.
        small_vector<double, 1000> merged_weights(
          weights.size());
        small_vector<Tensor<1, spacedim>>
        merged_directions(surrounding_points.size(), Point<spacedim>());
        small_vector<double> merged_distances(
          surrounding_points.size(), 0.0);

        unsigned int n_unique_directions = 0;
        for (unsigned int i = 0; i < surrounding_points.size(); ++i)
          {
            bool found_duplicate = false;

            // This inner loop is of $O(N^2)$ complexity, but
            // surrounding_points.size() is usually at most 8 points large.
            for (unsigned int j = 0; j < n_unique_directions; ++j)
              {
                const double squared_distance =
                  (directions[i] - directions[j]).norm_square();
                if (!found_duplicate && squared_distance < 1e-28)
                  {
                    found_duplicate = true;
                    for (unsigned int row = 0; row < weight_rows; ++row)
                      merged_weights[row * weight_columns + j] +=
                        weights[row * weight_columns + i];
                  }
              }

            if (found_duplicate == false)
              {
                merged_directions[n_unique_directions] = directions[i];
                merged_distances[n_unique_directions]  = distances[i];
                for (unsigned int row = 0; row < weight_rows; ++row)
                  merged_weights[row * weight_columns + n_unique_directions] =
                    weights[row * weight_columns + i];

                ++n_unique_directions;
              }
          }

        // Search for duplicate weight rows and merge them to minimize the cost of
        // the get_new_point function call below.
        small_vector<unsigned int> merged_weights_index(
          new_points.size(), numbers::invalid_unsigned_int);
        for (unsigned int row = 0; row < weight_rows; ++row)
          {
            for (unsigned int existing_row = 0; existing_row < row;
                 ++existing_row)
              {
                bool identical_weights = true;

                for (unsigned int weight_index = 0;
                     weight_index < n_unique_directions;
                     ++weight_index)
                  if (std::abs(
                        merged_weights[row * weight_columns + weight_index] -
                        merged_weights[existing_row * weight_columns +
                                       weight_index]) > tolerance)
                    {
                      identical_weights = false;
                      break;
                    }

                if (identical_weights)
                  {
                    merged_weights_index[row] = existing_row;
                    break;
                  }
              }
          }

        // Note that we only use the n_unique_directions first entries in the
        // ArrayView
        const ArrayView<const Tensor<1, spacedim>> array_merged_directions =
          make_array_view(merged_directions.begin(),
                          merged_directions.begin() + n_unique_directions);
        const ArrayView<const double> array_merged_distances =
          make_array_view(merged_distances.begin(),
                          merged_distances.begin() + n_unique_directions);

        for (unsigned int row = 0; row < weight_rows; ++row)
          if (!accurate_point_was_found[row])
            {
              if (merged_weights_index[row] == numbers::invalid_unsigned_int)
                {
                  const ArrayView<const double> array_merged_weights(
                    &merged_weights[row * weight_columns], n_unique_directions);
                  new_candidates[row].second =
                    internal::SphericalManifold::do_get_new_point(
                      array_merged_directions,
                      array_merged_distances,
                      array_merged_weights,
                      Point<spacedim>(new_candidates[row].second));
                }
              else
                new_candidates[row].second =
                  new_candidates[merged_weights_index[row]].second;

              new_points[row] =
                center + new_candidates[row].first * new_candidates[row].second;
            }
      }
  }



  template <int dim, int spacedim>
  std::pair<double, Tensor<1, spacedim>>
  SphericalManifold<dim, spacedim>::guess_new_point(
    const ArrayView<const Tensor<1, spacedim>> &directions,
    const ArrayView<const double>              &distances,
    const ArrayView<const double>              &weights) const
  {
    const double        tolerance = 1e-10;
    double              rho       = 0.;
    Tensor<1, spacedim> candidate;

    // Perform a simple average ...
    double total_weights = 0.;
    for (unsigned int i = 0; i < directions.size(); ++i)
      {
        // if one weight is one, return its direction
        if (std::abs(1 - weights[i]) < tolerance)
          return std::make_pair(distances[i], directions[i]);

        rho += distances[i] * weights[i];
        candidate += directions[i] * weights[i];
        total_weights += weights[i];
      }

    // ... and normalize if the candidate is different from the origin.
    const double norm = candidate.norm();
    if (norm == 0.)
      return std::make_pair(0.0, Point<spacedim>());
    candidate /= norm;
    rho /= total_weights;

    return std::make_pair(rho, candidate);
  }


  template class SphericalManifold<2>;
  template class SphericalManifold<3>;
}

#endif
