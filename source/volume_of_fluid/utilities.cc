/*
 Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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

#include <aspect/volume_of_fluid/utilities.h>

namespace aspect
{
  namespace VolumeOfFluid
  {
    namespace Utilities
    {
      using namespace dealii;


      double compute_fluid_fraction (const Tensor<1, 2> normal,
                                     const double d)
      {
        const int dim = 2;

        //Get 1-Norm
        double norm1 = 0.0;
        double max = 0.0;
        for (unsigned int i = 0; i < dim; ++i)
          {
            const double normal_component = std::abs (normal[i]);
            norm1 += normal_component;
            max = (max < normal_component) ? normal_component : max;
          }

        //Obtain volume
        if (d <= -0.5*norm1)
          {
            return 0.0;
          }
        if (d >= 0.5*norm1)
          {
            return 1.0;
          }
        const double dtest = d / norm1; // Threshold value for changes in computation behavior
        // Normalized parameter indicating the "slope" (positive and finite)
        // Chosen due to resulting in simple formulas
        // Equal to the absolute value of the smaller vector entry in the 1-norm normalized normal vector
        const double mpos = 1.0 - max/norm1;
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



      double compute_interface_location (const Tensor<1, 2> normal,
                                         const double vol)
      {
        const int dim = 2;

        //Get 1-Norm
        double norm1 = 0.0;
        double max = 0.0;
        for (unsigned int i = 0; i < dim; ++i)
          {
            double normal_component = std::abs (normal[i]);
            norm1 += normal_component;
            max = (max < normal_component) ? normal_component : max;
          }

        // Normalized parameter indicating the "slope" (positive and finite)
        // Chosen due to resulting in simple formulas
        // Equal to the absolute value of the smaller vector entry in the 1-norm normalized normal vector
        const double mpos = (norm1==0.0)?0.0:(1.0-max/norm1);
        norm1 = (norm1==0.0)? 1.0:norm1;

        // Obtain correct RHS term for interface position
        if (vol <= 0.0)
          {
            return -0.5 * norm1;
          }
        else if (vol >= 1.0)
          {
            return 0.5 * norm1;
          }
        else if (vol < 0.5 * mpos / (1 - mpos))
          {
            return norm1 * (-0.5 + sqrt (2.0*vol * mpos * (1 - mpos)));
          }
        else if (vol > 1.0 - 0.5 * mpos / (1 - mpos))
          {
            return norm1 * (0.5 - sqrt (2.0*(1.0 - vol) * mpos * (1 - mpos)));
          }
        else
          {
            return norm1 * (1 - mpos) * (vol - 0.5);
          }
      }



      double compute_fluid_fraction (const Tensor<1, 3> normal,
                                     const double d)
      {
        // Calculations done by Scardovelli and Zaleski in
        // doi:10.1006/jcph.2000.6567,
        // modified to fit chosen convention on normal interface
        const int dim = 3;

        Tensor<1, 3> nnormal;
        //Simplify calculation by reducing to case of d<=0
        if (d>0.0)
          {
            return 1.0-compute_fluid_fraction(-normal, -d);
          }

        //Get 1-Norm
        double norm1 = 0.0;
        for (unsigned int i = 0; i < dim; ++i)
          {
            double term = std::abs(normal[i]);
            norm1 += term;
          }

        //Return volume in simple cases
        if (d <= -0.5*norm1)
          {
            return 0.0;
          }
        if (norm1 == 0.0)
          {
            // Case should never occur, return reasonable output for edge circumstance
            return 0.5;
          }
        const double dtest = d / norm1;

        // sort normalized values for normal into nnormal in ascending order of absolute value
        double mprod = 1.0;
        for (unsigned int i = 0; i < dim; ++i)
          {
            nnormal[i]=std::abs(normal[i])/norm1;
            mprod *= nnormal[i];
          }
        for (unsigned int i =0; i < dim; ++i)
          {
            for (unsigned int j=i+1; j < dim; ++j)
              {
                if (nnormal[j]<nnormal[i])
                  {
                    //Swap
                    const double tmp = nnormal[j];
                    nnormal[j] = nnormal[i];
                    nnormal[i] = tmp;
                  }
              }
          }
        const double m12 = nnormal[0]+nnormal[1];
        const double mmin = (m12<nnormal[2])?m12:nnormal[2];
        const double eps = 1e-10;
        const double v1= (6.0*nnormal[1]*nnormal[2]>eps)
                         ?
                         nnormal[0]*nnormal[0]/(eps)
                         :
                         nnormal[0]*nnormal[0]/(6*nnormal[1]*nnormal[2]);

        // do computation for standard cases
        // Case 1 of Scardovelli and Zaleski (Tetrahedron)
        if (dtest<nnormal[0]-0.5)
          {
            return (dtest+0.5)*(dtest+0.5)*(dtest+0.5)/(6*mprod);
          }
        // Case 2 of Scardovelli and Zaleski ("Triangular prism-like")
        else if (dtest<nnormal[1]-0.5)
          {
            return (dtest+0.5)*(dtest+0.5-nnormal[0])/(2*nnormal[1]*nnormal[2])+v1;
          }
        // Case 3 of Scardovelli and Zaleski
        else if (dtest<mmin-0.5)
          {
            return (dtest+0.5)*(dtest+0.5)*(3*m12-dtest-0.5)/(6*mprod) +
                   nnormal[0]*nnormal[0]*(nnormal[0]-3*dtest-1.5)/(6*mprod) +
                   nnormal[1]*nnormal[1]*(nnormal[1]-3*dtest-1.5)/(6*mprod);
          }
        // Case 4 of Scardovelli and Zaleski
        else if ( nnormal[2]<m12)
          {
            return (dtest+0.5)*(dtest+0.5)*(3-2*dtest-1.0)/(6*mprod) +
                   nnormal[0]*nnormal[0]*(nnormal[0]-3*dtest-1.5)/(6*mprod) +
                   nnormal[1]*nnormal[1]*(nnormal[1]-3*dtest-1.5)/(6*mprod) +
                   nnormal[2]*nnormal[2]*(nnormal[2]-3*dtest-1.5)/(6*mprod);
          }

        // Case 5 of Scardovelli and Zaleski
        return 0.5*(2.0*dtest+1.0-m12)/nnormal[2];
      }



      double compute_interface_location (const Tensor<1, 3, double> normal,
                                         const double vol)
      {
        // Calculations done by Scardovelli and Zaleski in
        // doi:10.1006/jcph.2000.6567,
        // modified to fit chosen convention on normal interface
        const int dim = 3;
        // Simplify to vol<0.5 case
        if (vol>0.5)
          {
            return -compute_interface_location(-normal, 1.0-vol);
          }

        //Get 1-Norm
        double norm1 = 0.0;
        for (unsigned int i = 0; i < dim; ++i)
          {
            double term = std::abs(normal[i]);
            norm1 += term;
          }
        const double eps = 1e-10;
        double mprod = 1.0;
        Tensor<1, 3> nnormal;
        if (norm1<eps)
          {
            nnormal[0] = 0.0;
            nnormal[1] = 0.0;
            nnormal[2] = 1.0;
            norm1 = 1.0;
            mprod =0.0;
          }
        else
          {
            // sort normalized values for normal into nnormal in ascending order of absolute value
            for (unsigned int i = 0; i < dim; ++i)
              {
                nnormal[i]=std::abs(normal[i])/norm1;
                mprod *= nnormal[i];
              }
            for (unsigned int i =0; i < dim; ++i)
              {
                for (unsigned int j=i+1; j < dim; ++j)
                  {
                    if (nnormal[j]<nnormal[i])
                      {
                        //Swap
                        const double tmp = nnormal[j];
                        nnormal[j] = nnormal[i];
                        nnormal[i] = tmp;
                      }
                  }
              }
          }
        //Simple cases
        if (vol<=0.0)
          {
            return -0.5*norm1;
          }

        const double m1 = nnormal[0];
        const double m2 = nnormal[1];
        const double m3 = nnormal[2];
        const double m12 = m1+m2;

        // Case 1 of Scardovelli and Zaleski
        const double v1=(6.0*m2*m3>eps)
                        ?
                        m1*m1/(eps)
                        :
                        m1*m1/(6*m2*m3);

        if (vol<v1)
          {
            return -0.5+std::pow(6*mprod*vol, 1./3.);
          }

        // Case 2 of Scardovelli and Zaleski
        const double v2 = v1 + 0.5*(m2-m1)/m3;
        if (vol<v2)
          {
            return 0.5*(-1+m1+sqrt(m1*m1+8*m2*m3*(vol-v1)));
          }

        // Case 3 of Scardovelli and Zaleski
        double v3 = (m12<m3)
                    ?
                    m3*m3*(3*m12-m3)/(6.0*mprod)
                    + m1*m1*(m1-3*m3)/(6*mprod)
                    + m2*m2*(m2-3*m3)/(6*mprod)
                    :
                    0.5*m12/m3;

        if (vol<v3)
          {
            // Solve appropriate cubic
            const double a2 = -3.0*m12;
            const double a1 = 3.0*(m1*m1+m2*m2);
            const double a0 = 6*mprod*vol-m1*m1*m1-m2*m2*m2;
            const double np0 = a2*a2/9.0-a1/3.0;
            const double q0 = (a1*a2-3.0*a0)/6.0-a2*a2*a2/27.0;
            const double theta = acos(q0/sqrt(np0*np0*np0))/3.0;
            return sqrt(np0)*(sqrt(3.0)*sin(theta)-cos(theta))-a2/3.0;
          }

        // Case 4
        if (m3<m12)
          {
            // Solve appropriate cubic
            double a2 = -1.5;
            double a1 = 1.5*(m1*m1+m2*m2+m3*m3);
            double a0 = 6*mprod*vol-m1*m1*m1-m2*m2*m2-m3*m3*m3;
            double np0 = a2*a2/9.0-a1/3.0;
            double q0 = (a1*a2-3.0*a0)/6.0-a2*a2*a2/27.0;
            double theta = acos(q0/sqrt(np0*np0*np0))/3.0;
            return sqrt(np0)*(sqrt(3.0)*sin(theta)-cos(theta))-a2/3.0;
          }

        return -0.5+m1*vol+0.5*m12;
      }



      void xFEM_Heaviside(const unsigned int degree,
                          const Tensor<1, 2> normal,
                          const double d,
                          const std::vector<Point<2>> &points,
                          std::vector<double> &values)
      {
        const int basis_count=4;
        std::vector<double> coeffs(basis_count);

        const double n_xp = fabs(normal[0]), n_yp = fabs(normal[1]);
        const double sign_n_x = (((normal[0]) > 0) - ((normal[0]) < 0)),
                     sign_n_y = (((normal[1]) > 0) - ((normal[1]) < 0));

        const double norm1 = n_xp + n_yp;
        const double triangle_break = 0.5*fabs(n_xp-n_yp);

        const unsigned int max_degree = 1;

        AssertThrow(degree<=max_degree,
                    ExcMessage("Cannot generate xFEM polynomials. Only implemented for degrees<2."));

        // The formulas below calculate the correct coefficients for a given
        // basis in order to form a polynomial $f$ which will satisfy $\int
        // fpdx=\int pH(d-n\cdot x)dx$ for all polynomials $p$ less than or
        // equal to the given degree
        //
        // The functions for the correct values were calculated and exported using sympy
        if (d<-0.5*norm1)
          {
            for (unsigned int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d>0.5*norm1)
          {
            // Full cell
            coeffs[0] = 1.0;
            for (unsigned int i =1; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (norm1< 1e-7)
          {
            coeffs[0] = 0.5;
            for (unsigned int i =1; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d<=-triangle_break)
          {
            //Triangle
            const double d_n = d + 0.5*norm1;
            coeffs[0]=0.5*d_n*d_n/(n_xp*n_yp); // 1
            coeffs[1]=d_n*d_n*(d_n - 1.5*n_yp)/(n_xp*n_yp*n_yp)*sign_n_y; // 2*y - 1
            coeffs[2]=d_n*d_n*(d_n - 1.5*n_xp)/(n_xp*n_xp*n_yp)*sign_n_x; // 2*x - 1
            coeffs[3]=1.5*d_n*d_n*(d_n*d_n - 2.0*d_n*n_xp - 2.0*d_n*n_yp + 3.0*n_xp*n_yp)/
                      (n_xp*n_xp*n_yp*n_yp)*sign_n_x*sign_n_y; // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_xp<n_yp)
          {
            //Trapezoid X
            coeffs[0]=(d + 0.5*n_yp)/n_yp; // 1
            coeffs[1]=0.25*(12.0*d*d + n_xp*n_xp - 3.0*n_yp*n_yp)/(n_yp*n_yp)*sign_n_y; // 2*y - 1
            coeffs[2]=-0.5*n_xp/n_yp*sign_n_x; // 2*x - 1
            coeffs[3]=-3.0*d*n_xp/(n_yp*n_yp)*sign_n_x*sign_n_y; // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_yp<n_xp)
          {
            //Trapezoid Y
            coeffs[0]=(d + 0.5*n_xp)/n_xp; // 1
            coeffs[1]=-0.5*n_yp/n_xp*sign_n_y; // 2*y - 1
            coeffs[2]=0.25*(12.0*pow(d, 2) - 3.0*n_xp*n_xp + n_yp*n_yp)/(n_xp*n_xp)*sign_n_x; // 2*x - 1
            coeffs[3]=-3.0*d*n_yp/(n_xp*n_xp)*sign_n_x*sign_n_y; // (2*x - 1)*(2*y - 1)
          }
        else
          {
            //ITriangle
            const double d_nn = 0.5*norm1-d;
            coeffs[0]=1.0L-0.5L*d_nn*d_nn/(n_xp*n_yp); // 1
            coeffs[1]=0.5L*(d_nn*d_nn)*sign_n_y*(2*d_nn - 3*n_yp)/(n_xp*(n_yp*n_yp)); // 2*y - 1
            coeffs[2]=0.5L*(d_nn*d_nn)*sign_n_x*(2*d_nn - 3*n_xp)/((n_xp*n_xp)*n_yp); // 2*x - 1
            coeffs[3]=1.5L*(d_nn*d_nn)*sign_n_x*sign_n_y*(-(d_nn*d_nn) + 2*d_nn*n_xp + 2*d_nn*n_yp - 3*n_xp*n_yp)/((n_xp*n_xp)*(n_yp*n_yp)); // (2*x - 1)*(2*y - 1)
          }

        // Calculate the correct values at the provided quadrature points by
        // multiplying coefficients by the basis polynomials.
        for (unsigned int i = 0; i<points.size(); ++i)
          {
            const Point<2> point = points[i];
            const double x = point[0], y = point[1];
            values[i] = coeffs[0];
            if (degree>=1)
              {
                values[i] += coeffs[1]*(2.0*y-1.0) +
                             coeffs[2]*(2.0*x-1.0) +
                             coeffs[3]*(2.0*x - 1)*(2.0*y - 1.0);
              }
          }
      }



      void xFEM_Heaviside(const unsigned int /*degree*/,
                          const Tensor<1, 3> /*normal*/,
                          const double /*d*/,
                          const std::vector<Point<3>> &/*points*/,
                          std::vector<double> &/*values*/)
      {
        AssertThrow(false, ExcNotImplemented());
      }



      void xFEM_Heaviside_derivative_d(const unsigned int degree,
                                       const Tensor<1, 2> normal,
                                       const double d,
                                       const std::vector<Point<2>> &points,
                                       std::vector<double> &values)
      {
        const int basis_count=4;
        std::vector<double> coeffs(basis_count);

        const double n_xp = fabs(normal[0]), n_yp = fabs(normal[1]);
        const double sign_n_x = (((normal[0]) > 0) - ((normal[0]) < 0)),
                     sign_n_y = (((normal[1]) > 0) - ((normal[1]) < 0));

        const double norm1 = n_xp + n_yp;
        const double triangle_break = 0.5L*fabs(n_xp-n_yp);

        const int max_degree = 1;

        AssertThrow(degree<=max_degree,
                    ExcMessage("Cannot generate xFEM polynomials are only functional for degrees<2."));


        // The formulas below calculate the correct coefficients for a given
        // basis in order to form a polynomial $f$ which will satisfy $\int
        // fpdx=\int pH(d-n\cdot x)dx$ for all polynomials $p$ less than or
        // equal to the given degree
        //
        // The functions for the correct values were calculated and exported using sympy
        if (d<-0.5*norm1)
          {
            for (int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d>0.5*norm1)
          {
            // Full cell
            for (int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (norm1<1e-7)
          {
            for (int i =0; i < basis_count; ++i)
              coeffs[i] = 0.0;
          }
        else if (d<=-triangle_break)
          {
            //D Triangle
            const double d_n = d + 0.5*norm1;
            coeffs[0]=d_n/(n_xp*n_yp); // 1
            coeffs[1]=3*d_n*sign_n_y*(d_n - n_yp)/(n_xp*(n_yp*n_yp)); // 2*y - 1
            coeffs[2]=3*d_n*sign_n_x*(d_n - n_xp)/((n_xp*n_xp)*n_yp); // 2*x - 1
            coeffs[3]=3*d_n*sign_n_x*sign_n_y*(2*(d_n*d_n) - 3*d_n*n_xp - 3*d_n*n_yp + 3*n_xp*n_yp)/((n_xp*n_xp)*(n_yp*n_yp)); // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_xp<n_yp)
          {
            //D Trapezoid X
            coeffs[0]=1.0/n_yp; // 1
            coeffs[1]=6*d*sign_n_y/(n_yp*n_yp); // 2*y - 1
            coeffs[2]=0; // 2*x - 1
            coeffs[3]=-3*n_xp*sign_n_x*sign_n_y/(n_yp*n_yp); // (2*x - 1)*(2*y - 1)
          }
        else if (d<triangle_break && n_yp<n_xp)
          {
            //D Trapezoid Y
            coeffs[0]=1.0/n_xp; // 1
            coeffs[1]=0; // 2*y - 1
            coeffs[2]=6*d*sign_n_x/(n_xp*n_xp); // 2*x - 1
            coeffs[3]=-3*n_yp*sign_n_x*sign_n_y/(n_xp*n_xp); // (2*x - 1)*(2*y - 1)
          }
        else
          {
            //D ITriangle
            const double d_nn = 0.5*norm1-d;
            coeffs[0]=d_nn/(n_xp*n_yp); // 1
            coeffs[1]=3*d_nn*sign_n_y*(-d_nn + n_yp)/(n_xp*(n_yp*n_yp)); // 2*y - 1
            coeffs[2]=3*d_nn*sign_n_x*(-d_nn + n_xp)/((n_xp*n_xp)*n_yp); // 2*x - 1
            coeffs[3]=3*d_nn*sign_n_x*sign_n_y*(2*(d_nn*d_nn) - 3*d_nn*n_xp - 3*d_nn*n_yp + 3*n_xp*n_yp)/((n_xp*n_xp)*(n_yp*n_yp)); // (2*x - 1)*(2*y - 1)
          }


        // Calculate the correct values at the provided quadrature points by
        // multiplying coefficients by the basis polynomials.
        for (unsigned int i = 0; i<points.size(); ++i)
          {
            const Point<2> point = points[i];
            const double x = point[0], y = point[1];
            values[i] = coeffs[0];
            if (degree>=1)
              {
                values[i] += coeffs[1]*(2*y-1.0) +
                             coeffs[2]*(2*x-1.0) +
                             coeffs[3]*(2*x - 1)*(2*y - 1);
              }
          }
      }



      void xFEM_Heaviside_derivative_d(const unsigned int /*degree*/,
                                       const Tensor<1, 3> /*normal*/,
                                       const double /*d*/,
                                       const std::vector<Point<3>> &/*points*/,
                                       std::vector<double> &/*values*/)
      {
        AssertThrow(false, ExcNotImplemented());
      }



      template<int dim>
      double compute_interface_location_newton(const unsigned int degree,
                                               const Tensor<1, dim> normal,
                                               const double volume_fraction,
                                               const double vol,
                                               const double epsilon,
                                               const std::vector<Point<dim>> &points,
                                               const std::vector<double> &weights)
      {
        double norm1=0.0;
        for (int i=0; i<dim; ++i)
          norm1+=fabs(normal[i]);
        double d_l=-0.5L*norm1, d_h=0.5L*norm1;
        double f_l=0.0, f_h=1.0;
        double d_guess= d_l + (volume_fraction-f_l)*(d_h-d_l)/(f_h-f_l);

        std::vector<double> f_values(points.size());
        std::vector<double> df_values(points.size());

        for (int iter=0; iter<40; ++iter)
          {
            xFEM_Heaviside(degree, normal, d_guess, points, f_values);
            xFEM_Heaviside_derivative_d(degree, normal, d_guess, points, df_values);

            double f_guess=0.0;
            double df_guess=0.0;
            for (unsigned int i=0; i<points.size(); ++i)
              {
                const double factor = weights[i]/vol;
                f_guess  += f_values[i]*factor;
                df_guess += df_values[i]*factor;
              }

            // Break if within tolerance
            if (fabs(f_guess-volume_fraction)<epsilon)
              {
                break;
              }

            if (volume_fraction<f_guess)
              {
                d_h = d_guess;
                f_h = f_guess;
              }
            else
              {
                d_l = d_guess;
                f_l = f_guess;
              }

            if (fabs(df_guess)<epsilon)
              {
                d_guess = (volume_fraction-f_l)/(f_h-f_l)*(d_h-d_l);
              }
            else
              {
                d_guess += (volume_fraction-f_guess)/(df_guess);

                if (d_guess < d_l || d_guess > d_h)
                  {
                    d_guess = d_l + (volume_fraction-f_l)*(d_h-d_l)/(f_h-f_l);
                  }
              }
          }

        return d_guess;
      }

      template<int dim>
      double compute_fluid_volume(const unsigned int degree,
                                  const Tensor<1, dim> normal,
                                  const double d,
                                  const std::vector<Point<dim>> &points,
                                  const std::vector<double> &weights)
      {
        std::vector<double> f_values(points.size());

        xFEM_Heaviside(degree, normal, d, points, f_values);

        double fluid_volume=0.0;
        for (unsigned int i=0; i<points.size(); ++i)
          fluid_volume += f_values[i]*weights[i];

        return fluid_volume;
      }

      template<int dim>
      double calculate_volume_flux(const unsigned int dir,
                                   const double time_direction_derivative,
                                   const Tensor<1, dim> normal,
                                   const double d_face)
      {
        Tensor<1, dim> i_normal;
        double i_d;

        // Get d value at center of "aperture" (cell face cross timestep)
        i_d = d_face+0.5*time_direction_derivative;
        // Get normal vector on "aperture" by replacing the appropriate component with the time_direction_derivative
        i_normal = normal;
        i_normal[dir] = time_direction_derivative;

        return compute_fluid_fraction (i_normal, i_d);
      }
    }
  }
}

namespace aspect
{
  namespace VolumeOfFluid
  {
    namespace Utilities
    {
#define INSTANTIATE(dim) \
  template double calculate_volume_flux<dim>(const unsigned int dir, \
                                             const double time_direction_derivative, \
                                             const Tensor<1, dim> normal, \
                                             const double d); \
  template double compute_interface_location_newton<dim>(const unsigned int degree, \
                                                         const Tensor<1, dim> normal, \
                                                         const double volume_fraction, \
                                                         const double vol, \
                                                         const double epsilon, \
                                                         const std::vector<Point<dim>> &points, \
                                                         const std::vector<double> &weights); \
  template double compute_fluid_volume<dim>(const unsigned int degree, \
                                            const Tensor<1, dim> normal,\
                                            const double d,\
                                            const std::vector<Point<dim>> &points,\
                                            const std::vector<double> &weights);

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
