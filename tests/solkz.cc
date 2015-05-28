#include <aspect/material_model/simple.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  /**
   * This is the "Sol Kz" benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {T. Duretz and D. A. May and T. V. Gerya and P. J. Tackley},
   *    title =        {Discretization errors and free surface stabilization in the
   *                  finite difference and marker-in-cell method for applied
   *                  geodynamics: {A} numerical study},
   *    journal =      {Geochemistry Geophysics Geosystems},
   *    year =         2011,
   *    volume =       12,
   *    pages =        {Q07004/1--26}}
   * @endcode
   *
   * The results are published in Kronbichler, Heister and Bangerth paper.
   */
  namespace InclusionBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      // based on http://geodynamics.org/hg/cs/AMR/Discontinuous_Stokes with permission
      // The following code has been taken from http://www.underworldproject.org/,
      // release 1.7.0. As mentioned in the Underworld Manual, this code has been
      // released under the GNU General Public License (GPL).


      void _Velic_solKz(
        double pos[],
        double _sigma, /* density */
        double _km, int _n, /* wavelength in z, wavenumber in x */
        double _B, /* viscosity parameter */
        double vel[], double *presssure,
        double total_stress[], double strain_rate[] )
      {
        double Z;
        double u1,u2,u3,u4,u5,u6,SS;
        double sum1,sum2,sum3,sum4,sum5,sum6,mag,sum7,x,z;
        double sigma;
        int n;
        double kn;
        double _C1,_C2,_C3,_C4;
        double B, Rp, UU, VV;
        double rho,a,b,r,_aa,_bb,AA,BB,Rm,km;


        double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
        double t11,t12,t13,t14,t15,t16,t17,t18,t19,t20;
        double t21,t22,t23,t24,t25,t26,t27,t28,t29,t31;
        double t33,t34,t35,t37,t38,t40,t41,t42,t43,t45;
        double t47,t51,t52,t53,t54,t55,t56,t57,t58,t59;
        double t60,t61,t62,t64,t65,t66,t67,t68,t69,t70;
        double t71,t72,t73,t74,t75,t76,t77,t78,t79,t80;
        double t81,t82,t83,t84,t85,t86,t89,t90,t92,t94;
        double t96,t97,t98,t99,t100,t101,t103,t104,t105,t106;
        double t107,t108,t109,t110,t111,t112,t113,t114,t115,t116;
        double t117,t118,t119,t120,t121,t122,t123,t124,t125,t126;
        double t127,t130,t131,t132,t134,t135,t141,t144,t147,t148;
        double t150,t151,t152,t161,t171;



        /*************************************************************************/
        /*************************************************************************/
        /* rho = -sigma*sin(km*z)*cos(kn*x) */
        /* viscosity  Z= exp(2*B*z)  */
        B = _B; /* viscosity parameter must be non-zero*/
        km = _km; /* solution valid for km not zero -- should get trivial solution if km=0 */
        n = _n; /* solution valid for n not zero */
        sigma = _sigma;
        /*************************************************************************/
        /*************************************************************************/
        kn = (double) _n*M_PI;
        a = B*B + kn*kn;
        b = 2.0*kn*B;
        r = sqrt(a*a + b*b);
        Rp = sqrt( (r+a)/2.0 );
        Rm  = sqrt( (r-a)/2.0 );
        UU  = Rp - B;
        VV = Rp + B;


        x = pos[0];
        z = pos[1];

        sum1=0.0;
        sum2=0.0;
        sum3=0.0;
        sum4=0.0;
        sum5=0.0;
        sum6=0.0;
        sum7=0.0;




        /*******************************************/
        /*         calculate the constants         */
        /*******************************************/

        t3 = kn * kn;
        t4 = km * km;
        t6 = B * B;
        t8 = 0.4e1 * t3 * t6;
        t10 = 0.4e1 * t4 * t6;
        t13 = 0.8e1 * kn * t6 * km;
        t14 = t4 * t4;
        t16 = 0.2e1 * t3 * t4;
        t17 = t3 * t3;
        _aa = -0.4e1 * B * km * kn * (t3 + t4) / (t8 + t10 + t13 + t14 + t16 + t17) / (-t13 + t8 + t10 + t14 + t16 + t17);

        t1 = kn * kn;
        t2 = t1 * t1;
        t3 = B * B;
        t5 = 0.4e1 * t1 * t3;
        t6 = km * km;
        t7 = t6 * t6;
        t9 = 0.2e1 * t1 * t6;
        t11 = 0.4e1 * t3 * t6;
        t16 = 0.8e1 * kn * t3 * km;
        _bb = kn * (t2 + t5 + t7 + t9 - t11) / (t5 + t11 + t16 + t7 + t9 + t2) / (-t16 + t5 + t11 + t7 + t9 + t2);

        AA = _aa;
        BB = _bb;

        t1 = B * B;
        t2 = t1 * Rp;
        t4 = Rm * Rm;
        t5 = t4 * Rp;
        t7 = t4 * B;
        t8 = km * km;
        t12 = Rp * Rp;
        t13 = B * t12;
        t21 = 0.8e1 * t1 * km * BB * Rp;
        t23 = 0.2e1 * Rm;
        t24 = cos(t23);
        t26 = Rm * Rp;
        t38 = sin(t23);
        t51 = exp(-0.2e1 * Rp);
        t53 = B + Rp;
        t54 = Rm * t53;
        t55 = Rm * B;
        t57 = 0.2e1 * B * km;
        t58 = t55 + t57 - t26;
        t62 = 0.3e1 * t1;
        t64 = 0.2e1 * Rp * B;
        t65 = t62 + t64 + t4 - t8 - t12;
        t67 = t54 * t65 * BB;
        t69 = Rm - km;
        t70 = cos(t69);
        t72 = -t57 + t55 - t26;
        t77 = Rm + km;
        t78 = cos(t77);
        t81 = t54 * t65 * AA;
        t86 = sin(t77);
        t92 = sin(t69);
        t96 = exp(-t53);
        t98 = B - Rp;
        t99 = Rm * t98;
        t100 = t55 + t57 + t26;
        t104 = t62 - t64 + t4 - t8 - t12;
        t106 = t99 * t104 * BB;
        t109 = -t57 + t55 + t26;
        t116 = t99 * t104 * AA;
        t130 = exp(-0.3e1 * Rp - B);
        t135 = exp(-0.4e1 * Rp);
        t144 = t4 * t1;
        t150 = t4 * t12;
        _C1 = (((0.2e1 * Rp * (0.2e1 * t2 + 0.2e1 * t5 + t7 + B * t8 - 0.3e1 * t1 * B + t13) * AA + t21) * t24 + (-0.2e1 * t26 * (t4 - t8 - t12 + 0.5e1 * t1) * AA + 0.8e1 * B * BB * km * Rm * Rp) * t38 - 0.2e1 * B * (0.2e1 * t13 + t12 * Rp - 0.3e1 * t2 + t5 + 0.2e1 * t7 + t8 * Rp) * AA - t21) * t51 + ((0.2e1 * t54 * t58 * AA + t67) * t70 + (0.2e1 * t54 * t72 * AA - t67) * t78 + (t81 + 0.2e1 * t54 * t72 * BB) * t86 + (t81 - 0.2e1 * t54 * t58 * BB) * t92) * t96 + ((-0.2e1 * t99 * t100 * AA - t106) * t70 + (-0.2e1 * t99 * t109 * AA + t106) * t78 + (-t116 - 0.2e1 * t99 * t109 * BB) * t86 + (-t116 + 0.2e1 * t99 * t100 * BB) * t92) * t130 + 0.4e1 * t4 * t98 * t53 * AA * t135) / (((-0.8e1 * t4 - 0.8e1 * t1) * t12 * t24 + 0.8e1 * t144 + 0.8e1 * t12 * t1) * t51 + (0.4e1 * t150 - 0.4e1 * t144) * t135 + 0.4e1 * t150 - 0.4e1 * t144);

        t1 = Rm * Rp;
        t2 = Rm * Rm;
        t3 = km * km;
        t4 = Rp * Rp;
        t5 = B * B;
        t12 = km * Rm;
        t17 = 0.2e1 * Rm;
        t18 = cos(t17);
        t22 = t2 * Rp;
        t25 = B * t3;
        t26 = t5 * B;
        t33 = t5 * km;
        t38 = sin(t17);
        t40 = Rm * B;
        t41 = 0.3e1 * t5;
        t51 = exp(-0.2e1 * Rp);
        t53 = B + Rp;
        t54 = Rm * t53;
        t57 = t41 + 0.2e1 * Rp * B + t2 - t3 - t4;
        t59 = t54 * t57 * AA;
        t60 = B * km;
        t61 = 0.2e1 * t60;
        t62 = t40 + t61 - t1;
        t67 = Rm - km;
        t68 = cos(t67);
        t70 = -t61 + t40 - t1;
        t75 = Rm + km;
        t76 = cos(t75);
        t82 = t54 * t57 * BB;
        t84 = sin(t75);
        t90 = sin(t67);
        t94 = exp(-t53);
        t97 = 0.3e1 * Rm * t26;
        t98 = t2 * Rm;
        t99 = t98 * B;
        t100 = t3 * Rm;
        t101 = t100 * Rp;
        t103 = Rm * t4 * B;
        t104 = t4 * Rp;
        t105 = Rm * t104;
        t107 = 0.8e1 * t33 * Rp;
        t109 = 0.5e1 * t1 * t5;
        t110 = t98 * Rp;
        t111 = t100 * B;
        t112 = t97 + t99 - t101 + t103 - t105 + t107 + t109 + t110 - t111;
        t114 = t2 * t4;
        t116 = 0.2e1 * t60 * t1;
        t117 = t2 * t5;
        t119 = 0.3e1 * t26 * Rp;
        t120 = t104 * B;
        t121 = t4 * t5;
        t122 = 0.2e1 * t121;
        t123 = t22 * B;
        t125 = 0.2e1 * t33 * Rm;
        t126 = t25 * Rp;
        t127 = t114 + t116 + t117 - t119 + t120 + t122 + t123 + t125 + t126;
        t132 = -t107 + t103 - t105 - t101 + t97 - t111 + t110 + t109 + t99;
        t134 = t120 - t125 + t123 - t116 + t122 + t117 + t114 + t126 - t119;
        t152 = exp(-0.3e1 * Rp - B);
        t161 = exp(-0.4e1 * Rp);
        _C2 = (((0.2e1 * t1 * (t2 - t3 - t4 + 0.5e1 * t5) * AA - 0.8e1 * B * BB * t12 * Rp) * t18 + (0.2e1 * Rp * (0.2e1 * t5 * Rp + 0.2e1 * t22 + t2 * B + t25 - 0.3e1 * t26 + B * t4) * AA + 0.8e1 * t33 * BB * Rp) * t38 + 0.2e1 * t40 * (t41 + t4 + t2 - t3) * AA - 0.8e1 * t5 * BB * t12) * t51 + ((-t59 + 0.2e1 * t54 * t62 * BB) * t68 + (-t59 - 0.2e1 * t54 * t70 * BB) * t76 + (0.2e1 * t54 * t70 * AA - t82) * t84 + (0.2e1 * t54 * t62 * AA + t82) * t90) * t94 + ((t112 * AA - 0.2e1 * t127 * BB) * t68 + (t132 * AA + 0.2e1 * t134 * BB) * t76 + (-0.2e1 * t134 * AA + t132 * BB) * t84 + (-0.2e1 * t127 * AA - t112 * BB) * t90) * t152 + (-0.2e1 * t59 + 0.8e1 * t40 * km * t53 * BB) * t161) / (((-0.8e1 * t2 - 0.8e1 * t5) * t4 * t18 + 0.8e1 * t117 + 0.8e1 * t121) * t51 + (0.4e1 * t114 - 0.4e1 * t117) * t161 + 0.4e1 * t114 - 0.4e1 * t117);

        t1 = B * B;
        t2 = t1 * Rp;
        t4 = Rm * Rm;
        t5 = t4 * Rp;
        t7 = Rp * Rp;
        t8 = B * t7;
        t11 = km * km;
        t13 = t4 * B;
        t21 = 0.8e1 * t1 * km * BB * Rp;
        t23 = 0.2e1 * Rm;
        t24 = cos(t23);
        t26 = Rm * Rp;
        t38 = sin(t23);
        t51 = exp(-0.2e1 * Rp);
        t53 = B + Rp;
        t54 = Rm * t53;
        t55 = Rm * B;
        t57 = 0.2e1 * B * km;
        t58 = t55 + t57 - t26;
        t62 = 0.3e1 * t1;
        t64 = 0.2e1 * Rp * B;
        t65 = t62 + t64 + t4 - t11 - t7;
        t67 = t54 * t65 * BB;
        t69 = Rm - km;
        t70 = cos(t69);
        t72 = -t57 + t55 - t26;
        t77 = Rm + km;
        t78 = cos(t77);
        t81 = t54 * t65 * AA;
        t86 = sin(t77);
        t92 = sin(t69);
        t96 = exp(-t53);
        t98 = B - Rp;
        t99 = Rm * t98;
        t100 = t55 + t57 + t26;
        t104 = t62 - t64 + t4 - t11 - t7;
        t106 = t99 * t104 * BB;
        t109 = -t57 + t55 + t26;
        t116 = t99 * t104 * AA;
        t130 = exp(-0.3e1 * Rp - B);
        t141 = t4 * t1;
        t147 = t4 * t7;
        t151 = exp(-0.4e1 * Rp);
        _C3 = (((-0.2e1 * Rp * (-0.2e1 * t2 - 0.2e1 * t5 + t8 - 0.3e1 * t1 * B + B * t11 + t13) * AA - t21) * t24 + (0.2e1 * t26 * (t4 - t11 - t7 + 0.5e1 * t1) * AA - 0.8e1 * B * BB * km * Rm * Rp) * t38 - 0.2e1 * B * (0.2e1 * t8 + 0.2e1 * t13 + 0.3e1 * t2 - t7 * Rp - t5 - t11 * Rp) * AA + t21) * t51 + ((-0.2e1 * t54 * t58 * AA - t67) * t70 + (-0.2e1 * t54 * t72 * AA + t67) * t78 + (-t81 - 0.2e1 * t54 * t72 * BB) * t86 + (-t81 + 0.2e1 * t54 * t58 * BB) * t92) * t96 + ((0.2e1 * t99 * t100 * AA + t106) * t70 + (0.2e1 * t99 * t109 * AA - t106) * t78 + (t116 + 0.2e1 * t99 * t109 * BB) * t86 + (t116 - 0.2e1 * t99 * t100 * BB) * t92) * t130 + 0.4e1 * t4 * t98 * t53 * AA) / (((-0.8e1 * t4 - 0.8e1 * t1) * t7 * t24 + 0.8e1 * t141 + 0.8e1 * t7 * t1) * t51 + (0.4e1 * t147 - 0.4e1 * t141) * t151 + 0.4e1 * t147 - 0.4e1 * t141);

        t1 = Rm * Rp;
        t2 = Rm * Rm;
        t3 = km * km;
        t4 = Rp * Rp;
        t5 = B * B;
        t12 = km * Rm;
        t17 = 0.2e1 * Rm;
        t18 = cos(t17);
        t22 = t2 * Rp;
        t25 = t5 * B;
        t27 = B * t3;
        t33 = t5 * km;
        t38 = sin(t17);
        t40 = Rm * B;
        t41 = 0.3e1 * t5;
        t51 = exp(-0.2e1 * Rp);
        t53 = t2 * Rm;
        t54 = t53 * B;
        t56 = 0.5e1 * t1 * t5;
        t58 = Rm * t4 * B;
        t59 = t3 * Rm;
        t60 = t59 * Rp;
        t62 = 0.8e1 * t33 * Rp;
        t64 = 0.3e1 * Rm * t25;
        t65 = t53 * Rp;
        t66 = t59 * B;
        t67 = t4 * Rp;
        t68 = Rm * t67;
        t69 = t54 - t56 + t58 + t60 - t62 + t64 - t65 - t66 + t68;
        t71 = t2 * t4;
        t73 = 0.3e1 * t25 * Rp;
        t74 = t2 * t5;
        t75 = t27 * Rp;
        t76 = B * km;
        t78 = 0.2e1 * t76 * t1;
        t80 = 0.2e1 * t33 * Rm;
        t81 = t22 * B;
        t82 = t4 * t5;
        t83 = 0.2e1 * t82;
        t84 = t67 * B;
        t85 = t71 + t73 + t74 - t75 - t78 + t80 - t81 + t83 - t84;
        t89 = Rm - km;
        t90 = cos(t89);
        t92 = t60 - t66 - t65 + t58 + t54 - t56 + t62 + t68 + t64;
        t94 = t73 + t78 - t81 + t74 - t80 - t84 - t75 + t83 + t71;
        t98 = Rm + km;
        t99 = cos(t98);
        t105 = sin(t98);
        t111 = sin(t89);
        t115 = exp(-Rp - B);
        t117 = B - Rp;
        t118 = Rm * t117;
        t121 = t41 - 0.2e1 * Rp * B + t2 - t3 - t4;
        t123 = t118 * t121 * AA;
        t124 = 0.2e1 * t76;
        t125 = t40 + t124 + t1;
        t131 = -t124 + t40 + t1;
        t141 = t118 * t121 * BB;
        t152 = exp(-0.3e1 * Rp - B);
        t171 = exp(-0.4e1 * Rp);
        _C4 = (((-0.2e1 * t1 * (t2 - t3 - t4 + 0.5e1 * t5) * AA + 0.8e1 * B * BB * t12 * Rp) * t18 + (-0.2e1 * Rp * (-0.2e1 * t5 * Rp - 0.2e1 * t22 + t4 * B - 0.3e1 * t25 + t27 + t2 * B) * AA - 0.8e1 * t33 * BB * Rp) * t38 + 0.2e1 * t40 * (t41 + t4 + t2 - t3) * AA - 0.8e1 * t5 * BB * t12) * t51 + ((t69 * AA - 0.2e1 * t85 * BB) * t90 + (t92 * AA + 0.2e1 * t94 * BB) * t99 + (-0.2e1 * t94 * AA + t92 * BB) * t105 + (-0.2e1 * t85 * AA - t69 * BB) * t111) * t115 + ((-t123 + 0.2e1 * t118 * t125 * BB) * t90 + (-t123 - 0.2e1 * t118 * t131 * BB) * t99 + (0.2e1 * t118 * t131 * AA - t141) * t105 + (0.2e1 * t118 * t125 * AA + t141) * t111) * t152 - 0.2e1 * t123 + 0.8e1 * t40 * km * t117 * BB) / (((-0.8e1 * t2 - 0.8e1 * t5) * t4 * t18 + 0.8e1 * t74 + 0.8e1 * t82) * t51 + (0.4e1 * t71 - 0.4e1 * t74) * t171 + 0.4e1 * t71 - 0.4e1 * t74);

        /******************************************************************/
        /******************************************************************/

        /*******************************************/
        /*       calculate the velocities etc      */
        /*******************************************/

        t2 = exp(UU * z);
        t3 = Rm * z;
        t4 = cos(t3);
        t6 = sin(t3);
        t11 = exp(-VV * z);
        t18 = exp(-0.2e1 * z * B);
        t19 = km * z;
        t20 = cos(t19);
        t22 = sin(t19);
        u1 = kn * (t2 * (_C1 * t4 + _C2 * t6) + t11 * (_C3 * t4 + _C4 * t6) + t18 * (AA * t20 + BB * t22));

        t1 = Rm * z;
        t2 = cos(t1);
        t4 = sin(t1);
        t14 = exp(UU * z);
        t26 = exp(-VV * z);
        t28 = km * z;
        t29 = cos(t28);
        t31 = sin(t28);
        t43 = exp(-0.2e1 * z * B);
        u2 = (-UU * (_C1 * t2 + _C2 * t4) + _C1 * t4 * Rm - _C2 * t2 * Rm) * t14 + (VV * (_C3 * t2 + _C4 * t4) + _C3 * t4 * Rm - _C4 * t2 * Rm) * t26 + (0.2e1 * B * (AA * t29 + BB * t31) + AA * t31 * km - BB * t29 * km) * t43;

        t2 = 0.2e1 * z * B;
        t3 = exp(t2);
        t4 = t3 * kn;
        t5 = Rm * z;
        t6 = cos(t5);
        t8 = sin(t5);
        t18 = exp(UU * z);
        t31 = exp(-VV * z);
        t34 = km * z;
        t35 = cos(t34);
        t37 = sin(t34);
        t47 = exp(-t2);
        u3 = 0.2e1 * t4 * (UU * (_C1 * t6 + _C2 * t8) - _C1 * t8 * Rm + _C2 * t6 * Rm) * t18 + 0.2e1 * t4 * (-VV * (_C3 * t6 + _C4 * t8) - _C3 * t8 * Rm + _C4 * t6 * Rm) * t31 + 0.2e1 * t4 * (-0.2e1 * B * (AA * t35 + BB * t37) - AA * t37 * km + BB * t35 * km) * t47;

        t1 = Rm * Rm;
        t3 = UU * UU;
        t8 = kn * kn;
        t11 = Rm * z;
        t12 = sin(t11);
        t14 = cos(t11);
        t20 = t14 * Rm;
        t27 = 0.2e1 * z * B;
        t28 = exp(t27);
        t31 = exp(UU * z);
        t38 = VV * VV;
        t54 = exp(-VV * z);
        t56 = km * km;
        t59 = B * B;
        t66 = km * z;
        t67 = sin(t66);
        t69 = cos(t66);
        t83 = exp(-t27);
        u4 = ((_C2 * t1 - t3 * _C2 + 0.2e1 * UU * _C1 * Rm - _C2 * t8) * t12 + _C1 * t14 * t1 - t3 * _C1 * t14 - 0.2e1 * UU * _C2 * t20 - t8 * _C1 * t14) * t28 * t31 + ((-0.2e1 * VV * _C3 * Rm + _C4 * t1 - _C4 * t8 - t38 * _C4) * t12 + 0.2e1 * VV * _C4 * t20 + _C3 * t14 * t1 - t8 * _C3 * t14 - t38 * _C3 * t14) * t28 * t54 + ((BB * t56 - t8 * BB - 0.4e1 * t59 * BB - 0.4e1 * B * AA * km) * t67 + AA * t69 * t56 - t8 * AA * t69 - 0.4e1 * t59 * AA * t69 + 0.4e1 * B * BB * t69 * km) * t28 * t83;


        t1 = Rm * z;
        t2 = sin(t1);
        t3 = Rm * Rm;
        t4 = t3 * Rm;
        t5 = t2 * t4;
        t6 = UU * UU;
        t7 = t6 * UU;
        t8 = cos(t1);
        t15 = 0.2e1 * B * t8 * t3;
        t19 = B * UU;
        t20 = t2 * Rm;
        t23 = kn * kn;
        t24 = B * t23;
        t26 = 0.2e1 * t24 * t8;
        t27 = t23 * UU;
        t29 = B * t6;
        t33 = t23 * t2 * Rm;
        t35 = 0.1e1 / kn;
        t42 = 0.2e1 * B * t2 * t3;
        t43 = t8 * t4;
        t45 = 0.2e1 * t24 * t2;
        t52 = t23 * t8 * Rm;
        t53 = t8 * Rm;
        t64 = 0.2e1 * z * B;
        t65 = exp(t64);
        t68 = exp(UU * z);
        t70 = B * VV;
        t76 = t23 * VV;
        t78 = VV * VV;
        t79 = t78 * VV;
        t84 = B * t78;
        t108 = exp(-VV * z);
        t111 = km * z;
        t112 = sin(t111);
        t113 = km * km;
        t118 = cos(t111);
        t119 = t118 * km;
        t121 = B * B;
        t123 = t112 * km;
        t130 = t113 * km;
        t148 = exp(-t64);
        u5 = (-(-t5 - t7 * t8 + 0.3e1 * UU * t8 * t3 + t15 + 0.3e1 * t6 * t2 * Rm + 0.4e1 * t19 * t20 - t26 + t27 * t8 - 0.2e1 * t29 * t8 - t33) * t35 * _C1 - (-t7 * t2 + t27 * t2 + t42 + t43 - t45 + 0.3e1 * UU * t2 * t3 - 0.2e1 * t29 * t2 + t52 - 0.4e1 * t19 * t53 - 0.3e1 * t6 * t8 * Rm) * t35 * _C2) * t65 * t68 + (-(t15 - 0.4e1 * t70 * t20 - t33 - 0.3e1 * VV * t8 * t3 - t76 * t8 + t79 * t8 + 0.3e1 * t78 * t2 * Rm - 0.2e1 * t84 * t8 - t26 - t5) * t35 * _C3 - (t52 - 0.3e1 * VV * t2 * t3 + t79 * t2 + 0.4e1 * t70 * t53 - 0.3e1 * t78 * t8 * Rm - 0.2e1 * t84 * t2 + t43 - t76 * t2 + t42 - t45) * t35 * _C4) * t65 * t108 - t65 * (-0.4e1 * B * BB * t112 * t113 + t23 * BB * t119 + 0.4e1 * t121 * AA * t123 - 0.4e1 * t121 * BB * t119 + BB * t118 * t130 - AA * t112 * t130 - 0.4e1 * B * AA * t118 * t113 - t23 * AA * t123 - 0.4e1 * t24 * AA * t118 - 0.4e1 * t24 * BB * t112) * t35 * t148;



        t2 = 0.2e1 * z * B;
        t3 = exp(t2);
        t4 = t3 * kn;
        t5 = Rm * z;
        t6 = cos(t5);
        t8 = sin(t5);
        t18 = exp(UU * z);
        t31 = exp(-VV * z);
        t34 = km * z;
        t35 = cos(t34);
        t37 = sin(t34);
        t47 = exp(-t2);
        u6 = -0.2e1 * t4 * (UU * (_C1 * t6 + _C2 * t8) - _C1 * t8 * Rm + _C2 * t6 * Rm) * t18 - 0.2e1 * t4 * (-VV * (_C3 * t6 + _C4 * t8) - _C3 * t8 * Rm + _C4 * t6 * Rm) * t31 - 0.2e1 * t4 * (-0.2e1 * B * (AA * t35 + BB * t37) - AA * t37 * km + BB * t35 * km) * t47;




        /******************************************************************/
        /******************************************************************/



        sum5 += u5*cos(n*M_PI*x);  /* pressure */
        u6 -= u5; /* get total stress */
        sum6 += u6*cos(n*M_PI*x);  /* xx stress */

        u1 *= cos(n*M_PI*x); /* z velocity */
        sum1 += u1;
        u2 *= sin(n*M_PI*x); /* x velocity */
        sum2 += u2;
        u3 -= u5; /* get total stress */
        u3 *= cos(n*M_PI*x); /* zz stress */
        sum3 += u3;
        u4 *= sin(n*M_PI*x); /* zx stress */
        sum4 += u4;

        rho = -sigma*sin(km*z)*cos(kn*x); /* density */
        sum7 += rho;

        SS = exp(UU*z)*(_C1*cos(Rm*z)+_C2*sin(Rm*z)) +exp(-VV*z)*(_C3*cos(Rm*z)+_C4*sin(Rm*z)) + exp(-2*z*B)*(AA*cos(km*z)+BB*sin(km*z));
        SS *= sin(kn*x); /* stream function */

        mag=sqrt(u1*u1+u2*u2);
        /*printf("%0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f\n",x,z,sum1,sum2,sum3,sum4,sum5,sum6,mag,sum7,SS);*/


        /* Output */
        if ( vel != NULL )
          {
            vel[0] = sum2;
            vel[1] = sum1;
          }
        if ( presssure != NULL )
          {
            (*presssure) = sum5;
          }
        if ( total_stress != NULL )
          {
            total_stress[0] = sum6;
            total_stress[1] = sum3;
            total_stress[2] = sum4;
          }
        if ( strain_rate != NULL )
          {
            /* sigma = tau - p, tau = sigma + p, tau[] = 2*eta*strain_rate[] */
            Z = exp( 2.0 * B * z );
            strain_rate[0] = (sum6+sum5)/(2.0*Z);
            strain_rate[1] = (sum3+sum5)/(2.0*Z);
            strain_rate[2] = (sum4)/(2.0*Z);
          }
        /* Value checks, could be cleaned up if needed. Julian Giordani 9-Oct-2006*/
        if ( fabs( sum5 - ( -0.5*(sum6+sum3) ) ) > 1e-5 )
          {
            assert(0);
          }
      }



      /**
       * The exact solution for the SolKz benchmark.
       */
      template <int dim>
      class FunctionSolKz : public Function<dim>
      {
        public:
          FunctionSolKz () : Function<dim>() {}

          virtual void vector_value (const Point< dim >   &p,
                                     Vector< double >   &values) const
          {
            double pos[2]= {p(0),p(1)};
            double total_stress[3], strain_rate[3];
            static const double B = 0.5 * std::log(1e6);
            AnalyticSolutions::_Velic_solKz
            (pos,
             1.0, 2, 3,
             B,
             &values[0], &values[2], total_stress, strain_rate );
          }
      };


    }



    template <int dim>
    class SolKzMaterial : public MaterialModel::InterfaceCompatibility<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return true if the viscosity() function returns something that
         * may depend on the variable identifies by the argument.
         */
        virtual bool
        viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the density() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        density_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the compressibility() function returns something
         * that may depend on the variable identifies by the argument.
         *
         * This function must return false for all possible arguments if the
         * is_compressible() function returns false.
         */
        virtual bool
        compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the specific_heat() function returns something
         * that may depend on the variable identifies by the argument.
         */
        virtual bool
        specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the thermal_conductivity() function returns
         * something that may depend on the variable identifies by the
         * argument.
         */
        virtual bool
        thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the contuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

//TODO: should we make this a virtual function as well? where is it used?
        double reference_thermal_diffusivity () const;

        double reference_cp () const;
        /**
         * @}
         */
    };



    template <int dim>
    double
    SolKzMaterial<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      // defined as given in the Duretz et al. paper
      static const double B = 0.5 * std::log(1e6);
      return std::exp(2*B*p[1]);
    }


    template <int dim>
    double
    SolKzMaterial<dim>::
    reference_viscosity () const
    {
      return 1;
    }

    template <int dim>
    double
    SolKzMaterial<dim>::
    reference_density () const
    {
      return 0;
    }

    template <int dim>
    double
    SolKzMaterial<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 0;
    }

    template <int dim>
    double
    SolKzMaterial<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    SolKzMaterial<dim>::
    reference_cp () const
    {
      return 0;
    }

    template <int dim>
    double
    SolKzMaterial<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    SolKzMaterial<dim>::
    reference_thermal_diffusivity () const
    {
      return 0;
    }

    template <int dim>
    double
    SolKzMaterial<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &p) const
    {
      // defined as given in the paper
      return -std::sin(2*p[1])*std::cos(3*numbers::PI*p[0]);
    }


    template <int dim>
    double
    SolKzMaterial<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    SolKzMaterial<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }



    template <int dim>
    bool
    SolKzMaterial<dim>::
    viscosity_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    SolKzMaterial<dim>::
    density_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SolKzMaterial<dim>::
    compressibility_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SolKzMaterial<dim>::
    specific_heat_depends_on (const MaterialModel::NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    SolKzMaterial<dim>::
    thermal_conductivity_depends_on (const MaterialModel::NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    SolKzMaterial<dim>::
    is_compressible () const
    {
      return false;
    }






    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Duretz et al. reference above.
      */
    template <int dim>
    class SolKzPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };

    template <int dim>
    std::pair<std::string,std::string>
    SolKzPostprocessor<dim>::execute (TableHandler &statistics)
    {
      AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                  ExcNotImplemented());

      std_cxx1x::shared_ptr<Function<dim> > ref_func;
      if (dynamic_cast<const SolKzMaterial<dim> *>(&this->get_material_model()) != NULL)
        {
          const SolKzMaterial<dim> *
          material_model
            = dynamic_cast<const SolKzMaterial<dim> *>(&this->get_material_model());

          ref_func.reset (new AnalyticSolutions::FunctionSolKz<dim>());
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Postprocessor DuretzEtAl only works with the material model SolCx, SolKz, and Inclusion."));
        }

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                          dim+2);
      ComponentSelectFunction<dim> comp_p(dim, dim+2);

      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_u,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_p);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_ul2,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_pl2,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p);

      std::ostringstream os;
      os << std::scientific << cellwise_errors_u.l1_norm()
         << ", " << cellwise_errors_p.l1_norm()
         << ", " << cellwise_errors_ul2.l2_norm()
         << ", " << cellwise_errors_pl2.l2_norm();

      return std::make_pair("Errors u_L1, p_L1, u_L2, p_L2:", os.str());
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolKzMaterial,
                                   "SolKzMaterial",
                                   "A material model that corresponds to the 'SolKz' benchmark "
                                   "defined in Duretz et al., G-Cubed, 2011.")


    ASPECT_REGISTER_POSTPROCESSOR(SolKzPostprocessor,
                                  "SolKzPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Duretz et al., G-Cubed, 2011, paper with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")
  }
}
