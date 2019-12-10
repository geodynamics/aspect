/*
 * This module contains the various subroutines
Usage import """
#######################################################################################
*/
//#####################  IMPORT STANDARD MODULES   ######################################   
/*
import sys,os
import numpy as np #for numerical analysis
from scipy import integrate  #integration routine
from subprocess import call #for calling unix commands
from datetime import date  #to give a timestamp to output
import pdb	#for the debugger pdb.set_trace()
from matplotlib.pyplot import * # for the figures
from math import pi,exp,log,sqrt
from ConfigParser import SafeConfigParser
*/
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <string>
#include <map>
#include <random>
#include <assert.h>

#include <gsl/gsl_integration.h>

#include "refvalues_and_utilities.h"

using namespace std;



/*
 * ------------------------------------------------------------------
 * 
 * 		CONSTANTS
 * 
 * ------------------------------------------------------------------
 */


/*
 * ***********************************************************
 * 
 * A series of constants relevant for a specific modeling approach.
 * 
 * The modeling approach name is contained in the name of namespace.
 * 
 * The namespace allows for fetching the correct constant value
 * for the specified Model and Parameter.
 * 
 * !!!!!
 * 			The constants are simple copy and pasted here
 * 			from the constants_complete.ini file.
 * 			These values should not be hardcoded, and
 * 			should be changed later to read from file.
 * !!!!!
 * 
 ************************************************************
 */
namespace Rheology_Constants
{
	namespace Common
	{
		// Ideal gas constant [J K^-1 mol^-1]
		const double R = 8.3144621;
	}
	//###########################  Jackson & Faul 2010  ###########################
	//###########################   Extended Burgers    ###########################
	namespace JF10_eBurg
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 66.5e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 1.8;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 360000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1e-5;
		
		// Reference Temperature [K]
		const double TR = 1173;
		
		// Reference Pressure [Pa]
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 13.4e-6;
		
		// Anelastic grain size exponent
		const double ma = 1.31;
		
		// Frequency dependence
		const double alpha = 0.274;
		
		// Relaxation strength - Burgers
		const double Delta = 1.04;
		
		// Reference relaxation time - upper bound [s]
		const double tauHR = 1e7;
		
		// Reference relaxation time - lower bound [s]
		const double tauLR = 1e-3;
		
		// Reference relaxation time - Maxwell [s]
		const double tauMR = 3.02e+7;
		
		
		// "extended" dissipation peak parms
		// Peak width
		const double sigma = 4;
		
		// Relaxation strength - Peak
		const double DeltaP = 0.057;
		
		// Reference relaxation time - Peak [s]
		const double tauPR = 3.9811e-4;
		
		// Viscous grain size exponent
		const double mv = 3;
	}
	//###########################  Jackson & Faul 2010  ###########################
	//###########################        Andrade        ###########################
	namespace Andrade
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 62.2e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 1.8;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 303000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1e-5;
		
		// Reference Temperature [K]
		const double TR = 1173;
		
		// Reference Pressure [Pa]
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 3.1e-6;
		
		// grainsize dependence exponent ------------ ?????? who cares ??????
		const double m = 3;
		
		// frequency exponent
		const double n = 0.33;
		
		// Beta star = Beta / J_unrelaxed
		const double Bstar = 0.020;
		
		// Reference relaxation time - Maxwell [s]
		const double tauMR = 1.9953e+5;
	}
	//##############################################################################
	//###########################  McCarthy et al. 2011  ###########################
	namespace M11
	{
		// Shear modulus at TR, PR, Isaak, 1992 [Pa]
		const double GUR = 82e+9;
		
		// T-derivative of Shear modulus, Isaak, 1992 [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus, Isaak, 1992
		const double dGdP = 1.8;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 505000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1.2e-5;
		
		// Reference Temperature [K]
		const double TR = 1473;
		
		// Reference Pressure [Pa]
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 1e-3;
		
		// grainsize dependence exponent
		const double m = 3;
		
		// Reference viscosity [Pa s]
		const double eta0 = 6.6e+19;
		
		// fn polynomial fit parms (equation 26)
		const double a0 = 5.5097e-1;
		const double a1 = 5.4332e-2;
		const double a2 = -2.3615e-3;
		const double a3 = -5.7175e-5;
		const double a4 = 9.9473e-6;
		const double a5 = -3.4761e-7;
		const double a6 = 3.9461e-9;
	}
	//##########################################################################
	//###########################  Takei et al 2014  ###########################
	namespace Tak14
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 82e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.0136e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 1.8;
		
		// Reference Temperature [K]
		const double TR = 1473;
		
		// Reference Pressure (by ref to M11, not in paper??)  Pa
		const double PR = 0.2e+9;
		
		// Reference Grainsize [m]
		const double gsR = 1e-3;
		
		// Activation energy ("H" in the paper) [J mol^-1]
		const double E = 505000;
		
		// Activation volume [m^3 mol^-1]
		const double V = 1.2e-5;
		
		// Reference Viscosity [Pa s]
		const double eta0 = 6.6e+19;
		
		// grainsize dependence exponent
		const double m = 3;
		
		// Peak standard deviation (value for melt-free Ol)
		const double sigmap = 4;
		
		// Peak pre-exponent (value for melt-free Ol)
		const double Ap = 0.007;
	}
	//###################################################################################
	//###########################  Priestley & McKenzie 2013  ###########################
	namespace P_M13
	{
		// Shear modulus at TR, PR [Pa]
		const double GUR = 72.66e+9;
		
		// T-derivative of Shear modulus [Pa K^-1]
		const double dGdT = -0.00871e+9;
		
		// P-derivative of Shear modulus
		const double dGdP = 2.04;
		
		// Activation energy ("U" in the paper) [J mol^-1]
		const double E = 402900;
		
		// Activation volume - m^3 mol^-1]
		const double V = 7.81e-6;
		
		// Reference Temperature [K]
		const double TR = 1473;
		
		// Reference Pressure [Pa]
		const double PR = 1.5e+9;
		
		// Reference Viscosity [Pa s]
		const double eta0 = 2.3988e+22;
		
		// Reference Grainsize [m]
		const double gsR = 1e-3;
		
		// grainsize dependence exponent
		const double m = 3;
		
		// fn polynomial fit parms (equation 26)
		const double a0 = 5.5097e-1;
		const double a1 = 5.4332e-2;
		const double a2 = -2.3615e-3;
		const double a3 = -5.7175e-5;
		const double a4 = 9.9473e-6;
		const double a5 = -3.4761e-7;
		const double a6 = 3.9461e-9;
	}
	//##############################################################################
	//###########################  LOWER MANTLE GUESSES  ###########################
	namespace LOWM
	{
		// Activation volume - m^3 mol^-1
		const double V = 1e-6;
	}
}


/*
 * Forward declaration
 */
double getfn(double freq, double P, double T, double gs, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN());
double getJ1byJu(double freq, double P, double T, double gs, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN());
double getJ2byJu(double freq, double P, double T, double gs, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN());
double getJu(double T, double P, std::string model, double Vs0 = std::numeric_limits<double>::quiet_NaN(), double rho = std::numeric_limits<double>::quiet_NaN() );
double gettau(double tau0, double m, double P, double T, double gs, std::string model);
double intgrate_j1b(double limL, double limH, double alpha, double omega, int N=1000);
double mccarthyXn(double taun);
double intgrate_j2b(double limL, double limH, double alpha, double omega, int N=1000);
double gettauM_P_M13(double T, double P, double Ju);



/*
 * ------------------------------------------------------------------
 *
 *		Integrations
 *
 *-------------------------------------------------------------------
 */

struct intgrnd_jp_params { double omega; double tauP; double sigma; };

double intgrnd_j1p (double tau, void *parameters)
{
  struct intgrnd_jp_params *p = (struct intgrnd_jp_params *)parameters;
  return ( (1.0/tau) * (1.0 / (1.0 + std::pow(p->omega*tau, 2.0) ) ) * exp(-std::pow(log(tau/p->tauP), 2.0) / (2.0*std::pow(p->sigma, 2.0) ) ) );
}



double intgrnd_j2p (double tau, void *parameters)
{
  struct intgrnd_jp_params *p = (struct intgrnd_jp_params *)parameters;
  return ( (1.0/(1.0 + std::pow(p->omega*tau, 2.0) ) ) * exp(-std::pow(log(tau/p->tauP), 2.0) / (2.0*std::pow(p->sigma, 2.0) ) ) );
}



/*
 * ------------------------------------------------------------------
 * 
 *              Integrations for Takei 2014
 * 
 *-------------------------------------------------------------------
 */

struct intgrnd_jtak_params {double tau; double omega; double tauM; double Ap; double sigmap;};

double intgrnd_j1tak(double tau, void *parameters) //*
{
         struct intgrnd_jtak_params *p = (struct intgrnd_jtak_params *)parameters;
        // Integration part of J1byJu term (eqn 8 in Takei et al 2014)
        return ( (1.0 / tau)*(1.0 / (1.0 + std::pow(p->omega*tau, 2.0))) * (0.444*std::pow(tau/p->tauM, 0.38) + p->Ap*exp(-0.5*std::pow(((log(tau/p->tauM) + 8.1)/p->sigmap), 2.0))) );
}



double intgrnd_j2tak(double tau, void *parameters) //*
{
        struct intgrnd_jtak_params *p = (struct intgrnd_jtak_params *)parameters;
        // Integration part of J2byJu term (eqn 8 in Takei et al 2014)
        return ( (p->omega/(1.0 + std::pow(p->omega*tau, 2.0) ) ) * (0.444*std::pow(tau/p->tauM, 0.38) + p->Ap*exp(-0.5*std::pow(((log(tau/p->tauM) + 8.1) / p->sigmap), 2.0) ) ) );
}



double Vact_lowM(double P, double Vtop, double Vbot, double Ptop, double Pfold) //*
{
        return ( Vbot + (Vtop-Vbot) * exp( (Ptop-P) / Pfold ) );
}



std::pair<double, double> PTd_VsQ( double freq, double P, double T, double gs, double rho, std::string model, double highQapprox, double Vs0 = std::numeric_limits<double>::quiet_NaN() )
{
	/*
	 * This subroutine calculates the values of Vs and Q at a given
	 * PTd condition at a frequency f and density rho. If highQapprox is 1
	 * approximate values of Vs and Q can be calculated from
	 * e.g. McCarthy et al. (2011) eqn. 24, else use their more explicit form eqn. B6
	 */

	//---------[ NB this uses V_activation = sloping from 11 - 3 (e-6) in the LM ]---------


	double J1byJu, J2byJu, Ju, Vs, G;

	// Value of the real and complex part of compliance
	if ( std::isnan(Vs0) )
	{
		J1byJu = getJ1byJu(freq, P, T, gs, model);
		J2byJu = getJ2byJu(freq, P, T, gs, model);
		Ju = getJu(T, P, model);
	}
	else
	{
		J1byJu = getJ1byJu(freq, P, T, gs, model, Vs0, rho);
		J2byJu = getJ2byJu(freq, P, T, gs, model, Vs0, rho);
		Ju = getJu(T, P, model, Vs0, rho);
	}


	double J1 = J1byJu * Ju;
	double J2 = J2byJu * Ju;


	// EVERYONE AGREES
	double qs = J2 / J1;


	// VS DISAGREEMENT
	if (model == "M11") Vs = 1.0 / ( sqrt(rho * J1) );
	else if (model == "P_M13") Vs = 1.0/( sqrt(rho * J1) );
	else if (model == "JF10_eBurg")
	{
		G = std::pow((J1*J1 + J2*J2), -0.5);
		Vs = sqrt(G/rho);
	}
	else if (model == "Tak14")
	{
		// pdb.set_trace()
		G = std::pow(J1*J1 + J2*J2, -0.5);
		Vs = sqrt(G/rho);
	}


	// **** When Qs^-1 << 1 (J2/J1 << 1), highQapprox can be 1. If Qs^-1 is not small, Qs^-1 does not equal Qmu^-1 ****
	if (highQapprox == 0)
	{
		double factor = ( 1.0 + sqrt(1.0 + std::pow((J2/J1), 2.0) ) )/2.0;
		qs = qs/factor;
		Vs = Vs/sqrt(factor);
	}


	double Q = 1.0 / qs;


	// if (Q > 1e9) pdb.set_trace();

	return std::make_pair(Vs,Q);
}



double getJ1byJu(double freq, double P, double T, double gs, std::string model, double Vs0, double rho)
{
	// This subroutine gets the value of J1/Ju based on expressions in a paper
	
	
	double J1byJu, Ju, omega;
	double fn, val;
	double a0, a1, a2, a3, a4, a5, a6;
	
	
	// If model is McCarthy 2011
	if (model == "M11")
	{
		if ( std::isnan(Vs0) ) fn = getfn(freq, P, T, gs, model);
		else fn = getfn(freq, P, T, gs, model, Vs0, rho);
		
		
		// Getting J1 values from McCarthy et al., 2011
		a0 = Rheology_Constants::M11::a0;
		a1 = Rheology_Constants::M11::a1;
		a2 = Rheology_Constants::M11::a2;
		a3 = Rheology_Constants::M11::a3;
		a4 = Rheology_Constants::M11::a4;
		a5 = Rheology_Constants::M11::a5;
		a6 = Rheology_Constants::M11::a6;

		if (fn <= 1.0e+13)
		{
			val = 0.0;
			val = val + a0*(std::pow(log(fn), 0));
			val = val + a1*(std::pow(log(fn), 1));
			val = val + a2*(std::pow(log(fn), 2));
			val = val + a3*(std::pow(log(fn), 3));
			val = val + a4*(std::pow(log(fn), 4));
			val = val + a5*(std::pow(log(fn), 5));
			val = val + a6*(std::pow(log(fn), 6));
			J1byJu = 1.0 / val;
		}
		else J1byJu = 1.0;
	}
	else if (model == "P_M13")
	{
		fn = getfn(freq, P, T, gs, model, Vs0, rho);
		
		// Getting J1 values from McCarthy et al., 2011
		a0 = Rheology_Constants::P_M13::a0;
		a1 = Rheology_Constants::P_M13::a1;
		a2 = Rheology_Constants::P_M13::a2;
		a3 = Rheology_Constants::P_M13::a3;
		a4 = Rheology_Constants::P_M13::a4;
		a5 = Rheology_Constants::P_M13::a5;
		a6 = Rheology_Constants::P_M13::a6;
		
		
		if (fn <= 1.0e+13)
		{
			val = 0.0;
			val = val + a0*(std::pow(log(fn), 0));
			val = val + a1*(std::pow(log(fn), 1));
			val = val + a2*(std::pow(log(fn), 2));
			val = val + a3*(std::pow(log(fn), 3));
			val = val + a4*(std::pow(log(fn), 4));
			val = val + a5*(std::pow(log(fn), 5));
			val = val + a6*(std::pow(log(fn), 6));
			J1byJu = 1.0 / val;
		}
		else J1byJu = 1.0;
	}
	else if (model == "JF10_eBurg")
	{
		double alpha = Rheology_Constants::JF10_eBurg::alpha;
		double sigma = Rheology_Constants::JF10_eBurg::sigma;
		double Delta = Rheology_Constants::JF10_eBurg::Delta;
		double tauHR = Rheology_Constants::JF10_eBurg::tauHR;
		double tauLR = Rheology_Constants::JF10_eBurg::tauLR;
		double tauPR = Rheology_Constants::JF10_eBurg::tauPR;
		double DeltaP = Rheology_Constants::JF10_eBurg::DeltaP;
		double ma = Rheology_Constants::JF10_eBurg::ma;
		
		omega = 2.0 * UniversalConst::PI * freq;
		double tauH = gettau(tauHR, ma, P, T, gs, model);
		double tauL = gettau(tauLR, ma, P, T, gs, model);
		double tauP = gettau(tauPR, ma, P, T, gs, model);
		
		double j1b = (alpha*Delta)/(std::pow(tauH,alpha) - std::pow(tauL,alpha));
		double j1p = DeltaP/(sigma*sqrt(2*UniversalConst::PI));
		
		double i1b = intgrate_j1b(tauL,tauH,alpha,omega,1000);
                double i1p;

		gsl_integration_cquad_workspace *workspace = gsl_integration_cquad_workspace_alloc(100);

		double abserr;
		size_t nevals;

		struct intgrnd_jp_params params = {omega,tauP,sigma};
                gsl_function F;
		F.function = &intgrnd_j1p;
		F.params = &params;

		gsl_integration_cquad(&F,
		                      0,
		                      std::numeric_limits<double>::infinity(),
		                      1e-12,
		                      1e-12,
		                      workspace,
		                      &i1p,
		                      &abserr,
		                      &nevals);

		gsl_integration_cquad_workspace_free(workspace);
		
		J1byJu = 1.0 + (j1b*i1b) + (j1p*i1p);
	}
	else if (model == "Tak14")
	{
		double m = Rheology_Constants::Tak14::m; // grainsize exponent
		double eta0 = Rheology_Constants::Tak14::eta0; // ref viscosity
		double Ap = Rheology_Constants::Tak14::Ap; // pre-exponent for peak (normalisation factor)
		double sigmap = Rheology_Constants::Tak14::sigmap; // controls peak width - see Takei 2014 eqn 17
		
		// Value of the real and complex part of compliance			
		if ( std::isnan(Vs0) ) Ju = getJu(T, P, model);
		else Ju = getJu(T, P, model, Vs0, rho);
		
		omega = 2.0 * UniversalConst::PI * freq;
		double tauM = gettau(Ju*eta0, m, P, T, gs, model);
		double i1;

                gsl_integration_cquad_workspace *workspace = gsl_integration_cquad_workspace_alloc(100);

                double abserr;
                size_t nevals;

                struct intgrnd_jtak_params params = {omega, tauM, Ap, sigmap};
                gsl_function F;
                F.function = &intgrnd_j1tak;
                F.params = &params;

                gsl_integration_cquad(&F,
                                      0,
                                      std::numeric_limits<double>::infinity(),
                                      1e-12,
                                      1e-12,
                                      workspace,
                                      &i1,
                                      &abserr,
                                      &nevals);

                gsl_integration_cquad_workspace_free(workspace);
		
		J1byJu = 1.0 + i1;
	}
	
	
	return J1byJu;
}



double getJ2byJu(double freq, double P, double T, double gs, std::string model, double Vs0, double rho)
{
	// This subroutine gets the value of J2 based on expressions in a paper
	
	double fn, taun, Xn, J2byJu, omega, tauM;
	
	if (model == "M11")  // *** If model is McCarthy 2011
	{
		fn = getfn(freq, P, T, gs, model, Vs0, rho);
		taun = 1.0/(2.0 * UniversalConst::PI * fn);
		
		// Getting J2 values from McCarthy et al., 2011
		Xn = mccarthyXn(taun);
		J2byJu = (UniversalConst::PI / 2.0*Xn) + taun / (2.0 * UniversalConst::PI);
	}
	else if (model == "P_M13")  // *** If model is Priestley and McKenzie 2013
	{
		fn = getfn(freq, P, T, gs, model, Vs0, rho);
		taun = 1.0/(2.0 * UniversalConst::PI * fn);
		
		//#Getting J2 values from McCarthy et al., 2011
		Xn = mccarthyXn(taun);
		J2byJu = (UniversalConst::PI/2.0*Xn)+taun/(2.0*UniversalConst::PI);
	}
	else if (model == "JF10_eBurg")  // *** If model is extended Burgers from Jackson & Faul 2010
	{
		double alpha = Rheology_Constants::JF10_eBurg::alpha;
		double sigma = Rheology_Constants::JF10_eBurg::sigma;
		double Delta = Rheology_Constants::JF10_eBurg::Delta;
		double tauHR = Rheology_Constants::JF10_eBurg::tauHR;
		double tauLR = Rheology_Constants::JF10_eBurg::tauLR;
		double tauPR = Rheology_Constants::JF10_eBurg::tauPR;
		double tauMR = Rheology_Constants::JF10_eBurg::tauMR;
		double DeltaP = Rheology_Constants::JF10_eBurg::DeltaP;
		double ma = Rheology_Constants::JF10_eBurg::ma;
		double mv = Rheology_Constants::JF10_eBurg::mv;
		
		omega = 2.0 * UniversalConst::PI * freq;
		
		double tauH = gettau(tauHR, ma, P, T, gs, model);
		double tauL = gettau(tauLR, ma, P, T, gs, model);
		double tauP = gettau(tauPR, ma, P, T, gs, model);
		tauM = gettau(tauMR, mv, P, T, gs, model);
		
		double j2b = omega*(alpha*Delta)/(std::pow(tauH,alpha) - std::pow(tauL,alpha));
		double j2p = omega*DeltaP/(sigma*sqrt(2*UniversalConst::PI));
		
		std::vector<double> tau_arr = MyUtilities::linspace(tauL,tauH,100); //tau_arr=np.MyUtilities::linspace(tauL,tauH,100);
		double i2b = intgrate_j2b(tauL, tauH, alpha, omega, 1000);
                double i2p;

                gsl_integration_cquad_workspace *workspace = gsl_integration_cquad_workspace_alloc(100);

                double abserr;
                size_t nevals;

                struct intgrnd_jp_params params = {omega,tauP,sigma};
                gsl_function F;
                F.function = &intgrnd_j1p;
                F.params = &params;

                gsl_integration_cquad(&F,
                                      0,
                                      std::numeric_limits<double>::infinity(),
                                      1e-12,
                                      1e-12,
                                      workspace,
                                      &i2p,
                                      &abserr,
                                      &nevals);

                gsl_integration_cquad_workspace_free(workspace);
		
		J2byJu = (j2b*i2b)+(j2p*i2p)+(1.0/(omega*tauM));
	}
	else if (model == "Tak14")  // *** If model is Takei 2014
	{
		double m = Rheology_Constants::Tak14::m; // grainsize exponent
		double eta0 = Rheology_Constants::Tak14::eta0; // ref viscosity
		double Ap = Rheology_Constants::Tak14::Ap; // pre-exponent for peak (normalisation factor)
		double sigmap = Rheology_Constants::Tak14::sigmap; // controls peak width - see Takei 2014 eqn 17
		
		
		// Value of the real and complex part of compliance			
		double Ju;
		if ( std::isnan(Vs0) ) Ju = getJu(T, P, model);
		else Ju = 1.0 / (rho * std::pow(Vs0, 2.0));
		
		
		omega = 2.0 * UniversalConst::PI * freq;
		tauM = gettau(Ju*eta0, m, P, T, gs, model);
                double i2;

                gsl_integration_cquad_workspace *workspace = gsl_integration_cquad_workspace_alloc(100);

                double abserr;
                size_t nevals;

                struct intgrnd_jtak_params params = {omega, tauM, Ap, sigmap};
                gsl_function F;
                F.function = &intgrnd_j2tak;
                F.params = &params;

                gsl_integration_cquad(&F,
                                      0,
                                      std::numeric_limits<double>::infinity(),
                                      1e-12,
                                      1e-12,
                                      workspace,
                                      &i2,
                                      &abserr,
                                      &nevals);

                gsl_integration_cquad_workspace_free(workspace);
		
		J2byJu = i2 + (1.0 / (omega*tauM) );
	}
	
	
	return J2byJu;
}



double mccarthyXn(double taun)
{
	// This subroutine gets value of relaxation spectrum at a value of 
	// normalized time scale (taun) from McCarthy et al. (2011) eqn. 25
	
	double Xn;
	
	if (taun >= 1.0e-11) Xn = 0.32 * (std::pow(taun, (0.39-0.28/(1.0+2.6 * std::pow(taun, 0.1) ) ) ) );
	else Xn = 1853.0 * sqrt(taun);
	
	
	return Xn;
}



double getfn(double freq, double P, double T, double gs, std::string model, double Vs0, double rho)
{
	// Get the normalized frequency from McCarthy et al. (2011) eqn. 19 and 22
	
	
	// Check for an invalid model
	assert (model != "JF10_eBurg" || model != "Andrade" || model != "LOWM");
	
	
	double eta0, m;
	
	if (model == "M11")
	{
		eta0 = Rheology_Constants::M11::eta0;
		m = Rheology_Constants::M11::m;
	}
	else if (model == "Tak14")
	{
		eta0 = Rheology_Constants::Tak14::eta0;
		m = Rheology_Constants::Tak14::m;
	}
	else if (model == "P_M13")
	{
		eta0 = Rheology_Constants::P_M13::eta0;
		m = Rheology_Constants::P_M13::m;
	}
	
	
	double Ju;
	
	if ( std::isnan(Vs0) ) Ju = getJu(T,P,model);
	else Ju = getJu(T,P,model,Vs0,rho);
	
	
	double tau0 = Ju * eta0;
	double taur = gettau(tau0, m, P, T, gs, model);
	
	if (model == "P_M13") taur = gettauM_P_M13(T,P,Ju);
	
	double fn = freq * taur;
	
	return fn;
}



double getJu(double T, double P, std::string model, double Vs0, double rho)
{
	// This gets the unrelaxed modulus at T,P conditions based on constants from study
	
	double Ju;
	
	
	// Check for an invalid model
	assert (model != "LOWM");
	
	
	if ( std::isnan(Vs0) )
	{
		double GUR, dGdT, dGdP, TR, PR;
		
		if (model == "JF10_eBurg")
		{
			GUR = Rheology_Constants::JF10_eBurg::GUR;
			dGdT = Rheology_Constants::JF10_eBurg::dGdT; 
			dGdP = Rheology_Constants::JF10_eBurg::dGdP;
			TR = Rheology_Constants::JF10_eBurg::TR;
			PR = Rheology_Constants::JF10_eBurg::PR;
		}
		else if (model == "Andrade")
		{
			GUR = Rheology_Constants::Andrade::GUR;
			dGdT = Rheology_Constants::Andrade::dGdT; 
			dGdP = Rheology_Constants::Andrade::dGdP;
			TR = Rheology_Constants::Andrade::TR;
			PR = Rheology_Constants::Andrade::PR;
		}
		else if (model == "M11")
		{
			GUR = Rheology_Constants::M11::GUR;
			dGdT = Rheology_Constants::M11::dGdT;
			dGdP = Rheology_Constants::M11::dGdP;
			TR = Rheology_Constants::M11::TR;
			PR = Rheology_Constants::M11::PR;
		}
		else if (model == "Tak14")
		{
			GUR = Rheology_Constants::Tak14::GUR;
			dGdT = Rheology_Constants::Tak14::dGdT; 
			dGdP = Rheology_Constants::Tak14::dGdP;
			TR = Rheology_Constants::Tak14::TR;
			PR = Rheology_Constants::Tak14::PR;
		}
		else if (model == "P_M13")
		{
			GUR = Rheology_Constants::P_M13::GUR;
			dGdT = Rheology_Constants::P_M13::dGdT; 
			dGdP = Rheology_Constants::P_M13::dGdP;
			TR = Rheology_Constants::P_M13::TR;
			PR = Rheology_Constants::P_M13::PR;
		}
		
		
		Ju = 1.0*(GUR + dGdT*(T-TR) + dGdP*(P-PR));
	}
	else Ju = 1.0/(rho*std::pow(Vs0, 2.0));
	
	return Ju;
}



double gettau(double tau0, double m, double P, double T, double gs, std::string model)
{
	// calculate the relaxation time at P,T,gs, given model and reference(0) relaxation 
	// time and grainsize exponent
	
	
	// Check for an invalid model
	assert (model != "LOWM");
	
	
	double PR, TR, gsR, E, V;
	
	
	double R = Rheology_Constants::Common::R;
	
	
	if (model == "JF10_eBurg")
	{
		PR = Rheology_Constants::JF10_eBurg::PR;
		TR = Rheology_Constants::JF10_eBurg::TR;
		gsR = Rheology_Constants::JF10_eBurg::gsR;
		E = Rheology_Constants::JF10_eBurg::E;
		V = Rheology_Constants::JF10_eBurg::V;
	}
	else if (model == "Andrade")
	{
		PR = Rheology_Constants::Andrade::PR;
		TR = Rheology_Constants::Andrade::TR;
		gsR = Rheology_Constants::Andrade::gsR;
		E = Rheology_Constants::Andrade::E;
		V = Rheology_Constants::Andrade::V;
	}
	else if (model == "M11")
	{
		PR = Rheology_Constants::M11::PR;
		TR = Rheology_Constants::M11::TR;
		gsR = Rheology_Constants::M11::gsR;
		E = Rheology_Constants::M11::E;
		V = Rheology_Constants::M11::V;
	}
	else if (model == "Tak14")
	{
		PR = Rheology_Constants::Tak14::PR;
		TR = Rheology_Constants::Tak14::TR;
		gsR = Rheology_Constants::Tak14::gsR;
		E = Rheology_Constants::Tak14::E;
		V = Rheology_Constants::Tak14::V;
	}
	else if (model == "P_M13")
	{
		PR = Rheology_Constants::P_M13::PR;
		TR = Rheology_Constants::P_M13::TR;
		gsR = Rheology_Constants::P_M13::gsR;
		E = Rheology_Constants::P_M13::E;
		V = Rheology_Constants::P_M13::V;
	}
	
	
	if (P > 24.3e9)
	{
		V = 5e-6;	//# constant
//# 		V = -7.2727e-17*P + 1.1770e-05		# 10 to 2
//# 		V = -5.4545e-17*P + 1.13250e-05		# 10 to 4
//# 		V = -7.2727e-17*P + 1.3767e-05		# 12 to 4
//# 		Vact_lowM(P,Vtop,Vbot,Ptop,Pfold) where Pfold of 11e9 gives 1/2 LM die-out
//# 		V = Vact_lowM(P,15e-6,4e-6,24.3e9,11e9)
	}
	
	double tau = tau0 * std::pow((gs/gsR), m) * exp( (E/R)*(TR-T)/(TR*T) + (V/R)*(P*TR-PR*T)/(T*TR) );
	
	return tau;
}



double gettauM_P_M13(double T, double P, double Ju)
{
	// This gets the Maxwell relaxation time at T,P conditions based on constants
	// from study model='P_M13'
	
	double R = Rheology_Constants::Common::R;
	
	
	double TR = Rheology_Constants::P_M13::TR;
	double PR = Rheology_Constants::P_M13::PR;
	double E = Rheology_Constants::P_M13::E;
	double V = Rheology_Constants::P_M13::V;
	double eta0 = Rheology_Constants::P_M13::eta0;
	
	
	if (P > 24.3e9)
	{
		V = 5e-6;	//# constant
//# 		V = -7.2727e-17*P + 1.1770e-05		# 10 to 2
//# 		V = -5.4545e-17*P + 1.13250e-05		# 10 to 4
//# 		V = -7.2727e-17*P + 1.3767e-05		# 12 to 4
//# 		Vact_lowM(P,Vtop,Vbot,Ptop,Pfold) where Pfold of 11e9 gives 1/2 LM die-out
//# 		V = Vact_lowM(P,15e-6,4e-6,24.3e9,11e9)
	}
	
	
	double astar = exp((E - PR*V)/(R*TR))/exp((E-P*V)/(R*T));
	double eta = eta0/astar;
	double tauM = Ju*eta;
	
	return tauM;
}



/*
 * ------------------------------------------------------------------
 * 
 *		Integrations for JF10
 * 
 * ------------------------------------------------------------------
 */

/**
 * Function to calculate the function from J1B at a
 * specific coordinate xx.
 * 
 * Needs to be integrated over xx
 */
template <class T>
struct Function_J1B
{
	double alpha, omega;
	
	Function_J1B (double A, double w) : alpha(A), omega(w) {}
	double operator()(double xx)
	{
		return ( std::pow(xx, alpha-1.0) / (1.0 + std::pow(omega*xx, 2.0) ) );
	}
};
double intgrate_j1b(double limL, double limH, double alpha, double omega, int N) //*
{
	std::vector<double> xx = MyUtilities::logspace(log10(limL), log10(limH), N, 10.0);
	std::vector<double> yy(xx.size(),0.0);
	//for (unsigned int ii=0; ii<xx.size(); ++ii)
	//{
	//	yy[ii] = std::pow(xx[ii], alpha-1.0) / (1.0 + std::pow(omega*xx[ii], 2.0) );
	//}
	
	
	Function_J1B<double> J1Bfunc(alpha, omega);
	
	//yy[0] = 0.0; // !!!!! The vector is already initialized to zero, maybe this can go !!!!!
	
	
	// !!! SHIFT TO START FROM 1 TO XX.SIZE() AND FIRST ELEMENT IS 0 !!!
	for (unsigned int inc=1; inc<xx.size(); ++inc)
	{
		unsigned int inc_prev = inc - 1;
		double inc1 = xx[inc_prev]; //radius[irad];
		double inc2 = xx[inc]; //radius[irad+1];
		
		
		// Calculate the cumulative integral
		double current_increment = IntegrationRoutines::qromb(J1Bfunc,inc1,inc2, 1.0e-10);
		if (inc == 0)
		{
			yy[inc] = current_increment;
		}
		else
		{
			yy[inc] = yy[inc-1] + current_increment;
		}
	}
	
	// lines are from the python script. These can be deleted when everything else tested
	//I=integrate.cumtrapz(yy,xx)
	//I = I[len(I)-1]
	//return I;
	
	return yy[xx.size()-1]; // Only return the last element
}



/**
 * Function to calculate the function from J2B at a
 * specific coordinate xx.
 * 
 * Needs to be integrated over xx
 */
template <class T>
struct Function_J2B
{
	double alpha, omega;
	
	Function_J2B (double A, double w) : alpha(A), omega(w) {}
	double operator()(double xx)
	{
		return ( std::pow(xx, alpha) / (1.0 + std::pow(omega*xx, 2.0) ) );
	}
};
double intgrate_j2b(double limL, double limH, double alpha, double omega, int N) //*
{
	std::vector<double> xx = MyUtilities::logspace(log10(limL), log10(limH), N, 10.0);
	std::vector<double> yy(xx.size(),0.0);
	/*
	for (unsigned int ii=0; ii<xx.size(); ++ii)
	{
		yy[ii] = std::pow(xx[ii], alpha) / (1.0 + std::pow(omega*xx[ii], 2.0) );
	}
	*/
	
	Function_J2B<double> J2Bfunc(alpha, omega);
	
	
	// !!! SHIFT TO START FROM 1 TO XX.SIZE() AND FIRST ELEMENT IS 0 !!!
	for (unsigned int inc=1; inc<xx.size(); ++inc)
	{
		unsigned int inc_prev = inc - 1;
		double inc1 = xx[inc_prev];
		double inc2 = xx[inc];
		
		
		// Calculate the cumulative integral
		double current_increment = IntegrationRoutines::qromb(J2Bfunc,inc1,inc2, 1.0e-10);
		if (inc == 0)
		{
			yy[inc] = current_increment;
		}
		else
		{
			yy[inc] = yy[inc-1] + current_increment;
		}
	}
	/*
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I;
	*/
	return yy[xx.size()-1]; // Only return the last element
}



/**
 * Function to calculate the function from J1P at a
 * specific coordinate xx.
 * 
 * Needs to be integrated over xx
 */
template <class T>
struct Function_J1P
{
	double omega, tauP, sigma;
	
	Function_J1P (double w, double Tau, double Sig) : omega(w), tauP(Tau), sigma(Sig) {}
	double operator()(double xx)
	{
		return ( (1.0/xx) * (1.0/(1.0 + std::pow(omega*xx,2.0))) * exp(-std::pow(log(xx/tauP), 2.0) / (2.0*std::pow(sigma, 2.0))) );
	}
};
double intgrate_j1p(double limL, double limH, double omega, double tauP, double sigma, int N=1000) //*
{
	std::vector<double> xx = MyUtilities::logspace(log10(limL), log10(limH), N, 10.0);
	std::vector<double> yy(xx.size(),0.0);
	/*
	for (unsigned int ii=0; ii<xx.size(); ++ii)
	{
		yy[ii] = (1.0/xx[ii]) * (1.0/(1.0 + std::pow(omega*xx[ii],2.0))) * exp(-std::pow(log(xx[ii]/tauP), 2.0) / (2.0*std::pow(sigma, 2.0)));
	}
	*/
	
	Function_J1P<double> J1Pfunc(omega, tauP, sigma);
	
	
	// !!! SHIFT TO START FROM 1 TO XX.SIZE() AND FIRST ELEMENT IS 0 !!!
	for (unsigned int inc=1; inc<xx.size(); ++inc)
	{
		unsigned int inc_prev = inc - 1;
		double inc1 = xx[inc_prev];
		double inc2 = xx[inc];
		
		
		// Calculate the cumulative integral
		double current_increment = IntegrationRoutines::qromb(J1Pfunc,inc1,inc2, 1.0e-10);
		if (inc == 0)
		{
			yy[inc] = current_increment;
		}
		else
		{
			yy[inc] = yy[inc-1] + current_increment;
		}
	}
	
	/*
	I=integrate.cumtrapz(yy,xx)
	I = I[len(I)-1]
	return I;
	*/
	return yy[xx.size()-1]; // Only return the last element
}


int main ()
{
  std::pair<double, double> VsQ = PTd_VsQ(1.0,0.0,1000,1e-4, 3300, "M11", 1.0, std::numeric_limits<double>::quiet_NaN());
  std::cout << VsQ.first << " " << VsQ.second << std::endl;
  return 0;
}
