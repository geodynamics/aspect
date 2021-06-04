#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <assert.h>

#include <string>
#include <map>


using namespace std;


/*
 * ------------------------------------------------------------------
 * TODO:
 * ------------------------------------------------------------------
 */


namespace UniversalConst
{
	const double PI = 3.141592653589793;
    //const long double PI = 3.14159265358979323846264338328L;
	
	// Conversion from degree to radians
	const double DEG2RAD = PI/180.0;
	
	// Conversion from radians to degrees
	const double RAD2DEG = 180.0/PI;
	
	// Gravitational constant
	// G=(6.67384+-0.00080)x10^11 m3\kg s2
	// [Mohr et al., 2012], Williams et al. 2014
	const double G_gravconst = 6.67384e-11;
}



namespace EarthConst
{
	// Mean radius:
	// Various (2000). David R. Lide, ed. Handbook of Chemistry and Physics (81st ed.). CRC. ISBN 0-8493-0481-4. -- Wikipedia [6]
	const double mean_radius = 6.3710e6;  // meters [6]
	
	
	// Equatorial radius:
	//"Selected Astronomical Constants, 2011". The Astronomical Almanac. Archived from the original on 26 August 2013. Retrieved 25 February 2011. -- Wikipedia [7]
	// World Geodetic System (WGS-84). Available online from National Geospatial-Intelligence Agency. -- Wikipedia [8]
	const double equatorial_radius = 6.3781e6;  // meters [7][8]
	
	
	// Polar radius:
	// Cazenave, Anny (1995). "Geoid, Topography and Distribution of Landforms" (PDF). In Ahrens, Thomas J. Global Earth Physics: A Handbook of Physical Constants. Washington, DC: American Geophysical Union. ISBN 0-87590-851-9. Archived from the original (PDF) on 16 October 2006. Retrieved 3 August 2008. -- Wikipedia [9]
	const double polar_radius = 6.3568e6;  // meters  [9]
	
	
	// Flattening:
	// International Earth Rotation and Reference Systems Service (IERS) Working Group (2004). "General Definitions and Numerical Standards" (PDF). In McCarthy, Dennis D.; Petit, Gerard. IERS Conventions (2003) (PDF). IERS Technical Note No. 32. Frankfurt am Main: Verlag des Bundesamts fur Kartographie und Geodasie. p. 12. ISBN 3-89888-884-3. Retrieved 29 April 2016. -- Wikipedia [10]
	const double flattening = 0.0033528;  // [10]
	//earth_flat = 1.0/298.257222101 # (ETRS89)
	
	
	// Average Atmospheric Pressure on the surface of the Earth [Pa]
	// !!! CITATION NEED !!!
	const double surface_pressure = 101325.0;
	
	
	// PREM Outer Core Radius
	const double PREM_CoreRad = 3480000.000000;

}



namespace VectorUtilities
{


	/**
	 * Function to create and return a linearly
	 * increasing, evenly spaced vector.
	 * 
	 * Inputs:
	 * 			a = "start point" of type T
	 * 			b = "end point" of type T
	 * 			N = "number of increments" of type int
	 * 
	 * Output:	std::vector of type T
	 *
	 * Modified from https://gist.github.com/lorenzoriano/5414671
	 *	20181005
	 */
	template <typename T>
	std::vector<T> linspace(T a, T b, size_t N)
	{
		T h = (b - a) / static_cast<T>(N-1);
		std::vector<T> xs(N);
		typename std::vector<T>::iterator x;
		T val;
		for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) *x = val;
		
		
		return xs;
	}
	
	
	
	/**
	 * Function to create and return a linearly
	 * increasing, evenly spaced vector on log scale.
	 * 
	 * Inputs:
	 * 			a = "start point" of type T
	 * 			b = "end point" of type T
	 * 			N = "number of increments" of type int
	 * 			base = "the log base value" of type T
	 * 					the default value is 10.0
	 * 
	 * Output:	std::vector of type T
	 *
	 * Modified from https://gist.github.com/lorenzoriano/5414671
	 *	20181005
	 */
	template <typename T>
	std::vector<T> logspace(T a, T b, size_t N, T base=10)
	{
		T h = (b - a) / static_cast<T>(N-1);
		std::vector<T> xs(N);
		typename std::vector<T>::iterator x;
		T val;
		for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) *x = pow(base,val);
		
		
		return xs;
	}
	
}



namespace Utilities
{

	/**
	 * Function to calculate the mass of a spherical
	 * shell with a given density, where density is
	 * constant within the shell.
	 */
	template <class T>
	struct CalcCumulativeMass
	{
		double density;
		
		CalcCumulativeMass (double den) : density(den) {}
		double operator()(double radius)
		{
			return (4.0*UniversalConst::PI*density*radius*radius); //Needs to be integrated over radius
		}
	};
	
	
	
		
	/**
	 * Funciton to calculate pressure at radii
	 * with given density and gravity.
	 * 
	 */
	template <class T>
	struct ComputePressure
	{
		std::vector<double> density, gravity, depth_increasing;
		
		ComputePressure (std::vector<double> den, std::vector<double> grav, std::vector<double> dep) : density(den), gravity(grav), depth_increasing(dep) {}
		double operator()(unsigned int irad, unsigned int idep1, unsigned int idep2)
		{
			return (density[irad]*gravity[irad] * (depth_increasing[idep2] - depth_increasing[idep1]) );
		}
	};
	
	
	
	/**
	 * Compute the average density of a series of spherical
	 * shells from the surface down to a current depth.
	 */
	template <class T>
	struct ComputeAverageDensity
	{
		std::vector<double> density, radius;
		
		ComputeAverageDensity (std::vector<double> den, std::vector<double> rad) : density(den), radius(rad) {}
		double operator()(int r1, int r2)
		{
			double avgdensity = 0.0; //= ( radius[r1] - radius[r1-1] ) * density[r1];
			
			for (unsigned int i=r1; i<r2; ++i)
			{
				avgdensity += ( radius[i+1] - radius[i] ) * density[i];
			}
			/*
			for (unsigned int i=r1+1; i==r2; --i)
			{
				avgdensity += ( radius[i] - radius[i+1] ) * density[i];
				//avgdensity = density[i];
				//cout << i << "\n";
			}
			*/
			
			return ( avgdensity / (radius[r2] - radius[r1]) );
		}
	};
	
	
	
	/**
	 * Function to compute the average density of several material components.
	 * 
	 * The general scheme is to determine the total mass and total volume
	 * of the combined material components. This is done by determining
	 * the fraction of volume and mass for each component, and summing
	 * them up. The average density is then simply:
	 * 
	 * 	AVGDENSITY = TOTMASS / TOTVOL
	 * 
	 * This particular function is valid for non-interactive solids.
	 * In the case of interactive solids and liquids, it is possible
	 * that the combined volume is not the same as the sum of
	 * individual volumes. Thus, more information would be needed to
	 * calculate the average density.
	 */
	double avg_density_multmatcomp (std::vector<double> densities, std::vector<double> volumes, std::vector<double> fractions)
	{
		double total_mass, total_volume;
		double avg_density;
		
		total_mass = total_volume = 0.0;
		
		for (unsigned int i=0; i<densities.size(); ++i)
		{
			total_mass = total_mass + fractions[i]*densities[i]*volumes[i];
			total_volume = total_volume + fractions[i]*volumes[i];
		}
		avg_density = total_mass / total_volume;
		
		return avg_density;
	}
	
	
	
	/**
	 * The two following functions compute the Voigt averages of the
	 * elastic moduli, Bulk Modulus and Shear Modulus.
	 * 
	 * Voigt averages are the "upper" bound estimate of the true average.
	 * It can be interpreted as the ratio of average strain within a
	 * composite material. In other words, the Voigt average is found by
	 * assuming that the strain is uniform throughout the composite
	 * material (iso-strain), which yields an upper bound estimate.
	 */
	double BulkMod_voigt (std::vector<double> bulkmods, std::vector<double> fractions)
	{
		double avg_bulkmod = 0.0;
		
		for (unsigned int i=0; i<bulkmods.size(); ++i)
		{
			avg_bulkmod = avg_bulkmod + fractions[i]*bulkmods[i];
		}
		
		return avg_bulkmod;
	}
	double ShearMod_voigt (std::vector<double> shearmods, std::vector<double> fractions)
	{
		double avg_shearmod = 0.0;
		
		for (unsigned int i=0; i<shearmods.size(); ++i)
		{
			avg_shearmod = avg_shearmod + fractions[i]*shearmods[i];
		}
		
		return avg_shearmod;
	}
	
	
	
	/**
	 * The two following functions compute the Reuss averages of the
	 * elastic moduli, Bulk Modulus and Shear Modulus.
	 * 
	 * Reuss averages are the "lower" bound estimate of the true average.
	 * It can be interpreted as the ratio of average stress within a
	 * composite material. In other words, the Reuss average is found by
	 * assuming that the stress is uniform throughout the composite
	 * material (iso-stress), which yields a lower bound estimate.
	 */
	double BulkMod_reuss (std::vector<double> bulkmods, std::vector<double> fractions)
	{
		double avg_bulkmod = 0.0;
		
		for (unsigned int i=0; i<bulkmods.size(); ++i)
		{
			avg_bulkmod = avg_bulkmod + ( fractions[i] / bulkmods[i] );
		}
		
		avg_bulkmod = 1.0 / avg_bulkmod;
		
		return avg_bulkmod;
	}
	double ShearMod_reuss (std::vector<double> shearmods, std::vector<double> fractions)
	{
		double avg_shearmod = 0.0;
		
		for (unsigned int i=0; i<shearmods.size(); ++i)
		{
			avg_shearmod = avg_shearmod + ( fractions[i] / shearmods[i] );
		}
		
		avg_shearmod = 1.0 / avg_shearmod;
		
		return avg_shearmod;
	}
	
	
	
	/**
	 * The two following functions compute the Voigt-Reuss-Hill averages
	 * of the elastic moduli, Bulk Modulus and Shear Modulus, from the
	 * upper (Voigt) and lower (Reuss) bound estimates.
	 * 
	 * These funcions call the Voigt and Reuss averaging functions,
	 * then compute the average of the results to provide an estimate
	 * of the true average.
	 */
	double bulkmod_voigtreusshill (std::vector<double> bulkmods, std::vector<double> fractions)
	{
		double bulkmod_voigtavg = BulkMod_voigt(bulkmods, fractions);
		double bulkmod_reussavg = BulkMod_reuss(bulkmods, fractions);
		
		return ( (bulkmod_voigtavg + bulkmod_reussavg) / 2.0);
	}
	double shearmod_voigtreusshill (std::vector<double> shearmods, std::vector<double> fractions)
	{
		double shearmod_voigtavg = ShearMod_voigt(shearmods, fractions);
		double shearmod_reussavg = ShearMod_reuss(shearmods, fractions);
		
		return ( (shearmod_voigtavg + shearmod_reussavg) / 2.0);
	}
	
	
	
	/**
	 * The two following functions compute average seismic velocities
	 * given average density, bulk and shear moduli.
	 */
	double vp_avg (double BulkMod_avg, double ShearMod_avg, double Density_avg)
	{
		return sqrt((BulkMod_avg+(4.0/3.0)*ShearMod_avg)/Density_avg);
	}
	double vs_avg (double ShearMod_avg, double Density_avg)
	{
		return sqrt(ShearMod_avg/Density_avg);
	}
	
	
	
	/**
	 * The following functions calculate bulk sound speed, vphi,
	 * given the isotropic Vp, Vs, and density.
	 */
	double
	calc_vphi_from_vpvs (double vp, double vs)
	{
		return ( sqrt( (vp*vp) - ((4.0/3.0)*(vs*vs)) ) );
	}
	double
	calc_vphi_from_bulk_den (double bulkmod, double density)
	{
		return ( sqrt( bulkmod/density ) );
	}
	
	
	
	/**
	 * The following two functions calculate the
	 * isotropic bulk and shear moduli given the
	 * isotropic Vp, Vs, and density.
	 */
	double
	calc_bulkmod_from_vpvsden (double vp, double vs, double density)
	{
		return ( ( (vp*vp) - ( (4.0/3.0)*(vs*vs) ) ) * density );
	}
	double
	calc_shearmod_from_vsden (double vs, double density)
	{
		return ( ( (vs*vs) * density ) );
	}
	
	
	
	
	/**
    * The following functions convert the anisotropic Pv, Ph, Sv, and Sh values
    * to the Voigt average isotropic Vp, Vs.
    * 
    * Conversion equations are from Appendix A of:
    *      Panning, Mark, and Barbara Romanowicz.
    *      "A three-dimensional radially anisotropic model of shear velocity in the whole mantle."
    *      Geophysical Journal International 167.1 (2006): 361-379.
    *
    * Eqn (A8): Vs*Vs = (2*Vsv*Vsv + Vsh*Vsh) / 3
    * 
    * Eqn (A9): Vp*Vp = (Vpv*Vpv + 4*Vph*Vph) / 5
    * 
    * 
    * These are approximations of:
    * Eqn (A5): rho*Vp*Vp = (1/15) * (3C + (8 + 4*eta)*A + 8*(1 − eta)*L)
    * 
    * and
    * 
    * Eqn (A6): rho*Vs*Vs = (1/15) * (C + (1 + 2*eta)*A + (6 + 4*eta)*L + 5*N)
    * where eta = F / (A - 2*L), and A, C, L, N, F are the Love parameters (Love 1927).
    * 
    * For eta ~ 1, eqn's A5 and A6 simplify to A8 and A9, respectively.
    * See reference for more detail.
    * 
    * The eta value in the anisotropic PREM is equal to 1.0 except in the upper
    * mantle, where it is >0.9. A more accurate calculation would be to use eqn's
    * A5 and A6 for these depths, and should be revisited.
    * 
    * 
    * Written 21 Sep 2018. PMBremner
    **/
    double convert_PvPh_to_isoVp (double Vpv, double Vph)
    {
        double isoVp=0;
        
        // Calculate the Voigt average isotropic velocity
        
        // Eqn (A9):
        isoVp = sqrt( (Vpv*Vpv + 4.0*Vph*Vph) / 5.0 );
        
        return isoVp;
    }
    double convert_SvSh_to_isoVs (double Vsv, double Vsh)
    {
        double isoVs=0;
        
        // Calculate the Voigt average isotropic velocity
        
        // Eqn (A8)
        isoVs = sqrt( (2.0*Vsv*Vsv + Vsh*Vsh) / 3.0 );
        
        return isoVs;
    }
        

	
	
	
	/**
	 * Function to calculate and populate a vector of depth increments
	 * from a vector of radii. The depth vector is populated such that
	 * the first element is the surface, and succesive elements move
	 * deeper. In other words, depth increases with increasing radii
	 * elements.
	 * 
	 * It is assumed the input radii are increasing from some depth
	 * onward toward the surface with increasing elements, and that
	 * the final radius element is the surface.
	 */
	std::vector<double> make_increasing_depth(std::vector<double> radius_in)
	{
		std::vector<double> depth_increasing(radius_in.size());
		for (unsigned int irad=0, idepth=radius_in.size()-1; irad<radius_in.size(); ++irad, --idepth)
		{
			depth_increasing[idepth] = EarthConst::mean_radius - radius_in[irad];
		}
		return depth_increasing;
	}
	
	
}




namespace PREMUtilities
{
	
	class PREMLookup
	{
	public:
		std::vector<double> load_coord (std::string filename, bool IncludeCore);
		std::vector<double> load_avg_den_vp_vs_overdepthrange (std::string filename, double mindepth, double maxdepth);
		std::vector< std::vector<double> > load_den_vp_vs (std::string filename, bool IncludeCore);
		std::vector< std::vector<double> > load_coord_and_den_vp_vs (std::string filename, bool IncludeCore);
		std::vector< std::vector<double> > load_coord_and_all_data (std::string filename, bool IncludeCore);
	};
	
	
	
	std::vector<double>
	PREMLookup::load_coord (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<double> radii;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			
			// Read in radius
			double radius;
			line >> radius;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true) radii.emplace_back (radius);
			else if (radius > EarthConst::PREM_CoreRad) radii.emplace_back (radius);
			
		}
		
		infile.close();
		
		return radii;
	}
	
	
	
	std::vector<double>
	PREMLookup::load_avg_den_vp_vs_overdepthrange (std::string filename, double mindepth, double maxdepth)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<double> data;
		double density_avg, Vp_avg, Vs_avg;
		double density_sum=0.0, Vp_sum=0.0, Vs_sum=0.0;
		int depthcount = 0;
		
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (radius and data: Vp,Vs,Rho)
			double radius;
			double density, Vp, Vs;
			line >> radius >> density >> Vp >> Vs;
			
			double PREMdepth = EarthConst::mean_radius - radius;
			if (PREMdepth>=mindepth && PREMdepth<=maxdepth)
			{
				density_sum = density_sum + density;
				Vp_sum = Vp_sum + Vp;
				Vs_sum = Vs_sum + Vs;
				depthcount = ++depthcount;
			}
			
		}
		
		infile.close();
		
		density_avg = density_sum / depthcount;
		Vp_avg = Vp_sum / depthcount;
		Vs_avg = Vs_sum / depthcount;
		
		data.emplace_back (density_avg);
		data.emplace_back (Vp_avg);
		data.emplace_back (Vs_avg);
		
		
		return data;
	}
	
	
	
	std::vector< std::vector<double> >
	PREMLookup::load_den_vp_vs (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (radius and data: Vp,Vs,Rho)
			double radius;
			double density, Vp, Vs;
			line >> radius >> density >> Vp >> Vs;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true)
			{
				
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
			else if (radius > EarthConst::PREM_CoreRad)
			{
				
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
		}
		
		infile.close();
		
		return data;
	}
	
	
	std::vector< std::vector<double> >
	PREMLookup::load_coord_and_den_vp_vs (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data Point (radius and data: Vp,Vs,Rho)
			double radius;
			double density, Vp, Vs;
			line >> radius >> density >> Vp >> Vs;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
			else if (radius > EarthConst::PREM_CoreRad)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				
				data.emplace_back (values);
			}
		}
		
		infile.close();
		
		return data;
	}
	
	
	std::vector< std::vector<double> >
	PREMLookup::load_coord_and_all_data (std::string filename, bool IncludeCore=true)
	{
		// Open Input File
		std::ifstream infile;
		infile.open(filename);
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		// Skip the header lines
		for (unsigned int ihead=0; ihead<3; ++ihead)
		{
			std::getline(infile,temp); //throw these lines away
			std::cout << temp << "\n";
		}
		
		
		// Collect radius and data
		while (getline(infile,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double radius;
			double density, Vp, Vs, Q_kappa, Q_mu, eta;
			line >> radius >> density >> Vp >> Vs >> Q_kappa >> Q_mu >> eta;
			
			if (radius <= EarthConst::PREM_CoreRad && IncludeCore==true)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				values.emplace_back (Q_kappa);
				values.emplace_back (Q_mu);
				values.emplace_back (eta);
				
				
				data.emplace_back (values);
			}
			else if (radius > EarthConst::PREM_CoreRad)
			{
				values.emplace_back (radius);
				values.emplace_back (density);
				values.emplace_back (Vp);
				values.emplace_back (Vs);
				values.emplace_back (Q_kappa);
				values.emplace_back (Q_mu);
				values.emplace_back (eta);
				
				
				data.emplace_back (values);
			}
		}
		
		infile.close();
		
		return data;
	}
}



namespace GypsumUtilities
{
	
	//template <class T>
	class GyPSuMLookup
	{
	public:
		std::vector<double> load_minmaxdepth (std::string filename);
		std::vector< std::vector<double> > load_coord (std::string filename);
		std::vector< std::vector<double> > load_data_perturb (std::string filename);
		std::vector< std::vector<double> > load_data_absolute (std::string filename, std::string reffile);
		std::vector< std::vector<double> > load_coord_data_absolute (std::string filename, std::string reffile);
		double perturb_to_absolute (double pertval, double refval);
		std::vector<double> load_temperature_variation (std::string filename);
		std::vector<double> load_pressure_variation (std::string filename);
		void make_combined_data (std::string densityfile, std::string vpfile, std::string vsfile, std::string outfile);
	};
	
	
	
	std::vector<double>
	GyPSuMLookup::load_minmaxdepth (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<double> depths;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		depths.emplace_back (mindepth);
		depths.emplace_back (maxdepth);
		
		in_f.close();
		
		
		return depths;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_coord (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<std::vector <double> > coordinates;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			line >> lat >> lon;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			
			coordinates.emplace_back (values);
		}
		
		in_f.close();
		
		
		return coordinates;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_data_perturb (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<std::vector <double> > data;
		
		
		in_f.open(filename);
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double density, vp, vs;
			line >> lat >> lon >> density >> vp >> vs;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (density);
			values.emplace_back (vp);
			values.emplace_back (vs);
			
			data.emplace_back (values);
		}
		
		in_f.close();
		
		
		return data;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_data_absolute (std::string filename, std::string reffile)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		PREMUtilities::PREMLookup lookup;
		GyPSuMLookup pert_abs;
		std::string temp;
		std::vector<std::vector <double> > data;
		
		
		in_f.open(filename);
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		// Retrieve the avg PREM values over the considered depth range
		// *** BE SURE TO FEED IN MINDEPTH AND MAXDEPTH IN METERS ***
		std::vector<double> ref_den_vp_vs = lookup.load_avg_den_vp_vs_overdepthrange(reffile, mindepth*1000.0, maxdepth*1000.0);
		std::cout << ref_den_vp_vs[0] << " " << ref_den_vp_vs[1] << " " << ref_den_vp_vs[2] << "\n";
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double density, vp, vs;
			line >> lat >> lon >> density >> vp >> vs;
			
			// Covert to absolute values
			// *** BE SURE TO CONVERT dlnrho, dlnVp, dlnVs from PERCENT PERTURBATION TO DECIMAL PERTURBATION ***
			density = pert_abs.perturb_to_absolute(density/100.0, ref_den_vp_vs[0]);
			vp = pert_abs.perturb_to_absolute(vp/100.0, ref_den_vp_vs[1]);
			vs = pert_abs.perturb_to_absolute(vs/100.0, ref_den_vp_vs[2]);
			
			//values.emplace_back (lat);
			//values.emplace_back (lon);
			values.emplace_back (density);
			values.emplace_back (vp);
			values.emplace_back (vs);
			
			data.emplace_back (values);
		}
		
		in_f.close();
		
		
		return data;
	}
	
	
	
	std::vector< std::vector<double> >
	GyPSuMLookup::load_coord_data_absolute (std::string filename, std::string reffile)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		PREMUtilities::PREMLookup lookup;
		GyPSuMLookup pert_abs;
		std::string temp;
		std::vector<std::vector <double> > data;
		
		
		in_f.open(filename);
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		// Retrieve the avg PREM values over the considered depth range
		// *** BE SURE TO FEED IN MINDEPTH AND MAXDEPTH IN METERS ***
		std::vector<double> ref_den_vp_vs = lookup.load_avg_den_vp_vs_overdepthrange(reffile, mindepth*1000.0, maxdepth*1000.0);
		std::cout << ref_den_vp_vs[0] << " " << ref_den_vp_vs[1] << " " << ref_den_vp_vs[2] << "\n";
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double density, vp, vs;
			line >> lat >> lon >> density >> vp >> vs;
			
			// Covert to absolute values
			// *** BE SURE TO CONVERT dlnrho, dlnVp, dlnVs from PERCENT PERTURBATION TO DECIMAL PERTURBATION ***
			density = pert_abs.perturb_to_absolute(density/100.0, ref_den_vp_vs[0]);
			vp = pert_abs.perturb_to_absolute(vp/100.0, ref_den_vp_vs[1]);
			vs = pert_abs.perturb_to_absolute(vs/100.0, ref_den_vp_vs[2]);
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (density);
			values.emplace_back (vp);
			values.emplace_back (vs);
			
			data.emplace_back (values);
		}
		
		in_f.close();
		
		
		return data;
	}
	
	
	
	double
	GyPSuMLookup::perturb_to_absolute (double pertval, double refval)
	{
		return ( (pertval*refval) + refval);
	}
	
	
	
	std::vector<double>
	GyPSuMLookup::load_temperature_variation (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<double> temperature_v;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			
			// Read In Input Data
			double lat, lon, Tvar;
			line >> lat >> lon >> Tvar;
			
			temperature_v.emplace_back (Tvar);
		}
		
		in_f.close();
		
		
		return temperature_v;
	}
	
	
	
	std::vector<double>
	GyPSuMLookup::load_pressure_variation (std::string filename)
	{
		// Open Input File
		std::ifstream in_f;
		
		
		std::string temp;
		std::vector<double> pressure_v;
		
		
		in_f.open(filename);
		
		
		// Read in the min and max depths  (in km)
		double mindepth, maxdepth;
		std::getline(in_f,temp);
		std::istringstream line(temp);
		
		line >> mindepth >> maxdepth;
		
		
		// Collect coordinates and data
		while (getline(in_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon, Pvar;
			line >> lat >> lon;
			
			pressure_v.emplace_back (Pvar);
		}
		
		in_f.close();
		
		
		return pressure_v;
	}
	
	
	
	void
	GyPSuMLookup::make_combined_data (std::string densityfile, std::string vpfile, std::string vsfile, std::string outfile)
	{
		// Open Input File
		std::ifstream density_f, vp_f, vs_f;
		
		
		std::string temp;
		std::vector<std::vector <double> > densities, vpvals, vsvals;
		
		
		// -------------------------------------------
		//
		//  READ IN THE DENSITY FILE
		//
		// -------------------------------------------
		
		density_f.open(densityfile);
		
		// Read in the min and max depths  (in km)
		double density_mindepth, density_maxdepth;
		std::getline(density_f,temp);
		std::istringstream line(temp);
		
		line >> density_mindepth >> density_maxdepth;
		
		
		// Collect coordinates and data
		while (getline(density_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double value;
			line >> lat >> lon >> value;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (value);
			
			densities.emplace_back (values);
		}
		
		density_f.close();
		
		
		// -------------------------------------------
		//
		//  READ IN THE VP FILE
		//
		// -------------------------------------------
		
		vp_f.open(vpfile);
		
		// Read in the min and max depths  (in km)
		double vp_mindepth, vp_maxdepth;
		std::getline(vp_f,temp);
		std::istringstream line2(temp);
		
		line2 >> vp_mindepth >> vp_maxdepth;
		
		
		// Collect coordinates and data
		while (getline(vp_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double value;
			line >> lat >> lon >> value;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (value);
			
			vpvals.emplace_back (values);
		}
		
		vp_f.close();
		
		
		// -------------------------------------------
		//
		//  READ IN THE VS FILE
		//
		// -------------------------------------------
		
		vs_f.open(vsfile);
		
		// Read in the min and max depths  (in km)
		double vs_mindepth, vs_maxdepth;
		std::getline(vs_f,temp);
		std::istringstream line3(temp);
		
		line3 >> vs_mindepth >> vs_maxdepth;
		
		
		// Collect coordinates and data
		while (getline(vs_f,temp))
		{
			std::istringstream line(temp);
			std::vector<double> values;
			
			// Read In Input Data
			double lat, lon;
			double value;
			line >> lat >> lon >> value;
			
			values.emplace_back (lat);
			values.emplace_back (lon);
			values.emplace_back (value);
			
			vsvals.emplace_back (values);
		}
		
		vs_f.close();
		
		// -------------------------------------------
		//
		//  CHECK THE (DEPTH, LAT, LON) AND COMBINE
		//
		// -------------------------------------------
		
		// Check that all depths are the same
		if (density_mindepth != vp_mindepth)
		{
			std::cout << density_mindepth << " " << vp_mindepth << "\n";
			throw("Min depth ranges don't match, check input files for density and vp!");
		}
		if (density_maxdepth != vp_maxdepth)
		{
			std::cout << density_maxdepth << " " << vp_maxdepth << "\n";
			throw("Max depth ranges don't match, check input files for density and vp!");
		}
		
		if (vs_mindepth != vp_mindepth)
		{
			std::cout << vs_mindepth << " " << vp_mindepth << "\n";
			throw("Min depth ranges don't match, check input files for vp and vs!");
		}
		if (vs_maxdepth != vp_maxdepth)
		{
			std::cout << vs_maxdepth << " " << vp_maxdepth << "\n";
			throw("Max depth ranges don't match, check input files for vp and vs!");
		}
		
		
		// Check the number of data points
		if (densities.size() != vpvals.size()) throw("Data lengths don't match! Check input files for density and vp.");
		if (vsvals.size() != vpvals.size()) throw("Data lengths don't match! Check input files for vp and vs.");
		
		
		// Name and open the output file
		std::ofstream out_f;
		out_f.open(outfile);
		
		out_f << density_mindepth << " " << density_maxdepth << "\n";
		
		
		// Check coordinates...if good, paste the columns together
		for (unsigned int i=0; i<densities.size(); ++i)
		{
			if (densities[i][0]==vpvals[i][0] && densities[i][1]==vpvals[i][1])
			{
				if (vsvals[i][0]==vpvals[i][0] && vsvals[i][1]==vpvals[i][1])
				{
					out_f << densities[i][0] << " " << densities[i][1] << " " << densities[i][2] << " " << vpvals[i][2] << " " << vsvals[i][2] << "\n";
				}
			}
		}
		
		out_f.close();
		
		
		return;
	}
	
	
	
}


