#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include <string>
#include <map>

//#include "nr3.h"
//#include "hash_mod.h"
//#include "ran_mod.h"

using namespace std;


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


namespace MyUtilities
{


	/**
	 * Some templated functions from:
	 * Numerical Recipes. The Art of Scientific Computing. 3rd ed. Press et. al.
	 */
	template<class T>
	inline const T &MAX(const T &a, const T &b)
			{return b > a ? (b) : (a);}
	
	inline float MAX(const double &a, const float &b)
			{return b > a ? (b) : float(a);}
	
	inline float MAX(const float &a, const double &b)
			{return b > a ? float(b) : (a);}
	
	
	
	/**
	* Function to calculate the normalized difference
	*/
	double normdiff (double newval, double refval)
	{
		if (refval == 0.0) return (newval - refval);
		else return ((newval - refval) / refval);
	}
	
	
	
	/**
	 * Split a string by a delimiter into a vector of string parts.
	 */
	std::vector<std::string> split(const std::string& s, char delimiter)
	{
		// Vector of string to save tokens 
		std::vector<std::string> tokens;
		
		std::string token;
		
		// Turn input string into a stringstream
		std::istringstream tokenStream(s);
		
		// Tokenizing w.r.t. delimiter
		while (std::getline(tokenStream, token, delimiter))
		{
			tokens.push_back(token);
		}
		return tokens;
	}
	
	
	
	/**
	 * Simple function that takes in a vector as input and returns
	 * a sum over elements starting from element a to element b,
	 * including the endpoints elements a and b.
	 * 
	 * This function assumes that element a < element b,
	 * and that every element inbetween is summed and counted,
	 * (consecutive).
	 * 
	 * There are two versions of this function to deal with
	 * a simple 1-d vector, and 2-d vector.
	 * 
	 * Specified in the input arguments are:
	 * 		input vector
	 * 		element a
	 * 		element b
	 * 		vector row (for 2-d vectors only)
	 */
	double sumvec_consecutive (std::vector<double> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		
		for (unsigned int i=element_a, count=0; i<(element_b+1); ++i, ++count)
		{
			sum = sum + vector_in[i];
		}
		return sum;
	}
	int sumvec_consecutive (std::vector<int> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			sum = sum + vector_in[i];
		}
		return sum;
	}
	double sumvec_consecutive (std::vector<std::vector<double> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
		}
		return sum;
	}
	int sumvec_consecutive (std::vector<std::vector<int> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
		}
		return sum;
	}
	
	
	
	/**
	 * Simple function that takes in a vector as input,
	 * sums over elements starting from element a to element b,
	 * including the endpoints elements a and b, and then
	 * returns the arithmitic average.
	 * 
	 * This function assumes that element a < element b,
	 * and that every element inbetween is summed and counted,
	 * (consecutive).
	 * 
	 * There are two versions of this function to deal with
	 * a simple 1-d vector, and 2-d vector.
	 * 
	 * Specified in the input arguments are:
	 * 		input vector
	 * 		element a
	 * 		element b
	 * 		vector row (for 2-d vectors only)
	 */
	double avgvec_consecutive (std::vector<double> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			sum = sum + vector_in[i];
			++count;
		}
		return (sum/count);
	}
	int avgvec_consecutive (std::vector<int> vector_in, int element_a, int element_b)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			sum = sum + vector_in[i];
			++count;
		}
		return (sum/count);
	}
	double avgvec_consecutive (std::vector<std::vector<double> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
			++count;
		}
		return (sum/count);
	}
	double avgvec_consecutive (std::vector<std::vector<int> > vector_in, int element_a, int element_b, int sumaxis, int fixed_axis_element)
	{
		double sum = 0;
		int count = 0;
		
		for (unsigned int i=element_a; i<(element_b+1); ++i)
		{
			if (sumaxis==0) sum = sum + vector_in[i][fixed_axis_element];
			else if (sumaxis==1) sum = sum + vector_in[fixed_axis_element][i];
			++count;
		}
		return (sum/count);
	}
	
	
	
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
	 * !!!	UNFINISHED	!!!
	 * Funciton to calculate gravity at radii
	 * with given density.
	 * 
	 * !!!	MAY NOT NEED THIS FUNCTION AT ALL !!!
	 * !!!	THIS FUNCTION DOES NOTHING AT THE MOMENT	!!!
	 */
	double compute_gravity(std::vector <double> density, std::vector <double> radius, int irad, const double Gconst)
	{
		// Version based on Mark Panning's
		// Matlab code for Mars.
		
		const double PI = 3.141592653589793;
		return 4.0*PI;
		//return (4.0*PI*density*radius*radius);
		
		
		//gravity_grav[1:] = Gconst*mass_at_radii[1:]/(radii_grav[1:]**2.0)
		//gravity_grav[0] = 0.0
	}
	
	
	
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



namespace SearchUtilities
{
  // A place to put file search utilities
  
  
  /**
   * Functions to define Points and Boxes
   */
  //  --------------------------------------------------
  template <int dim>
  struct Point
  {
    // Simple structure to represent a point in dim dimensions
    
    double x[dim];		// the coordinates
    
    // Copy constructor
    Point(const Point &p)
    {
      for (int i=0; i<dim; i++) x[i] = p.x[i];
    }
    
    // Assignment operator
    Point& operator= (const Point &p)
    {
      for (int i=0; i<dim; i++) x[i] = p.x[i];
      return *this;
    }
    
    bool operator== (const Point &p) const
    {
      for (int i=0; i<dim; i++) if (x[i] != p.x[i]) return false;
      return true;
    }
    
    Point(double x0=0.0, double x1=0.0, double x2=0.0)
    {
      // Constructor by coordinate values. Arguments beyond the required
      // number are not used and can be omitted.
      
      x[0] = x0;
      if (dim > 1) x[1] = x1;
      if (dim > 2) x[2] = x2;
      if (dim > 3) throw("Point not implemented for dim > 3");
    }
  };
  
  
  template <int dim>
  double dist(const Point<dim> &p, const Point<dim> &q)
  {
    // Returns the distance between two points in dim dimensions
    
    double dd = 0.0;
    
    for (int j=0; j<dim; j++) dd += (q.x[j]-p.x[j])*(q.x[j]-p.x[j]);
    
    return sqrt(dd);
  }
  
  
  template <int dim>
  struct Box
  {
    // Structure to represent a Cartesian box in dim dimensions
    
    // Diagonally opposite corners (min of all coodinates and max of all coordinates) are stored as two points
    Point<dim> lo, hi;
    
    Box() {}
    Box(const Point<dim> &mylo, const Point<dim> &myhi) : lo(mylo), hi(myhi) {}
  };
  
  
  template <int dim>
  double dist(const Box<dim> &b, const Point<dim> &p)
  {
    // If point p lies outside box b, the distance to the nearest point on b is returned.
    // If p is inside b or on its surface, zero is returned.
    
    double dd = 0;
    
    for (int i=0; i<dim; i++)
    {
      if (p.x[i]<b.lo.x[i]) dd += (p.x[i]-b.lo.x[i])*(p.x[i]-b.lo.x[i]);
      if (p.x[i]>b.hi.x[i]) dd += (p.x[i]-b.lo.x[i])*(p.x[i]-b.hi.x[i]);
    }
    return sqrt(dd);
  }
  
  
  template <int dim>
  struct Boxnode : Box<dim>
  {
    // Node in a binary tree of boxes containing points. See text for details.
    
    int mom, dau1, dau2, ptlo, pthi;
    
    Boxnode() {}
    Boxnode(Point<dim> mylo, Point<dim> myhi, int mymom, int myd1, int myd2, int myptlo, int mypthi) : Box<dim>(mylo, myhi), mom(mymom), dau1(myd1), dau2(myd2), ptlo(myptlo), pthi(mypthi) {}
  };
}
  //  --------------------------------------------------


namespace InterpolationRoutines
{
	// A place to put general interpolation stuff
	
	struct Base_interp
	{
		int n, mm, jsav, cor, dj;
		const double *xx, *yy;
		Base_interp(std::vector<double> &x, const double *y, int m)
			: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
			dj = std::max(1,(int)pow((double)n,0.25));
		}

		double interp(double x) {
			int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo,x);
		}

		int locate(const double x);
		int hunt(const double x);
		
		double virtual rawinterp(int jlo, double x) = 0;

	};
	int Base_interp::locate(const double x)
	{
		int ju,jm,jl;
		if (n < 2 || mm < 2 || mm > n) throw("locate size error");
		bool ascnd=(xx[n-1] >= xx[0]);
		jl=0;
		ju=n-1;
		while (ju-jl > 1) {
			jm = (ju+jl) >> 1;
			if (x >= xx[jm] == ascnd)
				jl=jm;
			else
				ju=jm;
		}
		cor = abs(jl-jsav) > dj ? 0 : 1;
		jsav = jl;
		return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
	}
	int Base_interp::hunt(const double x)
	{
		int jl=jsav, jm, ju, inc=1;
		if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
		bool ascnd=(xx[n-1] >= xx[0]);
		if (jl < 0 || jl > n-1) {
			jl=0;
			ju=n-1;
		} else {
			if (x >= xx[jl] == ascnd) {
				for (;;) {
					ju = jl + inc;
					if (ju >= n-1) { ju = n-1; break;}
					else if (x < xx[ju] == ascnd) break;
					else {
						jl = ju;
						inc += inc;
					}
				}
			} else {
				ju = jl;
				for (;;) {
					jl = jl - inc;
					if (jl <= 0) { jl = 0; break;}
					else if (x >= xx[jl] == ascnd) break;
					else {
						ju = jl;
						inc += inc;
					}
				}
			}
		}
		while (ju-jl > 1) {
			jm = (ju+jl) >> 1;
			if (x >= xx[jm] == ascnd)
				jl=jm;
			else
				ju=jm;
		}
		cor = abs(jl-jsav) > dj ? 0 : 1;
		jsav = jl;
		return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
	}
	
	
	
	struct Linear_interp : Base_interp
	{
		Linear_interp(std::vector<double> &xv, std::vector<double> &yv)
			: Base_interp(xv,&yv[0],2)  {}
		double rawinterp(int j, double x) {
			if (xx[j]==xx[j+1]) return yy[j];
			else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
		}
	};
	
	
	
	struct Poly_interp : Base_interp
	{
		double dy;
		Poly_interp(std::vector<double> &xv, std::vector<double> &yv, int m)
			: Base_interp(xv,&yv[0],m), dy(0.) {}
		double rawinterp(int jl, double x);
	};

	double Poly_interp::rawinterp(int jl, double x)
	{
		int i,m,ns=0;
		double y,den,dif,dift,ho,hp,w;
		const double *xa = &xx[jl], *ya = &yy[jl];
		std::vector<double> c(mm),d(mm);
		dif=abs(x-xa[0]);
		for (i=0;i<mm;i++) {
			if ((dift=abs(x-xa[i])) < dif) {
				ns=i;
				dif=dift;
			}
			c[i]=ya[i];
			d[i]=ya[i];
		}
		y=ya[ns--];
		for (m=1;m<mm;m++) {
			for (i=0;i<mm-m;i++) {
				ho=xa[i]-x;
				hp=xa[i+m]-x;
				w=c[i+1]-d[i];
				if ((den=ho-hp) == 0.0) throw("Poly_interp error");
				den=w/den;
				d[i]=hp*den;
				c[i]=ho*den;
			}
			y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
		}
		return y;
	}
}

namespace IntegrationRoutines
{
	/**
	 * Polynominal interpolation routines from:
	 * Numerical Recipes. The Art of Scientific Computing. 3rd ed. Press et. al.
	 * 
	 * Each algorithm is marked with the section number it is found in.
	 * Comments are completely preserved.
	 */
	
	// Sec 3.1
	struct Base_interp
	{
		// Abstract base class used by all interpolation routines in this
		// chapter. Only the routine interp is called directly by the user.
		
		int n, mm, jsav, cor, dj;
		const double *xx, *yy;
		
		// Constructor: Set up for interpolating on a table of
		// x's and y's of length m. Normally called by a derived class,
		// not by the user.
		Base_interp (const std::vector <double> &x, const double *y, int m) : n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y)
		{
			dj = std::max(1, (int)pow((double)n, 0.25 ) );
		}
		
		// Given a value x, return an interpolated value, using data pointed
		// to by xx and yy.
		double interp(double x)
		{
			int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo, x);
		}
		
		int locate(const double x);  // see definitions below
		int hunt(const double x);
		
		// Derived classes provide this as the actual interpolation method.
		double virtual rawinterp(int jlo, double x) = 0;
		
	};
	
	// Sec 3.2
	struct Poly_interp : Base_interp
	{
		// Polynomial interpolation object. Construct with x and y vectors,
		// and the number M of points to be used locally (polynomial order
		// plus one), then call interp for interpolated values.
		
		double dy;
		
		Poly_interp (const std::vector <double> &xv, const std::vector <double> &yv, int m) : Base_interp (xv, &yv[0], m), dy(0.) {}
		double rawinterp (int jl, double x);
	};
	
	// Sec 3.2
	double Poly_interp::rawinterp(int jl, double x)
	{
		// Given a value x, and using pointers to data xx and yy, this routine
		// returns an interpolated value y, and stores an error estimate dy.
		// The returned value is obtained by mm-point polynomial interpolation
		// on the subrange xx[jl..jl+mm-1].
		
		int i, m, ns=0;
		double y, den, dif, dift, ho, hp, w;
		const double *xa = &xx[jl], *ya = &yy[jl];
		
		std::vector <double> c(mm), d(mm);
		
		dif = abs(x - xa[0]);
		
		// Here we find the index ns of the closest table entry,
		// and initialize the tableau of c's and d's.
		for (i=0; i<mm; i++)
		{
			if ( (dift = abs(x-xa[i]) ) < dif )
			{
				ns = i;
				dif = dift;
			}
			
			c[i] = ya[i];
			d[i] = ya[i];
		}
		
		// This is the initial approximation to y.
		y = ya[ns--];
		
		// For each column of the tableau,
		// we loop over the current c's and d's and update them.
		for (m=1; m<mm; m++)
		{
			for (i=0; i<mm-m; i++)
			{
				ho = xa[i] - x;
				hp = xa[i+m] - x;
				w = c[i+1] - d[i];
				
				// This error can occur only if two input xa's are
				// (to within roundoff) identical.
				if ( (den = ho - hp) == 0.0 ) throw("Poly_interp error");
				
				den = w / den;
				
				// Here the c's and d's are updated.
				d[i] = hp * den;
				c[i] = ho * den;
			}
			
			// After each column in the tableau is completed, we decide
			// which correction, c or d, we want to add to our accumulating
			// value of y, i.e., which path to take through the tableau
			// --- forking up or down. We do this in such a way as to take
			// the most "straight line" route through the tableau to its
			// apex, updating ns accordingly to keep track of where we are.
			// This route keeps the partial approximations centered
			// (insofar as possible) on the target x. The last dy added
			// is thus the error indication.
			y += ( dy = (2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]) );
		}
		
		return y;
	}
	
	
	
	
	/**
	 * Integration routines from:
	 * Numerical Recipes. The Art of Scientific Computing. 3rd ed. Press et. al.
	 * 
	 * Each algorithm is marked with the section number it is found in.
	 * Comments are completely preserved.
	 */
	
	// Sec 4.2
	struct Quadrature{
		// Abstract base class for elementary quadrature algorithms
		int n;  // Current level of refinement
		virtual double next() = 0;
		// Returns the value of the integral at the nth stage of refinement.
		// The function next() must be defined in the derived class.
	};
	
	
	
	// Sec 4.2
	template<class T>
	struct Trapzd : Quadrature
	{
		// Routine implementing the extended trapezoidal rule.
		double a, b, s;  // Limits of integration and current value of integral
		T &func;
		Trapzd() {};
		Trapzd(T &funcc, const double aa, const double bb) : func(funcc), a(aa), b(bb) {n=0;}
		// The constructor takes as inputs func, the function or functor
		// to be integrated between limits a and b, also input.
		double next()
		{
			// Returns the nth stage of refinement of the extended trapezoidal
			// rule. On the first call (n=1), the routine returns the crudest
			// estimate of INTERGRAL{f(x)dx} from a to h. Subsequent calls
			// set n=2,3,... and improve the accuracy by adding
			// 2^(n-2) additional interior points.
			
			double x, tnm, sum, del;
			int it, j;
			n++;
			
			if (n == 1)
			{
				return ( s = 0.5 * (b-a) * (func(a) + func(b)) );
			}
			else
			{
				for (it=1,j=1; j<n-1; j++) it <<= 1;
				tnm = it;
				del = (b-a) / tnm;  // This is the spacing of the points to be added
				x = a + 0.5 * del;
				for (sum=0.0,j=0; j<it; j++,x+=del) sum += func(x);
				s = 0.5 * (s + (b-a) * sum / tnm);  // This replaces s by its refined value
				return s;
			}
		}
	};
	
	
	// Sec 4.2
	template<class T>
	struct Trapzd_int : Quadrature
	{
		// Routine implementing the extended trapezoidal rule.
		double s;  // Limits of integration and current value of integral
		int a, b;
		T &func;
		Trapzd_int() {};
		Trapzd_int(T &funcc, const int aa, const int bb) : func(funcc), a(aa), b(bb) {n=0;}
		// The constructor takes as inputs func, the function or functor
		// to be integrated between limits a and b, also input.
		double next()
		{
			// Returns the nth stage of refinement of the extended trapezoidal
			// rule. On the first call (n=1), the routine returns the crudest
			// estimate of INTERGRAL{f(x)dx} from a to h. Subsequent calls
			// set n=2,3,... and improve the accuracy by adding
			// 2^(n-2) additional interior points.
			
			int x;
			double tnm, sum, del;
			int it, j;
			n++;
			cout << n << "\n";
			if (n == 1)
			{
				return ( s = 0.5 * (b-a) * (func(a) + func(b)) );
			}
			else
			{
				for (it=1,j=1; j<n-1; j++) it <<= 1;
				tnm = it;
				del = (b-a) / tnm;  // This is the spacing of the points to be added
				x = int(a + 0.5 * del);
				for (sum=0.0,j=0; j<it; j++,x+=int(del)) sum += func(x);
				s = 0.5 * (s + (b-a) * sum / tnm);  // This replaces s by its refined value
				cout << s << "\n";
				return s;
			}
		}
	};
	
	
	
	// Sec 4.3
	template <class T>
	double qromb (T &func, double a, double b, const double eps=1.0e-10)
	{
		// Returns the integral of the function or functor func from a to b.
		// Integration is performed by Romberg's method of order 2K, where
		// , e.g., K=2 is Simpson's rule.
		
		
		// Here EPS is the fractional accuracy desired, as determined by the
		// extrapolation error estimate; JMAX limits the total number of steps;
		// K is the number of points used in the extrapolation.
		const int JMAX=20, JMAXP=JMAX+1, K=5;
		
		
		// These store the successive trapezoidal approximations and their
		// relative stepsizes.
		std::vector <double> s(JMAX), h(JMAXP);
		
		Poly_interp polint(h, s, K);
		
		h[0] = 1.0;
		
		Trapzd<T> t(func, a, b);
		
		for (int j=1; j<=JMAX; j++)
		{
			s[j-1] = t.next();
			
			if (j >= K)
			{
				double ss = polint.rawinterp (j-K, 0.0);
				
				if ( abs(polint.dy) <= eps*abs(ss) ) return ss;
			}
			
			h[j] = 0.25 * h[j-1];
			// This is a key step: The factor is 0.25 even though the stepsize
			// is decreased by only 0.5. This makes the extrapolation a
			// polynomial in h^2 as allowed by NumericalRecipes_eq (4.2.1),
			// not just a polynomial in h.
		}
		
		throw("Too many steps in routine qromb");
	}
	
	
	// Sec 4.3 -- Overloaded Class -- Modified to accept integers a and b
	template <class T>
	double qromb (T &func, int a, int b, const double eps=1.0e-10)
	{
		// Returns the integral of the function or functor func from a to b.
		// Integration is performed by Romberg's method of order 2K, where
		// , e.g., K=2 is Simpson's rule.
		
		
		// Here EPS is the fractional accuracy desired, as determined by the
		// extrapolation error estimate; JMAX limits the total number of steps;
		// K is the number of points used in the extrapolation.
		const int JMAX=20, JMAXP=JMAX+1, K=5;
		
		
		// These store the successive trapezoidal approximations and their
		// relative stepsizes.
		std::vector <double> s(JMAX), h(JMAXP);
		
		Poly_interp polint(h, s, K);
		
		h[0] = 1.0;
		
		Trapzd_int<T> t(func, a, b);
		
		for (int j=1; j<=JMAX; j++)
		{
			s[j-1] = t.next();
			
			if (j >= K)
			{
				double ss = polint.rawinterp (j-K, 0.0);
				
				if ( abs(polint.dy) <= eps*abs(ss) ) return ss;
			}
			
			h[j] = 0.25 * h[j-1];
			// This is a key step: The factor is 0.25 even though the stepsize
			// is decreased by only 0.5. This makes the extrapolation a
			// polynomial in h^2 as allowed by NumericalRecipes_eq (4.2.1),
			// not just a polynomial in h.
		}
		
		throw("Too many steps in routine qromb");
	}
}



namespace PerpleXUtilities
{
	template <int querydim, int datadim>
	class PerpleXLookup
	{
		/** Things and steps this function should do:
		 * ------------------------------------------
		 * -- Load the necessary n-dim coordinate and m-dim data from a table
		 * -- Read in a target point that needs to be assessed
		 * -- Read in the interpolation type to find the corresponding
		 * -- Setup the KDTree from dealii
		 *    -- Load only coordinate points into a KDTree
		 *    -- Use the KDTree to return the n-number of closest points to a target point using:
		 *    -- KDTree< dim >::get_closest_points (const Point< dim > &target, const unsigned int n_points) const
		 *    data that goes with the target point
		 * -- Return data values for use
		 */
		
		
		public:
		std::vector<double> get_data (std::string filein, double T, double P); //(const MPI_Comm &comm); // Load the data table into memory
		//void closest_points; // Setup the KDTree and return the specified number of closest points
		//void data_interpolation; // Interpolate the data values based on data at the closest points
		std::string get_header (std::string filein);
		std::vector< SearchUtilities::Point<querydim> > load_coord (std::string filein);
		std::vector<std::vector <double> > load_data (std::string filein);
	};
	
	
	
	template <int querydim, int datadim>
	std::vector<double>
	PerpleXLookup<querydim,datadim>::get_data (std::string filein, double T, double P) { //(const MPI_Comm &comm) {
		
		
		std::vector<double> coords(querydim); // !!!	CREATE A POINT CLASS TO DEAL WITH COORDINATES	!!!
		std::vector<double> data;
		double dummy, rho, vp, vs;
		double mindist;
		int counter;
		
		std::vector<std::string> keys = {"T(K)","P(bar)","rho,kg/m3","alpha,1/K","cp,J/K/kg","vp,km/s","vs,km/s","h,J/kg"};
		//P(bar) T(K) rho,kg/m3 alpha,1/K beta,1/bar Ks,bar Gs,bar v0,km/s vp,km/s vs,km/s s,J/K/kg h,J/kg cp,J/K/kg V,J/bar/mol
		
		
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
			
		// Skip header
		for (unsigned int i=0; i<12; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
		
		double P_found, T_found;  // !!!	TEMPORARY STATEMENT FOR CHECKING CODE	!!!
		
		counter = 0;
		// Read in the full table
		while (getline(in,temp)) //(in.good() && !in.eof())   // ????
		{
			std::istringstream line(temp);
			
			// read the two components of the coordinates
			double P_in, T_in;
			line >> P_in >> T_in;
			
			// !!!	TEMPORARY CONVERSION -- BAR -> Pa	!!!
			P_in = P_in * 1.0e5;
			
			
			// Find the distance
			double Pdiff = P_in - P;
			double Tdiff = T_in - T;
			if (counter == 0)
			{
				mindist = std::sqrt( Pdiff*Pdiff + Tdiff*Tdiff );
				line >> rho >> dummy >> dummy >> dummy >> dummy >> dummy >> vp >> vs;
				//P_found = P_in;
				//T_found = T_in;
			}
			
			double currentdist = std::sqrt( Pdiff*Pdiff + Tdiff*Tdiff );
			
			if (currentdist < mindist)
			{
				mindist = currentdist;
				line >> rho >> dummy >> dummy >> dummy >> dummy >> dummy >> vp >> vs;
				//P_found = P_in;
				//T_found = T_in;
			}
			
			++counter;
			
		}
		
		// Populate the data vector
		data.emplace_back (rho);
		data.emplace_back (vp*1.0e3);
		data.emplace_back (vs*1.0e3);
		
		//std::cout << " T " << T << "  T_found " << T_found << "  P " << P << "  P_found " << P_found;
		//std::cout << "  RHO " << rho << "  VP " << vp << "  VS " << vs << "\n";
		
		in.close();
		
		return data;
	}
	
	
	template <int querydim, int datadim>
	std::string
	PerpleXLookup<querydim,datadim>::get_header (std::string filein) { //(const MPI_Comm &comm) {
		
				
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
		// Skip header
		for (unsigned int i=0; i<12; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
        
        
		// Read the header that describes the data columns
		std::getline(in,temp);
		
		in.close();
		
		
		return temp;
	}
	
	
	
	template <int querydim, int datadim>
	std::vector< SearchUtilities::Point<querydim> >
	PerpleXLookup<querydim,datadim>::load_coord (std::string filein) { //(const MPI_Comm &comm) {
		
		
		std::vector< SearchUtilities::Point<querydim> > coords;
		//std::vector<double> coords(querydim); // !!!	CREATE A POINT CLASS TO DEAL WITH COORDINATES	!!!
				
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
			
		// Skip header
		for (unsigned int i=0; i<13; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
		
		
		while (getline(in,temp)) //(in.good() && !in.eof())   // ????
		{
			std::istringstream line(temp);
			SearchUtilities::Point<querydim> coordinates;
			// read the two components of the coordinates
			double P_in, T_in;
			line >> P_in >> T_in;
			
			//std::get<0>(coordinates) = P_in;
			coordinates.x[0] = P_in;
			coordinates.x[1] = T_in;
			//std::cout << "P " << P_in << " T " << T_in << " Coord " << coordinates.x[0] << "," << coordinates.x[1] << "\n"; 
			////std::vector<double> coordinates(querydim);
			//line >> coordinates;
			coords.emplace_back(coordinates);	
			//coords.emplace_back(coordinates.x[0],coordinates.x[1]);
		}
		
		in.close();
		
		//std::cout << "A point: " << coords[2].x[0] << "," << coords[2].x[1] << "\n";
		
		return coords;
	}
	
	
	
	template <int querydim, int datadim>
	std::vector<std::vector <double> >
	PerpleXLookup<querydim,datadim>::load_data (std::string filein) { //(const MPI_Comm &comm) {
		
		
		std::vector<std::vector <double> > data;
		
		std::string temp;
		// Read data from disk and distribute among processes
		std::ifstream in;//(Utilities::read_and_distribute_file_content(filename,comm));
		in.open(filein);
		
		
			
		// Skip header
		for (unsigned int i=0; i<13; ++i)
		{
			std::getline(in,temp); // throw these lines away
		}
		
		
		while (getline(in,temp)) //(in.good() && !in.eof())   // ????
		{
			std::istringstream line(temp);
			
			// read the two components of the coordinates
			double P_in, T_in;
			line >> P_in >> T_in;
			
			std::vector<double> values;
			for (unsigned int i=0; i<datadim; ++i) // !!! NEED TO CHANGE THIS TO FIND END OF LINE INSTEAD OF HARD-CODED NUMBER OF COLUMNS !!!
			{
				double value;
				line >> value;

				values.emplace_back(value);
			}
			data.emplace_back (values);
		}
		
		in.close();
		
		return data;
	}
	
}



namespace SeismicUtilities
{
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


