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
