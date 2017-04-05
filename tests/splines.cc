/**
 * This tests the Utilities::tk::spline class for the linear case,
 * the spline case and the monotone spline case.
 */
#include <aspect/utilities.h>
#include <iostream>

int f()
{
  using namespace aspect::Utilities;
  tk::spline linear,spline,monotone_spline;

  std::vector<double> x_list(5,0.0);
  std::vector<double> y_list(5,0.0);

  x_list = {1,2,3,4,5};
  y_list = {0,5,15,2,6};

  linear.set_points(x_list,y_list,false,false);
  spline.set_points(x_list,y_list,true,false);
  monotone_spline.set_points(x_list,y_list,true,true);

  std::cout << "Linear interpolates: "
            << "0 = " << linear(0) << ", 0.25 = " << linear(0.25) << ", 0.5 = " << linear(0.5) << ", 0.75 = " << linear(0.75)
            << ", 1 = " << linear(1) << ", 1.25 = " << linear(1.25) << ", 1.5 = " << linear(1.5) << ", 1.75 = " << linear(1.75)
            << ", 2 = " << linear(2) << ", 2.25 = " << linear(2.25) << ", 2.5 = " << linear(2.5) << ", 2.75 = " << linear(2.75) << std::endl
            << "                     3 = " << linear(3) << ", 3.25 = " << linear(3.25) << ", 3.5 = " << linear(3.5) << ", 3.75 = " << linear(3.75)
            << ", 4 = " << linear(4) << ", 4.25 = " << linear(4.25) << ", 4.5 = " << linear(4.5) << ", 4,75 = " << linear(4.75)
            << ", 5 = " << linear(5) << ", 5.25 = " << linear(5.25) << ", 5.5 = " << linear(5.5) << ", 5.75 = " << linear(5.75)
            << ", 6 = " << linear(6) << std::endl;

  std::cout << "Spline interpolates: "
            << "0 = " << spline(0) << ", 0.25 = " << spline(0.25) << ", 0.5 = " << spline(0.5) << ", 0.75 = " << spline(0.75)
            << ", 1 = " << spline(1) << ", 1.25 = " << spline(1.25) << ", 1.5 = " << spline(1.5) << ", 1.75 = " << spline(1.75)
            << ", 2 = " << spline(2) << ", 2.25 = " << spline(2.25) << ", 2.5 = " << spline(2.5) << ", 2.75 = " << spline(2.75) << std::endl
            << "                     3 = " << spline(3) << ", 3.25 = " << spline(3.25) << ", 3.5 = " << spline(3.5) << ", 3.75 = " << spline(3.75)
            << ", 4 = " << spline(4) << ", 4.25 = " << spline(4.25) << ", 4.5 = " << spline(4.5) << ", 4,75 = " << spline(4.75)
            << ", 5 = " << spline(5) << ", 5.25 = " << spline(5.25) << ", 5.5 = " << spline(5.5) << ", 5.75 = " << spline(5.75)
            << ", 6 = " << spline(6) << std::endl;

  std::cout << "Monotone spline interpolates: "
            << "0 = " << monotone_spline(0) << ", 0.25 = " << monotone_spline(0.25) << ", 0.5 = " << monotone_spline(0.5) << ", 0.75 = " << monotone_spline(0.75)
            << ", 1 = " << monotone_spline(1) << ", 1.25 = " << monotone_spline(1.25) << ", 1.5 = " << monotone_spline(1.5) << ", 1.75 = " << monotone_spline(1.75)
            << ", 2 = " << monotone_spline(2) << ", 2.25 = " << monotone_spline(2.25) << ", 2.5 = " << monotone_spline(2.5) << ", 2.75 = " << monotone_spline(2.75) << std::endl
            << "                              3 = " << monotone_spline(3) << ", 3.25 = " << monotone_spline(3.25) << ", 3.5 = " << monotone_spline(3.5) << ", 3.75 = " << monotone_spline(3.75)
            << ", 4 = " << monotone_spline(4) << ", 4.25 = " << monotone_spline(4.25) << ", 4.5 = " << monotone_spline(4.5) << ", 4,75 = " << monotone_spline(4.75)
            << ", 5 = " << monotone_spline(5) << ", 5.25 = " << monotone_spline(5.25) << ", 5.5 = " << monotone_spline(5.5) << ", 5.75 = " << monotone_spline(5.75)
            << ", 6 = " << monotone_spline(6) << std::endl;

  exit(0);
  return 42;
}
// run this function by initializing a global variable by it
int i = f();

