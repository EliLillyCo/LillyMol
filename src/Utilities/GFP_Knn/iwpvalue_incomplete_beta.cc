#include <stdlib.h>
#include <iostream>

#include <boost/math/special_functions/beta.hpp>

double
niwpvalue (int d, double x)
{
/*std::cerr << "Computing value for d = " << d << " x = " << x << std::endl;
//return 1.0 - 0.5 * boost::math::ibeta(d / (d + x*x), static_cast<double>(d)/2.0, 0.5);
  std::cerr << (d / (d + x*x)) << "\n";
  for (int i = 0; i < 20; i++)
  {
    double y = i / 20.0;

    std::cerr << y << ' ' << boost::math::ibeta(d/2.0, 0.5, y) << "\n";
  }*/
  return 1.0 - 0.5 * boost::math::ibeta(static_cast<double>(d)/2.0, 0.5, d / (d+x*x));
}
