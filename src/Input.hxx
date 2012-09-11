#ifndef _c011f36f_e983_4b63_95d8_7318fe25da5f
#define _c011f36f_e983_4b63_95d8_7318fe25da5f

#include "Constants.hxx"

namespace Input
{
  // IO:
  // const double plot_time_freq = 100000000.0e-09;
  const double plot_time_freq = 1.0e-20;
  
  // EOS:
  const double gamma = 5.0/3.0;
  const double abar = 4.0;

  // GEOMETRY:
  const GeomType::E geom = GeomType::SPHERICAL;
  const unsigned ncells = 2000;

  // TIME:
  // const unsigned nmax = 10000000;
  const unsigned nmax = 1000;
  const double tmax = 500.0e-09;
  const double dtinit = 1.0e-11;
  const double growth = 1.05;
  
  // HYDRO:
  const double clin = 0.2;
  const double cquad = 1.0;
  const double cfl = 0.25;
  const HydroBc::E hydbc_left = HydroBc::FIXED;
  const HydroBc::E hydbc_right = HydroBc::FIXED;  

  const bool UseArtVisc = true;
}

#endif
