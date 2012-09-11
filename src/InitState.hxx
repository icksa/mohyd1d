#ifndef _f9ee293a_3c50_48b4_8234_09936cd5c62f
#define _f9ee293a_3c50_48b4_8234_09936cd5c62f

// Set the initial conditions:
template <typename E>
void init_state(E &env)
{
  const double XMAX = 100.0e-04;
  unsigned nverts = env.nverts;
  unsigned ncells = env.ncells;

  double *coords = env.coords;
  double *velx_data = env.velx_data;
  for(unsigned i = 0; i < nverts; i++) {
    coords[i] = i*(XMAX/ncells);
    velx_data[i] = 0.0;
  }

  double *dens_data = env.dens_data;
  double *tele_data = env.tele_data;
  for(unsigned i = 0; i < ncells; i++) {
    double xc = 0.5 * (coords[i] + coords[i+1]);
    if(xc < 50.0e-04) {
      dens_data[i] = 0.0018; 
    } else {
      dens_data[i] = 0.0006; 
    }
    tele_data[i] = 0.025;
  }

  eos_dens_temp(env);
}

#endif
