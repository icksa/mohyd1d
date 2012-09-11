#include "moab/Core.hpp"
#include "moab/ScdInterface.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include "InitState.hxx"
#include "Input.hxx"
#include "Constants.hxx"

#define MBERR(msg, rval) mberr(msg, rval, __LINE__, __FILE__)
#define ERR(msg) err(msg, __LINE__, __FILE__)

void mberr(const std::string& msg, moab::ErrorCode rval, 
	   unsigned line, const std::string &file)
{
  if(rval != moab::MB_SUCCESS) {
    std::cerr << file << ":" << line << ", MOAB Error: " << rval << "\n"
	      << msg << std::endl;
    throw(rval);
  }
}

void err(const std::string& msg, 
	 unsigned line, const std::string &file)
{
  std::cerr << file << ":" << line << ", ERROR\n"
	    << msg << std::endl;
  throw(msg);
}

class Mesh
{

public:
  Mesh(unsigned nc) :
    M_ncells(nc),
    M_mbcore(),
    M_scdint(NULL)
  {
    moab::ErrorCode rval;

    // Create the structured mesh interface:
    rval = this->mbint().query_interface(M_scdint);
    MBERR("mbcore.query_interface", rval);

    // Create the mesh:
    std::vector<double> coords(3*this->nverts());

    for (unsigned i = 0; i < coords.size(); i += 3) {
      coords[i+0] = double(i/3);
      coords[i+1] = 0.0;
      coords[i+2] = 0.0;
    }

    moab::ScdBox *scdbox = NULL;
    
    rval = M_scdint->construct_box(moab::HomCoord(0,0,0,1), 
				   moab::HomCoord(this->ncells(),0,0,1),
				   coords.data(), this->nverts(), scdbox);
    MBERR("scdint->construct_box", rval);

    // Mesh has been constructed, create ranges which represent the
    // vertexes and cells:
    M_verts = moab::Range(scdbox->start_vertex(), scdbox->start_vertex() + this->nverts()-1);
    M_cells = moab::Range(scdbox->start_element(), scdbox->start_element() + this->ncells()-1);
  }

  const unsigned ncells() const {
    return M_ncells;
  }

  const unsigned nverts() const {
    return M_ncells+1;
  }

  moab::Tag create_tag(const std::string& name)
  {
    moab::ErrorCode rval;
    double zero = 0.0;
    moab::Tag ret;
    rval = this->mbint().tag_get_handle(name.c_str(), 1, moab::MB_TYPE_DOUBLE, ret, 
					moab::MB_TAG_DENSE | moab::MB_TAG_CREAT, 
					&zero);
    MBERR("mbcore.tag_get_handle("+name+")", rval);
    return ret;
  }

  double * tag_data(moab::Tag tag, moab::EntityType etype)
  {
    moab::ErrorCode rval;
    void *ret;
    int count;
    moab::Range *prange = NULL;

    switch(etype) {
    case moab::MBVERTEX:
      prange = &this->M_verts;
      break;
    case moab::MBEDGE:
      prange = &this->M_cells;
      break;
    default :
      ERR("Bad entity type");
    }

    rval = this->mbint().tag_iterate(tag, prange->begin(), prange->end(), count, ret);
    MBERR("mbcore.tag_iterate", rval);

    if(unsigned(count) != prange->size()) {
      ERR("Bad count from tag_iterate");
    }

    return static_cast<double*>(ret);
  }

  /// Create tag and return data
  double * tag_data(const std::string& name, moab::EntityType etype) {
    return this->tag_data(this->create_tag(name), etype);
  }

  const double * tag_data(moab::Tag tag, moab::EntityType etype) const 
  {
    return const_cast<Mesh*>(this)->tag_data(tag, etype);
  }

  double *coords()
  {
    int count;
    double *xvals;
    double *yvals; 
    double *zvals;
    this->mbint().coords_iterate(M_verts.begin(), M_verts.end(),
				  xvals, yvals, zvals, count);
    return xvals;
  }

  const double *coords() const { return const_cast<Mesh*>(this)->coords(); }

  void write_file(const std::string& name) {
    this->mbint().write_file(name.c_str());
  }

protected:
  unsigned M_ncells;
  moab::Core M_mbcore;
  moab::ScdInterface *M_scdint;
  moab::Range M_verts;
  moab::Range M_cells;

  const moab::Interface& mbint() const {
    return M_mbcore;
  }

  moab::Interface& mbint() {
    return M_mbcore;
  }
};


template <bool UseArtVisc> struct Env_ArtVisc_Base; 

template <> struct Env_ArtVisc_Base<false> {}; 

template <> struct Env_ArtVisc_Base<true>
{
  Env_ArtVisc_Base(Mesh &mesh) :
    pvis_data( mesh.tag_data("dens", moab::MBEDGE) )
  {
  }

  double * const pvis_data;
};

typedef Env_ArtVisc_Base<Input::UseArtVisc> Env_ArtVisc;


struct Env_Top : Env_ArtVisc {

  Env_Top(Mesh &mesh) :
    Env_ArtVisc(mesh),
    ncells(mesh.ncells()),
    nverts(mesh.nverts()),
    coords( mesh.coords() ), 
    velx_data( mesh.tag_data("velx", moab::MBVERTEX) ), 
    dens_data( mesh.tag_data("dens", moab::MBEDGE) ), 
    tele_data( mesh.tag_data("tele", moab::MBEDGE) ), 
    eint_data( mesh.tag_data("eint", moab::MBEDGE) ), 
    pres_data( mesh.tag_data("pres", moab::MBEDGE) ), 
    cmas_data( mesh.tag_data("cmas", moab::MBEDGE) ), 
    vmas_data( mesh.tag_data("vmas", moab::MBVERTEX) ), 
    cspd_data( mesh.tag_data("cspd", moab::MBEDGE) ) ,
    pvis_data( mesh.tag_data("pvis", moab::MBEDGE) )
  {
  }

  const unsigned ncells;
  const unsigned nverts;

  double * const coords;
  double * const velx_data;
  double * const dens_data;
  double * const tele_data;
  double * const eint_data;
  double * const pres_data;
  double * const cmas_data;
  double * const vmas_data;
  double * const cspd_data;
  double * const pvis_data;

  static const double plot_time_freq;

  // Configure time parameters:
  static const GeomType::E geom = Input::geom;

  static const double gamma;
  static const double abar;
  static const double growth;
  static const double cfl;
  static const double clin;
  static const double cquad;

  static const HydroBc::E hydbc_left = Input::hydbc_left;
  static const HydroBc::E hydbc_right = Input::hydbc_right;
};

const double Env_Top::plot_time_freq = Input::plot_time_freq;

const double Env_Top::gamma  = Input::gamma;
const double Env_Top::abar   = Input::abar;
const double Env_Top::growth = Input::growth;
const double Env_Top::cfl    = Input::cfl;
const double Env_Top::clin   = Input::clin;
const double Env_Top::cquad  = Input::cquad;

// ***************************
// *                         *
// *     LOCAL FUNCTIONS     *
// *                         *
// ***************************

template <typename E>
void eos_dens_temp(E &env) {
  unsigned ncells = env.ncells;

  for(unsigned i = 0; i < ncells; i++) {
    double dens = env.dens_data[i];
    double tele = env.tele_data[i];

    double nion = NAVO*dens/env.abar;
    double pres = nion * KBOL * tele;
    double eint = (1/(env.gamma-1)) * pres / dens;

    env.pres_data[i] = pres;
    env.eint_data[i] = eint;
    env.cspd_data[i] = std::sqrt(env.gamma*pres/dens);
  }
}

template <typename E>
void eos_dens_eint(E &env)
{
  const unsigned ncells = env.ncells;
  const double gamma = env.gamma;
  const double *dens_data = env.dens_data;
  const double *eint_data = env.eint_data;
  double *pres_data = env.pres_data;
  double *tele_data = env.tele_data;
  double *cspd_data = env.cspd_data;


  for(unsigned i = 0; i < ncells; i++) {
    double dens = dens_data[i];
    double uint = eint_data[i] * dens;
    double nion = NAVO*dens/env.abar;
    double pres = (gamma-1) * uint;
    double tele = (gamma-1) * uint / (KBOL * nion);

    pres_data[i] = pres;
    tele_data[i] = tele;
    cspd_data[i] = std::sqrt(gamma*pres/dens);
  }
}

template <typename E>
double area(const E& env, double x)
{
  switch (env.geom) {
  case GeomType::CARTESIAN:
    return 1.0;
  case GeomType::SPHERICAL:
    return 4*PI*x*x;
  default:
    ERR("Shouldn't be here");
  }
}


template <typename E>
double volume(const E& env, double xm, double xp)
{
  switch (env.geom) {
  case GeomType::CARTESIAN:
    return xp-xm;
  case GeomType::SPHERICAL:
    return 4*PI/3 * (xp*xp*xp - xm*xm*xm);
  default:
    ERR("Shouldn't be here");
  }
}

template <typename E>
void update_state(E &env, double dt)
{
  double * dens_data = env.dens_data;
  double * eint_data = env.eint_data;
  const double * pres_data = env.pres_data;
  const double * cmas_data = env.cmas_data;
  const double * pvis_data = env.pvis_data;
  const double * coords = env.coords;

  unsigned ncells = env.ncells;
  for(unsigned i = 0; i < ncells; i++) {
    const double dens_old = dens_data[i];
    // Compute new density, needs to be modified for alternate
    // cordinate systems...
    dens_data[i] = cmas_data[i] / volume(env,coords[i],coords[i+1]);
    
    // Update internal energy:
    eint_data[i] += - (pres_data[i]+pvis_data[i]) * (1/dens_data[i] - 1/dens_old);
  }

  eos_dens_eint(env);
}


template <typename E>
void update_coords(E &env, double dt)
{
  const double *velx_data = env.velx_data;
  double *coords = env.coords;

  unsigned nverts = env.nverts;
  for(unsigned i = 0; i < nverts; i++)
    coords[i] += dt * velx_data[i];
}


template <typename E>
void update_velocity(E &env, double dt)
{
  unsigned nverts = env.nverts;
  unsigned ncells = env.ncells;
  double *velx_data = env.velx_data;
  const double *vmas_data = env.vmas_data;
  const double *pres_data = env.pres_data;
  const double *pvis_data = env.pvis_data;
  const double *coords = env.coords;

  // i-min boundary:
  switch(env.hydbc_left) {
  case HydroBc::FIXED:
    velx_data[0] = 0;
    break;
  case HydroBc::FREE:
    velx_data[0] -= dt/vmas_data[0] * area(env,coords[0]) * 
      (pres_data[0] + pvis_data[0]);
    break;    
  default:
    ERR("Shouldn't be here");
  }
  
  // Physical Cells:
  for(unsigned i = 1; i < nverts-1; i++) {
    double pm = pres_data[i-1] + pvis_data[i-1];
    double pp = pres_data[i] + pvis_data[i];
    velx_data[i] -= dt/vmas_data[i] * area(env,coords[i]) * (pp - pm);
  }
	
  // i-max boundary:
  switch(env.hydbc_right) {
  case HydroBc::FIXED:
    velx_data[nverts-1] = 0;
    break;
  case HydroBc::FREE:
    velx_data[nverts-1] += dt/vmas_data[nverts-1] * area(env,coords[nverts-1]) * 
      (pres_data[ncells-1] + pvis_data[ncells-1]);
    break;    
  default:
    ERR("Shouldn't be here");
  }
}

template <typename E>
double compute_dt(double dtold, const E& env)
{
  double dt;
  const double growth = env.growth;
  const double *coords = env.coords;
  const double *cspd_data = env.cspd_data;

  unsigned ncells = env.ncells;

  dt = dtold *growth;
  for(unsigned i = 0; i < ncells; i++) {
    dt = std::min(dt, env.cfl * (coords[i+1]-coords[i])/cspd_data[i]);
  }

  return dt;
}

template <typename E>
void compute_masses(const E& env)
{
  unsigned ncells = env.ncells;
  unsigned nverts = env.nverts;
  double *coords = env.coords;

  // Set cell/vertex masses:
  double *cmas_data = env.cmas_data;
  double *vmas_data = env.vmas_data;
  double *dens_data = env.dens_data;

  double xm, xc, xp;

  for(unsigned i = 0; i < ncells; i++) {
    xm = coords[i];
    xp = coords[i+1];
    cmas_data[i] = dens_data[i] * volume(env, xm,xp);
  }

  xm = coords[0];
  xp = coords[1];
  xc = 0.5*(xm+xp);
  vmas_data[0] = dens_data[0] * volume(env, xm,xc);

  for(unsigned i = 1; i < nverts-1; i++) {
    xm = coords[i-1];
    xp = coords[i];
    xc = 0.5 * (xm+xp);
    vmas_data[i] = dens_data[i-1] * volume(env,xc,xp);
    
    xm = coords[i];
    xp = coords[i+1];
    xc = 0.5 * (xm+xp);
    vmas_data[i] += dens_data[i] * volume(env,xm,xc);
  }
  
  xm = coords[nverts-2];
  xp = coords[nverts-1];
  xc = 0.5*(xm+xp);
  vmas_data[nverts-1] = dens_data[ncells-1] * volume(env, xc,xp);
}


template <typename E>
void compute_pvis(E& env)
{
  double *pvis_data = env.pvis_data;
  const double *dens_data = env.dens_data;
  const double *velx_data = env.velx_data;
  const double *cspd_data = env.cspd_data;
  
  unsigned ncells = env.ncells;
  for(unsigned i = 0; i < ncells; i++) {
    double dv = velx_data[i+1] - velx_data[i];
    if(dv < 0.0) {
      double dvmag = std::abs(dv);
      double dens = dens_data[i];
      double cspd = cspd_data[i];
      double pvis = dens * dvmag * (env.clin * cspd + env.cquad*dvmag);
      pvis_data[i] = pvis;
    } else {
      pvis_data[i] = 0.0;
    }
  }
}

// void plot(unsigned cycle, Mesh& mesh) {
//       std::cout << "Plotting...\n";
//       sstr << "plot-" 
// 	   << std::setw(5) 
// 	   << std::setfill('0')
// 	   << step << ".vtk";
//       mesh.write_file(sstr.str().c_str());
//       sstr.str("");
// }

// ****************
// *              *
// *     MAIN     *
// *              *
// ****************
int main()
{

  // ******************
  // *   INITIALIZE   *
  // ******************

  // Make object to store mesh and all simulation data:
  // SimData sim(100, GeomType::SPHERICAL);

  Mesh mesh(Input::ncells);
  Env_Top env(mesh);

  // Set initial condition
  init_state(env);
  compute_masses(env);

  // ************
  // *   LOOP   *
  // ************

  // Do the simulation...
  unsigned cycle = 0;
  double dt = Input::dtinit;
  double time = 0.0;
  double plot_last_time = - std::numeric_limits<double>::max();
  std::stringstream sstr;

  while (time < Input::tmax && cycle < Input::nmax) {
    std::cout << std::setw(10) << cycle 
    	      << std::setw(15) << dt
    	      << std::setw(15) << time
    	      << std::endl;

    if(time - plot_last_time > env.plot_time_freq) {
      std::cout << "Plotting at cycle " << cycle << ", time " << time << std::endl;
      plot_last_time = time;
      sstr << "plot-" 
	   << std::setw(5) 
	   << std::setfill('0')
	   << cycle << ".h5m";
      mesh.write_file(sstr.str().c_str());
      sstr.str("");
    }

    compute_pvis(env);
    update_velocity(env, dt);
    update_coords(env, dt);
    update_state(env,dt);
    
    // Compute new dt and new time...
    time += dt;
    dt = compute_dt(dt,env);
    cycle += 1;
  }

  std::cout << "Plotting at cycle " << cycle << ", time " << time << std::endl;
  sstr << "plot-" 
       << std::setw(5) 
       << std::setfill('0')
       << cycle << ".h5m";
  mesh.write_file(sstr.str().c_str());
  sstr.str("");

  return 0;
}
