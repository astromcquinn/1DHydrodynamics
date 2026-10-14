////////////////////////////////////////////////////////////////////////
/*!
  1D cooling cloud / cooling medium with a small sinusoidal perturbation.

  Units: length = 1 kpc, velocity = 100 km/s (time unit = 9.78 Myr),
         density = m_p (code density = rho / m_p  [cm^-3]).

  Two configurations, selected by L_CLOUD:

  L_CLOUD > 0 : "tophat" cloud of length L_CLOUD centered in the box, with
                density CHI*RHO_BG and temperature T_BG/CHI (pressure
                equilibrium with the hot background).  The edges are
                smoothed with a tanh of width EDGE_WIDTH (0 = sharp).
  L_CLOUD <= 0: "infinite cloud": the whole (periodic) box is filled with
                gas of density CHI*RHO_BG and temperature T_BG/CHI.

  Inside the cloud a sinusoidal perturbation of type PERT_TYPE is added, with
  c(x) = cos(2 pi (x - x_c)/LAMBDA + PHASE), s(x) = sin(...):
    0  density, isobaric:     delta rho/rho = AMP c,  delta p = 0
    1  temperature, isochoric: delta T/T = delta p/p = AMP c,  delta rho = 0
    2  standing sound wave:   delta rho/rho = AMP c,  delta p/p = gamma AMP c
    3  velocity:              v = AMP c_s s  (c_s of the unperturbed cloud gas)
    4  metallicity:           delta Z/Z = AMP c at uniform rho and p
  Tracer 2 carries the metallicity Z (solar units, ZMET in the cloud); it is
  used by radiat.cpp when the 2D table stores primordial and metal rates
  separately.

  Thermal model of the background (two options, can be combined):
   COOL_MASK = 1 : cooling is multiplied by the passive tracer, so only cloud
                   material cools (cooling "off outside the cloud").
   HEAT_CGS > 0  : see below.
  Cooling model (COOL_MODE, see radiat.cpp): 0 = 1D table Lambda(T) with a
  floor at T_FLOOR; 1 = 2D table of cooling and heating rates vs (n_H, T),
  e.g. photoionisation equilibrium in a UV background.
  Heating: a uniform, time-independent volumetric heating rate HEAT_CGS
  (erg cm^-3 s^-1) is added in radiat.cpp.  Set it to the cooling rate of
  the hot background so that the background is in thermal balance, or to 0
  for pure cooling.

  The passive tracer is 1 in the cloud and 0 in the background.
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

static void SetUnits(RunConfig &run_config)
{
  double *p = run_config.inputParam;
  run_config.unit_density  = CONST_mp;
  run_config.unit_length   = 1.e3*CONST_pc;
  run_config.unit_velocity = 1.e7;
  run_config.minCoolingTemp = p[T_FLOOR];
  run_config.maxCoolingRate = p[MAX_COOL_RATE];
  run_config.small_prs = 1.e-12;
}

////////////////////////////////////////////////////////////////////////
void Init (double *v, double x1, double x2, double x3, RunConfig &run_config)
{
  double *p = run_config.inputParam;
  SetUnits(run_config);

  double mu  = MeanMolecularWeight(v);
  double xc  = 0.5*(run_config.domBeg[IDIR] + run_config.domEnd[IDIR]);
  double rho_bg = p[RHO_BG];
  double prs = rho_bg*p[T_BG]/(Kelvin(run_config)*mu);   // uniform pressure

  // cloud profile f in [0,1]
  double f;
  if (p[L_CLOUD] <= 0.0) {
    f = 1.0;
  } else {
    double d = fabs(x1 - xc) - 0.5*p[L_CLOUD];
    if (p[EDGE_WIDTH] > 0.0) f = 0.5*(1.0 - tanh(d/p[EDGE_WIDTH]));
    else                     f = (d < 0.0 ? 1.0 : 0.0);
  }

  double phase = 2.0*CONST_PI*(x1 - xc)/p[LAMBDA] + p[PHASE];
  double c = cos(phase), s = sin(phase);
  double A = p[AMP];
  int    type = (int)(p[PERT_TYPE] + 0.5);
  double g = run_config.gamma;

  double drho = 0.0, dprs = 0.0, vel = 0.0, dZ = 0.0;
  if      (type == 0) drho = A*c;
  else if (type == 1) dprs = A*c;
  else if (type == 2) { drho = A*c; dprs = g*A*c; }
  else if (type == 3) vel = A*sqrt(g*prs/(p[CHI]*rho_bg))*s;
  else if (type == 4) dZ = A*c;

  double rho_c = p[CHI]*rho_bg*(1.0 + drho);

  v[RHO] = rho_bg + f*(rho_c - rho_bg);
  v[VX1] = f*vel;
  v[VX2] = v[VX3] = 0.0;
  v[PRS] = prs*(1.0 + f*dprs);
  v[TRC]     = f;
  v[TRC + 1] = p[ZMET]*(1.0 + f*dZ);
}

////////////////////////////////////////////////////////////////////////
void InitDomain (Data *data, Grid *grid)
{
}

////////////////////////////////////////////////////////////////////////
void Analysis (Data *data, Grid *grid)
{
}

////////////////////////////////////////////////////////////////////////
void UserDefBoundary (const Data *data, RBox *box, int side, Grid *grid)
{
}
