/* ///////////////////////////////////////////////////////////////////// */
/*!
  1D cooling cloud / cooling medium with a small sinusoidal perturbation.
  PLUTO 4.4 port of ../init.cpp (gPLUTO).  UNTESTED: written without
  access to a PLUTO 4.4 installation -- compile and compare one run
  against the gPLUTO version before using it for science.

  Units: length 1 kpc, velocity 100 km/s (time unit 9.78 Myr),
  density m_p (code density = rho/m_p in cm^-3).
  See ../init.cpp for a description of the parameters.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/* ********************************************************************* */
{
  double mu, xc, rho_bg, prs, f, d, pert, rho_c;

  g_minCoolingTemp = g_inputParam[T_FLOOR];
  g_maxCoolingRate = g_inputParam[MAX_COOL_RATE];

  mu     = MeanMolecularWeight(v);
  xc     = 0.5*(g_domBeg[IDIR] + g_domEnd[IDIR]);
  rho_bg = g_inputParam[RHO_BG];
  prs    = rho_bg*g_inputParam[T_BG]/(KELVIN*mu);

  if (g_inputParam[L_CLOUD] <= 0.0) f = 1.0;
  else {
    d = fabs(x1 - xc) - 0.5*g_inputParam[L_CLOUD];
    if (g_inputParam[EDGE_WIDTH] > 0.0) f = 0.5*(1.0 - tanh(d/g_inputParam[EDGE_WIDTH]));
    else                                f = (d < 0.0 ? 1.0 : 0.0);
  }

  pert  = g_inputParam[AMP]*cos(2.0*CONST_PI*(x1 - xc)/g_inputParam[LAMBDA]
                                + g_inputParam[PHASE]);
  rho_c = g_inputParam[CHI]*rho_bg*(1.0 + pert);

  v[RHO] = rho_bg + f*(rho_c - rho_bg);
  v[VX1] = v[VX2] = v[VX3] = 0.0;
  v[PRS] = prs;
  v[TRC] = f;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
{
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
{
}
