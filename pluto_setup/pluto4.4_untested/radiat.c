/* ///////////////////////////////////////////////////////////////////// */
/*!
  Local override of Src/Cooling/Tabulated/radiat.c for PLUTO 4.4.
  UNTESTED port of ../radiat.cpp (gPLUTO).  Same physics:
   - fully ionized n_e = n_H (1 + 2 f_He + A_Z f_Z/2);
   - Lambda(T) read from cooltable.dat (T, Lambda[n_e n_H norm.]),
     linear interpolation in T;
   - optional cooling mask by the passive tracer (COOL_MASK = 1);
   - uniform volumetric heating HEAT_CGS [erg cm^-3 s^-1];
   - smooth cutoff below g_minCoolingTemp for both terms.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

void Radiat (double *v, double *rhs)
{
  static int ntab;
  static double *L_tab, *T_tab, E_cost;
  int    klo, khi, kmid;
  double mu, T, prs, dT, lam, nH, ne, cutoff, trc;
  FILE *fcool;

  if (T_tab == NULL){
    print (" > Reading table from disk...\n");
    fcool = fopen("cooltable.dat","r");
    if (fcool == NULL){
      print ("! Radiat: cooltable.dat could not be found.\n");
      QUIT_PLUTO(1);
    }
    L_tab = ARRAY_1D(20000, double);
    T_tab = ARRAY_1D(20000, double);
    ntab = 0;
    while (fscanf(fcool, "%lf  %lf\n", T_tab + ntab, L_tab + ntab) != EOF) ntab++;
    fclose(fcool);
    E_cost = UNIT_LENGTH/UNIT_DENSITY/pow(UNIT_VELOCITY, 3.0);
  }

  prs = v[RHOE]*(g_gamma - 1.0);
  if (prs < 0.0) {
    prs     = g_smallPressure;
    v[RHOE] = prs/(g_gamma - 1.0);
  }
  mu = MeanMolecularWeight(v);
  T  = prs/v[RHO]*KELVIN*mu;

  lam = 0.0;
  if (T >= T_tab[0] && T <= T_tab[ntab-1]) {
    klo = 0; khi = ntab - 1;
    while (klo != (khi - 1)) {
      kmid = (klo + khi)/2;
      if (T <= T_tab[kmid]) khi = kmid;
      else                  klo = kmid;
    }
    dT  = T_tab[khi] - T_tab[klo];
    lam = (L_tab[klo]*(T_tab[khi] - T) + L_tab[khi]*(T - T_tab[klo]))/dT;
  }
  if (g_inputParam[COOL_MASK] > 0.5) {
    trc = MIN(MAX(v[TRC], 0.0), 1.0);
    lam *= trc;
  }

  nH = UNIT_DENSITY/CONST_amu*H_MASS_FRAC/CONST_AH*v[RHO];
  ne = nH*(1.0 + 2.0*FRAC_He + 0.5*CONST_AZ*FRAC_Z);

  cutoff    = 1.0 - 1.0/cosh(pow(T/g_minCoolingTemp, 12));
  rhs[RHOE] = (g_inputParam[HEAT_CGS] - ne*nH*lam)*cutoff*E_cost;
}
