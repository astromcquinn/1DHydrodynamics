////////////////////////////////////////////////////////////////////////
/*!
  Local override of Src/Cooling/Tabulated/radiat.cpp.

  Two cooling modes, selected by inputParam[COOL_MODE]:

  COOL_MODE = 0 : 1D table cooltable.dat (T [K], Lambda [erg cm^3 s^-1]
                  normalised to n_e n_H), read by PLUTO's read_tabfile.cpp,
                  linear interpolation in T.  Electron density for a fully
                  ionised H+He+Z plasma, n_e = n_H (1 + 2 f_He + A_Z f_Z/2)
                  (the distributed file omits the helium electrons).
  COOL_MODE = 1 : 2D table cooltable2d.dat of cooling AND heating rates,
                  Lambda/n_H^2 [erg cm^3 s^-1] as functions of (log n_H, log T)
                  (e.g. Cloudy photoionisation equilibrium in a UV background,
                  written by python/cooling_models.py).  Bilinear interpolation
                  of log10(rate) in (log n_H, log T); net loss
                  n_H^2 (Lambda_cool - Lambda_heat).
                  File format:  nn nT ncomp / nn values of log n_H / nT values
                  of log T / ncomp blocks of nn rows of log10(rate) (-99 for 0).
                  ncomp = 2: total cooling, total heating (fixed Z);
                  ncomp = 4: primordial cooling, solar-metal cooling,
                  primordial heating, solar-metal heating, combined with the
                  local metallicity Z carried by tracer 2:
                  rate = primordial + Z * metal.

  Common to both:
   - a uniform constant volumetric heating inputParam[HEAT_CGS]
     [erg cm^-3 s^-1] is added;
   - COOL_MASK = 1: the tabulated net rate is multiplied by the passive tracer
     (1 = cloud material, 0 = background), i.e. only the cloud cools/is heated
     ("cooling turned off outside the cloud"); mixed cells at the cloud edge
     cool partially.
     COOL_MASK = 2: sharp version, cool only cells with tracer > 0.5, which
     limits cooling of hot gas numerically mixed into the edge cells;
   - everything is switched off below run_conf.minCoolingTemp with PLUTO's
     smooth cutoff (set T_FLOOR well below 1e4 K when using a UV background).
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

static int     g2_nn = 0, g2_nT = 0, g2_nc = 2;
static double *g2_logn = NULL, *g2_logT = NULL, *g2_tab[4] = {NULL, NULL, NULL, NULL};

static void Read2DTable(void)
{
  char line[256];
  FILE *fp = fopen("cooltable2d.dat", "r");
  if (fp == NULL) QuitPLUTO("! Radiat(): cooltable2d.dat could not be found", 1);
  if (fgets(line, sizeof(line), fp) == NULL) QuitPLUTO("! Radiat(): empty cooltable2d.dat", 1);
  int nread = sscanf(line, "%d %d %d", &g2_nn, &g2_nT, &g2_nc);
  if (nread < 2) QuitPLUTO("! Radiat(): bad cooltable2d.dat header", 1);
  if (nread == 2) g2_nc = 2;
  if (g2_nc != 2 && g2_nc != 4) QuitPLUTO("! Radiat(): ncomp must be 2 or 4", 1);
  g2_logn = (double *) malloc(g2_nn*sizeof(double));
  g2_logT = (double *) malloc(g2_nT*sizeof(double));
  int ok = 1;
  for (int i = 0; i < g2_nn; i++) ok &= fscanf(fp, "%lf", g2_logn + i) == 1;
  for (int j = 0; j < g2_nT; j++) ok &= fscanf(fp, "%lf", g2_logT + j) == 1;
  for (int c = 0; c < g2_nc; c++) {
    g2_tab[c] = (double *) malloc(g2_nn*g2_nT*sizeof(double));
    for (int k = 0; k < g2_nn*g2_nT; k++) ok &= fscanf(fp, "%lf", g2_tab[c] + k) == 1;
  }
  fclose(fp);
  if (!ok) QuitPLUTO("! Radiat(): cooltable2d.dat is truncated", 1);
  printLog (" > Read 2D cooling table (%d x %d, %d components)\n", g2_nn, g2_nT, g2_nc);
}

static int Locate(const double *a, int n, double x)
{
  int lo = 0, hi = n - 1;
  if (x <= a[0])   return 0;
  if (x >= a[n-1]) return n - 2;
  while (hi - lo > 1) {
    int mid = (lo + hi)/2;
    if (x < a[mid]) hi = mid; else lo = mid;
  }
  return lo;
}

static double Bilinear(const double *tab, double x, double y)
{
  int i = Locate(g2_logn, g2_nn, x);
  int j = Locate(g2_logT, g2_nT, y);
  double tx = (x - g2_logn[i])/(g2_logn[i+1] - g2_logn[i]);
  double ty = (y - g2_logT[j])/(g2_logT[j+1] - g2_logT[j]);
  tx = MIN(MAX(tx, 0.0), 1.0);
  ty = MIN(MAX(ty, 0.0), 1.0);
  return   (1.0 - tx)*(1.0 - ty)*tab[i*g2_nT + j]     + tx*(1.0 - ty)*tab[(i+1)*g2_nT + j]
         + (1.0 - tx)*ty        *tab[i*g2_nT + j + 1] + tx*ty        *tab[(i+1)*g2_nT + j + 1];
}

void Radiat (double *v, double *rhs, RunConfig& run_conf)
{
  int    klo, khi, kmid;
  int    ntab = run_conf.ntab;
  double unit_density  = run_conf.unit_density;
  double unit_length   = run_conf.unit_length;
  double unit_velocity = run_conf.unit_velocity;
  double g_gamma = run_conf.gamma;

  double nH = unit_density/CONST_amu*H_MASS_FRAC/CONST_AH*v[RHO];
  double ne = nH*(1.0 + 2.0*FRAC_He + 0.5*CONST_AZ*FRAC_Z);

  double E_cost = unit_length/unit_density/POW3(unit_velocity);
  double prs = v[RHOE]*(g_gamma - 1.0);
  if (prs < 0.0) {
    v[RHOE] = run_conf.small_prs/(g_gamma - 1.0);
    prs = run_conf.small_prs;
  }

  double mu = (CONST_AH + FRAC_He*CONST_AHe + FRAC_Z*CONST_AZ) /
              (2.0 + 3.0*FRAC_He + FRAC_Z*(1.0 + CONST_AZ*0.5));
  double T  = prs/v[RHO]*Kelvin(run_conf)*mu;

  double loss = 0.0;   /* tabulated net loss per unit volume [erg cm^-3 s^-1] */

  if (run_conf.inputParam[COOL_MODE] < 0.5) {
    double lam = 0.0;
    if (T >= run_conf.T_tab[0] && T <= run_conf.T_tab[ntab-1]) {
      klo = 0;
      khi = ntab - 1;
      while (klo != (khi - 1)) {
        kmid = (klo + khi)/2;
        T <= run_conf.T_tab[kmid] ? (khi = kmid) : (klo = kmid);
      }
      double dT = run_conf.T_tab[khi] - run_conf.T_tab[klo];
      lam = (  run_conf.L_tab[klo]*(run_conf.T_tab[khi] - T)
             + run_conf.L_tab[khi]*(T - run_conf.T_tab[klo]))/dT;
    }
    loss = ne*nH*lam;
  } else {
    if (g2_tab[0] == NULL) Read2DTable();
    double x = log10(nH), y = log10(T);
    if (g2_nc == 2) {
      loss = nH*nH*(pow(10.0, Bilinear(g2_tab[0], x, y)) - pow(10.0, Bilinear(g2_tab[1], x, y)));
    } else {
      double Z = MAX(v[TRC + 1], 0.0);
      double cool = pow(10.0, Bilinear(g2_tab[0], x, y)) + Z*pow(10.0, Bilinear(g2_tab[1], x, y));
      double heat = pow(10.0, Bilinear(g2_tab[2], x, y)) + Z*pow(10.0, Bilinear(g2_tab[3], x, y));
      loss = nH*nH*(cool - heat);
    }
  }

  if (run_conf.inputParam[COOL_MASK] > 1.5) {          /* sharp mask */
    if (v[TRC] < 0.5) loss = 0.0;
  } else if (run_conf.inputParam[COOL_MASK] > 0.5) {   /* proportional mask */
    double trc = MIN(MAX(v[TRC], 0.0), 1.0);
    loss *= trc;
  }

  double cutoff = 1.0 - 1.0/cosh(pow(T/run_conf.minCoolingTemp, 12));
  double heat   = run_conf.inputParam[HEAT_CGS];
  rhs[RHOE] = (heat - loss)*cutoff*E_cost;
}
