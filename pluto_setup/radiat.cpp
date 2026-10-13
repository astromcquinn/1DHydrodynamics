////////////////////////////////////////////////////////////////////////
/*!
  Local override of Src/Cooling/Tabulated/radiat.cpp.

  Differences with the distributed version:
   - electron density for a fully ionized H+He+Z plasma
       n_e = n_H (1 + 2 f_He + A_Z f_Z / 2)
     (the distributed file omits the helium electrons);
   - a uniform, constant volumetric heating rate inputParam[HEAT_CGS]
     (erg cm^-3 s^-1) is added to the right hand side;
   - if inputParam[COOL_MASK] = 1 the cooling rate is multiplied by the
     passive tracer (1 = cloud material, 0 = background), i.e. only the
     cloud cools ("cooling turned off outside the cloud");
   - cooling is switched off below run_conf.minCoolingTemp with the
     same smooth cutoff used by PLUTO; heating is switched off below
     the floor too, so the floor is a fixed point.
  The cooling table cooltable.dat (T [K], Lambda [erg cm^3 s^-1],
  normalised to n_e n_H) is read by Src/Cooling/Tabulated/read_tabfile.cpp
*/
////////////////////////////////////////////////////////////////////////
#include "gpluto.hpp"

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

  double cutoff = 1.0 - 1.0/cosh(pow(T/run_conf.minCoolingTemp, 12));
  double heat   = run_conf.inputParam[HEAT_CGS];
  if (run_conf.inputParam[COOL_MASK] > 0.5) {
    double trc = MIN(MAX(v[TRC], 0.0), 1.0);
    lam *= trc;
  }
  rhs[RHOE] = (heat - ne*nH*lam)*cutoff*E_cost;
}
