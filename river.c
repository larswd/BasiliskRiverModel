#include "grid/multigrid.h"
#include "saint-venant.h"
#include "discharge.h"

/*
Plotting utility
*/
#if LAYERED
#else
scalar w;
#endif

#include "libs/output_vts.h"

static int level;

/*
  Simulation parameters,
  Q1, Q2 are flux in river 1 and 2 respectively
  slope is the steepness of the decline (in meters per 10 meter)
  xslopefac is a linear xslope-term for steepness of decline in positive x direction
  a4,a2 are coefficients of fourth and second order term
*/
static double Q1        = 2;
static double Q2        = 1;
static double slope     = 1;
static double xslopefac = 1e-9; 
static double a4        = 0.05;
static double a2        = 0.9;


int main(int argc, char * argv[]){
  // Numerical resolution given as command line argument
  if (argc <=1){
    level = 7;
  } else {
    level = atoi(argv[1]);
  }
    
  size(10);
  origin(-L0/2., -L0/2.);
  G = 9.81;
  N = 1 << level;
  run();
  return 0;
}

// We initialize two empty rivers. River scalar acts as identifier. 
scalar river[];
event init(i = 0){
  DT = 1e-2;
  foreach(serial){
    zb[] = a4*pow(x,4) - a2*x*x + 2. + slope*(1 + xslopefac*x)*(y + Y0);
    river[] = x < 0 ? 1 : 2;
    u.x[]  = 0.0;
  }
  u.n[top] = neumann(0);
  u.t[top] = dirichlet(0);

  u.n[bottom] = neumann(0);
}

double eta1, eta2;
event inflow(i++){
  eta1 = eta_b(Q1, top, river, 1);
  eta2 = eta_b(Q2, top, river, 2);
  h[top] = max ((river[] == 1. ? eta1 : eta2) - zb[], 0.);
  eta[top] = max ((river[] == 1. ? eta1 : eta2) - zb[], 0.) + zb[];
}


event volume (i += 10) {
  double volume1 = 0, volume2 = 0;
  foreach(reduction(+:volume1) reduction(+:volume2)) {
    double dv = h[]*sq(Delta);
    if (x < 0) volume1 += dv;
    else volume2 += dv;
  }
  fprintf (stderr, "%g %g %g %g %g\n",
	   t, volume1, volume2, eta1, eta2);
}

event animation (t <= 1; i += 10) {
  static int j = 0;
  char name[100];
  printf("Making plotfile\n");
  sprintf(name, "plots/field_%.6i.vts", j++);
  FILE* fp = fopen(name, "w");
  output_vts_ascii_all_layers(fp, {eta,zb}, N);
  fclose(fp);
}