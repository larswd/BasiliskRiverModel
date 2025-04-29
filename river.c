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
    zb[] = 0.05*pow(x,4) - x*x + 2. + 0.2*(y + Y0);
    river[] = x < 0 ? 1 : 2;
  }
  u.n[top] = neumann(0);
  u.t[top] = dirichlet(0);

  u.n[bottom] = neumann(0);
}

double eta1, eta2;
event inflow(i++){
  eta1 = eta_b(4, top, river, 1);
  eta2 = eta_b(2, top, river, 2);
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

event animation (t <= 1.2; i += 10) {
  static int j = 0;
  char name[100];
  printf("Making plotfile\n");
  sprintf(name, "plots/field_%.6i.vts", j++);
  FILE* fp = fopen(name, "w");
  output_vts_ascii_all_layers(fp, {eta,zb}, N);
  fclose(fp);
}