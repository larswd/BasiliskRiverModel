#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "libs/output_vts.h"
#include "libs/discharge_ml.h"

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
  river = new scalar[nl];
  foreach(){
    zb[] = 0.05*pow(x,4) - x*x + 2. + 0.2*(y + Y0);
    river = river[] = x < 0 ? 1 : 2;
  }
  u.n[top] = neumann(0);
  u.t[top] = dirichlet(0)

  u.n[bottom] = neumann(0);
}

double eta1, eta2;
event inflow(i++){
  eta1 = eta_b(4, top, river, 1);
  eta2 = eta_b(2, top, river, 2);
  h[top] = max ((river[] == 1. ? eta1 : eta2) - zb[], 0.);
  eta[top] = max ((river[] == 1. ? eta1 : eta2) - zb[], 0.) + zb[];
}

event cleanup(t=end){
  delete((scalar*){river});
}