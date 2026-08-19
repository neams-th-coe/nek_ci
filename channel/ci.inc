#include <math.h>

static int ciMode = 0;

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("FLUID VISCOSITY", std::string("0.1"));
  options.setArgs("FLUID DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.1"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("FLUID VELOCITY SOLVER", "CG+BLOCK");
  options.setArgs("BDF ORDER", "3");
  options.setArgs("ADVECTION", "FALSE");
  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", std::string("1e-08"));
  options.setArgs("FLUID STRESSFORMULATION", "TRUE");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE");

  if (ciMode == 1) {
    platform->par->set("casedata", "p_rot", 0.0);
  }
  if (ciMode == 2) {
    platform->par->set("casedata", "p_rot", 0.7853981633974483);
  }
}

