#include <math.h>
#include <vector>

static int ciMode = 0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", std::string("0"));
  options.setArgs("FLUID VISCOSITY", std::string("1e-2"));
  options.setArgs("FLUID DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.05"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("FLUID VELOCITY SOLVER", "CG+BLOCK");
  options.setArgs("BDF ORDER", std::string("3"));
  options.setArgs("EXT ORDER", std::string("3"));
  options.setArgs("ADVECTION", "FALSE");
  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", std::string("1e-08"));
  options.setArgs("ADVECTION TYPE", "CONVECTIVE");

  if (ciMode == 1) {
    platform->par->set("casedata", "p_rot", std::string("0.0"));
  }
  if (ciMode == 2) {
    platform->par->set("casedata", "p_rot", std::string("0.7853981633974483"));
  }
}

