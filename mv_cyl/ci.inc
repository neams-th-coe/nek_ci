#include <cmath>

static int ciMode = 0;

#define EPS 5e-2

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("NUMBER OF SCALARS", std::string("1"));
  options.setArgs("END TIME", std::string("0.1"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("FLUID PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("BDF ORDER", "3");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-09"));
  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", std::string("1e-10"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));

  options.setArgs("GEOM SOLVER", "NONE");
  options.setArgs("MOVING MESH", "TRUE");

  platform->par->set("casedata", "p_rot", 0.0);

  if (ciMode == 2) {
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  if (ciMode == 3) {
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("GEOM INITIAL GUESS", "PROJECTION");
    options.setArgs("GEOM RESIDUAL PROJECTION VECTORS", std::string("1"));
    options.setArgs("GEOM SOLVER", "PCG+BLOCK");
    options.setArgs("GEOM SOLVER TOLERANCE", "1e-6");
  }
  if (ciMode == 4) {
    nekrsAbort(comm, EXIT_FAILURE, "%s\n", "ciMode=4 not defined");
  }
  if (ciMode == 5) {
    platform->par->set("casedata", "p_rot", 0.5);
  }
  if (ciMode == 6 || ciMode == 8) {
    options.setArgs("GEOM SOLVER", "PCG+BLOCK");
    options.setArgs("GEOM SOLVER TOLERANCE", "1e-6");
    options.setArgs("GEOM INITIAL GUESS", "PROJECTION");
    options.setArgs("GEOM RESIDUAL PROJECTION VECTORS", std::string("1"));
    platform->par->set("casedata", "p_rot", 0.5);
  }
  if (ciMode == 7 || ciMode == 8) {
    options.setArgs("CVODE", "TRUE");
    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-08");
    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-06");
    options.setArgs("SCALAR00 SOLVER", "CVODE");
  }
}

void ciTestErrors()
{
  if(!nrs->lastStep) 
    return;

  auto err = nek::ptr<double>("errors");

  constexpr dfloat expectedP0thError = 0.2973762E-05;
  constexpr dfloat expectedDpdtError = 0.5820622E-06;
  constexpr int expectedVIter = 2;
  constexpr int expectedPIter = 5;

  auto errVal = std::abs(err[2] - expectedP0thError) / expectedP0thError;
  CiEvalTest("p0th (LinfErr: " + to_string_f(errVal) + ")", errVal < EPS);

  errVal = std::abs(err[3] - expectedDpdtError) / expectedDpdtError;
  CiEvalTest("dpdt (LinfErr: " + to_string_f(errVal) + ")", errVal < EPS);

  errVal = std::abs(expectedVIter - nrs->fluid->ellipticSolver[0]->Niter());
  CiEvalTest("velIter (LinfErr: " + to_string_f(errVal) + ")", errVal < 2);

  errVal = std::abs(expectedPIter - nrs->fluid->ellipticSolverP->Niter());
  CiEvalTest("pIter (LinfErr: " + to_string_f(errVal) + ")", errVal < 2);
}
