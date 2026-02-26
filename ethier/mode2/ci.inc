#include <cstdlib>
#include <math.h>
#include <limits>
#include <vector>
#include <numeric>
#include <array>
#include <fstream>
#include <random>

#include "pointInterpolation.hpp"
#include "randomVector.hpp"
#include "LVector.hpp"

double numCVODERHS = 0;

static int ciMode = 0;

#define EPS 1e-1

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("9"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("VISCOSITY", std::string("0.01"));
  options.setArgs("DENSITY", std::string("1"));
  options.setArgs("NUMBER OF SCALARS", std::string("2"));
  options.setArgs("SCALAR00 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR00 DENSITY", std::string("1"));
  options.setArgs("SCALAR01 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR01 DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.06"));
  options.setArgs("DT", std::string("2e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("VELOCITY BLOCK SOLVER", "FALSE");
  options.setArgs("SCALAR INITIAL GUESS", "PREVIOUS");
  options.setArgs("VELOCITY INITIAL GUESS", "PREVIOUS");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");

  options.setArgs("SCALAR00 SOLVER", "PCG+COMBINED");
  options.setArgs("SCALAR01 SOLVER", "PCG");

  if (ciMode == 2 || ciMode == 19) {
    options.setArgs("VELOCITY BLOCK SOLVER", "TRUE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }

  options.setArgs("BDF ORDER", "3");
  options.setArgs("VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  const auto pTol = (sizeof(dfloat) == sizeof(double)) 
                    ? std::string("1e-10") : std::string("1e-05"); 
  options.setArgs("PRESSURE SOLVER TOLERANCE", pTol);
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("VARIABLEPROPERTIES", "FALSE");
}

void ciTestErrors(nrs_t *nrs,
                  double time,
                  int tstep,
                  lpm_t *particles,
                  occa::kernel exactUVWPKernel,
                  occa::kernel RKKernel)
{
  const int rank = platform->comm.mpiRank;
  if (tstep == 1 && ciMode != 7 && ciMode != 13) {
    int NiterP = nrs->pSolver->Niter;

    int expectedNiterP = 7;

    const int pIterErr = abs(NiterP - expectedNiterP);
    if (pIterErr >= 2) {
      if (rank == 0) {
        printf("Expected %d P iteration in first step, but got %d!\n", expectedNiterP, NiterP);
      }
      CIFAIL;
#if 0
      nrsFinalize(nrs);
      exit(platform->exitValue);
#endif
    }
  }

  if(nrs->cds->cvode){
    numCVODERHS += nrs->cds->cvode->numRHSEvals();
    numCVODERHS += nrs->cds->cvode->numLinIters();
  }

  if (!nrs->lastStep) {
    return;
  }

  nek::ocopyToNek(time, tstep);
  nek::userchk();

  double *err = (double *)nek::scPtr(1);

  double vxErr, prErr;
  double s1Err, s2Err;

  int pIterErr;
  int velIterErr;

  int NiterU;

  if (platform->options.compareArgs("VELOCITY BLOCK SOLVER", "TRUE")) {
    NiterU = nrs->uvwSolver->Niter;
  } else if (nrs->uSolver) {
    NiterU = nrs->uSolver->Niter;
  } else {
    NiterU = -1;
  }

  int NiterP = nrs->pSolver ? nrs->pSolver->Niter : -1;

  int NiterS01 = -1;
  if (nrs->cds->solver[0]) {
    NiterS01 = nrs->cds->solver[0]->Niter;
  }
  int NiterS02 = -1;
  if (nrs->cds->solver[1]) {
    NiterS02 = nrs->cds->solver[1]->Niter;
  }
  int s01IterErr, s02IterErr;

  s01IterErr = abs(NiterS01 - 2);
  s02IterErr = abs(NiterS02 - 2);

  velIterErr = abs(NiterU - 10);
  s1Err = abs((err[2] - 6.67E-12) / err[2]);
  s2Err = abs((err[3] - 7.49E-12) / err[3]);
  pIterErr = abs(NiterP - 4);
  vxErr = abs((err[0] - 2.77E-10) / err[0]);
  prErr = abs((err[1] - 6.98E-10) / err[1]);
  // on ci modes 12, 13, confirm that the correct solvers are present

  if (rank == 0) {
    printf("relative error to target: vx=%g pr=%g s1=%g s2=%g velIter=%d pIter=%d s01Iter=%d s02Iter=%d\n",
           vxErr,
           prErr,
           s1Err,
           s2Err,
           velIterErr,
           pIterErr,
           s01IterErr,
           s02IterErr);
  }

  std::vector<double> relErrors = {vxErr, prErr, s1Err, s2Err};
  std::vector<double> absErrors = {abs(err[0]), abs(err[1]), abs(err[2]), abs(err[3])};
  bool passTest = true;
  const double absTol = 5e-11; // values below this threshold are ignored
  for (size_t i = 0; i < relErrors.size(); ++i) {

    const double absError = absErrors[i];
    const double relError = relErrors[i];
    if (absError > absTol) {
      passTest &= relError < EPS;
    }
  }

  passTest &= velIterErr <= 1;
  passTest &= pIterErr <= 2;
  passTest &= s01IterErr <= 1;
  passTest &= s02IterErr <= 1;

  passTest ? (CIPASS) : (CIFAIL);
}
