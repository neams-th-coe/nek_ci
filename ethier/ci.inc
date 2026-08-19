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

std::fexcept_t flag;

std::vector<dfloat> xp0;
std::vector<dfloat> yp0;
std::vector<dfloat> zp0;

// initial coordinates, based on global particle id (duplicated on each proc)
std::vector<dfloat> xpg0;
std::vector<dfloat> ypg0;
std::vector<dfloat> zpg0;

std::vector<dfloat> pid0;


// solve ODE system dx_i/dt = u_i using 4th order Runge-Kutta
occa::memory solveODE(nrs_t *nrs,
                      dfloat t0,
                      dfloat tf,
                      dfloat dt,
                      occa::memory o_y0,
                      occa::kernel exactUVWPKernel,
                      occa::kernel RKKernel)
{
  auto rhs = [exactUVWPKernel](dfloat t, occa::memory o_y, occa::memory o_ydot) {
    auto n = o_y.length();
    n /= 3;

    if (n == 0) {
      return;
    }

    auto o_xcoord = o_y + 0 * n;
    auto o_ycoord = o_y + 1 * n;
    auto o_zcoord = o_y + 2 * n;

    static occa::memory o_Pscratch;
    if (o_Pscratch.length() != n) {
      o_Pscratch = platform->device.malloc<dfloat>(n);
    }

    exactUVWPKernel(n, t, o_xcoord, o_ycoord, o_zcoord, n, o_Pscratch, o_ydot);
  };

  std::array<dfloat, 4> rk = {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
  auto o_rk = platform->device.malloc<dfloat>(rk.size(), rk.data());
  auto o_k1 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_k2 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_k3 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_k4 = platform->device.malloc<dfloat>(o_y0.length());
  auto o_y = platform->device.malloc<dfloat>(o_y0.length());
  auto o_tmp = platform->device.malloc<dfloat>(o_y0.length());

  const auto n = o_y0.length();

  if (n == 0) {
    return o_y;
  }

  o_y.copyFrom(o_y0);

  dfloat tn = t0;
  while (tn < tf) {
    rhs(tn, o_y, o_k1);

    platform->linAlg->axpbyz(n, 0.5 * dt, o_k1, 1.0, o_y, o_tmp);
    rhs(tn + 0.5 * dt, o_tmp, o_k2);

    platform->linAlg->axpbyz(n, 0.5 * dt, o_k2, 1.0, o_y, o_tmp);
    rhs(tn + 0.5 * dt, o_tmp, o_k3);

    platform->linAlg->axpbyz(n, dt, o_k3, 1.0, o_y, o_tmp);
    rhs(tn + dt, o_tmp, o_k4);

    RKKernel(n, dt, o_rk, o_k1, o_k2, o_k3, o_k4, o_y);

    if ((tn + dt) > tf) {
      dt = tf - tn;
    }
    tn += dt;
  }

  o_k1.free();
  o_k2.free();
  o_k3.free();
  o_k4.free();
  o_tmp.free();
  return o_y;
}

static int ciMode = 0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  std::feclearexcept(FE_ALL_EXCEPT);
  std::fegetexceptflag(&flag, FE_ALL_EXCEPT);
  std::fesetexceptflag(&flag, FE_ALL_EXCEPT);  

  options.setArgs("POLYNOMIAL DEGREE", std::string("9"));
  if (platform->options.getArgs("CASENAME") == "ethierRefine") {
    options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  }

  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("FLUID VISCOSITY", std::string("0.01"));
  options.setArgs("FLUID DENSITY", std::string("1"));
  options.setArgs("NUMBER OF SCALARS", std::string("2"));
  options.setArgs("SCALAR00 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR00 DENSITY", std::string("1"));
  options.setArgs("SCALAR01 DIFFUSIVITY", std::string("0.01"));
  options.setArgs("SCALAR01 DENSITY", std::string("1"));
  options.setArgs("END TIME", std::string("0.06"));
  options.setArgs("DT", std::string("2e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  options.setArgs("FLUID PRESSURE INITIAL GUESS", "PREVIOUS");

  options.setArgs("SCALAR INITIAL GUESS", "PREVIOUS");
  options.setArgs("FLUID VELOCITY INITIAL GUESS", "PREVIOUS");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");

  options.setArgs("FLUID VELOCITY SOLVER", "PCG");
  options.setArgs("SCALAR00 SOLVER", "PCG+COMBINED");
  options.setArgs("SCALAR01 SOLVER", "PCG");

  options.setArgs("BDF ORDER", "3");
  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  const auto pTol = (sizeof(dfloat) == sizeof(double)) 
                    ? std::string("1e-10") : std::string("1e-06"); 
  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", pTol);
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-12"));

  if (ciMode == 2 || ciMode == 19 || ciMode == 22) {
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");

    options.setArgs("SUBCYCLING STEPS", std::string("1"));

    if (ciMode == 22) {
      options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-10"));
      options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-10"));
      options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-10"));

      nekrsCheck(sizeof(dfloat) != sizeof(float),
                 platform->comm.mpiComm(),
                 EXIT_FAILURE,
                 "%s\n",
                 "CI test requires sizeof(dfloat) == sizeof(float)");
    } 
  }
  if (ciMode == 3) {
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("FLUID VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("FLUID PRESSURE PRECONDITIONER", "MULTIGRID+SEMFEM");
    options.setArgs("FLUID PRESSURE MULTIGRID COARSE SOLVER LOCATION", "DEVICE");
    if (platform->serial()) { 
      options.setArgs("FLUID PRESSURE MULTIGRID COARSE SOLVER LOCATION", "CPU");
    } 
  }
  if (ciMode == 4) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("FLUID VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 5) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");
    options.setArgs("SUBCYCLING STEPS", std::string("0"));
    std::string casename = platform->options.getArgs("CASENAME");
    nekrsCheck(casename != "mv_ethier",
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "CI mode %d requires mv_ethier.par\n",
             ciMode);
  }
  if (ciMode == 6) {
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    std::string casename = platform->options.getArgs("CASENAME");
    nekrsCheck(casename != "mv_ethier",
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "CI mode %d requires mv_ethier.par\n",
             ciMode);
  }
  if (ciMode == 7) {
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("FLUID VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("FLUID PRESSURE SOLVER MAXIMUM ITERATIONS", "1000");
    options.setArgs("FLUID PRESSURE PRECONDITIONER", "JACOBI");
    options.setArgs("END TIME", std::string("0.012"));
  }
  if (ciMode == 8) {
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
//    options.setArgs("FLUID VELOCITY INITIAL GUESS", "PROJECTION-ACONJ");
    options.setArgs("DT", std::string("1e-5"));
    options.setArgs("END TIME", std::string("0.2"));
    options.setArgs("VARIABLE DT", "TRUE");
    options.setArgs("TARGET CFL", "0.5");

#if 0
    options.setArgs("DT", std::string("2e-3"));
    options.setArgs("TARGET CFL", "0.274");
#endif

  }
  if (ciMode == 9) {
    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 10) {
    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 11) {
    options.setArgs("FLUID PRESSURE MULTIGRID SMOOTHER", "DAMPEDJACOBI,CHEBYSHEV");
    options.setArgs("FLUID PRESSURE MULTIGRID CHEBYSHEV MIN EIGENVALUE BOUND FACTOR", "0.1");
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }
  if (ciMode == 12) {
    options.setArgs("SCALAR00 SOLVER", "NONE");
  }
  if (ciMode == 13) {
    options.setArgs("SCALAR00 SOLVER", "NONE");
    options.setArgs("FLUID VELOCITY SOLVER", "NONE");

    std::string casename = platform->options.getArgs("CASENAME");
    nekrsCheck(casename != "ethierScalar",
             platform->comm.mpiComm(),
             EXIT_FAILURE,
             "CI mode %d requires ethierScalar.par\n",
             ciMode);
  }
  if (ciMode == 14) {
    options.setArgs("FLUID PRESSURE MULTIGRID SMOOTHER", "ASM");
    options.setArgs("FLUID PRESSURE MGSOLVER CYCLE", "VCYCLE+ADDITIVE+OVERLAPCRS");
    options.setArgs("FLUID VELOCITY SOLVER", options.getArgs("FLUID VELOCITY SOLVER") + "+BLOCK");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
    options.setArgs("FLUID PRESSURE INITIAL GUESS", "PROJECTION-ACONJ");
  }

  // ci mode requiring no timestepping
  if (ciMode == 15) {
    options.setArgs("END TIME", std::string("0.0"));
    options.setArgs("NUMBER TIMESTEPS", std::string("0"));
    options.setArgs("SCALAR00 SOLVER", "NONE");
    options.setArgs("SCALAR01 SOLVER", "NONE");
    options.setArgs("FLUID", "FALSE");
  }

  // CVODE, both scalars
  if (ciMode == 16) {

    options.setArgs("CVODE", "TRUE");
    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("SCALAR01 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-10");

    options.setArgs("CVODE GMRES BASIS VECTORS", "20");

    options.setArgs("SCALAR00 SOLVER", "CVODE");
    options.setArgs("SCALAR01 SOLVER", "CVODE");
  }

  // CVODE, just for S01
  if (ciMode == 17) {

    options.setArgs("CVODE", "TRUE");

    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("SCALAR01 CVODE ABSOLUTE TOLERANCE", "1e-10");

    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-10");
    options.setArgs("CVODE GMRES BASIS VECTORS", "20");

    options.setArgs("SCALAR01 SOLVER", "CVODE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }

  // CVODE, both scalars, no dealiasing
  if (ciMode == 18) {

    options.setArgs("CVODE", "TRUE");

    options.setArgs("SCALAR00 CVODE ABSOLUTE TOLERANCE", "1e-10");
    options.setArgs("SCALAR01 CVODE ABSOLUTE TOLERANCE", "1e-10");

    options.setArgs("ADVECTION TYPE", "CONVECTIVE");
    options.setArgs("CVODE ADVECTION TYPE", "CONVECTIVE");

    options.setArgs("SCALAR00 REGULARIZATION METHOD", "HPFRT");
    options.setArgs("SCALAR01 REGULARIZATION METHOD", "HPFRT");
    options.setArgs("SCALAR00 HPFRT MODES", "1");
    options.setArgs("SCALAR01 HPFRT MODES", "1");
    options.setArgs("SCALAR00 HPFRT STRENGTH", "10.0");
    options.setArgs("SCALAR01 HPFRT STRENGTH", "10.0");

    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-10");
    options.setArgs("CVODE GMRES BASIS VECTORS", "20");

    options.setArgs("SCALAR00 SOLVER", "CVODE");
    options.setArgs("SCALAR01 SOLVER", "CVODE");
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  if (ciMode == 23) {
    options.setArgs("NUMBER TIMESTEPS", std::string("0"));
    options.setArgs("END TIME", std::string("0"));
    options.setArgs("CHECKPOINT INTERVAL", std::string("0"));
  }

  if (ciMode == 24 || ciMode == 25 || ciMode == 26 || ciMode == 27 || ciMode == 28) {
    std::string casename = platform->options.getArgs("CASENAME");
    nekrsCheck(casename != "ethierRefine",
               platform->comm.mpiComm(),
               EXIT_FAILURE,
               "Unexpected input file for CI mode %d\n",
               ciMode);
  }

  auto checkCheckpointingFile = [&](std::string fname, int ciPrev) {
    nekrsCheck(!fs::exists(fname),
               platform->comm.mpiComm(),
               EXIT_FAILURE,
               "Missing required checkpoint file %s. Run cimode=%d first!\n",
               fname.c_str(), ciPrev);
  };

  // refine = 2 
  if (ciMode == 24) {
    options.setArgs("END TIME", std::string("0.01"));
    options.setArgs("MESH HREFINEMENT SCHEDULE", "2");
    options.setArgs("DT", std::string("1e-3"));
    options.setArgs("CHECKPOINT INTERVAL", std::string("0"));
  }
  // refine = 2,2
  if (ciMode == 25) {
    options.setArgs("END TIME", std::string("0.005"));
    options.setArgs("MESH HREFINEMENT SCHEDULE", "2,2");
    options.setArgs("DT", std::string("5e-4"));
    options.setArgs("CHECKPOINT INTERVAL", std::string("0"));
  }
  // refine = 2 + checkpoint
  if (ciMode == 26) {
    options.setArgs("MESH HREFINEMENT SCHEDULE", "2");
    options.setArgs("NUMBER TIMESTEPS", std::string("0"));
    options.setArgs("END TIME", std::string("0"));
    options.setArgs("CHECKPOINT INTERVAL", std::string("0"));

    checkCheckpointingFile("ci23.fld", 23);
    checkCheckpointingFile("ci23.bp", 23);
  }
  // refine = 2,2 + checkpoint
  if (ciMode == 27) {
    options.setArgs("MESH HREFINEMENT SCHEDULE", "2,2");
    options.setArgs("NUMBER TIMESTEPS", std::string("0"));
    options.setArgs("END TIME", std::string("0"));
    options.setArgs("CHECKPOINT INTERVAL", std::string("0"));
    
    checkCheckpointingFile("ci23.fld", 23);
    checkCheckpointingFile("ci23.bp", 23);
    checkCheckpointingFile("ci26.fld", 26);
    checkCheckpointingFile("ci26.bp", 26);
  }
  // refine = 2,3 + checkpoint
  if (ciMode == 28) {
    options.setArgs("MESH HREFINEMENT SCHEDULE", "2,3");
    options.setArgs("NUMBER TIMESTEPS", std::string("0"));
    options.setArgs("END TIME", std::string("0"));
    options.setArgs("CHECKPOINT INTERVAL", std::string("0"));
    checkCheckpointingFile("ci23.fld", 23);
    checkCheckpointingFile("ci23.bp", 23);
    checkCheckpointingFile("ci26.fld", 26);
    checkCheckpointingFile("ci26.bp", 26);
  }
}

auto generatePoints(int nPoints, double R)
{
  auto pointsX = randomVector<dfloat>(nPoints, 0, 1, true);
  auto pointsY = randomVector<dfloat>(nPoints, 0, 1, true);
  auto pointsZ = randomVector<dfloat>(nPoints, 0, 1, true);

  const auto third = 1.0 / 3.0;

  for (int i = 0; i < nPoints; ++i) {
    auto x = pointsX[i];
    auto y = pointsY[i];
    auto z = pointsZ[i];

    // map [0,1] to [-1,1]
    x = 2 * x - 1;
    y = 2 * y - 1;
    z = 2 * z - 1;

    // map unit-cube [-1,1]^3 to sphere
    // adapted from https://stackoverflow.com/questions/2656899
    const auto sx = x * std::sqrt(1.0 - 0.5 * y * y - 0.5 * z * z + third * y * y * z * z);
    const auto sy = y * std::sqrt(1.0 - 0.5 * z * z - 0.5 * x * x + third * z * z * x * x);
    const auto sz = z * std::sqrt(1.0 - 0.5 * x * x - 0.5 * y * y + third * x * x * y * y);

    // rescale to desired radius
    auto Rp = randomVector<dfloat>(1, 0, R)[0];
    const auto r = std::sqrt(sx * sx + sy * sy + sz * sz);
    const auto scale = Rp / r;
    pointsX[i] = sx * scale;
    pointsY[i] = sy * scale;
    pointsZ[i] = sz * scale;
  }
  return std::make_tuple(pointsX, pointsY, pointsZ);
}

void ciTestParticles(nrs_t *nrs, lpm_t *particles, occa::kernel exactUVWPKernel, occa::kernel RKKernel)
{
  double endTime;
  platform->options.getArgs("END TIME", endTime);

  // gid of particles currently on rank
  std::vector<int> gids(particles->numParticles());
  std::vector<dfloat> xp0_g(particles->numParticles()); // based on current rank
  std::vector<dfloat> yp0_g(particles->numParticles());
  std::vector<dfloat> zp0_g(particles->numParticles());
  auto o_gidFloat = particles->getProp("gid");

  std::vector<dfloat> gidFloat(o_gidFloat.length());
  o_gidFloat.copyTo(gidFloat.data());

  for (int i = 0; i < particles->numParticles(); ++i) {
    int gid = gidFloat[i];
    gids[i] = gid;
    xp0_g[i] = xpg0[gid];
    yp0_g[i] = ypg0[gid];
    zp0_g[i] = zpg0[gid];
  }

  // compare against RK4 solution
  std::vector<dfloat> xExpected(xp0_g.size(), 0.0);
  std::vector<dfloat> yExpected(yp0_g.size(), 0.0);
  std::vector<dfloat> zExpected(zp0_g.size(), 0.0);

  auto o_y0 = platform->device.malloc<dfloat>(nrs->meshV->dim * xp0_g.size());
  if (particles->numParticles() > 0) {
    auto o_x0coord = o_y0 + 0 * xp0_g.size();
    auto o_y0coord = o_y0 + 1 * yp0_g.size();
    auto o_z0coord = o_y0 + 2 * zp0_g.size();
    o_x0coord.copyFrom(xp0_g.data(), xp0_g.size());
    o_y0coord.copyFrom(yp0_g.data(), yp0_g.size());
    o_z0coord.copyFrom(zp0_g.data(), zp0_g.size());
  }

  auto o_y = solveODE(nrs, 0.0, endTime, 5e-3 * nrs->dt[0], o_y0, exactUVWPKernel, RKKernel);
  if (particles->numParticles() > 0) {
    auto o_xcoord = o_y + 0 * xp0_g.size();
    auto o_ycoord = o_y + 1 * xp0_g.size();
    auto o_zcoord = o_y + 2 * xp0_g.size();

    o_xcoord.copyTo(xExpected.data(), xExpected.size());
    o_ycoord.copyTo(yExpected.data(), yExpected.size());
    o_zcoord.copyTo(zExpected.data(), zExpected.size());
  }

  CiEvalTest("number of particles", xExpected.size() == particles->numParticles()); 

  auto relTol = 1e-7;
  auto floorTol = 1e-12;

  auto LinfErrX = 0.0;
  auto LinfErrY = 0.0;
  auto LinfErrZ = 0.0;

  auto o_xP = particles->getDOF("x");
  auto o_yP = particles->getDOF("y");
  auto o_zP = particles->getDOF("z");
  
  std::vector<dfloat> xP(o_xP.length());
  std::vector<dfloat> yP(o_yP.length());
  std::vector<dfloat> zP(o_xP.length());
  o_xP.copyTo(xP.data());
  o_yP.copyTo(yP.data());
  o_zP.copyTo(zP.data());

  bool passed = true;
  for (int pid = 0; pid < particles->numParticles(); ++pid) {
    auto expectedX = xExpected.at(pid);
    auto expectedY = yExpected.at(pid);
    auto expectedZ = zExpected.at(pid);

    auto x = xP.at(pid);
    auto y = yP.at(pid);
    auto z = zP.at(pid);

    double errX = std::abs(expectedX - x);
    double errY = std::abs(expectedY - y);
    double errZ = std::abs(expectedZ - z);
    if (std::abs(expectedX) > floorTol) {
      errX /= std::abs(expectedX);
    }
    if (std::abs(expectedY) > floorTol) {
      errY /= std::abs(expectedY);
    }
    if (std::abs(expectedZ) > floorTol) {
      errZ /= std::abs(expectedZ);
    }

    LinfErrX = std::max(LinfErrX, errX);
    LinfErrY = std::max(LinfErrY, errY);
    LinfErrZ = std::max(LinfErrZ, errZ);
  }

  std::array<double, 3> LinfErrs = {LinfErrX, LinfErrY, LinfErrZ};
  MPI_Allreduce(MPI_IN_PLACE, LinfErrs.data(), 3, MPI_DOUBLE, MPI_MAX, platform->comm.mpiComm());
  LinfErrX = LinfErrs[0];
  LinfErrY = LinfErrs[1];
  LinfErrZ = LinfErrs[2];

  CiEvalTest("Particle integration (LinfErr: " + to_string_f(LinfErrX) + ", " + to_string_f(LinfErrY) + ", " + to_string_f(LinfErrZ) + ")",
             LinfErrX < relTol && LinfErrY < relTol && LinfErrZ < relTol
            );
}

// Construct multiple distinct pointInterpolation_t objects corresponding to different interpolation points.
// The call-order to the pointInterpolation_t objects should not matter,
// however, pointInterpolation_t itself is not thread-safe, so the two objects cannot be used concurrently.
void ciTestPointInterpolation(nrs_t *nrs, double time, occa::kernel exactUVWPKernel)
{
  auto mesh = nrs->meshV;

  auto o_Ue = platform->device.malloc<dfloat>(nrs->meshV->dim * nrs->fieldOffset);
  auto o_Pe = platform->device.malloc<dfloat>(nrs->fieldOffset);
  exactUVWPKernel(mesh->Nlocal, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->fieldOffset, o_Pe, o_Ue);

  const auto verbosity = pointInterpolation_t::VerbosityLevel::Detailed;
  const auto threshold = (sizeof(dfloat) == sizeof(double)) ? 3e-8 : 1e-5;

  const auto computeRelErrors = [&](dlong n, dlong offset, occa::memory o_Uexact, occa::memory o_Uinterp) {
    const auto absTol = 1e-6;

    auto o_UexactX = o_Uexact.slice(0 * offset, n);
    auto o_UexactY = o_Uexact.slice(1 * offset, n);
    auto o_UexactZ = o_Uexact.slice(2 * offset, n);

    auto o_UinterpX = o_Uinterp.slice(0 * offset, n);
    auto o_UinterpY = o_Uinterp.slice(1 * offset, n);
    auto o_UinterpZ = o_Uinterp.slice(2 * offset, n);

    const auto LinfUx = 
      platform->linAlg->maxRelativeError(n, 1, 0, absTol, o_UinterpX, o_UexactX, platform->comm.mpiComm()); 
    const auto LinfUy = 
      platform->linAlg->maxRelativeError(n, 1, 0, absTol, o_UinterpY, o_UexactY, platform->comm.mpiComm()); 
    const auto LinfUz = 
      platform->linAlg->maxRelativeError(n, 1, 0, absTol, o_UinterpZ, o_UexactZ, platform->comm.mpiComm()); 

    return std::make_tuple(LinfUx, LinfUy, LinfUz);
  };

  // Call find/eval on 0 points with a nullptr
  occa::memory o_NULL;
  auto interpNoop = pointInterpolation_t(mesh, platform->comm.mpiComm());
  interpNoop.setPoints(o_NULL, o_NULL, o_NULL);
  interpNoop.find(verbosity);
  interpNoop.eval(nrs->meshV->dim, nrs->fieldOffset, o_Ue, 0, o_NULL);

  // Interpolate onto points that coincide with the mesh, call findpts
  auto interpFluid = pointInterpolation_t(mesh, platform->comm.mpiComm());
  interpFluid.setPoints(mesh->o_x, mesh->o_y, mesh->o_z);
  interpFluid.find(verbosity);

  auto o_UintFluid = platform->device.malloc<dfloat>(nrs->meshV->dim * nrs->fieldOffset);
  interpFluid.eval(nrs->meshV->dim, nrs->fieldOffset, o_Ue, nrs->fieldOffset, o_UintFluid);

  {
    auto [linfErrUx, linfErrUy, linfErrUz] = computeRelErrors(mesh->Nlocal, nrs->fieldOffset, o_Ue, o_UintFluid);
    CiEvalTest("point interp (fluid) (RelLinfErr: " + 
               to_string_f(linfErrUx) + ", " + to_string_f(linfErrUy) + ", " + to_string_f(linfErrUz) + ")",
               linfErrUx < threshold && linfErrUy < threshold && linfErrUz < threshold
              );
  }

  // Interpolate onto a random number of a random set of points inside the unit ball
  // Repeat for nRandomInterp sets of points
  const double R = 0.75;

  constexpr int nRandomInterp = 3;
  std::array<std::vector<dfloat>, nRandomInterp> pointsX;
  std::array<std::vector<dfloat>, nRandomInterp> pointsY;
  std::array<std::vector<dfloat>, nRandomInterp> pointsZ;
  std::array<pointInterpolation_t*, nRandomInterp> interpRand;

  for (int n = 0; n < nRandomInterp; ++n) {
    std::tie(pointsX[n], pointsY[n], pointsZ[n]) = generatePoints(mesh->Nlocal, R);

    interpRand[n] = new pointInterpolation_t(mesh, platform->comm.mpiComm());
    interpRand[n]->setPoints(pointsX[n], pointsY[n], pointsZ[n]);
    interpRand[n]->find(verbosity);
  }

  std::array<occa::memory, nRandomInterp> o_UintRand;
  for (int n = 0; n < nRandomInterp; ++n) {
    o_UintRand[n] = platform->device.malloc<dfloat>(nrs->meshV->dim * pointsX[n].size());
    if (n == nRandomInterp - 1) {
      interpRand[n]->eval(nrs->meshV->dim, nrs->fieldOffset, o_Ue, pointsX[n].size(), o_UintRand[n]);

      std::tie(pointsX[n], pointsY[n], pointsZ[n]) = generatePoints(mesh->Nlocal / 4, R);
      interpRand[n]->setPoints(pointsX[n], pointsY[n], pointsZ[n]);
      interpRand[n]->find(verbosity);

      interpRand[n]->eval(nrs->meshV->dim, nrs->fieldOffset, o_Ue, pointsX[n].size(), o_UintRand[n]);
    } else {
      interpRand[n]->eval(nrs->meshV->dim, nrs->fieldOffset, o_Ue, pointsX[n].size(), o_UintRand[n]);
    }
  }

  std::array<occa::memory, nRandomInterp> o_UeRand;
  for (int n = 0; n < nRandomInterp; ++n) {
    const dlong nPoints = pointsX[n].size();

    o_UeRand[n] = platform->device.malloc<dfloat>(nrs->meshV->dim * nPoints);

    auto o_Pe = platform->device.malloc<dfloat>(nPoints);

    auto o_x = platform->device.malloc<dfloat>(nPoints);
    o_x.copyFrom(pointsX[n].data());
    auto o_y = platform->device.malloc<dfloat>(nPoints);
    o_y.copyFrom(pointsY[n].data());
    auto o_z = platform->device.malloc<dfloat>(nPoints);
    o_z.copyFrom(pointsZ[n].data());

    exactUVWPKernel(nPoints, time, o_x, o_y, o_z, nPoints, o_Pe, o_UeRand[n]);
  }

  for (int n = 0; n < nRandomInterp; ++n) {
    auto [linfErrUxRand, linfErrUyRand, linfErrUzRand] =
      computeRelErrors(pointsX[n].size(), pointsX[n].size(), o_UeRand[n], o_UintRand[n]);

    CiEvalTest("point interp (random group" + std::to_string(n) + ") (RelLinfErr: " +
               to_string_f(linfErrUxRand) + ", " + to_string_f(linfErrUyRand) + ", " + to_string_f(linfErrUzRand) + ")",
               linfErrUxRand < threshold && linfErrUyRand < threshold && linfErrUzRand < threshold
              );
  }
}

void ciTestLVector(nrs_t *nrs, double time, occa::kernel exactUVWPKernel)
{
  auto mesh = nrs->meshV;

  auto o_U = platform->device.malloc<dfloat>(nrs->meshV->dim * nrs->fieldOffset);
  auto o_P = platform->device.malloc<dfloat>(nrs->fieldOffset);
  exactUVWPKernel(mesh->Nlocal, time, mesh->o_x, mesh->o_y, mesh->o_z, nrs->fieldOffset, o_P, o_U);

  
  auto o_Usave = platform->deviceMemoryPool.reserve<dfloat>(nrs->meshV->dim * nrs->fieldOffset);
  o_U.copyTo(o_Usave, nrs->meshV->dim * nrs->fieldOffset);

  std::vector<mesh_t*> meshes = {nrs->meshV, nrs->meshV, nrs->meshV};
  LVector_t<dfloat> UL(meshes, true);
  LVector_t<dfloat> VL(meshes, false); // create secondary vector, not managed

  UL.copyFromE(nrs->fieldOffset, o_U);
  platform->linAlg->fill(nrs->meshV->dim * nrs->fieldOffset, 0.0, o_U);
  UL.copyToE(nrs->fieldOffset, o_U);

  bool failed = false;

  // sanity check: is o_U the same as before?
  auto checkErr = [&](bool subtract = true){
    auto o_Uerr = platform->deviceMemoryPool.reserve<dfloat>(nrs->meshV->dim * nrs->fieldOffset);
    if(subtract){
      platform->linAlg->axpbyzMany(nrs->meshV->Nlocal,
                                  nrs->meshV->dim,
                                  nrs->fieldOffset,
                                  1.0,
                                  o_Usave,
                                  -1.0,
                                  o_U,
                                  o_Uerr);
    }
    platform->linAlg->aydxMany(nrs->meshV->Nlocal,
                              nrs->meshV->dim,
                              nrs->fieldOffset,
                              1,
                              1.0,
                              o_Usave,
                              o_Uerr);
  
    auto o_Uxerr = o_Uerr + 0 * nrs->fieldOffset;
    auto o_Uyerr = o_Uerr + 1 * nrs->fieldOffset;
    auto o_Uzerr = o_Uerr + 2 * nrs->fieldOffset;
    
    platform->linAlg->abs(nrs->meshV->Nlocal, o_Uxerr);
    platform->linAlg->abs(nrs->meshV->Nlocal, o_Uyerr);
    platform->linAlg->abs(nrs->meshV->Nlocal, o_Uzerr);

    auto LinfUx = platform->linAlg->max(nrs->meshV->Nlocal, o_Uxerr, platform->comm.mpiComm());
    auto LinfUy = platform->linAlg->max(nrs->meshV->Nlocal, o_Uyerr, platform->comm.mpiComm());
    auto LinfUz = platform->linAlg->max(nrs->meshV->Nlocal, o_Uzerr, platform->comm.mpiComm());
    if (platform->comm.mpiRank() == 0) {
      printf("relLinf Ux/y/z error: %g/%g/%g\n", LinfUx, LinfUy, LinfUz);
    }

    const auto tol = 1e2 * std::numeric_limits<dfloat>::epsilon();
    if (LinfUx > tol || LinfUy > tol || LinfUz > tol) {
      failed = true;
    }
  };
  
  checkErr();

  // restore o_U
  o_U.copyFrom(o_Usave, nrs->meshV->dim * nrs->fieldOffset);
  
  // test packing compactly, as in CVODE
  std::vector<dlong> lengths_(nrs->meshV->dim, 0);
  dlong nEq = 0;
  for(int is = 0; is < nrs->meshV->dim; ++is){
    const auto Nlocal = UL.Nlocal(is);
    nEq += Nlocal;
    lengths_[is] = Nlocal;
  }

  auto o_UL = platform->deviceMemoryPool.reserve<dfloat>(nEq);
  UL.fieldOffsets(lengths_);
  UL.optr(o_UL);
  
  UL.copyFromE(nrs->fieldOffset, o_U);
  platform->linAlg->fill(nrs->meshV->dim * nrs->fieldOffset, 0.0, o_U);
  UL.copyToE(nrs->fieldOffset, o_U);
  
  checkErr();
  
  // restore o_U
  o_U.copyFrom(o_Usave, nrs->meshV->dim * nrs->fieldOffset);

  // test same packing, but with un-allocated version
  VL.fieldOffsets(lengths_);
  VL.optr(o_UL);
  VL.copyFromE(nrs->fieldOffset, o_U);
  platform->linAlg->fill(nrs->meshV->dim * nrs->fieldOffset, 0.0, o_U);
  VL.copyToE(nrs->fieldOffset, o_U);
  
  checkErr();

  // restore o_U
  o_U.copyFrom(o_Usave, nrs->meshV->dim * nrs->fieldOffset);

  // different test: fill VL with a constant
  const auto val = 2.8;
  platform->linAlg->fill(nEq, val, VL.optr());

  platform->linAlg->fill(nrs->meshV->dim * nrs->fieldOffset, 0.0, o_U);
  VL.copyToE(nrs->fieldOffset, o_U);

  // subtract val
  auto o_Ux = o_U + 0 * nrs->fieldOffset;
  auto o_Uy = o_U + 1 * nrs->fieldOffset;
  auto o_Uz = o_U + 2 * nrs->fieldOffset;
  platform->linAlg->add(nrs->meshV->Nlocal, -val, o_Ux);
  platform->linAlg->add(nrs->meshV->Nlocal, -val, o_Uy);
  platform->linAlg->add(nrs->meshV->Nlocal, -val, o_Uz);
  
  checkErr(false);


  if(failed)
    CiFailTest();
  
  // Test failed -- exit
  nekrsCheck(platform->exitValue, platform->comm.mpiComm(), EXIT_FAILURE, "%s", "");
}


void ciCheckpointing(nrs_t *nrs)
{
  auto run = [&](const std::string &tag, const std::string& engine, const std::vector<std::string>& fileList, bool readOnly = false)
  {
    auto FP64 = (!readOnly) && platform->options.compareArgs("CHECKPOINT PRECISION", "FP64");
    const auto tol = (FP64) ? 10 * std::numeric_limits<double>::epsilon() : 10 * std::numeric_limits<float>::epsilon();
    const double time = 0.123; // just for testing, actual solution is based on t==0

    if (!readOnly) {
      platform->options.setArgs("CHECKPOINT ENGINE", engine);
      nrs->checkpointWriter = iofldFactory::create();
      nrs->writeCheckpoint(time);
    }

    // wipe out solutions after write
    platform->linAlg->fill(nrs->fluid->fieldOffsetSum, 0.0, nrs->fluid->o_U);
    if (nrs->Nscalar > 0) {
      platform->linAlg->fill(nrs->scalar->fieldOffsetSum, 0.0, nrs->scalar->o_S);
    }

    nrs->restartFromFiles(fileList);

    double startTime;
    platform->options.getArgs("START TIME", startTime);

    nrs->copyToNek(0.0, 0, true);
    nek::userchk();

    auto err = nek::ptr<double>("errors");
    CiEvalTest(tag + " velocity", err[0] < tol);
    CiEvalTest(tag + " pressure", err[1] < tol);
    CiEvalTest(tag + " scalar00", err[2] < tol);
    CiEvalTest(tag + " scalar01", err[3] < tol);
    CiEvalTest(tag + " time", std::abs(startTime - time) < tol);

    nrs->checkpointWriter.reset();
  };

  auto breakMesh = [&]() {
    auto mesh = nrs->meshT;
    platform->linAlg->fill(mesh->Nlocal, -0.12, mesh->o_x);
    platform->linAlg->fill(mesh->Nlocal, -0.12, mesh->o_y);
    platform->linAlg->fill(mesh->Nlocal, -0.12, mesh->o_z);
    auto [x, y, z] = mesh->xyzHost(); // device to host
    for (int i = 0; i < mesh->Nlocal; i++) {
      nekData.xm1[i] = x[i];
      nekData.ym1[i] = y[i];
      nekData.zm1[i] = z[i];
    }
  };

  if (ciMode==23) {
    platform->options.setArgs("CHECKPOINT PRECISION", "FP64");
    platform->options.setArgs("CHECKPOINT OUTPUT MESH", "TRUE");
    run("1", "NEK", {"ethier0.f00000"}); // FP64, with mesh
    platform->options.setArgs("CHECKPOINT PRECISION", "FP32");

    run("2", "NEK", {"ethier0.f00001"}); // FP32, with mesh
    run("3", "NEK", {"ethier0.f00001+U+P+S"}, true);

    breakMesh(); // force loading mesh from file
    run("4", "NEK", {"ethier0.f00000+X", "ethier0.f00001+U+P+S"}, true);

    run("5", "ADIOS", {"ethier.bp"});
    run("6", "ADIOS", {"ethier.bp+U+P+S"}, true);
    run("7", "ADIOS", {"ethier.bp+int"}, true);

    MPI_Barrier(platform->comm.mpiComm());
    if (platform->comm.mpiRank() == 0) { // for cimode = 26, 27
      fs::rename("ethier0.f00000", "ci23.fld");

      if (fs::exists("ci23.bp")) fs::remove_all("ci23.bp");
      fs::rename("ethier.bp", "ci23.bp");
    }
    MPI_Barrier(platform->comm.mpiComm());
  }
  if (ciMode==26) { // ref = 2
    run("1", "NEK", {"ethierRefine0.f00000"});
    run("2", "NEK", {"ci23.fld+U+P+S"}, true);
    breakMesh();
    run("3", "NEK", {"ethierRefine0.f00000+X", "ci23.fld+U+P+S"}, true);
    run("4", "NEK", {"ethierRefine0.f00001"});
    run("5", "ADIOS", {"ethierRefine.bp"});

    run("6", "ADIOS", {"ci23.bp+U+P+S"}, true); // test if read flds match orig mesh
    breakMesh();
    run("7", "ADIOS", {"ci23.bp"}, true);

    MPI_Barrier(platform->comm.mpiComm());
    if (platform->comm.mpiRank() == 0) { // for cimode = 27 and 28
      std::rename("ethierRefine0.f00000", "ci26.fld");

      if (fs::exists("ci26.bp")) fs::remove_all("ci26.bp");
      fs::rename("ethierRefine.bp", "ci26.bp");
    }
    MPI_Barrier(platform->comm.mpiComm());
  }
  if (ciMode==27 || ciMode==28) { // refine = 2,2 or 2,3
    run("1", "NEK", {"ethierRefine0.f00000"});
    run("2", "NEK", {"ci26.fld+U+P+S"}, true);
    run("3", "NEK", {"ci23.fld+U+P+S"}, true);
    breakMesh();
    run("4", "NEK", {"ethierRefine0.f00000+X", "ci23.fld+U", "ci26.fld+P+S"}, true);
    run("5", "NEK", {"ethierRefine0.f00001"});
    run("6", "ADIOS", {"ethierRefine.bp"});

    run("7", "ADIOS", {"ci26.bp+U+P+S"}, true);
    breakMesh();
    run("8", "ADIOS", {"ci26.bp"}, true);
    run("9", "ADIOS", {"ci23.bp+U+P+S"}, true);
    breakMesh();
    run("10", "ADIOS", {"ci23.bp"}, true);

    breakMesh();
    run("11", "ADIOS", {"ethierRefine.bp+X", "ci23.bp+U", "ci26.bp+P+S"}, true);
  }
}

void ciTestDistance(nrs_t *nrs)
{
  auto mesh = nrs->meshV;

  int nbID = 1;
  int bID = 1;
  occa::memory o_bID = platform->deviceMemoryPool.reserve<int>(nbID);
  o_bID.copyFrom(&bID);

  auto o_dist = mesh->minDistance(nbID, o_bID, "cheap_dist");

  std::vector<dfloat> dist(mesh->Nlocal);
  o_dist.copyTo(dist.data(), mesh->Nlocal);

  std::vector<dfloat> x(mesh->Nlocal);
  std::vector<dfloat> y(mesh->Nlocal);
  std::vector<dfloat> z(mesh->Nlocal);

  mesh->o_x.copyTo(x.data());
  mesh->o_y.copyTo(y.data());
  mesh->o_z.copyTo(z.data());

  // distance from wall is 1-r
  dfloat maxRelErr = 0.0;
  std::vector<dfloat> refVal(mesh->Nlocal, 0.0);
  for(int n = 0; n < mesh->Nlocal; ++n){
    const auto r = std::sqrt(x[n]*x[n] + y[n]*y[n] + z[n]*z[n]);

    const auto ref = std::abs(1. - r);
    refVal[n] = ref;
  }

  auto o_refVal = platform->deviceMemoryPool.reserve<dfloat>(nrs->fieldOffset);
  o_refVal.copyFrom(refVal.data(), mesh->Nlocal);

  platform->linAlg->axpby(mesh->Nlocal, 1.0, o_refVal, -1.0, o_dist);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_Jw, o_dist);
  platform->linAlg->axmy(mesh->Nlocal, 1.0, mesh->o_Jw, o_refVal);
  auto intErr = platform->linAlg->sum(mesh->Nlocal, o_dist, platform->comm.mpiComm());
  auto intRef = platform->linAlg->sum(mesh->Nlocal, o_refVal, platform->comm.mpiComm());
  auto relErr = std::abs(intErr / intRef);
  const auto relTol = 2e-2;
  if(relErr > relTol){
    if(platform->comm.mpiRank() == 0){
      std::cout << "ciTestDistance: relative error " << relErr << " is too large!\n";
    }
    CiFailTest();
  }
}

dfloat sumRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat sumr = 0.0;
  for (int i = 0; i < size; ++i) {
    sumr += array[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &sumr, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sumr;
}

dfloat maxRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat maxr = -9e30;
  for (int i = 0; i < size; ++i) {
    maxr = (array[i] > maxr) ? array[i] : maxr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &maxr, 1, MPI_DFLOAT, MPI_MAX, comm);
  return maxr;
}

dfloat minRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat minr = 9e30;
  for (int i = 0; i < size; ++i) {
    minr = (array[i] < minr) ? array[i] : minr;
  }
  MPI_Allreduce(MPI_IN_PLACE, &minr, 1, MPI_DFLOAT, MPI_MIN, comm);
  return minr;
}

dfloat norm2Ref(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm += array[i] * array[i];
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return sqrt(norm);
}

dfloat norm1Ref(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm += fabs(array[i]);
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return norm;
}

dfloat amaxRef(dfloat const *const array, const int size, MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < size; ++i) {
    norm = std::max(norm, std::abs(array[i]));
  }
  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_MAX, comm);
  return norm;
}

dfloat weightedNorm2Many(dfloat const *const array,
                         dfloat const *const weight,
                         const int n,
                         const int Nfields,
                         const int fieldOffset,
                         MPI_Comm comm)
{
  dfloat norm = 0.0;
  for (int i = 0; i < n; ++i) {
    const auto wt = weight[i];
    for (int field = 0; field < Nfields; ++field) {
      const auto value = array[i + field * fieldOffset];
      norm += wt * value * value;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DFLOAT, MPI_SUM, comm);
  return std::sqrt(norm);
}

void ciTestLinAlg(nrs_t *nrs, const int N)
{
  linAlg_t *linAlg = platform->linAlg;
  MPI_Comm comm = platform->comm.mpiComm();
  int rank = platform->comm.mpiRank();
  auto x = randomVector<dfloat>(N, 0, 1, true);
  occa::memory o_x = platform->device.malloc<dfloat>(N, x.data());

  const dfloat testTol = 50. * std::numeric_limits<dfloat>::epsilon();

  auto evaluateMethod = [&](dfloat value, std::function<dfloat(dfloat const *const array, const int size, MPI_Comm comm)> referenceMethod, std::string id) 
  {
    const dfloat reference = referenceMethod(x.data(), N, comm);
    const dfloat err = abs((value - reference) / reference);
    CiEvalTest("linAlg " + id, err < testTol);
  };

  evaluateMethod(linAlg->sum(N, o_x, comm), &sumRef, "sum");
  evaluateMethod(linAlg->max(N, o_x, comm), &maxRef, "max");
  evaluateMethod(linAlg->min(N, o_x, comm), &minRef, "min");
  evaluateMethod(linAlg->norm2(N, o_x, comm), &norm2Ref, "norm2");
  evaluateMethod(linAlg->norm1(N, o_x, comm), &norm1Ref, "norm1");
  evaluateMethod(linAlg->amax(N, o_x, comm), &amaxRef, "amax");
  o_x.free();

  int fieldOffset = alignStride<dfloat>(N);

  constexpr int maxFields = 10;
  for (int Nfields = 1; Nfields <= maxFields; ++Nfields) {
    auto x = randomVector<dfloat>(Nfields * fieldOffset, 0, 1, true);
    auto wt = randomVector<dfloat>(N, 0, 1, true);
    occa::memory o_x = platform->device.malloc<dfloat>(x.size(), x.data());
    occa::memory o_wt = platform->device.malloc<dfloat>(wt.size(), wt.data());

    const auto reference = weightedNorm2Many(x.data(), wt.data(), N, Nfields, fieldOffset, comm);
    const auto value = linAlg->weightedNorm2Many(N, Nfields, fieldOffset, o_wt, o_x, comm);
    const dfloat err = abs((value - reference) / reference);
#if 0
    if (rank == 0) {
      printf("linAlg errs (N=%d, Nfields=%d): %s=%g\n", N, Nfields, "weightedNorm2Many", err);
    }
#endif
    CiEvalTest("linAlg weightedNorm2Many (N=" + std::to_string (N) + ", Nfields=" + std::to_string(Nfields) + ")", err < testTol);

    o_x.free();
    o_wt.free();
  }
}

