#include <math.h>

static int ciMode = 0;

dfloat hill_step(dfloat x, dfloat w, dfloat h);
dfloat hill_height(dfloat x, dfloat Lx, dfloat w, dfloat h);

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("CONSTANT FLOW RATE", "TRUE");
  options.setArgs("CONSTANT FLOW DIRECTION", "X");
  if (ciMode == 1) {
    options.setArgs("RESTART FILE NAME", std::string("sst.fld+time=0"));
  }
  if (ciMode == 2) {
    options.setArgs("RESTART FILE NAME", std::string("ddes.fld+time=0"));
  }
  if (ciMode == 3) {
    options.setArgs("RESTART FILE NAME", std::string("iddes.fld+time=0"));
  }
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("FLUID VISCOSITY", std::to_string(1.0 / 10595.0));
  options.setArgs("FLUID DENSITY", std::string("1.0"));
  options.setArgs("END TIME", std::string("0.5"));
  options.setArgs("ADVECTION TYPE", "CUBATURE+CONVECTIVE");

  options.setArgs("DT", std::string("5e-3"));
  options.setArgs("VARIABLE DT", std::string("TRUE"));
  options.setArgs("TARGET CFL", std::string("1.5"));
  options.setArgs("SUBCYCLING STEPS", std::string("1"));

  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", std::string("1e-4"));

  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-6"));

  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-6"));

  options.setArgs("SCALAR01 SOLVER TOLERANCE", std::string("1e-6"));
}

void computeCfError(nrs_t *nrs, double time, int tstep, const occa::memory& o_avgU)
{
  if (!nrs->lastStep) {
    return;
  }

  auto mesh = nrs->meshV;

  auto o_Sij = nrs->strainRate(nrs->fluid->fieldOffset, o_avgU);

  std::vector<int> bidWall = {1};
  auto o_bidWall = platform->device.malloc<int>(bidWall.size(), bidWall.data());

  auto o_tauw = nrs->viscousShearStress(o_bidWall, o_Sij);

  // https://turbmodels.larc.nasa.gov/Other_LES_Data/2Dhill_periodic/hill_LES_cf_digitized.dat
  std::ifstream file("refLes.dat");
  std::string line;
  std::vector<dfloat> lesData;
  while (file >> line) {
    lesData.push_back(std::stod(line));
  }
  file.close();

  deviceMemory<dfloat> o_lesData(lesData);

  auto o_cferr = nrs->scalar->o_solution("cferr");
  platform->linAlg->fill(mesh->Nlocal, 0.0, o_cferr);

  auto xmin = platform->linAlg->min(mesh->Nlocal, mesh->o_x, platform->comm.mpiComm());
  auto xmax = platform->linAlg->max(mesh->Nlocal, mesh->o_x, platform->comm.mpiComm());
  cfErrKernel(mesh->Nelements,
              bidWall.size(),
              o_bidWall,
              mesh->o_vmapM,
              mesh->o_EToB,
              mesh->o_x,
              o_tauw,
              lesData.size(),
              xmin,
              xmax,
              o_lesData,
              o_cferr);
}
