#include <math.h>
#include <vector>

static int ciMode = 0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("RESTART FILE NAME", std::string("bfs.fld+time=0"));  
  options.setArgs("END TIME", std::string("0.005"));
  options.setArgs("DT", std::string("1e-04"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
}

void extractLine(occa::memory& o_tauw)
{
  const auto np = (platform->comm.mpiRank() == 0) ? 200 : 0;
  const auto offset = np;

  static pointInterpolation_t *interpolator = nullptr;
  static std::vector<dfloat> xp, yp, zp;
  static deviceMemory<dfloat> o_tauwp;

  if (!interpolator) {
    auto mesh = nrs->meshV;

    if (np) {
      const auto x0 = 0.0;
      const auto y0 = -0.25;
      const auto z0 = 0.0;

      xp.push_back(x0);
      yp.push_back(y0);
      zp.push_back(z0);

      const auto dx = 10.0 / (np - 1);
      for (int i = 1; i < np - 1; i++) {
        xp.push_back(x0 + i * dx);
        yp.push_back(y0);
        zp.push_back(z0);
      }

      xp.push_back(10.0);
      yp.push_back(y0);
      zp.push_back(z0);
      o_tauwp.resize(offset);
    }

    interpolator = new pointInterpolation_t(mesh, platform->comm.mpiComm());
    interpolator->setPoints(xp, yp, zp);
    interpolator->find();
  }

  interpolator->eval(1, nrs->fluid->fieldOffset, o_tauw, offset, o_tauwp);

  if (platform->comm.mpiRank() == 0) {
    std::vector<dfloat> tauwp(np);
    o_tauwp.copyTo(tauwp);

    std::ofstream f("profile.dat");
    for (int i = 0; i < np; i++) {
      f << std::scientific << " " << xp[i] << " " << yp[i] << " " << zp[i] << " "
        << 2.0 * tauwp[i + 0 * offset] << std::endl;
    }
    f.close();
  }
}

void getCfErr(occa::memory &o_err)
{
  auto mesh = nrs->meshV;

  std::vector<int> bidWall = {4};
  auto o_bidWall = platform->device.malloc<int>(bidWall.size(), bidWall.data());

  auto o_Sij = nrs->strainRate();
  auto o_tau = nrs->viscousShearStress(o_bidWall, o_Sij);

  // extractLine(o_tau);

  cfErrKernel(mesh->Nelements,
              bidWall.size(),
              o_bidWall,
              mesh->o_vmapM,
              mesh->o_EToB,
              mesh->o_x,
              o_tau,
              o_err);
}

