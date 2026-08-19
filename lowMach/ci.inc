#include <math.h>

static int ciMode = 0;
static int NELGX;
static int NELGY;
static int NELGZ;

double numCVODERHS = 0;

void ciSetup(MPI_Comm comm, setupAide &options)
{
  options.setArgs("POLYNOMIAL DEGREE", std::string("7"));
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("END TIME", std::string("0.3"));
  options.setArgs("DT", std::string("1e-3"));
  options.setArgs("SUBCYCLING STEPS", std::string("0"));
  if (ciMode == 2) {
    options.setArgs("SUBCYCLING STEPS", std::string("1"));
  }
  options.setArgs("BDF ORDER", "2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", std::string("1e-12"));
  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", std::string("1e-08"));
  options.setArgs("SCALAR00 SOLVER TOLERANCE", std::string("1e-12"));

  // set up CVODE solver details
  if (ciMode == 3) {
    options.setArgs("CVODE", "TRUE");
    options.setArgs("CVODE INTEGRATOR", "BDF");
    options.setArgs("CVODE MAX TIMESTEPPER ORDER", "3");
    options.setArgs("CVODE EPS LIN", "0.05");
    options.setArgs("CVODE STOP TIME", "TRUE");

    const auto hmaxRatio = 1;
    options.setArgs("CVODE HMAX RATIO", to_string_f(hmaxRatio));

    const auto maxSteps = 10000;
    options.setArgs("CVODE MAX STEPS", std::to_string(maxSteps));

    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-08");
    options.setArgs("CVODE ABSOLUTE TOLERANCE", "1e-12");
    options.setArgs("CVODE GMR VECTORS", "20");

    options.setArgs("SCALAR00 SOLVER", "CVODE");
  }
}
