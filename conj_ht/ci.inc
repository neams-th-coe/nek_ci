#include <cmath>
#include <vector>

static int ciMode = 0;

void ciSetup(MPI_Comm comm, setupAide& options)
{
  options.setArgs("POLYNOMIAL DEGREE", "7");
  options.setArgs("RESTART FILE NAME", "restart.fld+time=0");
  options.setArgs("SOLUTION OUTPUT INTERVAL", "0");
  options.setArgs("END TIME", "1");
  options.setArgs("DT", "1.0e-2");
  options.setArgs("SUBCYCLING STEPS", "0");

  if (ciMode == 2) {
    options.setArgs("SUBCYCLING STEPS", "1");
  }

  options.setArgs("BDF ORDER", "2");
  options.setArgs("ADVECTION TYPE", "CONVECTIVE+CUBATURE");
  options.setArgs("FLUID VELOCITY SOLVER TOLERANCE", "1e-12");
  options.setArgs("FLUID PRESSURE SOLVER TOLERANCE", "1e-10");
  options.setArgs("SCALAR00 SOLVER TOLERANCE", "1e-12");
  options.setArgs("VARIABLEPROPERTIES", "TRUE");

  if (ciMode == 3) {
    options.setArgs("CVODE", "TRUE");
    options.setArgs("CVODE INTEGRATOR", "BDF");
    options.setArgs("CVODE MAX TIMESTEPPER ORDER", "3");
    options.setArgs("CVODE EPS LIN", "0.05");
    options.setArgs("CVODE STOP TIME", "TRUE");
    options.setArgs("CVODE HMAX RATIO", "1");
    options.setArgs("CVODE MAX STEPS", "10000");
    options.setArgs("CVODE RELATIVE TOLERANCE", "1e-06");
    options.setArgs("CVODE ABSOLUTE TOLERANCE", "1e-08");
    options.setArgs("CVODE GMR VECTORS", "20");
    options.setArgs("SCALAR00 SOLVER", "CVODE");
  }
}

