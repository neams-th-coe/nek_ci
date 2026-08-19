[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'mv_cyl'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  csv = true
  show = 'pass'
  execute_on = final
  console = false
  file_base = 'nek_out'
[]

# -----------------------------------------------------------------------------
# Physical parameters
# -----------------------------------------------------------------------------

XLEN = 1.0
YLEN = 1.5
ZLEN = 0.1

P_GAMMA = 1.4
P_OMEGA = 3.141592653589793
P_AMP = 1.5707963267948966
P_ROT = 0.0

pex0 = 1.0

CI_EPS = 5.0E-2

REF_VOL_ERROR     = 0.2465620E-06
REF_YPISTON_ERROR = 0.2465620E-05
REF_TEMP_ERROR    = 0.4929804E-06

TERMV_ABS_TOL = 1.0E-11

# -----------------------------------------------------------------------------
# Exact solution
# -----------------------------------------------------------------------------

[Functions]
  [unitFunction]
    type = ParsedFunction
    expression = '1.0'
  []

  [areap]
    type = ParsedFunction
    expression = '${XLEN} * ${ZLEN}'
  []

  [volex0]
    type = ParsedFunction
    expression = 'areap * ${YLEN}'
    symbol_names = 'areap'
    symbol_values = 'areap'
  []

  [vpex]
    type = ParsedFunction
    expression = '${P_AMP} * sin(${P_OMEGA} * t)'
  []

  [volex]
    type = ParsedFunction
    expression = 'volex0 + areap * ${P_AMP} * (cos(${P_OMEGA} * t) - 1.0) / ${P_OMEGA}'
    symbol_names = 'volex0 areap'
    symbol_values = 'volex0 areap'
  []

  [pex]
    type = ParsedFunction
    expression = '${pex0} * (volex0 / volex)^${P_GAMMA}'
    symbol_names = 'volex0 volex'
    symbol_values = 'volex0 volex'
  []

  [dpdtex]
    type = ParsedFunction
    expression = '${P_GAMMA} * ${pex0} * volex0^${P_GAMMA} * areap * vpex / volex^(${P_GAMMA} + 1.0)'
    symbol_names = 'volex0 areap vpex volex'
    symbol_values = 'volex0 areap vpex volex'
  []

  [qtlex]
    type = ParsedFunction
    expression = '((${P_GAMMA} - 1.0) / ${P_GAMMA} - 1.0) * dpdtex / pex'
    symbol_names = 'dpdtex pex'
    symbol_values = 'dpdtex pex'
  []

  [ypex]
    type = ParsedFunction
    expression = '-0.5 * (1.0 + cos(${P_OMEGA} * t)) * cos(${P_ROT})'
  []

  [tex]
    type = ParsedFunction
    expression = 'pex^((${P_GAMMA} - 1.0) / ${P_GAMMA})'
    symbol_names = 'pex'
    symbol_values = 'pex'
  []

  [termVex]
    type = ParsedFunction
    expression = 'vpex * ${XLEN} * ${ZLEN}'
    symbol_names = 'vpex'
    symbol_values = 'vpex'
  []
[]

# -----------------------------------------------------------------------------
# Postprocessors
# -----------------------------------------------------------------------------

[Postprocessors]
  # ---------------------------------------------------------------------------
  # Numerical NekRS quantities
  # ---------------------------------------------------------------------------

  [temperature_average]
    type = NekVolumeAverage
    field = temperature
    execute_on = final
  []

  [volume]
    type = NekVolumeIntegral
    field = unity
    function = unitFunction
    execute_on = final
  []

  # For this geometry, the fixed upper boundary is at y = 0.5 and:
  #
  #   volume = area * (0.5 - y_piston)
  #
  # Therefore:
  #
  #   y_piston = 0.5 - volume / area
  #
  # This reconstructs the numerical piston location without requiring direct
  # access to the Nek5000 ypist variable.
  [ypiston]
    type = ParsedPostprocessor
    expression = '0.5 - volume / (${XLEN} * ${ZLEN})'
    pp_names = 'volume'
    execute_on = final
  []

  # Additional Cardinal-only boundary-flux diagnostic.
  [termV_x]
    type = NekSideIntegral
    field = velocity_x
    boundary = '1'
    execute_on = final
  []

  [termV_y]
    type = NekSideIntegral
    field = velocity_y
    boundary = '1'
    execute_on = final
  []

  [termV]
    type = ParsedPostprocessor
    expression = '-termV_x * sin(${P_ROT}) + termV_y * cos(${P_ROT})'
    pp_names = 'termV_x termV_y'
    execute_on = final
  []

  # ---------------------------------------------------------------------------
  # Exact values
  # ---------------------------------------------------------------------------

  [temperature_exact]
    type = FunctionValuePostprocessor
    function = tex
    execute_on = final
  []

  [volume_exact]
    type = FunctionValuePostprocessor
    function = volex
    execute_on = final
  []

  [ypiston_exact]
    type = FunctionValuePostprocessor
    function = ypex
    execute_on = final
  []

  [termV_exact]
    type = FunctionValuePostprocessor
    function = termVex
    execute_on = final
  []

  # These exact values are retained for output or future tests. Cardinal does
  # not presently expose the corresponding numerical p0th and dp0thdt values
  # through standard Nek postprocessors.
  [p0th_exact]
    type = FunctionValuePostprocessor
    function = pex
    execute_on = final
  []

  [dpdt_exact]
    type = FunctionValuePostprocessor
    function = dpdtex
    execute_on = final
  []

  # ---------------------------------------------------------------------------
  # Physical numerical errors
  # ---------------------------------------------------------------------------

  [temperature_error]
    type = ParsedPostprocessor
    expression = 'abs(temperature_average - temperature_exact)'
    pp_names = 'temperature_average temperature_exact'
    execute_on = final
  []

  [vol_error]
    type = ParsedPostprocessor
    expression = 'abs(volume_exact - volume)'
    pp_names = 'volume_exact volume'
    execute_on = final
  []

  [ypiston_error]
    type = ParsedPostprocessor
    expression = 'abs(ypiston_exact - ypiston)'
    pp_names = 'ypiston_exact ypiston'
    execute_on = final
  []

  [termV_error]
    type = ParsedPostprocessor
    expression = 'abs(termV_exact - termV)'
    pp_names = 'termV_exact termV'
    execute_on = final
  []

  # ---------------------------------------------------------------------------
  # NekRS ciMode = 1 normalized regression differences
  #
  # This reproduces:
  #
  #   abs(error - expected_error) / expected_error
  #
  # used by ciTestErrors().
  # ---------------------------------------------------------------------------

  [vol_regression_error]
    type = ParsedPostprocessor
    expression = 'abs(vol_error - ${REF_VOL_ERROR}) / ${REF_VOL_ERROR}'
    pp_names = 'vol_error'
    execute_on = final
  []

  [ypiston_regression_error]
    type = ParsedPostprocessor
    expression = 'abs(ypiston_error - ${REF_YPISTON_ERROR}) / ${REF_YPISTON_ERROR}'
    pp_names = 'ypiston_error'
    execute_on = final
  []

  [temperature_regression_error]
    type = ParsedPostprocessor
    expression = 'abs(temperature_error - ${REF_TEMP_ERROR}) / ${REF_TEMP_ERROR}'
    pp_names = 'temperature_error'
    execute_on = final
  []

  # ---------------------------------------------------------------------------
  # Individual pass indicators
  # ---------------------------------------------------------------------------

  [volume_pass]
    type = ParsedPostprocessor
    expression = 'if(vol_regression_error < ${CI_EPS}, 1, 0)'
    pp_names = 'vol_regression_error'
    execute_on = final
  []

  [ypiston_pass]
    type = ParsedPostprocessor
    expression = 'if(ypiston_regression_error < ${CI_EPS}, 1, 0)'
    pp_names = 'ypiston_regression_error'
    execute_on = final
  []

  [temperature_pass]
    type = ParsedPostprocessor
    expression = 'if(temperature_regression_error < ${CI_EPS}, 1, 0)'
    pp_names = 'temperature_regression_error'
    execute_on = final
  []

  # Additional Cardinal-only diagnostic. It is intentionally not part of
  # `pass`, because ciMode = 1 does not test termV.
  [termV_pass]
    type = ParsedPostprocessor
    expression = 'if(termV_error < ${TERMV_ABS_TOL}, 1, 0)'
    pp_names = 'termV_error'
    execute_on = final
  []

  # ---------------------------------------------------------------------------
  # Overall NekRS-equivalent result
  #
  # Included ciMode = 1 tests:
  #   - volume
  #   - piston position
  #   - temperature
  #
  # Not included because standard Cardinal postprocessors do not expose them:
  #   - thermodynamic pressure p0th
  #   - dp0thdt
  #   - final velocity iteration count
  #   - final pressure iteration count
  # ---------------------------------------------------------------------------

  [pass]
    type = ParsedPostprocessor
    expression = 'if(volume_pass & ypiston_pass & temperature_pass, 1, 0)'
    pp_names = 'volume_pass ypiston_pass temperature_pass'
    execute_on = final
  []
[]

