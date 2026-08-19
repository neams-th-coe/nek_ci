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

  #[console]
  #  type = Console
  #  time_step_interval = 100
  #[]
[]

XLEN = 1.0
YLEN = 1.5
ZLEN = 0.1

P_GAMMA = 1.4
P_OMEGA = 3.141592653589793
P_AMP = 1.5707963267948966
P_ROT = 0.0

pex0 = 1.0

# NekRS ciMode = 2 uses a 5% relative tolerance around the
# stored reference errors.
CI_EPS = 5.0E-2

# ciMode = 2 reference errors
REF_VOL_ERROR     = 0.2465620E-06
REF_YPISTON_ERROR = 0.2465620E-05
REF_TEMP_ERROR    = 0.3000000E-05

# Additional Cardinal-only test; not part of NekRS ciMode = 2.
TERMV_ABS_TOL = 1.0E-11

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

  # Reconstruct the piston location from the computed volume:
  #
  #   volume = XLEN * ZLEN * (0.5 - ypiston)
  #
  # Therefore:
  #
  #   ypiston = 0.5 - volume / (XLEN * ZLEN)
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

  # ---------------------------------------------------------------------------
  # Physical errors
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
  # NekRS ciMode = 2 regression errors
  #
  # Match ciTestErrors:
  #
  #   abs(current_error - reference_error) / reference_error
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

  # Additional Cardinal-only check. This is deliberately excluded
  # from the ciMode = 2-equivalent overall pass result.
  [termV_pass]
    type = ParsedPostprocessor
    expression = 'if(termV_error < ${TERMV_ABS_TOL}, 1, 0)'
    pp_names = 'termV_error'
    execute_on = final
  []

  # ---------------------------------------------------------------------------
  # Overall result
  #
  # Included ciMode = 2 checks available through current Cardinal objects:
  #   - volume
  #   - piston position
  #   - temperature
  #
  # Not available without NekInfoPostprocessor:
  #   - p0th
  #   - dp0thdt
  #   - velocity iteration count
  #   - pressure iteration count
  # ---------------------------------------------------------------------------

  [pass]
    type = ParsedPostprocessor
    expression = 'if(volume_pass &
                     ypiston_pass &
                     temperature_pass,
                     1,
                     0)'
    pp_names = 'volume_pass
                ypiston_pass
                temperature_pass'
    execute_on = final
  []
[]
