[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'lowMach'
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
  #  time_step_interval = 300
  #[]
[]

P_DELTA = 0.2
PRESSURE_SHIFT = -1.4919313668E-01

TOL_V  = 1.59E-05
TOL_P  = 3.19E-03
TOL_T  = 3.07E-07
TOL_REL    = 1.00E-2

[Functions]
  [uexact]
    type = ParsedFunction
    expression = '0.5 * (3.0 + tanh(xd))'
    symbol_names = 'xd'
    symbol_values = 'xd'
  []
  [xd]
    type = ParsedFunction
    expression = 'x / ${P_DELTA}'
  []
  [aa]
    type = ParsedFunction
    expression = '3./2. - (tanh(1.0) - tanh(-1.0))/3.0'
  []
  [qtle]
    type = ParsedFunction
    expression = '0.5/${P_DELTA} * (1.0 - (tanh(xd) * tanh(xd)))'
    symbol_names = 'xd'
    symbol_values = 'xd'
  []
  [pexact]
    type = ParsedFunction
    expression = '4./3.0 * qtle - uexact + aa + ${PRESSURE_SHIFT}'
    symbol_names = 'qtle uexact aa'
    symbol_values = 'qtle uexact aa'
  []
[]

[Postprocessors]
  # Calculate L_infinity errors
  [uxerr]
    type = NekVolumeNorm
    field = velocity_x
    function = uexact
    N = infinity
    execute_on = final
  []
  [terr]
    type = NekVolumeNorm
    field = temperature
    function = uexact
    N = infinity
    execute_on = final
  []
  [perr]
    type = NekVolumeNorm
    field = pressure
    function = pexact
    N = infinity
    execute_on = final
  []
  
  # Check if all test passed
  [pass]
    type = ParsedPostprocessor
    expression = 'if(
      (abs((uxerr - ${TOL_V}) / uxerr) < ${TOL_REL}) &
      (abs((terr - ${TOL_T}) / terr) < ${TOL_REL}) &
      (abs((perr - ${TOL_P}) / perr) < ${TOL_REL}),
      1,
      0)'
    pp_names = 'uxerr terr perr'
    execute_on = final
  []
[]
