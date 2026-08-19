[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'channel'
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
  #  time_step_interval = 50
  #[]
[]

TOL_VX = 6.00E-15
TOL_VY = 3.89E-15
TOL    = 1.00E-11

[Functions]
  [h_tilde]
    type = ParsedFunction
    expression = 'y * cos(${P_ROT}) - x * sin(${P_ROT})'
  []
  [uexact]
    type = ParsedFunction
    expression = '1.5 * (1.0 - h_tilde * h_tilde)'
    symbol_names = 'h_tilde'
    symbol_values = 'h_tilde'
  []
  [vexact]
    type = ParsedFunction
    expression = '0.0'
  []
  [urotated]
    type = ParsedFunction
    expression = 'uexact * cos(${P_ROT}) - vexact * sin(${P_ROT})'
    symbol_names = 'uexact vexact'
    symbol_values = 'uexact vexact'
  []
  [vrotated]
    type = ParsedFunction
    expression = 'uexact * sin(${P_ROT}) + vexact * cos(${P_ROT})'
    symbol_names = 'uexact vexact'
    symbol_values = 'uexact vexact'
  []
[]

[Postprocessors]
  [uerrl2]
    type = NekVolumeNorm
    field = velocity_x
    function = urotated
    execute_on = final
  []
  [verrl2]
    type = NekVolumeNorm
    field = velocity_y
    function = vrotated
    execute_on = final
  []
  
  [pass]
    type = ParsedPostprocessor
    expression = 'if((uerrl2 < ${TOL_VX} | uerrl2 < ${TOL}) &
                     (verrl2 < ${TOL_VY} | verrl2 < ${TOL}), 1, 0)'
    pp_names = 'uerrl2 verrl2'
    execute_on = final
  []
[]
