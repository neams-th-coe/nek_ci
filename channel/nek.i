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
  #  time_step_interval = 100
  #[]
[]

P_U0 = 1.0
P_PHX = ${fparse pi}/2.0
P_PHY = ${fparse pi}/2.0

TOL_VX = 9.27E-07
TOL_VY = 9.27E-07
TOL    = 1.00E-11

[Functions]
  [uexact]
    type = ParsedFunction
    expression = '${P_U0} * cos(${fparse pi}*xrot + ${P_PHX}) * sin(2.0*${fparse pi}*yrot + ${P_PHY})'
    symbol_names = 'xrot yrot'
    symbol_values = 'xrot yrot'
  []
  [vexact]
    type = ParsedFunction
    expression = '-0.5 * ${P_U0} * sin(${fparse pi}*xrot + ${P_PHX}) * cos(2.0*${fparse pi}*yrot + ${P_PHY})'
    symbol_names = 'xrot yrot'
    symbol_values = 'xrot yrot'
  []
  [xrot]
    type = ParsedFunction
    expression = 'x * cos(${P_ROT}) + y * sin(${P_ROT})'
  []
  [yrot]
    type = ParsedFunction
    expression = '-x * sin(${P_ROT}) + y * cos(${P_ROT})'
  []
  [ur]
    type = ParsedFunction
    expression = 'uexact * cos(${P_ROT}) - vexact * sin(${P_ROT})'
    symbol_names = 'uexact vexact'
    symbol_values = 'uexact vexact'
  []
  [vr]
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
    function = ur
    execute_on = final
  []
  [verrl2]
    type = NekVolumeNorm
    field = velocity_y
    function = vr
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
