[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'turbPipe'
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
  #[csv]
  #  type = CSV
  #  time_step_interval = 1
  #[]
[]

# Reference value from:
# https://www.lstm.tf.fau.de/database/simulation-database/
P_utauRef = 5.7901316173625254E-02
P_EPS = 2.00E-02

[Postprocessors]
  [drag]
    type = NekViscousSurfaceForce
    boundary = '1'
    mesh = fluid
    component = z
  []
  [area]
    type = NekSideIntegral
    field = unity
    boundary = '1'
  []
  [utau]
    type = ParsedPostprocessor
    expression = 'sqrt(abs(drag) / area)'
    pp_names = 'drag area'
  []
  [err]
    type = ParsedPostprocessor
    expression = 'abs(utau - ${P_utauRef}) / ${P_utauRef}'
    pp_names = 'utau'
  []
  [pass]
    type = ParsedPostprocessor
    expression = 'if (err < ${P_EPS}, 1, 0)'
    pp_names = 'err'
  []
[]
