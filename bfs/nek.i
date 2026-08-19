[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'bfs'
  initialize_usrwrk = false
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
  console = true
  file_base = 'nek_out'
[]

P_EPS = 1.38E-03

[Postprocessors]
  [cferrInt]
    type = NekSideIntegral
    field = scalar03
    boundary = '4'
    execute_on = final
  []
  [area]
    type = NekSideIntegral
    field = unity
    boundary = '4'
    execute_on = final
  []
  [cferr]
    type = ParsedPostprocessor
    expression = 'cferrInt / area'
    pp_names = 'cferrInt area'
    execute_on = final
  []
  [pass]
    type = ParsedPostprocessor
    expression = 'if (cferr < ${P_EPS}, 1, 0)'
    pp_names = 'cferr'
    execute_on = final
  []
[]
