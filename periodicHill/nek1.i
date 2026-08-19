[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'periodicHill'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

TOL_CF = 1.04e-2

[Postprocessors]
  [cferr_integral]
    type = NekSideIntegral
    field = scalar03
    boundary = '1'
    execute_on = final
  []

  [surface_area]
    type = NekSideIntegral
    field = unity
    boundary = '1'
    execute_on = final
  []

  [cferr]
    type = ParsedPostprocessor
    expression = 'cferr_integral / surface_area'
    pp_names = 'cferr_integral surface_area'
    execute_on = final
  []

  [pass]
    type = ParsedPostprocessor
    expression = 'if(cferr < ${TOL_CF}, 1, 0)'
    pp_names = 'cferr'
    execute_on = final
  []
[]

[Outputs]
  csv = true
  show = 'pass'
  execute_on = final
  console = true
  file_base = 'nek_out'
[]
