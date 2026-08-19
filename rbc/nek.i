[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = rbc
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
[]

[Postprocessors]
  [wT_int]
    type = NekVolumeIntegral
    field = scalar01
    execute_on = final
  []

  [volume]
    type = NekVolumeIntegral
    field = unity
    execute_on = final
  []

  [wT_avg]
    type = ParsedPostprocessor
    pp_names = 'wT_int volume'
    expression = 'wT_int / volume'
    execute_on = final
  []

  [Nu_v]
    type = ParsedPostprocessor
    pp_names = 'wT_avg'
    expression = '1.0 + sqrt(1.7e5 * 0.7) * wT_avg'
    execute_on = final
  []

  [Nu_v_error]
    type = ParsedPostprocessor
    pp_names = 'Nu_v'
    expression = 'abs(Nu_v - 5.0)'
    execute_on = final
  []

  [pass]
    type = ParsedPostprocessor
    pp_names = 'Nu_v_error'
    expression = 'Nu_v_error < 9.8e-2'
    execute_on = final
  []
[]

