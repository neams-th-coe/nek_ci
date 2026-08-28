[Mesh]
  type = NekRSMesh
  boundary = '1'
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
  execute_on = final
  hide = 'drag area utau rel_err'
  file_base = 'nek1_out'
[]

[Postprocessors]
  [drag]
    type = NekViscousSurfaceForce
    boundary = '1'
    mesh = fluid
    component = x
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
  [rel_err]
    type = ParsedPostprocessor
    expression = 'abs(utau - 4.58794e-2) / 4.58794e-2'
    pp_names = 'utau'
  []
  [pass]
    type = ParsedPostprocessor
    expression = 'if (rel_err < 4.0e-3, 1, 0)'
    pp_names = 'rel_err'
  []
[]
