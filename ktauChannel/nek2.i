rhoRef = 998.0
lenScale = 0.01
Re = 43500.0
velScale = ${fparse Re * 1e-3 / (rhoRef * lenScale)}
utauRef = ${fparse 4.58794e-2 * velScale}

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
  file_base = 'nek2_out'
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
    mesh = fluid
  []

  [utau]
    type = ParsedPostprocessor
    expression = 'sqrt(abs(drag) / (${rhoRef} * area))'
    pp_names = 'drag area'
  []

  [rel_err]
    type = ParsedPostprocessor
    expression = 'abs(utau - ${utauRef}) / ${utauRef}'
    pp_names = 'utau'
  []

  [pass]
    type = ParsedPostprocessor
    #TODO: Update post v26.1
    expression = 'if(rel_err < 4.0e-1, 1, 0)'
    pp_names = 'rel_err'
  []
[]
