mesh_width = 0.0193145

[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'cavity'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]


[Postprocessors]
  [heat_flux]
    type = NekHeatFluxIntegral
    boundary = '1'
  []
  [nusselt]
    type = ParsedPostprocessor
    expression = '-heat_flux / width_y'
    pp_names = 'heat_flux'
    constant_names = 'width_y'
    constant_expressions = '${mesh_width}'
  []
[]

[Outputs]
  csv=true
  execute_on=final
[]

