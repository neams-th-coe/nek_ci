[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'conj_ht'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

# Analytical-solution parameters; these must match conj_ht.par.
Pe   = 1000.0
q    = 1.0
k    = 10.0
HP_H = 0.5

c1 = ${fparse -Pe * q / k}
c2 = ${fparse -Pe * q * HP_H}
c3 = ${fparse 2.0 * q * HP_H}
c4 = ${fparse -c2 * 17.0 / 70.0}

TOL_U = 7.0e-8
TOL_T = 1.2e-5

[Functions]
  [uexact]
    type = ParsedFunction
    expression = 'if(y < 0.0, 0.0, if(y > 1.0, 0.0, 6.0*y*(1.0-y)))'
  []

  [t_solid_superior]
    type = ParsedFunction
    expression = '${c1}*(0.5*y^2-y*(1.0+${HP_H})+(0.5+${HP_H}))+${c4}+${c3}*x'
  []

  [t_solid_inferior]
    type = ParsedFunction
    expression = '${c1}*(0.5*y^2+y*${HP_H})+${c4}+${c3}*x'
  []

  [t_fluid]
    type = ParsedFunction
    expression = '${c2}*(y^4-2.0*y^3+y)+${c4}+${c3}*x'
  []

  [texact]
    type = ParsedFunction
    expression = 'if(y > 1.0, t_upper, if(y < 0.0, t_lower, t_center))'
    symbol_names = 't_upper t_lower t_center'
    symbol_values = 't_solid_superior t_solid_inferior t_fluid'
  []
[]

[Postprocessors]
  [uxerrl2]
    type = NekVolumeNorm
    field = velocity_x
    function = uexact
    execute_on = final
  []

  [terrl2]
    type = NekVolumeNorm
    field = temperature
    function = texact
    execute_on = final
  []

  [pass]
    type = ParsedPostprocessor
    expression = 'if(uxerrl2 < ${TOL_U} & terrl2 < ${TOL_T}, 1, 0)'
    pp_names = 'uxerrl2 terrl2'
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
