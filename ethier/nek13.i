[Mesh]
  type = NekRSMesh
  volume = true
[]

[Problem]
  type = NekRSProblem
  casename = 'ethierScalar'
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

P_U0 = 0.5
P_V0 = 0.1
P_W0 = 0.2
P_A0 = 0.025
P_D0 = 0.5

a = ${fparse pi * P_A0}
nu = ${fparse 1.0 / 100.0}
d = ${fparse pi * P_D0}

TOL_S = 1.53E-05

[Functions]
  [xx]
    type = ParsedFunction
    expression = 'x - ${P_U0} * t'
  []

  [yy]
    type = ParsedFunction
    expression = 'y - ${P_V0} * t'
  []

  [zz]
    type = ParsedFunction
    expression = 'z - ${P_W0} * t'
  []

  [ex]
    type = ParsedFunction
    expression = 'exp(${a} * xx)'
    symbol_names = 'xx'
    symbol_values = 'xx'
  []

  [ez]
    type = ParsedFunction
    expression = 'exp(${a} * zz)'
    symbol_names = 'zz'
    symbol_values = 'zz'
  []

  [e2t]
    type = ParsedFunction
    expression = 'exp(-${nu} * ${d} * ${d} * t)'
  []

  [syz]
    type = ParsedFunction
    expression = 'sin(${a} * yy + ${d} * zz)'
    symbol_names = 'yy zz'
    symbol_values = 'yy zz'
  []

  [cxy]
    type = ParsedFunction
    expression = 'cos(${a} * xx + ${d} * yy)'
    symbol_names = 'xx yy'
    symbol_values = 'xx yy'
  []

  [uexact]
    type = ParsedFunction
    expression = '-${a} * (ex * syz + ez * cxy) * e2t + ${P_U0}'
    symbol_names = 'ex syz ez cxy e2t'
    symbol_values = 'ex syz ez cxy e2t'
  []
[]

[Postprocessors]
  [serrl2]
    type = NekVolumeNorm
    field = scalar02
    function = uexact
    execute_on = final
  []

  [pass]
    type = ParsedPostprocessor
    expression = 'if(serrl2 < ${TOL_S}, 1, 0)'
    pp_names = 'serrl2'
    execute_on = final
  []
[]
