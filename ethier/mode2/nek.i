[Mesh]
  type = NekRSMesh
  volume = true
  exact = true
[]

[Problem]
  type = NekRSStandaloneProblem
  casename = 'ethier'
  output = 'velocity'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = final
[]

# reference solution is placed into auxvariables
[AuxVariables]
  [u]
  []
[]

[AuxKernels]
  [u]
    type = FunctionAux
    variable = u
    function = u
  []
[]

P_U0 = 0.5
P_V0 = 0.1
P_W0 = 0.2
a = ${fparse pi/4/100.}
nu = ${fparse 1/100}
d = ${fparse pi/2.}

[Functions]
  [u]
    type = ParsedFunction
    value = '-${a} * (ex * syz + ez * cxy) * exp(-${nu} * ${d} * ${d} * t) + ${P_U0}'
    symbol_names = 'ex syz ez cxy'
    symbol_values = 'ex syz ez cxy'
  []

  # we create some intermediate functions to make the formulas above a bit more
  # human-readable (and to match the intermediate variables used in the original
  # NekRS test)
  [xx]
    type = ParsedFunction
    value = 'x - ${P_U0} * t'
  []
  [yy]
    type = ParsedFunction
    value = 'y - ${P_V0} * t'
  []
  [zz]
    type = ParsedFunction
    value = 'z - ${P_W0} * t'
  []
  [ex]
    type = ParsedFunction
    value = 'exp(${a} * xx)'
    symbol_names = 'xx'
    symbol_values = 'xx'
  []
  [ey]
    type = ParsedFunction
    value = 'exp(${a} * yy)'
    symbol_names = 'yy'
    symbol_values = 'yy'
  []
  [ez]
    type = ParsedFunction
    value = 'exp(${a} * zz)'
    symbol_names = 'zz'
    symbol_values = 'zz'
  []
  [sxy]
    type = ParsedFunction
    value = 'sin(${a} * xx + ${d} * yy )'
    symbol_names = 'xx yy'
    symbol_values = 'xx yy'
  []
  [syz]
    type = ParsedFunction
    value = 'sin(${a} * yy + ${d} * zz )'
    symbol_names = 'zz yy'
    symbol_values = 'zz yy'
  []
  [szx]
    type = ParsedFunction
    value = 'sin(${a} * zz + ${d} * xx )'
    symbol_names = 'xx zz'
    symbol_values = 'xx zz'
  []
  [cxy]
    type = ParsedFunction
    value = 'cos(${a} * xx + ${d} * yy )'
    symbol_names = 'xx yy'
    symbol_values = 'xx yy'
  []
  [cyz]
    type = ParsedFunction
    value = 'cos(${a} * yy + ${d} * zz )'
    symbol_names = 'zz yy'
    symbol_values = 'zz yy'
  []
  [czx]
    type = ParsedFunction
    value = 'cos(${a} * zz + ${d} * xx )'
    symbol_names = 'xx zz'
    symbol_values = 'xx zz'
  []
[]

[Postprocessors]
  [uxerrl2]
    type = NodalL2Error
    variable = vel_x
    function = u
  []
[]
